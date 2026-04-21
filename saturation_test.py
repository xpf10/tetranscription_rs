#!/usr/bin/env python3
"""
Saturation test: run TEcount (count/local) at different downsampling fractions and plot results.
Directly calls the Rust extension (no subprocess) for speed.
"""

import os
import sys
import subprocess
import tempfile
import shutil
import argparse

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description="Saturation test for TEcount")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
    parser.add_argument("--GTF", required=True, help="Gene GTF")
    parser.add_argument("--TE", required=True, help="TE GTF")
    parser.add_argument("--fractions", nargs="+", type=float,
                        default=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
                        help="Sampling fractions")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--outdir", default="saturation_results", help="Output directory")
    parser.add_argument("-p", "--prefix", default="sat", help="Output prefix")
    parser.add_argument("--mode", choices=["count", "local", "both"], default="count",
                        help="TEcount mode: 'count' (family-level), 'local' (locus-level), or 'both'")
    args = parser.parse_args()

    args.bam = os.path.abspath(args.bam)
    args.GTF = os.path.abspath(args.GTF)
    args.TE = os.path.abspath(args.TE)
    args.outdir = os.path.abspath(args.outdir)
    os.makedirs(args.outdir, exist_ok=True)

    run_count = args.mode in ("count", "both")
    run_local = args.mode in ("local", "both")

    # Import Rust extension
    from tetranscripts._core import GeneIndex, TEIndex, count_transcript_abundance

    # Build indices once
    print("Building gene index...")
    gene_idx = GeneIndex(args.GTF, "no", "exon", "gene_id")
    print("Building TE index...")
    te_idx = TEIndex(args.TE)
    print(f"  Gene features: {len(gene_idx.features)}, TE instances: {te_idx.num_instances()}")

    if run_local:
        locus_names = te_idx.get_locus_names()
        print(f"  TE loci for local mode: {len(locus_names)}")

    tmpdir = tempfile.mkdtemp(prefix="tetranscripts_sat_")
    print(f"Temp dir: {tmpdir}")

    results_count = []
    results_local = []

    for frac in sorted(args.fractions):
        frac_pct = int(frac * 100)
        print(f"\n{'='*50}")
        print(f"Fraction: {frac_pct}%")
        print(f"{'='*50}")

        if frac >= 1.0:
            bam = args.bam
        else:
            # Downsample
            bam = os.path.join(tmpdir, f"downsampled_{frac_pct}pct.bam")
            samtools_frac = args.seed + frac
            cmd = ["samtools", "view", "-b", "-s", str(samtools_frac), "-o", bam, args.bam]
            print(f"  Downsampling: {' '.join(cmd)}")
            subprocess.run(cmd, check=True, capture_output=True)

        # Count reads
        n_reads = int(subprocess.run(
            ["samtools", "view", "-c", bam],
            capture_output=True, text=True
        ).stdout.strip())
        print(f"  Reads: {n_reads:,}")

        # Run TEcount directly via Rust extension
        print("  Running count_transcript_abundance...")
        result = count_transcript_abundance(
            bam, gene_idx, te_idx,
            "no", "multi", False, 100, 0, 500,
        )

        if run_count:
            total_gene = int(sum(result.gene_counts.values()))
            total_te = int(sum(result.te_element_counts.values()))
            total = total_gene + total_te

            # Save cntTable (family-level)
            cntfile = os.path.join(args.outdir, f"{args.prefix}_{frac_pct}pct.cntTable")
            _write_cnttable(cntfile, bam, result.gene_counts, result.te_element_counts)

            results_count.append({
                "fraction": frac,
                "label": f"{frac_pct}%",
                "reads": n_reads,
                "gene_counts": total_gene,
                "te_counts": total_te,
                "total_counts": total,
                "annotated": result.total_annotated,
                "nonunique": result.total_nonunique,
                "unannotated": result.total_unannotated,
            })
            print(f"  [count] Gene: {total_gene:,}  TE: {total_te:,}  Total: {total:,}")

        if run_local:
            total_gene = int(sum(result.gene_counts.values()))
            total_te = int(sum(result.te_instance_counts))
            total = total_gene + total_te

            # Build locus-level TE map
            te_locus_counts = dict(zip(locus_names, result.te_instance_counts))

            # Save cntTable (locus-level)
            cntfile = os.path.join(args.outdir, f"{args.prefix}_local_{frac_pct}pct.cntTable")
            _write_cnttable(cntfile, bam, result.gene_counts, te_locus_counts)

            # Count loci with non-zero expression
            n_detected = sum(1 for v in result.te_instance_counts if v > 0)

            results_local.append({
                "fraction": frac,
                "label": f"{frac_pct}%",
                "reads": n_reads,
                "gene_counts": total_gene,
                "te_counts": total_te,
                "total_counts": total,
                "n_loci_detected": n_detected,
                "n_loci_total": len(locus_names),
                "annotated": result.total_annotated,
                "nonunique": result.total_nonunique,
                "unannotated": result.total_unannotated,
            })
            print(f"  [local] Gene: {total_gene:,}  TE(loci): {total_te:,}  "
                  f"Detected: {n_detected:,}/{len(locus_names):,}")

    shutil.rmtree(tmpdir, ignore_errors=True)

    if run_count:
        _save_tsv(args.outdir, args.prefix, results_count)
        _plot_saturation(args.outdir, args.prefix, args.bam, results_count, "count")

    if run_local:
        _save_tsv(args.outdir, f"{args.prefix}_local", results_local,
                  extra_cols=["n_loci_detected", "n_loci_total"])
        _plot_saturation(args.outdir, f"{args.prefix}_local", args.bam, results_local,
                         "local", locus_info=results_local)

    if run_count and run_local:
        _plot_comparison(args.outdir, args.prefix, args.bam, results_count, results_local)


def _write_cnttable(path, bam, gene_counts, te_counts):
    all_keys = sorted(set(gene_counts.keys()) | set(te_counts.keys()))
    with open(path, "w") as f:
        f.write(f"gene/TE\t{bam}\n")
        for key in all_keys:
            val = int(gene_counts.get(key, te_counts.get(key, 0)))
            f.write(f"{key}\t{val}\n")


def _save_tsv(outdir, prefix, results, extra_cols=None):
    if not results:
        return
    tsv_path = os.path.join(outdir, f"{prefix}_saturation.tsv")
    base_cols = ["fraction", "reads", "gene_counts", "te_counts", "total_counts",
                 "annotated", "nonunique", "unannotated"]
    cols = base_cols + (extra_cols or [])
    with open(tsv_path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in results:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"\nResults saved to {tsv_path}")


def _plot_saturation(outdir, prefix, bam_path, results, mode, locus_info=None):
    if not results:
        return

    reads_m = [r["reads"] / 1e6 for r in results]
    gene_m = [r["gene_counts"] / 1e6 for r in results]
    te_10k = [r["te_counts"] / 1e4 for r in results]
    total_m = [r["total_counts"] / 1e6 for r in results]

    if locus_info:
        n_panels = 4
        fig, axes = plt.subplots(1, n_panels, figsize=(6 * n_panels, 5))

        n_detected = [r["n_loci_detected"] for r in locus_info]
        n_total = locus_info[0]["n_loci_total"]
        pct_detected = [d / n_total * 100 for d in n_detected]

        axes[3].plot(reads_m, pct_detected, "^-", color="#7570b3", linewidth=2, markersize=8)
        axes[3].set_xlabel("Total Reads (M)", fontsize=12)
        axes[3].set_ylabel("Detected Loci (%)", fontsize=12)
        axes[3].set_title("Locus Detection Rate", fontsize=13, fontweight="bold")
        axes[3].axhline(80, ls="--", color="gray", alpha=0.5, label="80%")
        axes[3].axhline(90, ls="--", color="gray", alpha=0.5, label="90%")
        axes[3].set_ylim(0, 105)
        axes[3].legend(fontsize=9)
        axes[3].grid(True, alpha=0.3)
    else:
        n_panels = 3
        fig, axes = plt.subplots(1, n_panels, figsize=(6 * n_panels, 5))

    axes[0].plot(reads_m, total_m, "o-", color="#2c7fb8", linewidth=2, markersize=8)
    axes[0].set_xlabel("Total Reads (M)", fontsize=12)
    axes[0].set_ylabel("Annotated Reads (M)", fontsize=12)
    axes[0].set_title("Total Annotated", fontsize=13, fontweight="bold")
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(reads_m, gene_m, "s-", color="#31a354", linewidth=2, markersize=8)
    axes[1].set_xlabel("Total Reads (M)", fontsize=12)
    axes[1].set_ylabel("Gene Counts (M)", fontsize=12)
    axes[1].set_title("Gene Expression", fontsize=13, fontweight="bold")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(reads_m, te_10k, "D-", color="#e34a33", linewidth=2, markersize=8)
    axes[2].set_xlabel("Total Reads (M)", fontsize=12)
    axes[2].set_ylabel("TE Counts (x10k)", fontsize=12)
    axes[2].set_title(f"TE Expression ({mode})", fontsize=13, fontweight="bold")
    axes[2].grid(True, alpha=0.3)

    mode_label = "TElocal (locus)" if mode == "local" else "TEcount (family)"
    plt.suptitle(f"Saturation Test ({mode_label}) — {os.path.basename(bam_path)}",
                 fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    for ext in ["pdf", "png"]:
        path = os.path.join(outdir, f"{prefix}_saturation.{ext}")
        plt.savefig(path, bbox_inches="tight", dpi=150)
        print(f"Plot saved to {path}")
    plt.close()


def _plot_comparison(outdir, prefix, bam_path, results_count, results_local):
    """Overlay count vs local saturation curves."""
    reads_m = [r["reads"] / 1e6 for r in results_count]
    te_count = [r["te_counts"] / 1e4 for r in results_count]
    te_local = [r["te_counts"] / 1e4 for r in results_local]

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.plot(reads_m, te_count, "D-", color="#e34a33", linewidth=2, markersize=8, label="TEcount (family)")
    ax.plot(reads_m, te_local, "D--", color="#7570b3", linewidth=2, markersize=8, label="TElocal (locus)")
    ax.set_xlabel("Total Reads (M)", fontsize=12)
    ax.set_ylabel("TE Counts (x10k)", fontsize=12)
    ax.set_title("TE Saturation: count vs local", fontsize=13, fontweight="bold")
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    plt.suptitle(f"{os.path.basename(bam_path)}", fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    for ext in ["pdf", "png"]:
        path = os.path.join(outdir, f"{prefix}_comparison.{ext}")
        plt.savefig(path, bbox_inches="tight", dpi=150)
        print(f"Comparison plot saved to {path}")
    plt.close()


if __name__ == "__main__":
    main()
