#!/usr/bin/env python3
"""
RNA-seq saturation analysis — properly compute detected features at each depth.
Parses existing .cntTable files from saturation_test.py output.
"""

import os
import sys
import argparse
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


def parse_cnttable(path):
    """
    Parse a .cntTable file.
    Returns: n_detected_genes, n_detected_tes, total_gene_counts, total_te_counts
    A gene is any row without ':' in the name; a TE has ':' (e.g. "L1MdA:LINE:LINE").
    """
    n_genes = 0
    n_tes = 0
    total_gene = 0
    total_te = 0
    with open(path) as f:
        f.readline()  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            name = parts[0]
            count = float(parts[1])
            if ":" in name:
                total_te += count
                if count > 0:
                    n_tes += 1
            else:
                total_gene += count
                if count > 0:
                    n_genes += 1
    return n_genes, n_tes, total_gene, total_te


def main():
    parser = argparse.ArgumentParser(description="Plot RNA-seq saturation curves")
    parser.add_argument("--indir", default="saturation_results", help="Directory with cntTable files")
    parser.add_argument("--prefix", default="M1", help="Prefix used in saturation test")
    parser.add_argument("--tsv", default=None, help="TSV file from saturation_test.py (optional)")
    parser.add_argument("--outdir", default=None, help="Output directory (default: same as indir)")
    args = parser.parse_args()

    args.indir = os.path.abspath(args.indir)
    if args.outdir is None:
        args.outdir = args.indir
    else:
        args.outdir = os.path.abspath(args.outdir)
        os.makedirs(args.outdir, exist_ok=True)

    # Find cntTable files and parse
    fractions = []
    data = []

    for fname in sorted(os.listdir(args.indir)):
        if not fname.startswith(args.prefix) or not fname.endswith(".cntTable"):
            continue
        path = os.path.join(args.indir, fname)
        # Extract fraction from filename: M1_10pct.cntTable → 10
        m = re.search(r"(\d+)pct", fname)
        if not m:
            continue
        pct = int(m.group(1))
        frac = pct / 100.0

        n_genes, n_tes, total_gene, total_te = parse_cnttable(path)
        data.append({
            "frac": frac,
            "pct": pct,
            "n_genes": n_genes,
            "n_tes": n_tes,
            "total_gene": total_gene,
            "total_te": total_te,
        })
        fractions.append(frac)

    # Also read reads count from TSV if available
    tsv_path = args.tsv or os.path.join(args.indir, f"{args.prefix}_saturation.tsv")
    reads_map = {}
    if os.path.isfile(tsv_path):
        with open(tsv_path) as f:
            f.readline()  # header
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    frac = float(parts[0])
                    reads = int(parts[1])
                    reads_map[frac] = reads

    data.sort(key=lambda d: d["frac"])

    if not data:
        print("No cntTable files found.", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(data)} data points")
    print(f"{'Pct':>5} {'Reads':>12} {'Genes':>8} {'TEs':>8} {'Gene%':>8} {'TE%':>8}")
    print("-" * 60)
    for d in data:
        reads = reads_map.get(d["frac"], 0)
        total = d["total_gene"] + d["total_te"]
        gp = d["total_gene"] / total * 100 if total > 0 else 0
        tp = d["total_te"] / total * 100 if total > 0 else 0
        print(f"{d['pct']:>5}% {reads:>12,} {d['n_genes']:>8,} {d['n_tes']:>8,} {gp:>7.1f}% {tp:>7.1f}%")

    # =========================================================================
    # Plot
    # =========================================================================
    reads_m = np.array([reads_map.get(d["frac"], 0) / 1e6 for d in data])
    n_genes = np.array([d["n_genes"] for d in data])
    n_tes = np.array([d["n_tes"] for d in data])
    total_gene_m = np.array([d["total_gene"] / 1e6 for d in data])
    total_te_m = np.array([d["total_te"] / 1e6 for d in data])
    total_m = total_gene_m + total_te_m

    # Marginal gain: new features per additional million reads
    marg_genes = np.diff(n_genes) / np.diff(reads_m)
    marg_tes = np.diff(n_tes) / np.diff(reads_m)
    marg_reads_mid = (reads_m[:-1] + reads_m[1:]) / 2

    # Gene/TE proportion
    prop_gene = total_gene_m / total_m * 100
    prop_te = total_te_m / total_m * 100

    fig = plt.figure(figsize=(18, 14))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.35)

    # ---- (A) Detected genes vs reads ----
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(reads_m, n_genes / 1000, "o-", color="#2c7fb8", linewidth=2, markersize=7)
    ax1.set_xlabel("Total Reads (M)", fontsize=11)
    ax1.set_ylabel("Detected Genes (x1000)", fontsize=11)
    ax1.set_title("A. Gene Detection Saturation", fontsize=12, fontweight="bold")
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(bottom=0)

    # ---- (B) Detected TEs vs reads ----
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.plot(reads_m, n_tes, "D-", color="#e34a33", linewidth=2, markersize=7)
    ax2.set_xlabel("Total Reads (M)", fontsize=11)
    ax2.set_ylabel("Detected TE Elements", fontsize=11)
    ax2.set_title("B. TE Detection Saturation", fontsize=12, fontweight="bold")
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(bottom=0)

    # ---- (C) Marginal gain (new features per M reads) ----
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.plot(marg_reads_mid, marg_genes, "s-", color="#2c7fb8", linewidth=2, markersize=6, label="Genes / M reads")
    ax3.plot(marg_reads_mid, marg_tes, "^-", color="#e34a33", linewidth=2, markersize=6, label="TEs / M reads")
    ax3.set_xlabel("Total Reads (M)", fontsize=11)
    ax3.set_ylabel("Newly Detected Features per M Reads", fontsize=11)
    ax3.set_title("C. Marginal Detection Rate", fontsize=12, fontweight="bold")
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.3)
    ax3.set_ylim(bottom=0)

    # ---- (D) Saturation percentage ----
    max_genes = n_genes[-1]
    max_tes = n_tes[-1]
    sat_genes = n_genes / max_genes * 100
    sat_tes = n_tes / max_tes * 100

    ax4 = fig.add_subplot(gs[1, 0])
    ax4.plot(reads_m, sat_genes, "o-", color="#2c7fb8", linewidth=2, markersize=7, label="Genes")
    ax4.plot(reads_m, sat_tes, "D-", color="#e34a33", linewidth=2, markersize=7, label="TEs")
    ax4.axhline(y=80, color="gray", linestyle="--", alpha=0.5, linewidth=1)
    ax4.text(reads_m[0] + 0.5, 81, "80%", color="gray", fontsize=9)
    ax4.axhline(y=90, color="gray", linestyle="--", alpha=0.5, linewidth=1)
    ax4.text(reads_m[0] + 0.5, 91, "90%", color="gray", fontsize=9)
    ax4.set_xlabel("Total Reads (M)", fontsize=11)
    ax4.set_ylabel("Saturation (%)", fontsize=11)
    ax4.set_title("D. Detection Saturation %", fontsize=12, fontweight="bold")
    ax4.legend(fontsize=9)
    ax4.set_ylim(0, 105)
    ax4.grid(True, alpha=0.3)

    # ---- (E) Gene vs TE proportion ----
    ax5 = fig.add_subplot(gs[1, 1])
    ax5.fill_between(reads_m, 0, prop_gene, alpha=0.4, color="#2c7fb8", label="Gene reads")
    ax5.fill_between(reads_m, prop_gene, prop_gene + prop_te, alpha=0.4, color="#e34a33", label="TE reads")
    ax5.set_xlabel("Total Reads (M)", fontsize=11)
    ax5.set_ylabel("Proportion (%)", fontsize=11)
    ax5.set_title("E. Gene vs TE Read Proportion", fontsize=12, fontweight="bold")
    ax5.legend(fontsize=9)
    ax5.set_ylim(0, 100)
    ax5.grid(True, alpha=0.3)

    # ---- (F) Total annotated counts ----
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.plot(reads_m, total_gene_m, "s-", color="#2c7fb8", linewidth=2, markersize=6, label="Gene counts")
    ax6.plot(reads_m, total_te_m, "^-", color="#e34a33", linewidth=2, markersize=6, label="TE counts")
    ax6.plot(reads_m, total_m, "o-", color="#333333", linewidth=2, markersize=6, label="Total", alpha=0.7)
    ax6.set_xlabel("Total Reads (M)", fontsize=11)
    ax6.set_ylabel("Annotated Counts (M)", fontsize=11)
    ax6.set_title("F. Annotated Read Counts", fontsize=12, fontweight="bold")
    ax6.legend(fontsize=9)
    ax6.grid(True, alpha=0.3)
    ax6.set_ylim(bottom=0)

    plt.suptitle(
        f"RNA-seq Saturation Analysis — {args.prefix} ({os.path.basename(args.indir)})",
        fontsize=14, fontweight="bold", y=1.01
    )

    for ext in ["pdf", "png"]:
        path = os.path.join(args.outdir, f"{args.prefix}_saturation.{ext}")
        plt.savefig(path, bbox_inches="tight", dpi=150)
        print(f"Saved: {path}")
    plt.close()


if __name__ == "__main__":
    main()
