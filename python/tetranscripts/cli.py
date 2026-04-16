#!/usr/bin/env python3
"""
TEcount - Measuring TE expression per-sample.
Core computation in Rust via PyO3, CLI powered by typer.
"""

import sys
import os
import logging
from typing import Optional

import typer


app = typer.Typer(
    name="TEcount",
    help="Measuring TE expression per-sample.",
    no_args_is_help=True,
    rich_markup_mode="rich",
)


@app.command()
def count(
    bam: str = typer.Option(..., "-b", "--BAM", help="An RNAseq BAM file."),
    gtffile: str = typer.Option(..., "--GTF", help="GTF file for gene annotations."),
    tefile: str = typer.Option(..., "--TE", help="GTF file for TE annotations."),
    fileformat: str = typer.Option(
        "BAM", "--format", help="Input file format: BAM or SAM."
    ),
    stranded: str = typer.Option(
        "no",
        "--stranded",
        help="Stranded library? (no, forward, or reverse).",
    ),
    te_mode: str = typer.Option(
        "multi", "--mode", help="How to count TE: uniq or multi."
    ),
    prefix: str = typer.Option(
        "TEcount_out", "--project", help="Name of this project."
    ),
    outdir: Optional[str] = typer.Option(
        None, "--outdir", help="Directory for output files. Default: current directory."
    ),
    sort_by_pos: bool = typer.Option(
        False,
        "--sortByPos",
        help="Alignment file is sorted by chromosome position.",
    ),
    iteration: int = typer.Option(
        100, "-i", "--iteration", help="Number of iterations for optimization."
    ),
    max_length: int = typer.Option(
        500, "--maxL", help="Maximum fragment length."
    ),
    min_length: int = typer.Option(
        0, "--minL", help="Minimum fragment length."
    ),
    frag_length: int = typer.Option(
        0, "-L", "--fragmentLength",
        help="Average fragment length for single-end reads. 0 = auto-detect.",
    ),
    verbose: int = typer.Option(
        2, "--verbose", help="Verbose level (0-3). Default: 2."
    ),
    version: bool = typer.Option(
        False, "--version", help="Show version and exit."
    ),
):
    """Count TE and gene expression from an RNA-seq BAM file."""
    if version:
        typer.echo("TEcount 0.1.0")
        raise typer.Exit()

    # --- validation ---
    if not os.path.isfile(bam):
        typer.echo(f"Error: No such file: {bam}", err=True)
        raise typer.Exit(1)
    if not os.path.isfile(gtffile):
        typer.echo(f"Error: No such file: {gtffile}", err=True)
        raise typer.Exit(1)
    if not os.path.isfile(tefile):
        typer.echo(f"Error: No such file: {tefile}", err=True)
        raise typer.Exit(1)

    if fileformat not in ("BAM", "SAM"):
        typer.echo(f"Error: Invalid format: {fileformat}", err=True)
        raise typer.Exit(1)
    if stranded not in ("no", "forward", "reverse"):
        typer.echo(f"Error: Invalid stranded value: {stranded}", err=True)
        raise typer.Exit(1)
    if te_mode not in ("uniq", "multi"):
        typer.echo(f"Error: Invalid TE mode: {te_mode}", err=True)
        raise typer.Exit(1)
    if iteration < 0:
        iteration = 0
    if frag_length < 0:
        typer.echo("Error: Fragment length cannot be negative.", err=True)
        raise typer.Exit(1)
    if min_length < 0:
        typer.echo("Error: Minimum fragment length cannot be negative.", err=True)
        raise typer.Exit(1)
    if max_length < 0:
        typer.echo("Error: Maximum fragment length cannot be negative.", err=True)
        raise typer.Exit(1)

    outdir_str = outdir if outdir else None
    if outdir_str is not None and not os.path.isdir(outdir_str):
        typer.echo(f"Error: Output directory ({outdir_str}) does not exist.", err=True)
        raise typer.Exit(1)

    # --- logging ---
    log_level = (4 - min(verbose, 3)) * 10
    logging.basicConfig(
        level=log_level,
        format="%(levelname)-5s @ %(asctime)s: %(message)s ",
        datefmt="%a, %d %b %Y %H:%M:%S",
        stream=sys.stderr,
    )
    info = logging.info

    argtxt = "\n".join((
        "# ARGUMENTS LIST:",
        f"# name = {prefix}",
        f"# BAM file = {bam}",
        f"# GTF file = {gtffile}",
        f"# TE file = {tefile}",
        f"# multi-mapper mode = {te_mode}",
        f"# stranded = {stranded}",
        f"# number of iterations = {iteration}",
        f"# Alignments grouped by read ID = {not sort_by_pos}",
    ))
    info("\n" + argtxt + "\n")

    # --- import Rust extension ---
    try:
        from tetranscripts._core import (
            GeneIndex,
            TEIndex,
            count_transcript_abundance as count_fn,
        )
    except ImportError:
        typer.echo(
            "Error: Rust extension module not found. "
            "Please run 'maturin develop' or 'pip install .' first.",
            err=True,
        )
        raise typer.Exit(1)

    # --- build indices ---
    info("Building gene index...")
    gene_idx = GeneIndex(gtffile, stranded, "exon", "gene_id")
    info("Done building gene index.")

    info("Building TE index...")
    te_idx = TEIndex(tefile)
    info("Done building TE index.")

    # --- count ---
    info("\nReading sample file...")
    result = count_fn(
        bam,
        gene_idx,
        te_idx,
        stranded,
        te_mode,
        sort_by_pos,
        iteration,
        frag_length,
        max_length,
    )
    info("Finished processing sample file.")

    # --- output ---
    if outdir_str is not None:
        os.chdir(outdir_str)

    _output_count_tbl(result, bam, prefix)
    info(f"Output written to {prefix}.cntTable")


def _output_count_tbl(result, bam_path: str, prefix: str):
    """Write the count table to a .cntTable file."""
    fname = f"{prefix}.cntTable"
    try:
        f = open(fname, "w")
    except IOError:
        sys.stderr.write(f"Cannot create report file {fname}!\n")
        sys.exit(1)

    f.write(f"gene/TE\t{bam_path}\n")

    all_keys = sorted(set(result.gene_counts.keys()) | set(result.te_element_counts.keys()))
    for key in all_keys:
        val = 0
        if key in result.gene_counts:
            val = int(result.gene_counts[key])
        elif key in result.te_element_counts:
            val = int(result.te_element_counts[key])
        f.write(f"{key}\t{val}\n")

    f.close()


if __name__ == "__main__":
    try:
        app()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!\n")
        sys.exit(0)
