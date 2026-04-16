#!/usr/bin/env python3
"""
TEcount - Measuring TE expression per-sample.
Ported from TEtranscripts bin/TEcount, with core computation in Rust.

Usage:
    TEcount -b RNAseq.bam --GTF gene_annotation.gtf --TE TE_annotation.gtf --sortByPos --mode multi
"""

import sys
import os
import argparse
import logging


def prepare_parser():
    desc = "Measuring TE expression per-sample."
    exmp = "Example: TEcount -b RNAseq.bam --GTF gene_annotation.gtf --TE TE_annotation.gtf --sortByPos --mode multi"

    parser = argparse.ArgumentParser(prog='TEcount', description=desc, epilog=exmp)

    parser.add_argument('-b', '--BAM', metavar='RNAseq.bam', dest='bam', required=True,
                        help='An RNAseq BAM file.')
    parser.add_argument('--GTF', metavar='genic-GTF-file', dest='gtffile', type=str, required=True,
                        help='GTF file for gene annotations')
    parser.add_argument('--TE', metavar='TE-GTF-file', dest='tefile', type=str, required=True,
                        help='GTF file for transposable element annotations')
    parser.add_argument('--format', metavar='input file format', dest='fileformat', type=str,
                        default='BAM', choices=['BAM', 'SAM'],
                        help='Input file format: BAM or SAM. DEFAULT: BAM')
    parser.add_argument('--stranded', metavar='option', dest='stranded', type=str, default="no",
                        choices=['no', 'forward', 'reverse'],
                        help='Is this a stranded library? (no, forward, or reverse). DEFAULT: no.')
    parser.add_argument('--mode', metavar='TE counting mode', dest='te_mode', type=str, default='multi',
                        choices=['uniq', 'multi'],
                        help='How to count TE: uniq or multi. DEFAULT: multi')
    parser.add_argument('--project', metavar='name', dest='prefix', default='TEcount_out',
                        help='Name of this project. DEFAULT: TEcount_out')
    parser.add_argument('--outdir', metavar='directory', dest='outdir', nargs='?', default='NULL',
                        help='Directory for output files. DEFAULT: current directory')
    parser.add_argument('--sortByPos', dest='sortByPos', action="store_true",
                        help='Alignment file is sorted by chromosome position.')
    parser.add_argument('-i', '--iteration', metavar='iteration', dest='numItr', type=int, default=100,
                        help='Number of iterations for optimization. DEFAULT: 100')
    parser.add_argument('--maxL', metavar='maxL', dest='maxL', type=int, default=500,
                        help='Maximum fragment length. DEFAULT: 500')
    parser.add_argument('--minL', metavar='minL', dest='minL', type=int, default=0,
                        help='Minimum fragment length. DEFAULT: 0')
    parser.add_argument('-L', '--fragmentLength', metavar='fragLength', dest='fragLength', type=int, default=0,
                        help='Average fragment length for single end reads. DEFAULT: 0 (auto-detect)')
    parser.add_argument('--verbose', metavar='verbose', dest='verbose', type=int, nargs='?', default=2,
                        const=3,
                        help='Set verbose level (0-3). DEFAULT: 2')
    parser.add_argument('--version', action='version', version='%(prog)s 0.1.0')

    return parser


def parse_args(parser):
    args = parser.parse_args()

    # Validate input files
    if not os.path.isfile(args.bam):
        logging.error("No such file: %s", args.bam)
        sys.exit(1)
    if not os.path.isfile(args.gtffile):
        logging.error("No such file: %s", args.gtffile)
        sys.exit(1)
    if not os.path.isfile(args.tefile):
        logging.error("No such file: %s", args.tefile)
        sys.exit(1)

    # Validate stranded
    if args.stranded not in ['forward', 'no', 'reverse']:
        logging.error("Invalid stranded value: %s", args.stranded)
        sys.exit(1)

    # Validate TE mode
    if args.te_mode not in ['uniq', 'multi']:
        logging.error("Invalid TE mode: %s", args.te_mode)
        sys.exit(1)

    # Validate numeric parameters
    if args.numItr < 0:
        args.numItr = 0
    if args.fragLength < 0:
        logging.error("Fragment length cannot be negative.")
        sys.exit(1)
    if args.minL < 0:
        logging.error("Minimum fragment length cannot be negative.")
        sys.exit(1)
    if args.maxL < 0:
        logging.error("Maximum fragment length cannot be negative.")
        sys.exit(1)

    # Validate output directory
    if args.outdir != "NULL":
        if not os.path.isdir(args.outdir):
            logging.error("Output directory (%s) does not exist.", args.outdir)
            sys.exit(1)

    # Setup logging
    logging.basicConfig(
        level=(4 - args.verbose) * 10,
        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
        datefmt='%a, %d %b %Y %H:%M:%S',
        stream=sys.stderr,
        filemode="w",
    )

    args.sortByPos = bool(args.sortByPos)

    args.argtxt = "\n".join((
        "# ARGUMENTS LIST:",
        "# name = %s" % args.prefix,
        "# BAM file = %s" % args.bam,
        "# GTF file = %s" % args.gtffile,
        "# TE file = %s" % args.tefile,
        "# multi-mapper mode = %s" % args.te_mode,
        "# stranded = %s" % args.stranded,
        "# number of iterations = %d" % args.numItr,
        "# Alignments grouped by read ID = %s" % (not args.sortByPos),
    ))

    return args


def output_count_tbl(result, bam_path, prefix):
    """Output count table to file."""
    fname = "{}.cntTable".format(prefix)
    try:
        f = open(fname, 'w')
    except IOError:
        sys.stderr.write("Cannot create report file {}!\n".format(fname))
        sys.exit(1)

    header = "gene/TE\t{}".format(bam_path)
    f.write(header + "\n")

    # Merge gene and TE counts
    all_keys = set(result.gene_counts.keys()) | set(result.te_element_counts.keys())

    for gene in sorted(all_keys):
        val = 0
        if gene in result.gene_counts:
            val = int(result.gene_counts[gene])
        elif gene in result.te_element_counts:
            val = int(result.te_element_counts[gene])
        f.write("{}\t{}\n".format(gene, val))

    f.close()


def main():
    """Main entry point for TEcount."""
    parser = prepare_parser()
    args = parse_args(parser)

    info = logging.info
    info("\n" + args.argtxt + "\n")

    # Import Rust extension
    try:
        from tetranscripts._core import GeneIndex, TEIndex, count_transcript_abundance as count_fn
    except ImportError:
        sys.stderr.write(
            "Error: Rust extension module not found. "
            "Please run 'maturin develop' or 'pip install .' first.\n"
        )
        sys.exit(1)

    # Build gene index
    info("Building gene index...")
    gene_idx = GeneIndex(args.gtffile, args.stranded, "exon", "gene_id")
    info("Done building gene index.")

    # Build TE index
    info("Building TE index...")
    te_idx = TEIndex(args.tefile)
    info("Done building TE index.")

    # Count transcript abundance
    info("\nReading sample file...")
    result = count_fn(
        args.bam,
        gene_idx,
        te_idx,
        args.stranded,
        args.te_mode,
        args.sortByPos,
        args.numItr,
        args.fragLength,
        args.maxL,
    )
    info("Finished processing sample file.")

    # Change to output directory if specified
    if args.outdir != "NULL":
        os.chdir(args.outdir)

    output_count_tbl(result, args.bam, args.prefix)
    info("Output written to {}.cntTable".format(args.prefix))


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!\n")
        sys.exit(0)
