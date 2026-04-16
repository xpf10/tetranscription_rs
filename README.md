# TEtranscripts (Rust + Python)

Rust-accelerated TEtranscripts — quantifying transposable element (TE) expression from RNA-seq data.

The original [TEtranscripts](https://github.com/TEtranscripts/TEtranscripts) is implemented in pure Python. This version rewrites the performance-critical components (BAM parsing, interval overlap queries, EM iterative optimization) in Rust, compiled to a Python extension via [PyO3](https://pyo3.rs/) and [maturin](https://www.maturin.rs/). The command-line interface uses [typer](https://typer.tiangolo.com/).

## Performance

Compared to the original Python version, typical speedups are **5–10x** on the core counting loop, thanks to:

- **Rust BAM parsing** with `noodles` (pure Rust, no C dependencies)
- **BTreeMap-based TE index** with bin-based bucketing for O(log n) overlap queries
- **Center-based interval tree** for gene annotation lookups
- **Rayon parallelism** on overlap annotation across read groups

## Installation

**Prerequisites:** Rust toolchain (`rustup`), Python >= 3.8.

```bash
# Create a virtual environment (recommended)
python -m venv .venv
source .venv/bin/activate

# Install with maturin
pip install maturin
maturin develop --release
```

Or install as a package:

```bash
pip install .
```

## Usage

```
TEcount count [OPTIONS]
```

### Required Options

| Option | Description |
|--------|-------------|
| `-b`, `--BAM` | RNA-seq BAM file |
| `--GTF` | Gene annotation GTF file |
| `--TE` | TE annotation GTF file |

### Optional Options

| Option | Default | Description |
|--------|---------|-------------|
| `--format` | `BAM` | Input format: `BAM` or `SAM` |
| `--stranded` | `no` | Library strandedness: `no`, `forward`, or `reverse` |
| `--mode` | `multi` | TE counting mode: `uniq` or `multi` |
| `--project` | `TEcount_out` | Output file prefix |
| `--outdir` | current dir | Output directory |
| `--sortByPos` | off | BAM is sorted by position (will auto-sort by read name) |
| `-i`, `--iteration` | `100` | EM optimization iterations |
| `--maxL` | `500` | Maximum fragment length |
| `--minL` | `0` | Minimum fragment length |
| `-L`, `--fragmentLength` | `0` | Fragment length for single-end (0 = auto-detect) |
| `--verbose` | `2` | Verbose level (0–3) |

### Examples

```bash
# Basic usage
TEcount count -b sample.bam --GTF genes.gtf --TE te.gtf

# Stranded library, sorted by position
TEcount count -b sample.bam --GTF genes.gtf --TE te.gtf \
    --stranded reverse --sortByPos

# Custom output
TEcount count -b sample.bam --GTF genes.gtf --TE te.gtf \
    --project my_sample --outdir results/
```

### Output

Produces a tab-delimited `.cntTable` file:

```
gene/TE  sample.bam
BRCA1    142
TP53     87
LINE1:L1:LINE    523
ALU:Alu:SINE     1042
```

## Project Structure

```
tetranscripts/
├── Cargo.toml              # Rust crate config
├── pyproject.toml           # Python project + maturin config
├── src/                     # Rust source
│   ├── lib.rs               #   PyO3 module entry
│   ├── types.rs             #   Shared types (ExonInterval, Strand, CigarElement)
│   ├── gtf_parser.rs        #   GTF file parser
│   ├── interval_tree.rs     #   Center-based interval tree (gene queries)
│   ├── gene_index.rs        #   Gene annotation index
│   ├── te_index.rs          #   TE index with BTreeMap + bin bucketing
│   ├── annotation.rs        #   Annotation parsing & ambiguity resolution
│   ├── bam_reader.rs        #   BAM parsing + core counting loop (Rayon)
│   └── em_algorithm.rs      #   SQUAREM-accelerated EM optimization
├── python/tetranscripts/
│   ├── __init__.py
│   └── cli.py               # Typer CLI
└── tests/
```

## Compatibility

The command-line interface, input formats, and output format are compatible with the original TEtranscripts `TEcount` command. Existing gene GTF and TE GTF files can be used without modification.

## License

Artistic License (same as original TEtranscripts).
