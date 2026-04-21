# TEtranscripts (Rust + Python)

Rust-accelerated TEtranscripts — quantifying transposable element (TE) expression from RNA-seq data.

The original [TEtranscripts](https://github.com/TEtranscripts/TEtranscripts) is implemented in pure Python. This version rewrites the performance-critical components (BAM parsing, interval overlap queries, EM iterative optimization) in Rust, compiled to a Python extension via [PyO3](https://pyo3.rs/) and [maturin](https://www.maturin.rs/). The command-line interface uses [typer](https://typer.tiangolo.com/).

提供两个子命令：

- **`TEcount count`** — 按 TE family/subfamily 聚合定量（对应原版 TEtranscripts）
- **`TEcount local`** — 按 TE locus（单个插入位点）定量（对应原版 TElocal）

## Performance

相比原版 Python 实现，核心计数循环加速 **5–10x**，主要得益于：

- **Rust BAM 解析** — 使用 `noodles`（纯 Rust，无 C 依赖）
- **扁平排序数组 + end-max 增广** — TE overlap 查询 O(log n + k)
- **中心化区间树** — 基因注释查询
- **Rayon 并行** — 多线程 overlap annotation
- **流式处理** — 分批读取 BAM（50K groups/batch），控制内存占用

## Installation

### 1. 创建 Conda 环境

```bash
conda create -n tetranscripts python=3.11 -y
conda activate tetranscripts
```

### 2. 安装依赖

```bash
# Rust 工具链（如果尚未安装）
conda install rust -c conda-forge

# maturin（Rust/Python 构建工具）
pip install maturin

# CLI 依赖
pip install typer

# 可选：饱和度测试所需
pip install numpy matplotlib

# samtools（用于 BAM 下采样和按 read name 排序）
conda install samtools -c bioconda
```

### 3. 编译 Rust 扩展

```bash
git clone https://github.com/xpf10/tetranscription_rs.git
cd tetranscription_rs
maturin develop --release
```

如果遇到 `Both VIRTUAL_ENV and CONDA_PREFIX are set` 错误：

```bash
unset CONDA_PREFIX
maturin develop --release
```

### 4. 验证安装

```bash
TEcount --help
```

应看到两个子命令 `count` 和 `local`。

## Usage

### TEcount count — Family-level TE 定量

按 TE family/subfamily 聚合报告表达量，与原版 TEtranscripts 兼容。

```bash
TEcount count -b sample.bam --GTF genes.gtf --TE te.gtf
```

#### 输出格式

生成 tab 分隔的 `.cntTable` 文件：

```
gene/TE                                 sample.bam
0610009B22Rik                           744
TP53                                    87
LINE1:L1:LINE                           523
ALU:Alu:SINE                            1042
```

其中 TE 行格式为 `gene_id:family_id:class_id`。

### TEcount local — Locus-level TE 定量

按单个 TE 插入位点报告表达量，与原版 TElocal 兼容。

```bash
TEcount local -b sample.bam --GTF genes.gtf --TE te.gtf
```

#### 输出格式

生成 tab 分隔的 `.cntTable` 文件，TE 行格式为 `chrom:start-end:transcript_id:gene_id:family_id:class_id:strand`：

```
gene/TE                                                          sample.bam
0610009B22Rik                                                    744
TP53                                                             87
chr10:10000098-10001958:L1_Mus1_dup17558:L1_Mus1:L1:LINE:-      12
chr1:5000-6500:AluY_dup1234:AluY:Alu:SINE:+                     3
```

### 命令行参数

`count` 和 `local` 共享相同的参数：

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-b`, `--BAM` | （必需） | RNA-seq BAM 文件 |
| `--GTF` | （必需） | 基因注释 GTF 文件 |
| `--TE` | （必需） | TE 注释 GTF 文件 |
| `--format` | `BAM` | 输入格式：`BAM` 或 `SAM` |
| `--stranded` | `no` | 链特异性：`no`、`forward` 或 `reverse` |
| `--mode` | `multi` | TE 计数模式：`uniq`（仅唯一比对）或 `multi`（含多比对，EM 分配） |
| `--project` | `TEcount_out` / `TElocal_out` | 输出文件前缀 |
| `--outdir` | 当前目录 | 输出目录 |
| `--sortByPos` | off | BAM 按坐标排序（会自动用 samtools 按 read name 重排） |
| `-i`, `--iteration` | `100` | EM 优化迭代次数 |
| `--maxL` | `500` | 最大 fragment 长度 |
| `--minL` | `0` | 最小 fragment 长度 |
| `-L`, `--fragmentLength` | `0` | 单端 read 的 fragment 长度（0 = 自动检测） |
| `--verbose` | `2` | 日志级别（0-3，数字越小越详细） |

### Examples

```bash
# Family-level 定量
TEcount count -b sample.bam --GTF genes.gtf --TE te.gtf

# Locus-level 定量
TEcount local -b sample.bam --GTF genes.gtf --TE te.gtf

# 链特异性库 + BAM 按坐标排序
TEcount count -b sample.bam --GTF genes.gtf --TE te.gtf \
    --stranded reverse --sortByPos

# 自定义输出路径
TEcount local -b sample.bam --GTF genes.gtf --TE te.gtf \
    --project my_sample --outdir results/

# 仅统计唯一比对的 TE
TEcount count -b sample.bam --GTF genes.gtf --TE te.gtf --mode uniq
```

### 饱和度测试

项目包含一个饱和度测试脚本，通过不同 downsampling 比例评估测序深度对检测的影响：

```bash
# TEcount 饱和度（family-level）
python saturation_test.py \
    -b test_data/sample.bam.sort \
    --GTF test_data/mm10.ncbiRefSeq_fix.gtf \
    --TE test_data/GRCm38_GENCODE_rmsk_TE.gtf \
    --mode count \
    --prefix M1

# TElocal 饱和度（locus-level）
python saturation_test.py \
    -b test_data/sample.bam.sort \
    --GTF test_data/mm10.ncbiRefSeq_fix.gtf \
    --TE test_data/GRCm38_GENCODE_rmsk_TE.gtf \
    --mode local \
    --prefix M1_local

# 同时运行两种模式并生成对比图
python saturation_test.py \
    -b test_data/sample.bam.sort \
    --GTF test_data/mm10.ncbiRefSeq_fix.gtf \
    --TE test_data/GRCm38_GENCODE_rmsk_TE.gtf \
    --mode both \
    --prefix M1
```

可选参数：

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--fractions` | `0.1 0.2 ... 1.0` | 采样比例列表 |
| `--seed` | `42` | 随机种子 |
| `--outdir` | `saturation_results` | 输出目录 |
| `-p`, `--prefix` | `sat` | 输出文件前缀 |
| `--mode` | `count` | `count`、`local` 或 `both` |

输出文件：
- `{prefix}_{N}pct.cntTable` — 每个比例的计数表
- `{prefix}_saturation.tsv` — 汇总统计表
- `{prefix}_saturation.pdf/png` — 饱和度曲线图

## Architecture

### 数据流

```
GTF 文件 → GTF 解析 → 索引构建 → 区间树查询
                                    ↓
BAM 文件 → 流式读取 → 按 read name 分组 → Rayon 并行 overlap annotation → 计数聚合
                                                                              ↓
                                                                    EM/SQUAREM 多比对分配
                                                                              ↓
                                                                        .cntTable 输出
```

### Rust 模块

| 模块 | 功能 |
|------|------|
| `lib.rs` | PyO3 模块入口，导出 `GeneIndex`、`TEIndex`、`count_transcript_abundance` |
| `types.rs` | 共享类型定义：`ExonInterval`、`Strand`、`CigarElement` |
| `gtf_parser.rs` | GTF 文件解析器（基因 + TE） |
| `interval_tree.rs` | 中心化区间树（基因注释查询） |
| `gene_index.rs` | 基因注释索引，每条染色体 3 棵树（正链/负链/无链） |
| `te_index.rs` | TE 注释索引，扁平排序数组 + end-max 增广 |
| `annotation.rs` | 注释解析与歧义消解 |
| `bam_reader.rs` | BAM 解析 + 流式分组 + Rayon 并行 + 计数循环 |
| `em_algorithm.rs` | SQUAREM 加速 EM 迭代优化 |

### Python 接口

```
python/tetranscripts/
├── __init__.py    # 导出 Rust 类和函数
└── cli.py         # Typer CLI（count + local 子命令）
```

PyO3 桥接暴露以下接口：

- `GeneIndex(gtf_path, stranded, feature_type, id_attribute)` — 基因索引
- `TEIndex(te_gtf_path)` — TE 索引（含 `get_locus_names()` 方法）
- `count_transcript_abundance(bam, gene_idx, te_idx, ...)` — 计数函数
- 返回 `PyCountResult`：`gene_counts`、`te_instance_counts`（locus-level）、`te_element_counts`（family-level）

### 项目结构

```
tetranscripts/
├── Cargo.toml              # Rust crate 配置
├── pyproject.toml           # Python 项目 + maturin 配置
├── src/                     # Rust 源码
│   ├── lib.rs
│   ├── types.rs
│   ├── gtf_parser.rs
│   ├── interval_tree.rs
│   ├── gene_index.rs
│   ├── te_index.rs
│   ├── annotation.rs
│   ├── bam_reader.rs
│   └── em_algorithm.rs
├── python/tetranscripts/
│   ├── __init__.py
│   └── cli.py               # CLI 入口
├── saturation_test.py       # 饱和度测试脚本
├── plot_saturation.py       # 饱和度绘图脚本
└── test_data/               # 测试数据（.gitignore）
```

## Input Files

### Gene GTF

标准 GTF 格式，需包含 `gene_id` 属性。例如使用 GENCODE 或 RefSeq 注释：

```
chr1    hg38_refSeq  exon  1000  2000  .  +  .  gene_id "0610009B22Rik";
```

### TE GTF

GTF 格式，需包含以下属性：`gene_id`、`transcript_id`、`family_id`、`class_id`。例如使用 RepeatMasker 注释：

```
chr1    rmsk  exon  5000  6500  .  +  .  gene_id "AluY"; transcript_id "AluY_dup1234"; family_id "Alu"; class_id "SINE";
```

BAM 文件需按 read name 排序。如果按坐标排序，使用 `--sortByPos` 参数，程序会自动调用 `samtools sort -n` 进行重排。

## Compatibility

- 命令行接口、输入格式与原版 TEtranscripts/TElocal 兼容
- `count` 输出的 `.cntTable` 与原版 TEtranscripts 格式一致
- `local` 输出的 TE 行为 locus-level 标识，与原版 TElocal 功能对应
- 现有的 gene GTF 和 TE GTF 文件可直接使用，无需修改

## License

Artistic License（与原版 TEtranscripts 相同）
