# RNA-seq 差异表达与功能富集分析（拟南芥 DNRR）

## 项目概述

本项目使用拟南芥（*Arabidopsis thaliana*）离体叶片新生根再生（De novo Root Regeneration, DNRR）的公开 RNA-seq 数据，分析再生早期（0h → 12h）的转录组变化。通过 Snakemake 流程完成质控（fastp）、转录本定量（Salmon）、差异表达分析（edgeR）以及 GO 功能富集分析，揭示再生启动阶段的关键通路和核心转录特征。

## 数据来源

- **公开数据**：NCBI SRA (BioProject PRJNA613702, GEO GSE148292)
- **物种**：拟南芥 Col-0
- **组织**：离体叶片
- **样本**：
  - 对照组 (0h)：SRR11359743
  - 处理组 (12h)：SRR11359746
- **测序平台**：Illumina NovaSeq 6000 (双端 150 bp)

## 分析流程

原始数据（手动下载） → 质控修剪 (fastp) → 转录本定量 (Salmon) → 差异分析 (edgeR) → GO 富集分析

整个流程由 Snakemake 自动化管理，一键运行：

```bash
# 激活环境
micromamba activate rnaseq   # 或 conda activate rnaseq（命令兼容）

# 运行完整流程（自动处理所有依赖）
snakemake --cores 1 --resources mem=2000
```

各步骤详细说明：

- **fastp**：去除接头、poly‑G、低质量碱基（临时文件自动删除）。
- **salmon index**：基于 TAIR10 cDNA 构建索引（使用独立 conda 环境）。
- **salmon quant**：对每个样本进行定量，输出 `quant.sf`（使用独立 conda 环境）。
- **diff_expr**：使用 edgeR 无重复模式（固定 BCV=0.1）进行差异表达分析，生成火山图（颜色仅表示表达方向）、热图（Top 50 |logFC| 基因，log2 CPM，添加样本分组色条）。
- **GO enrichment**：基于绝对 logFC 最大的前 100 个上调/下调基因，利用超几何检验进行 GO 富集（仅 Biological Process），过滤 Count ≥3 的条目，按富集倍数排序，生成气泡图。

## 主要结果

- **差异基因数量**（基于 |log₂FC| > 1 筛选，无重复故 FDR 仅作参考）：
  上调基因：2,137 个，下调基因：1,139 个
- **火山图**：`results/diff_expr/volcano_plot.png`（Up/Down 仅表示表达方向，不表示统计显著性）
- **热图**：`results/diff_expr/heatmap_top50.png`
- **GO 富集气泡图**：`results/diff_expr/GO_up_bubble.png` / `GO_down_bubble.png`

**生物学结论**：上调基因富集于创伤响应、茉莉酸信号、防御反应等通路；下调基因富集于生长素信号、光响应、向性等过程，反映了离体叶片在再生早期从光合作用向防御/重编程转变的典型转录特征。

## 软件版本

- Snakemake 7.32.4
- fastp 1.3.3
- Salmon 1.10.2
- edgeR (无重复，固定 BCV=0.1)
- R 4.3.3
- 完整依赖见 `environment.yml`

## 技术难点与解决方案

1. **样本无生物学重复**
   仅 1 vs 1 样本，无法使用 DESeq2。
   **解决**：采用 edgeR 的 `exactTest`，手动设定生物学变异系数 `BCV = 0.1`，并在结果中如实说明局限性。
2. **缺少 tx2gene 映射文件**
   **解决**：从 TAIR10 cDNA FASTA 头部解析转录本 ID，去除版本号后作为基因 ID，自动构建映射表。
3. **计算资源有限（3GB 内存，23GB 硬盘）**
   **解决**：使用内存友好的 Salmon 定量；在 Snakefile 中将输出标记为 `temp()`，定量后自动删除；GO 富集仅使用前 100 个 |logFC| 基因，过滤 Count<3 的条目，避免内存溢出和假阳性。
4. **GO 富集出现宽泛或错误术语**
   **解决**：手动过滤 `nucleus`、`cytoplasm`、`response to stress` 等宽泛条目及动物特有术语（如 `egg-laying behavior`），并仅保留生物学过程（BP）。

## 方法局限性

- **缺乏生物学重复**：差异表达和富集分析的统计显著性不可靠，结果仅用于展示表达趋势，不可作为严格推断。
- **GO 富集手写实现**：未使用 clusterProfiler，但逻辑与标准方法一致，且额外增加了 Count 过滤和错误术语剔除。

## 如何复现

### 1. 克隆仓库

```bash
git clone https://github.com/zxitcast/rnaseq_heat_project_snakemake
cd snakemake_rnaseq
```

### 2. 安装依赖（需要 `micromamba` 或 `conda`）

```bash
# 创建主环境（包含 Snakemake, R, fastp 等）
micromamba env create -f environment.yml -n rnaseq
# 创建 salmon 专用环境（Snakemake 会自动切换）
micromamba env create -f envs/salmon.yaml -n salmon_only

# 激活主环境
micromamba activate rnaseq
```

> 注：如果你没有 `micromamba`，可以使用 `conda` 代替（命令完全相同）。或者先安装 `micromamba`（轻量、快速）。

### 3. 准备数据

- **安装下载工具** `aria2c`（如果尚未安装）：

  ```bash
  sudo apt install aria2   # Ubuntu/Debian
  # 或使用 conda: conda install -c conda-forge aria2
  ```

- **下载原始 FASTQ 数据**（约 11 GB）：

  ```bash
  bash scripts/download_data.sh
  ```

- **下载参考转录本、基因组和注释文件**（约 200 MB）：

  ```bash
  bash scripts/download_reference.sh
  ```

  你也可以手动将 FASTQ 文件软链接到 `data/raw/`，将参考转录本放入 `resources/genome/`。

### 4. 运行分析

```bash
snakemake --cores 1 --resources mem=2000
```

### 5. 查看结果

所有输出文件位于 `results/diff_expr/`。

## 目录结构

```text
snakemake_rnaseq/
├── LICENSE
├── README.md
├── environment.yml                 # 主环境（不含 salmon）
├── ath_kegg_pathways.csv           # KEGG 通路映射（未使用）
├── config/
│   └── config.yaml                 # 样本、路径、参数配置
├── data/
│   └── raw/                        # 原始 FASTQ（软链接或下载）
├── resources/
│   ├── genome/                     # 参考文件（软链接或下载）
│   │   ├── Arabidopsis_thaliana.TAIR10.cdna.all.fa
│   │   ├── Araport11_GTF_genes_transposons.20241001.gtf
│   │   └── TAIR10_chr_all.fasta
│   └── salmon_index/               # Salmon 索引（生成后存在）
├── results/
│   ├── diff_expr/                  # 差异表达与富集结果（CSV、图片）
│   └── quant/                      # Salmon 定量结果（含 quant.sf）
├── logs/                           # 日志文件（含 fastp HTML/JSON 报告）
├── scripts/
│   ├── download_data.sh            # 下载原始 FASTQ
│   ├── download_reference.sh       # 下载参考文件
│   ├── diff_expr.R                 # 差异表达和可视化
│   ├── 03_enrichment.R             # GO 富集分析
│   └── 04_kegg_enrichment.R        # KEGG 富集（可选）
├── envs/
│   └── salmon.yaml                 # 独立 salmon 环境文件
├── workflow/
│   └── Snakefile                   # Snakemake 主流程
└── Snakefile -> workflow/Snakefile # 软链接
```

## 许可证

MIT © 2026 zxitcast

## 参考文献

- Howard et al. (2013). PLOS ONE. (原始数据)
- Love, M.I., et al. (2014). Genome Biology. (DESeq2)
- Robinson, M.D., et al. (2010). Bioinformatics. (edgeR)
- Patro, R., et al. (2017). Nature Methods. (Salmon)
