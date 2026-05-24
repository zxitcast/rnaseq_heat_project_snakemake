#!/bin/bash
# =============================================
# download_reference.sh
# 功能: 下载拟南芥 TAIR10 参考基因组、注释和转录本
# 用法: bash scripts/download_reference.sh（可在任意目录执行，会自动定位项目根目录）
# 依赖: aria2c
# =============================================

set -e

# 自动切换到项目根目录
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_DIR"

mkdir -p resources/genome && cd resources/genome

echo "========================================="
echo "  下载拟南芥 TAIR10 参考数据"
echo "========================================="

# 1. 参考转录本 (cDNA) – 必须
echo "[1/3] 下载参考转录本 (cDNA)..."
aria2c -x 5 -s 5 -c --retry-wait=5 --max-tries=5 \
    -o Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz \
    "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-56/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz"
gunzip -f Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz

# 2. 参考基因组 (DNA) – 重命名为 TAIR10_chr_all.fasta
echo "[2/3] 下载参考基因组 (DNA)..."
aria2c -x 5 -s 5 -c --retry-wait=5 --max-tries=5 \
    -o TAIR10_chr_all.fasta.gz \
    "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-56/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
gunzip -f TAIR10_chr_all.fasta.gz

# 3. 注释文件 (GTF) – 重命名为 Araport11_GTF_genes_transposons.20241001.gtf
echo "[3/3] 下载注释文件 (GTF)..."
aria2c -x 5 -s 5 -c --retry-wait=5 --max-tries=5 \
    -o Araport11_GTF_genes_transposons.20241001.gtf.gz \
    "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-56/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.56.gtf.gz"
gunzip -f Araport11_GTF_genes_transposons.20241001.gtf.gz

# 删除压缩包（可选）
rm -f *.gz

echo "========================================="
echo "  所有参考文件已准备就绪！"
ls -lh
echo "========================================="
