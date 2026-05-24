#!/bin/bash
# =============================================
# download_data.sh
# 功能: 从 ENA 下载拟南芥 RNA-seq 原始 FASTQ 文件
# 用法: bash scripts/download_data.sh（可在任意目录执行，会自动定位项目根目录）
# 依赖: aria2c
# =============================================

set -e

# 自动切换到项目根目录（假设脚本位于 scripts/ 下）
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_DIR"

mkdir -p data/raw && cd data/raw

echo "===== 下载 SRR11359743 (0h) ====="
aria2c -x 4 -s 4 -c https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR113/043/SRR11359743/SRR11359743_1.fastq.gz
aria2c -x 4 -s 4 -c https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR113/043/SRR11359743/SRR11359743_2.fastq.gz

echo "===== 下载 SRR11359746 (12h) ====="
aria2c -x 4 -s 4 -c https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR113/046/SRR11359746/SRR11359746_1.fastq.gz
aria2c -x 4 -s 4 -c https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR113/046/SRR11359746/SRR11359746_2.fastq.gz

echo "===== 下载完成！文件位于 data/raw/ ====="
