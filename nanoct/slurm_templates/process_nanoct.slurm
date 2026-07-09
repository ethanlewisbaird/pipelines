#!/bin/bash
#SBATCH --job-name=nanoct_process
#SBATCH --output=logs/nanoct_process_%j.out
#SBATCH --error=logs/nanoct_process_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --partition=long

source /data/ebaird/miniconda3/etc/profile.d/conda.sh
conda activate nanoCT

export SCIT_PATH=/data/ebaird/scentinel/nanoCT/20260522.nanoCT/scit_src
export NANOCT_DATA_DIR=/data/ebaird/scentinel/nanoCT/20260522.nanoCT

python /data/ebaird/pipelines/nanoct/processing/process_nanoct.py \
    --fragments ${NANOCT_DATA_DIR}/H3K27ac/fragments.tsv.gz \
    --peaks ${NANOCT_DATA_DIR}/H3K27ac/peaks.bed \
    --output ${NANOCT_DATA_DIR}/analysis_05.26/nanoct_processed.h5ad \
    --mark H3K27ac
