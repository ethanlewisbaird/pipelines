#!/bin/bash
#SBATCH --job-name=nanoct_wnn
#SBATCH --output=logs/nanoct_wnn_%j.out
#SBATCH --error=logs/nanoct_wnn_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --partition=long

source /data/ebaird/miniconda3/etc/profile.d/conda.sh
conda activate nanoCT_R

Rscript /data/ebaird/pipelines/nanoct/processing/nanoCT_wnn_analysis.R
