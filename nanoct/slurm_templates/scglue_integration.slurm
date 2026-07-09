#!/bin/bash
#SBATCH --job-name=nanoct_scglue
#SBATCH --output=logs/nanoct_scglue_%j.out
#SBATCH --error=logs/nanoct_scglue_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --partition=long

source /data/ebaird/miniconda3/etc/profile.d/conda.sh
conda activate nanoCT

python /data/ebaird/pipelines/nanoct/integration/scglue_integration.py
