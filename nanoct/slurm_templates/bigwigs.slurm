#!/bin/bash
#SBATCH --job-name=nanoct_bigwigs
#SBATCH --output=logs/nanoct_bigwigs_%j.out
#SBATCH --error=logs/nanoct_bigwigs_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --partition=long

source /data/ebaird/miniconda3/etc/profile.d/conda.sh
conda activate deeptools

python /data/ebaird/pipelines/nanoct/downstream/create_cluster_bigwigs.py
