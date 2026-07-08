#!/bin/bash
#SBATCH --job-name=miRNA_corr
#SBATCH --output=miRNA_corr.out
#SBATCH --error=miRNA_corr.err
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=128G

source /data/ebaird/miniconda3/etc/profile.d/conda.sh
conda activate R_process7

Rscript /data/ebaird/scRNAseq/SCENTINELsep24/code/miRNA_corr.r

conda deactivate