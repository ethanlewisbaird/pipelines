#!/bin/bash
#SBATCH --job-name=soupx
#SBATCH --output=soupx.out
#SBATCH --error=soupx.err
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=72:00:00
#SBATCH --mem=128G

source /data/ebaird/miniconda3/etc/profile.d/conda.sh
conda activate soupx1

SAMPLES=("flp_10d","flp_12d","gal_10d","gal_12d")
OUTPUT_DIR="/data/ebaird/scRNAseq/SCENTINELsep24/soupx/"
INPUT_DIR="/data/ebaird/scRNAseq/SCENTINELsep24/cellranger/"

Rscript /data/ebaird/scRNAseq/SCENTINELsep24/code/soupx.r $SAMPLES $OUTPUT_DIR $INPUT_DIR

conda deactivate