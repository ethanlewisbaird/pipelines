#!/bin/bash
#SBATCH --job-name=sctenifold
#SBATCH --output=sctenifold.out
#SBATCH --error=sctenifold.err
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=128G

source /data/ebaird/miniconda3/etc/profile.d/conda.sh
conda activate sctenifold

python /data/ebaird/scRNAseq/pipeline_code/sctenifold.py

conda deactivate