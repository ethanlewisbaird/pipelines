#!/bin/bash
#SBATCH --job-name=osconet
#SBATCH --output=osconet.out
#SBATCH --error=osconet.err
#SBATCH --time=500:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=120G

source /data/ebaird/miniconda3/etc/profile.d/conda.sh
conda activate fullosconet38

vanilla Rscript /data/ebaird/scRNAseq/SCENTINELsep24/FullOscoNet/Code/1-pre-osconet.r
vanilla Rscript /data/ebaird/scRNAseq/SCENTINELsep24/FullOscoNet/Code/2-NormFilter.R
sbatch /data/ebaird/scRNAseq/SCENTINELsep24/FullOscoNet/Code/3-Run_All.sh
vanilla Rscript /data/ebaird/scRNAseq/SCENTINELsep24/FullOscoNet/Code/4-CommunityExtract.R
vanilla Rscript /data/ebaird/scRNAseq/SCENTINELsep24/FullOscoNet/Code/5-visualisation_osconet.R

conda deactivate