#!/bin/bash
#SBATCH --job-name=rnahybrid
#SBATCH --output=rnahybrid.out
#SBATCH --error=rnahybrid.err
#SBATCH --time=500:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=120G

# Load conda environment
source /data/ebaird/miniconda3/etc/profile.d/conda.sh
conda activate rnahybrid

cd /data/ebaird/scRNAseq/SCENTINELsep24/miRNA

# Input files
MIRNA=mir2279.fa
UTR=dmel-all-three_prime_UTR-r6.54.fasta
RAW_OUT=mir2279_results.txt
FILTERED_OUT=filtered_mir2279_hits.txt

echo "Running RNAhybrid..."
RNAhybrid -s 3utr_fly -t $UTR -q $MIRNA -b 1 -e -20 > "$RAW_OUT" -p 0.05 -g jpg

echo "Filtering strong binders (ΔG < -22 kcal/mol)..."
awk '$3 < -22' "$RAW_OUT" > "$FILTERED_OUT"

conda deactivate
