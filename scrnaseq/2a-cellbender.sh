#!/bin/bash
#SBATCH --job-name=cellbender
#SBATCH --output=cellbender.out
#SBATCH --error=cellbender.err
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=72:00:00
#SBATCH --mem=128G

conda activate cellbender

OUTPUT_DIR="/data/ebaird/scRNAseq/ProsRivsG4wRi/cellbender"
INPUT_DIR="/data/ebaird/scRNAseq/SCENTINELsep25/cellranger/scentinelsep25/outs/per_sample_outs"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

samples=("flp_10d" "flp_12d" "gal_10d" "gal_12d")

for sample in "${samples[@]}"; do
    cellbender remove-background \
        --input "$INPUT_DIR/$sample/outs/raw_feature_bc_matrix.h5" \
        --output "${sample}_CB_output.h5"
done

conda deactivate