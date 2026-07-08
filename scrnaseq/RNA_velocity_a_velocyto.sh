#!/bin/bash
#SBATCH --job-name=velocyto_run
#SBATCH --output=velocyto_run_.out
#SBATCH --error=velocyto_run.err
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --mem=240G

source /data/ebaird/miniconda3/etc/profile.d/conda.sh
conda activate velocyto

# Define input directory
INPUT_DIR="/data/ebaird/scRNAseq/SCENTINELsep25/cellranger/scentinelsep25/outs/per_sample_outs/"
GTF_FILE="/data/ebaird/scRNAseq/refs/genes.gtf"

# velocyto run -b "$INPUT_DIR/N+BreRi/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" -o "$INPUT_DIR/N+BreRi/velocyto" -m "/data/ebaird/scRNAseq/refs/repeat_mask/dm6.repeatmasker.sorted.gtf" "$INPUT_DIR/N+BreRi/outs/possorted_genome_bam.bam" "$GTF_FILE" --samtools-memory 2048 --samtools-threads 16
velocyto run -b "$INPUT_DIR/N+ProsRi/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" -o "$INPUT_DIR/N+ProsRi"/velocyto -m "/data/ebaird/scRNAseq/refs/repeat_mask/dm6.repeatmasker.sorted.gtf" "$INPUT_DIR/N+ProsRi/outs/possorted_genome_bam.bam" "$GTF_FILE" --samtools-memory 2048 --samtools-threads 16
# velocyto run -b "$INPUT_DIR/N+wRi/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" -o "$INPUT_DIR/N+wRi/velocyto" -m "/data/ebaird/scRNAseq/refs/repeat_mask/dm6.repeatmasker.sorted.gtf" "$INPUT_DIR/N+wRi/outs/possorted_genome_bam.bam" "$GTF_FILE" --samtools-memory 2048 --samtools-threads 16
velocyto run -b "$INPUT_DIR/ProsRi/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" -o "$INPUT_DIR/ProsRi"/velocyto -m "/data/ebaird/scRNAseq/refs/repeat_mask/dm6.repeatmasker.sorted.gtf" "$INPUT_DIR/ProsRi/outs/possorted_genome_bam.bam" "$GTF_FILE" --samtools-memory 2048 --samtools-threads 16

conda deactivate