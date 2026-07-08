#!/bin/bash
#SBATCH --job-name=demultiplex
#SBATCH --output=demultiplex.out
#SBATCH --error=demultiplex.err
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=72:00:00
#SBATCH --mem=128G

OUTPUT_DIR="/data/ebaird/scRNAseq/SCENTINELmar26/cellranger"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

### Non multiplexed samples
# /data/ebaird/miniconda3/envs/cellranger/opt/cellranger-9.0.1/cellranger count --id gal_12d --transcriptome /opt/ref/Drosophila_melanogaster/EnsemblBDGP6.46/dm6_BDGP_46_113_GFP_RFP --sample=2196_2024 --fastqs=/data/vtheodorou/SCENTINEL/SCENTINELSep24_10xGEX.IRB/Fastq_2223TF2NX/2196 --create-bam=true --localcores=40
# /data/ebaird/miniconda3/envs/cellranger/opt/cellranger-9.0.1/cellranger count --id flp_10d --transcriptome /opt/ref/Drosophila_melanogaster/EnsemblBDGP6.46/dm6_BDGP_46_113_GFP_RFP --sample=2197_2024 --fastqs=/data/vtheodorou/SCENTINEL/SCENTINELSep24_10xGEX.IRB/Fastq_2223TF2NX/2197 --create-bam=true --localcores=40
# /data/ebaird/miniconda3/envs/cellranger/opt/cellranger-9.0.1/cellranger count --id gal_10d --transcriptome /opt/ref/Drosophila_melanogaster/EnsemblBDGP6.46/dm6_BDGP_46_113_GFP_RFP --sample=2198_2024 --fastqs=/data/vtheodorou/SCENTINEL/SCENTINELSep24_10xGEX.IRB/Fastq_2223TF2NX/2198 --create-bam=true --localcores=40
# /data/ebaird/miniconda3/envs/cellranger/opt/cellranger-9.0.1/cellranger count --id flp_12d --transcriptome /opt/ref/Drosophila_melanogaster/EnsemblBDGP6.46/dm6_BDGP_46_113_GFP_RFP --sample=2199_2024 --fastqs=/data/vtheodorou/SCENTINEL/SCENTINELSep24_10xGEX.IRB/Fastq_2223TF2NX/2199 --create-bam=true --localcores=40

### Multiplexed samples CG000768_GEM-X_Universal_3_v4_4-plex_OCM kit
/data/ebaird/miniconda3/envs/cellranger/opt/cellranger-9.0.1/cellranger multi --id=SCentinelmar26 --csv=/data/ebaird/scRNAseq/SCENTINELmar26/raw/multi_config.csv
