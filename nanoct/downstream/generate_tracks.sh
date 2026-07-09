#!/bin/bash
#SBATCH --job-name=bigwigs_5kb
#SBATCH --output=/data/ebaird/scRNAseq/20260522.nanoCT/browser_tracks/bigwigs_%j.out
#SBATCH --error=/data/ebaird/scRNAseq/20260522.nanoCT/browser_tracks/bigwigs_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --partition=long

source /data/ebaird/miniconda3/etc/profile.d/conda.sh
conda activate deeptools

BASE="/data/ebaird/scRNAseq/20260522.nanoCT/SU.analysis.2026.05.22/Vasso_nanoCT_nanoscope"
OUT="/data/ebaird/scRNAseq/20260522.nanoCT/browser_tracks"

for MARK in H3K27ac H3K27me3; do
    BAM="${BASE}/pseudobulk/${MARK}/${MARK}_all_cells.bam"
    
    echo "=== ${MARK} ==="
    
    # 1. Copy existing 10bp bigwig for reference (peak-based approach)
    echo "Copying existing 10bp track..."
    cp "${BASE}/pseudobulk/${MARK}/${MARK}_all_cells.bw" "${OUT}/${MARK}_peaks_10bp.bw"
    
    # 2. Generate 5kb bin bigwig (smoother, matches scit binning approach)
    echo "Generating 5kb bin track..."
    bamCoverage \
        -b "$BAM" \
        -o "${OUT}/${MARK}_5kb_bins.bw" \
        --binSize 5000 \
        --normalizeUsing RPKM \
        --numberOfProcessors 4 \
        --minMappingQuality 30 \
        --ignoreDuplicates
    
    echo "Done with ${MARK}"
done

echo "=== All tracks generated ==="
ls -lh "${OUT}/"