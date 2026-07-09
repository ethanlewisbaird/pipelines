#!/bin/bash
#SBATCH --job-name=macs3_peaks
#SBATCH --output=/data/ebaird/scRNAseq/20260522.nanoCT/R_analysis_peaks/output/macs3_peaks_%j.out
#SBATCH --error=/data/ebaird/scRNAseq/20260522.nanoCT/R_analysis_peaks/output/macs3_peaks_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=long

BASE="/data/ebaird/scRNAseq/20260522.nanoCT/SU.analysis.2026.05.22/Vasso_nanoCT_nanoscope"
OUT="/data/ebaird/scRNAseq/20260522.nanoCT/R_analysis_peaks/macs3_peaks"

mkdir -p $OUT

echo "=== MACS3 peak calling ==="
echo "Start: $(date)"

# H3K27ac -- sharp peaks (default)
echo ""
echo "--- H3K27ac (sharp peaks) ---"
macs3 callpeak \
    -t "$BASE/H3K27ac/possorted_bam.bam" \
    -f BAM \
    -g dm \
    --outdir "$OUT" \
    -n H3K27ac_macs3 \
    --keep-dup all \
    --nomodel \
    --extsize 200 \
    --shift -100 \
    2>&1 | tail -20

H3K27ac_COUNT=$(grep -v "^#" "$OUT/H3K27ac_macs3_peaks.narrowPeak" 2>/dev/null | wc -l)
echo "H3K27ac MACS3 peaks (before filtering): $H3K27ac_COUNT"

# H3K27me3 -- broad peaks
echo ""
echo "--- H3K27me3 (broad peaks) ---"
macs3 callpeak \
    -t "$BASE/H3K27me3/possorted_bam.bam" \
    -f BAM \
    -g dm \
    --outdir "$OUT" \
    -n H3K27me3_macs3 \
    --keep-dup all \
    --broad \
    --broad-cutoff 0.1 \
    --nomodel \
    --extsize 300 \
    2>&1 | tail -20

H3K27me3_COUNT=$(grep -v "^#" "$OUT/H3K27me3_macs3_peaks.broadPeak" 2>/dev/null | wc -l)
echo "H3K27me3 MACS3 peaks (before filtering): $H3K27me3_COUNT"

# Convert narrowPeak/broadPeak to simple BED (chr start end) for Signac
echo ""
echo "--- Converting to BED format ---"
if [ -f "$OUT/H3K27ac_macs3_peaks.narrowPeak" ]; then
    awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' "$OUT/H3K27ac_macs3_peaks.narrowPeak" | \
        sort -k1,1 -k2,2n | \
        awk 'BEGIN{OFS="\t"} !seen[$1"_"$2"_"$3]++' \
        > "$OUT/H3K27ac_macs3_peaks.bed"
    echo "H3K27ac BED: $(wc -l < $OUT/H3K27ac_macs3_peaks.bed) peaks"
fi

if [ -f "$OUT/H3K27me3_macs3_peaks.broadPeak" ]; then
    awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' "$OUT/H3K27me3_macs3_peaks.broadPeak" | \
        sort -k1,1 -k2,2n | \
        awk 'BEGIN{OFS="\t"} !seen[$1"_"$2"_"$3]++' \
        > "$OUT/H3K27me3_macs3_peaks.bed"
    echo "H3K27me3 BED: $(wc -l < $OUT/H3K27me3_macs3_peaks.bed) peaks"
fi

# Also create a combined BED file for comparison with CellRanger peaks
echo ""
echo "=== Summary ==="
echo "CellRanger H3K27ac peaks: $(grep -v "^#" "$BASE/H3K27ac/peaks.bed" | wc -l)"
echo "CellRanger H3K27me3 peaks: $(grep -v "^#" "$BASE/H3K27me3/peaks.bed" | wc -l)"
echo "MACS3 H3K27ac peaks: $H3K27ac_COUNT"
echo "MACS3 H3K27me3 peaks: $H3K27me3_COUNT"

echo ""
echo "Done: $(date)"
