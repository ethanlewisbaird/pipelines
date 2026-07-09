#!/bin/bash

MARK=$1 # H3K27ac or H3K27me3
BASE_DIR="/data/ebaird/scRNAseq/20260522.nanoCT/SU.analysis.2026.05.22/Vasso_nanoCT_nanoscope"
BAM_FILE="${BASE_DIR}/${MARK}/possorted_bam.bam"
OUT_DIR="${BASE_DIR}/pseudobulk/${MARK}"
DEEPTOOLS_ENV="deeptools"

echo "Processing $MARK..."

for barcode_file in ${OUT_DIR}/${MARK}_cluster_*.txt; do
    cluster_id=$(basename $barcode_file .txt | sed "s/${MARK}_cluster_//")
    echo "  Cluster $cluster_id..."
    
    output_bam="${OUT_DIR}/${MARK}_cluster_${cluster_id}.bam"
    output_bw="${OUT_DIR}/${MARK}_cluster_${cluster_id}.bw"
    
    # 1. Extract cluster-specific reads
    if [ ! -f "$output_bam" ]; then
        echo "    Extracting BAM..."
        samtools view -@ 8 -D CB:${barcode_file} -b -o ${output_bam} ${BAM_FILE}
        samtools index ${output_bam}
    else
        echo "    BAM already exists, skipping extraction."
    fi
    
    # 2. Convert to BigWig
    if [ ! -f "$output_bw" ]; then
        echo "    Generating BigWig..."
        conda run -n ${DEEPTOOLS_ENV} bamCoverage -b ${output_bam} \
            -o ${output_bw} \
            --binSize 10 \
            --normalizeUsing RPKM \
            --numberOfProcessors 8 \
            --minMappingQuality 30 \
            --ignoreDuplicates
    else
        echo "    BigWig already exists, skipping conversion."
    fi
done

echo "Done with $MARK."
