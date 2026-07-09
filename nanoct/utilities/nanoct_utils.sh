#!/bin/bash
# nanoCT Utilities
# Format conversion and helper scripts
#
# Usage:
#   # Convert Seurat RDS to h5ad
#   bash nanoct_utils.sh rds2h5ad input.rds output.h5ad [assay]
#
#   # Export barcodes
#   bash nanoct_utils.sh export_barcodes clusters.csv output_dir
#
#   # Generate pseudobulk
#   bash nanoct_utils.sh pseudobulk mark output_dir

set -e

COMMAND=$1
shift

case $COMMAND in
    rds2h5ad)
        INPUT=$1
        OUTPUT=$2
        ASSAY=${3:-"RNA"}
        
        echo "Converting RDS to h5ad..."
        echo "  Input: $INPUT"
        echo "  Output: $OUTPUT"
        echo "  Assay: $ASSAY"
        
        # Create temp directory
        TMPDIR=$(mktemp -d)
        
        # Step 1: Dump Seurat to binary
        echo "  Dumping Seurat..."
        Rscript /data/ebaird/pipelines/nanoct/utilities/dump_seurat.R "$INPUT" "$TMPDIR" "$ASSAY"
        
        # Step 2: Reconstruct h5ad
        echo "  Reconstructing h5ad..."
        python /data/ebaird/pipelines/nanoct/utilities/reconstruct_h5ad.py "$TMPDIR" "$OUTPUT"
        
        # Cleanup
        rm -rf "$TMPDIR"
        
        echo "  Done: $OUTPUT"
        ;;
        
    export_barcodes)
        INPUT=$1
        OUTPUT_DIR=$2
        
        echo "Exporting barcodes..."
        echo "  Input: $INPUT"
        echo "  Output: $OUTPUT_DIR"
        
        mkdir -p "$OUTPUT_DIR"
        
        python -c "
import pandas as pd
import os

df = pd.read_csv('$INPUT')
for cluster in df['Cluster'].unique():
    bcs = df[df['Cluster'] == cluster]['Barcode']
    out = os.path.join('$OUTPUT_DIR', f'cluster_{cluster}.txt')
    bcs.to_csv(out, index=False, header=False)
    print(f'  Cluster {cluster}: {len(bcs)} barcodes')
"
        ;;
        
    pseudobulk)
        MARK=$1
        OUTPUT_DIR=$2
        BASE_DIR=${3:-"/data/ebaird/scentinel/nanoCT/20260522.nanoCT/SU.analysis.2026.05.22/Vasso_nanoCT_nanoscope"}
        
        echo "Generating pseudobulk for $MARK..."
        
        for BARCODE_FILE in "$OUTPUT_DIR"/cluster_*.txt; do
            CLUSTER=$(basename "$BARCODE_FILE" .txt | sed 's/cluster_//')
            echo "  Cluster $CLUSTER..."
            
            BAM="$BASE_DIR/$MARK/possorted_bam.bam"
            OUT_BAM="$OUTPUT_DIR/${MARK}_cluster_${CLUSTER}.bam"
            OUT_BW="$OUTPUT_DIR/${MARK}_cluster_${CLUSTER}.bw"
            
            # Extract cluster-specific reads
            if [ ! -f "$OUT_BAM" ]; then
                samtools view -@ 8 -D CB:"$BARCODE_FILE" -b -o "$OUT_BAM" "$BAM"
                samtools index "$OUT_BAM"
            fi
            
            # Convert to BigWig
            if [ ! -f "$OUT_BW" ]; then
                bamCoverage -b "$OUT_BAM" \
                    -o "$OUT_BW" \
                    --binSize 10 \
                    --normalizeUsing RPKM \
                    --numberOfProcessors 8 \
                    --minMappingQuality 30 \
                    --ignoreDuplicates
            fi
        done
        
        echo "  Done"
        ;;
        
    *)
        echo "Usage: nanoct_utils.sh <command> [args...]"
        echo ""
        echo "Commands:"
        echo "  rds2h5ad <input.rds> <output.h5ad> [assay]"
        echo "  export_barcodes <clusters.csv> <output_dir>"
        echo "  pseudobulk <mark> <output_dir> [base_dir]"
        exit 1
        ;;
esac
