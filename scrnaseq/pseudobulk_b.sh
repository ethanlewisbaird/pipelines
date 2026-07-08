#!/bin/bash
#SBATCH --job-name=bam_pseudobulk
#SBATCH --output=bam_pseudobulk.out
#SBATCH --error=bam_pseudobulk.err
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=72:00:00
#SBATCH --mem=128G

MAIN_DIR="/data/ebaird/scRNAseq/ProsRivsG4wRi"
CR_DATA_DIR="/data/ebaird/scRNAseq/SCENTINELsep25/cellranger/scentinelsep25/outs/per_sample_outs"
OUT_DIR=$MAIN_DIR"/pseudobulk"
mkdir -p $OUT_DIR

GENOTYPES=("ProsRi" "N+wRi" "N+ProsRi")
SAMPLES=("ProsRi_10d" "N+wRi_13d" "N+ProsRi_10d.13d.pld")

for geno in "${GENOTYPES[@]}"; do
  SAMPLE_FILES=$(ls $OUT_DIR/${geno}_*.txt 2>/dev/null)
  
  if [ -z "$SAMPLE_FILES" ]; then
    echo "No samples found for $geno"
    continue
  fi

  for barcode_file in $SAMPLE_FILES; do
    samp_name=$(basename $barcode_file .txt | cut -d_ -f2-)
    sample_bam=$CR_DATA_DIR"/${geno}/outs/possorted_genome_bam.bam"
    output_bam="$OUT_DIR/${samp_name}.bam"
    # Extract by CB tag using -D CB:
    samtools view -@ 4 -D CB:$barcode_file -b -o $output_bam $sample_bam
    samtools index $output_bam
  done

  ### Uncomment to Merge and index samples of same genotype or timepoint
  # merged_bam="$OUT_DIR/merged_${geno}.bam"
  # samtools merge -@ 8 $merged_bam $OUT_DIR/${geno}_*.bam
  # samtools index $merged_bam
done

### Generate BigWig files from BAM files

### Merged by genotype
# for geno in "${GENOTYPES[@]}"; do
#   bam="/data/ebaird/scRNAseqreports/res/Gal10d_Gal12d_Flp10d_Flp12d_070525/pseudobulk_output/merged_${geno}.bam" 
#   /data/ebaird/miniconda3/envs/deeptools/bin/bamCoverage -b $bam \
#     -o $OUT_DIR/${geno}.bw \
#     --binSize 50 \
#     --normalizeUsing RPKM \
#     --numberOfProcessors 8
# done

### Individual samples
for sample in "${SAMPLES[@]}"; do
  bam=$OUT_DIR"/*${sample}.bam" 
  /data/ebaird/miniconda3/envs/deeptools/bin/bamCoverage -b $bam \
    -o $OUT_DIR/${sample}.bw \
    --binSize 50 \
    --normalizeUsing RPKM \
    --numberOfProcessors 8

done

