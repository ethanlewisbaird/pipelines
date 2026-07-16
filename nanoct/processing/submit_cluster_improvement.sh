#!/bin/bash
#SBATCH --job-name=nanoct-cluster-improvement
#SBATCH --output=nanoct_cluster_%j.log
#SBATCH --error=nanoct_cluster_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=64G

# nanoCT Clustering Improvement Pipeline
# ATAC+DiffBind consensus peaks, dispersion-based variable peak selection,
# LSI+CCA dimensionality reduction, 3-strategy clustering comparison

MAIN_DIR="/data/ebaird/scentinel/nanoCT/20260522.nanoCT"
SCRIPT="/data/ebaird/pipelines/nanoct/processing/nanoct_combined.py"

source /data/ebaird/miniconda3/etc/profile.d/conda.sh
conda activate nanoCT

OUTPUT_DIR="${MAIN_DIR}/analysis_05.26"
mkdir -p "${OUTPUT_DIR}"
cp "$0" "${OUTPUT_DIR}/$(basename $0)"

cd "${MAIN_DIR}"
python "${SCRIPT}"

conda deactivate
echo "Done. Outputs written to: ${OUTPUT_DIR}"
