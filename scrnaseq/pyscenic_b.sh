#!/bin/bash
#SBATCH --job-name=pyscenic
#SBATCH --output=pyscenic.out
#SBATCH --error=pyscenic.err
#SBATCH --time=500:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=120G

source /data/ebaird/miniconda3/etc/profile.d/conda.sh
conda activate pyscenic
MAIN_DIR="/data/ebaird/scRNAseq/SCENTINELsep24/pyscenic"
OUTPUT_DIR="${MAIN_DIR}/pyscenic_res_5000_3"
mkdir -p "$OUTPUT_DIR"

# Set arguments
LOOM_FILE="/data/ebaird/scRNAseqreports/res/Gal10d_Gal12d_Flp10d_Flp12d_070525/t5000_Gal10d_Gal12d_Flp10d_Flp12d.loom"
TF_LIST="/data/ebaird/scRNAseq/SCENTINELsep24/refs/allTFs_dmel.txt"
DB_GLOB="/data/ebaird/scRNAseq/SCENTINELsep24/refs/dm6_v10_clust.genes_vs_motifs.rankings.feather"
NTASKS=24
ANNOTATIONS="/data/ebaird/scRNAseq/SCENTINELsep24/refs/motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl"

# Run the Python script with the arguments
python /data/ebaird/scRNAseq/SCENTINELsep24/code/pyscenic.py --output_dir "$OUTPUT_DIR" --loom_file "$LOOM_FILE" --tf_list "$TF_LIST" --db_glob "$DB_GLOB" --n_tasks "$NTASKS" --annotations "$ANNOTATIONS"

conda deactivate