# nanoCT Pipeline

Simplified pipeline for single-cell nanoCT (nano-scale Chromatin CUT&Tag) analysis.

## Pipeline Structure

```
nanoct/
├── processing/
│   ├── nanoct_pipeline.py    # Main Python pipeline (scit)
│   └── nanoct_analysis.R     # Main R pipeline (Signac/Seurat)
│
├── integration/
│   └── scglue_integration.py # nanoCT + scRNA-seq integration
│
├── downstream/
│   └── nanoct_downstream.py  # Tracks, markers, differential
│
├── utilities/
│   ├── nanoct_utils.sh       # Format conversion (rds2h5ad, barcodes, pseudobulk)
│   ├── dump_seurat.R         # Seurat to binary dump
│   └── reconstruct_h5ad.py   # Binary to h5ad
│
└── README.md
```

## Quick Start

### 1. Process nanoCT data

```bash
# Python processing with scit
baird jobs submit --id nanoct-process --project scentinel/nanoct \
  --context '{"pipeline_repo": "https://github.com/ethanlewisbaird/pipelines.git", "script_path": "nanoct/processing/nanoct_pipeline.py"}' \
  --command "cd /data/ebaird/scentinel/nanoCT/20260522.nanoCT && conda run -n nanoCT python /data/ebaird/pipelines/nanoct/processing/nanoct_pipeline.py --fragments H3K27ac/fragments.tsv.gz --peaks H3K27ac/peaks.bed --output output.h5ad"

# R analysis with Signac/Seurat
baird jobs submit --id nanoct-r --project scentinel/nanoct \
  --context '{"pipeline_repo": "https://github.com/ethanlewisbaird/pipelines.git", "script_path": "nanoct/processing/nanoct_analysis.R"}' \
  --command "cd /data/ebaird/scentinel/nanoCT/20260522.nanoCT && conda run -n nanoCT_R Rscript /data/ebaird/pipelines/nanoct/processing/nanoct_analysis.R --mode peak_analysis"
```

### 2. Integrate with scRNA-seq

```bash
baird jobs submit --id nanoct-scglue --project scentinel/nanoct \
  --context '{"pipeline_repo": "https://github.com/ethanlewisbaird/pipelines.git", "script_path": "nanoct/integration/scglue_integration.py"}' \
  --command "cd /data/ebaird/scentinel/nanoCT/20260522.nanoCT && conda run -n nanoCT python /data/ebaird/pipelines/nanoct/integration/scglue_integration.py --chrom chromatin.h5ad --rna rna.h5ad --output scglue_output"
```

### 3. Downstream analysis

```bash
# Find markers
baird jobs submit --id nanoct-markers --project scentinel/nanoct \
  --context '{"pipeline_repo": "https://github.com/ethanlewisbaird/pipelines.git", "script_path": "nanoct/downstream/nanoct_downstream.py"}' \
  --command "cd /data/ebaird/scentinel/nanoCT/20260522.nanoCT && conda run -n nanoCT python /data/ebaird/pipelines/nanoct/downstream/nanoct_downstream.py --h5ad output.h5ad --mode markers --output markers"

# Differential analysis
baird jobs submit --id nanoct-diff --project scentinel/nanoct \
  --context '{"pipeline_repo": "https://github.com/ethanlewisbaird/pipelines.git", "script_path": "nanoct/downstream/nanoct_downstream.py"}' \
  --command "cd /data/ebaird/scentinel/nanoCT/20260522.nanoCT && conda run -n nanoCT python /data/ebaird/pipelines/nanoct/downstream/nanoct_downstream.py --h5ad output.h5ad --mode differential --cluster1 0 --cluster2 1 --output diff"
```

### 4. Utilities

```bash
# Convert Seurat RDS to h5ad
baird jobs submit --id nanoct-convert --project scentinel/nanoct \
  --context '{"pipeline_repo": "https://github.com/ethanlewisbaird/pipelines.git", "script_path": "nanoct/utilities/nanoct_utils.sh"}' \
  --command "bash /data/ebaird/pipelines/nanoct/utilities/nanoct_utils.sh rds2h5ad input.rds output.h5ad RNA"

# Export barcodes
baird jobs submit --id nanoct-barcodes --project scentinel/nanoct \
  --context '{"pipeline_repo": "https://github.com/ethanlewisbaird/pipelines.git", "script_path": "nanoct/utilities/nanoct_utils.sh"}' \
  --command "bash /data/ebaird/pipelines/nanoct/utilities/nanoct_utils.sh export_barcodes clusters.csv output_dir"
```

## Script Options

### nanoct_pipeline.py

```
--fragments    Fragments TSV.GZ file
--peaks        Peaks BED file
--output       Output h5ad file
--mode         Pipeline mode: full, qc, process, finalize (default: full)
--remove-pc    PC to remove for finalize mode
--resolution   Clustering resolution (default: 0.8)
```

### nanoct_analysis.R

```
--mode         Analysis mode: peak_analysis, wnn, svd_tune (default: peak_analysis)
```

### nanoct_downstream.py

```
--h5ad         Input h5ad file
--output       Output directory
--mode         Analysis mode: tracks, markers, differential
--cluster1     First cluster for differential
--cluster2     Second cluster for differential
```

## Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `SCIT_PATH` | Path to scit library | `/data/ebaird/scentinel/nanoCT/20260522.nanoCT/scit_src` |
| `NANOCT_DATA_DIR` | Base data directory | `/data/ebaird/scentinel/nanoCT/20260522.nanoCT` |
| `NANOCT_OUTPUT_DIR` | Output directory | `analysis_R_output` |

## Conda Environments

| Environment | Usage |
|-------------|-------|
| `nanoCT` | Python processing |
| `nanoCT_R` | R analysis |
| `deeptools` | BigWig generation |

## Typical Workflow

```bash
# 1. Process with scit
python nanoct_pipeline.py --fragments frags.tsv.gz --peaks peaks.bed --output processed.h5ad

# 2. R analysis
Rscript nanoct_analysis.R --mode peak_analysis

# 3. Integrate with scRNA-seq
python scglue_integration.py --chrom processed.h5ad --rna rna.h5ad --output integrated/

# 4. Find markers
python nanoct_downstream.py --h5ad processed.h5ad --mode markers --output markers/
```
