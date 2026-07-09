# nanoCT Pipeline

Single-cell nanoCT (nano-scale Chromatin CUT&Tag) analysis pipeline for Drosophila melanogaster.

## Overview

nanoCT is a method for dual-mark chromatin accessibility profiling using H3K27ac and H3K27me3 histone marks.

## Pipeline Structure

```
nanoct/
├── qc/                    # Quality control scripts
│   └── functions_scCT.R  # Helper functions
├── processing/            # Main processing scripts
│   ├── process_nanoct.py # Python processing with scit
│   └── nanoCT_analysis.R # R analysis with Signac/Seurat
├── integration/           # Multi-modal integration
├── pseudobulk/            # Pseudobulk generation
│   ├── export_barcodes.py
│   └── generate_pseudobulk.sh
└── downstream/            # Downstream analysis
```

## Usage

### Python Processing (scit)

```bash
# Set environment variables
export NANOCT_FRAGMENTS=/path/to/fragments.tsv.gz
export NANOCT_PEAKS=/path/to/peaks.bed
export NANOCT_OUTPUT=/path/to/output.h5ad
export NANOCT_MARK=H3K27ac
export SCIT_PATH=/path/to/scit_src

# Run via BAIRD
baird jobs submit --id nanoct-process --project scentinel/nanoct \
  --context '{"pipeline_repo": "https://github.com/ethanlewisbaird/pipelines.git", "script_path": "nanoct/processing/process_nanoct.py"}' \
  --command "cd /data/ebaird/scentinel/nanoCT/20260522.nanoCT && conda run -n nanoCT python /data/ebaird/pipelines/nanoct/processing/process_nanoct.py --fragments fragments.tsv.gz --peaks peaks.bed --output output.h5ad"
```

### R Analysis (Signac/Seurat)

```bash
# Set environment variables
export NANOCT_DATA_DIR=/data/ebaird/scentinel/nanoCT/20260522.nanoCT
export NANOCT_MARK=H3K27ac
export NANOCT_OUTPUT_DIR=/path/to/output
export NANOCT_SCRIPT_DIR=/data/ebaird/pipelines/nanoct/qc

# Run via BAIRD
baird jobs submit --id nanoct-analysis --project scentinel/nanoct \
  --context '{"pipeline_repo": "https://github.com/ethanlewisbaird/pipelines.git", "script_path": "nanoct/processing/nanoCT_analysis.R"}' \
  --command "cd /data/ebaird/scentinel/nanoCT/20260522.nanoCT && conda run -n nanoCT Rscript /data/ebaird/pipelines/nanoct/processing/nanoCT_analysis.R"
```

### Pseudobulk Generation

```bash
# Export barcodes
python /data/ebaird/pipelines/nanoct/pseudobulk/export_barcodes.py

# Generate pseudobulk
bash /data/ebaird/pipelines/nanoct/pseudobulk/generate_pseudobulk.sh H3K27ac
bash /data/ebaird/pipelines/nanoct/pseudobulk/generate_pseudobulk.sh H3K27me3
```

## Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `NANOCT_DATA_DIR` | Base data directory | `/data/ebaird/scentinel/nanoCT/20260522.nanoCT` |
| `NANOCT_MARK` | Histone mark | `H3K27ac` |
| `NANOCT_FRAGMENTS` | Fragments file path | - |
| `NANOCT_PEAKS` | Peaks BED file | - |
| `NANOCT_H5AD` | Processed h5ad file | - |
| `NANOCT_OUTPUT_DIR` | Output directory | `analysis_R_output` |
| `NANOCT_RESOLUTION` | Clustering resolution | `0.8` |
| `SCIT_PATH` | Path to scit library | - |

## Dependencies

- **Python**: scanpy, anndata, scit, matplotlib, pandas, scipy
- **R**: Signac, Seurat, GenomicRanges, biomaRt, BSgenome.Dmelanogaster.UCSC.dm6
- **Conda environment**: `nanoCT` on hibu

## Data Layout

```
/data/ebaird/scentinel/nanoCT/20260522.nanoCT/
├── H3K27ac/
│   ├── fragments.tsv.gz
│   ├── possorted_bam.bam
│   └── analysis/
├── H3K27me3/
│   ├── fragments.tsv.gz
│   ├── possorted_bam.bam
│   └── analysis/
├── analysis_05.26/
│   ├── nanoCT_workshop.ipynb
│   ├── combined_*.h5ad
│   └── ...
└── SU.analysis.2026.05.22/
    └── Vasso_nanoCT_nanoscope/
        ├── H3K27ac/
        ├── H3K27me3/
        └── pseudobulk/
```
