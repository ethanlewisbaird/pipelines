# nanoCT Pipeline

Complete pipeline for single-cell nanoCT (nano-scale Chromatin CUT&Tag) analysis in Drosophila melanogaster.

## Overview

nanoCT profiles dual histone marks (H3K27ac + H3K27me3) at single-cell resolution. This pipeline handles the full analysis from raw fragments to publication-ready figures.

## Pipeline Structure

```
nanoct/
├── preprocessing/           # Data preparation and peak calling
│   ├── macs3_callpeak.sh   # MACS3 peak calling
│   └── fragment_qc.py      # Fragment quality control
│
├── processing/              # Core analysis
│   ├── process_nanoct.py   # Python: scit-based processing
│   ├── nanoCT_analysis.R   # R: Signac/Seurat analysis
│   ├── nanoCT_peak_analysis.R      # Peak-based analysis
│   ├── nanoCT_wnn_analysis.R       # WNN multi-modal integration
│   ├── nanoCT_svd_tuning.R         # SVD dimension tuning
│   ├── nanoCT_multires_clustering.R # Multi-resolution clustering
│   ├── run_peaks_pipeline_proper.py # Full scit pipeline
│   └── run_to_qc.py                # QC processing
│
├── integration/             # Multi-modal integration
│   ├── scglue_integration.py    # scGLUE: nanoCT + scRNA-seq
│   ├── scglue_peaks.py          # scGLUE with peaks
│   ├── harmony_patch.py         # Harmony batch correction
│   └── genelevel_integration.py # Gene-level integration
│
├── downstream/              # Downstream analysis
│   ├── create_cluster_bigwigs.py    # Per-cluster BigWig tracks
│   ├── create_merged_rpkm_bigwigs.py # Merged RPKM BigWigs
│   ├── generate_tracks.sh           # Browser track generation
│   ├── deeptools_heatmaps.py        # DeepTools heatmaps
│   ├── chromatin_marker_differential.py # Marker + differential
│   ├── generate_markers.py          # Marker gene analysis
│   └── plot_markers_differential.py # Marker visualization
│
├── visualization/           # Plotting and exploration
│   ├── explore_clustering.py    # Cluster exploration
│   ├── make_qc_plots.py         # QC visualization
│   ├── plot_per_cluster_umap.py # Per-cluster UMAP
│   └── resolution_sweep.py      # Resolution parameter sweep
│
├── utilities/               # Helper scripts
│   ├── dump_seurat.R       # Seurat to binary dump
│   ├── reconstruct_h5ad.py # Binary dump to h5ad
│   ├── rds2h5ad.sh         # RDS to h5ad conversion
│   ├── export_barcodes.py  # Barcode export
│   └── generate_pseudobulk.sh # Pseudobulk generation
│
├── qc/                      # Quality control functions
│   └── functions_scCT.R    # Helper functions
│
└── README.md               # This file
```

## Quick Start

### 1. Preprocessing: Peak Calling

```bash
# Call peaks with MACS3
baird jobs submit --id nanoct-peaks --project scentinel/nanoct \
  --context '{"pipeline_repo": "https://github.com/ethanlewisbaird/pipelines.git", "script_path": "nanoct/preprocessing/macs3_callpeak.sh"}' \
  --command "cd /data/ebaird/scentinel/nanoCT/20260522.nanoCT && bash /data/ebaird/pipelines/nanoct/preprocessing/macs3_callpeak.sh"
```

### 2. Processing: Core Analysis

```bash
# R-based analysis with Signac/Seurat
baird jobs submit --id nanoct-process --project scentinel/nanoct \
  --context '{"pipeline_repo": "https://github.com/ethanlewisbaird/pipelines.git", "script_path": "nanoct/processing/nanoCT_peak_analysis.R"}' \
  --command "cd /data/ebaird/scentinel/nanoCT/20260522.nanoCT && conda run -n nanoCT_R Rscript /data/ebaird/pipelines/nanoct/processing/nanoCT_peak_analysis.R"

# WNN multi-modal integration
baird jobs submit --id nanoct-wnn --project scentinel/nanoct \
  --context '{"pipeline_repo": "https://github.com/ethanlewisbaird/pipelines.git", "script_path": "nanoct/processing/nanoCT_wnn_analysis.R"}' \
  --command "cd /data/ebaird/scentinel/nanoCT/20260522.nanoCT && conda run -n nanoCT_R Rscript /data/ebaird/pipelines/nanoct/processing/nanoCT_wnn_analysis.R"
```

### 3. Integration: scGLUE

```bash
# Integrate nanoCT with scRNA-seq
baird jobs submit --id nanoct-scglue --project scentinel/nanoct \
  --context '{"pipeline_repo": "https://github.com/ethanlewisbaird/pipelines.git", "script_path": "nanoct/integration/scglue_integration.py"}' \
  --command "cd /data/ebaird/scentinel/nanoCT/20260522.nanoCT && conda run -n nanoCT python /data/ebaird/pipelines/nanoct/integration/scglue_integration.py"
```

### 4. Downstream: Tracks and Markers

```bash
# Generate BigWig tracks
baird jobs submit --id nanoct-bigwigs --project scentinel/nanoct \
  --context '{"pipeline_repo": "https://github.com/ethanlewisbaird/pipelines.git", "script_path": "nanoct/downstream/create_cluster_bigwigs.py"}' \
  --command "cd /data/ebaird/scentinel/nanoCT/20260522.nanoCT && conda run -n nanoCT python /data/ebaird/pipelines/nanoct/downstream/create_cluster_bigwigs.py"

# Find markers and differential peaks
baird jobs submit --id nanoct-markers --project scentinel/nanoct \
  --context '{"pipeline_repo": "https://github.com/ethanlewisbaird/pipelines.git", "script_path": "nanoct/downstream/chromatin_marker_differential.py"}' \
  --command "cd /data/ebaird/scentinel/nanoCT/20260522.nanoCT && conda run -n nanoCT python /data/ebaird/pipelines/nanoct/downstream/chromatin_marker_differential.py"
```

### 5. Utilities: Format Conversion

```bash
# Convert Seurat RDS to h5ad
baird jobs submit --id nanoct-rds2h5ad --project scentinel/nanoct \
  --context '{"pipeline_repo": "https://github.com/ethanlewisbaird/pipelines.git", "script_path": "nanoct/utilities/rds2h5ad.sh"}' \
  --command "cd /data/ebaird/scentinel/nanoCT/20260522.nanoCT && bash /data/ebaird/pipelines/nanoct/utilities/rds2h5ad.sh input.rds output.h5ad"
```

## Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `NANOCT_DATA_DIR` | Base data directory | `/data/ebaird/scentinel/nanoCT/20260522.nanoCT` |
| `NANOCT_MARK` | Histone mark | `H3K27ac` |
| `NANOCT_OUTPUT_DIR` | Output directory | `analysis_R_output` |
| `NANOCT_RESOLUTION` | Clustering resolution | `0.8` |
| `NANOCT_N_PCS` | Number of PCs | `50` |
| `NANOCT_MIN_CELLS` | Min cells per gene | `10` |
| `NANOCT_MIN_GENES` | Min genes per cell | `200` |
| `SCIT_PATH` | Path to scit library | `/data/ebaird/scentinel/nanoCT/20260522.nanoCT/scit_src` |

## Conda Environments

| Environment | Usage | Key Packages |
|-------------|-------|--------------|
| `nanoCT` | Python processing | scanpy, anndata, scit, scglue |
| `nanoCT_R` | R analysis | Signac, Seurat, GenomicRanges |
| `deeptools` | BigWig generation | deepTools, pyBigWig |

## Data Layout

```
/data/ebaird/scentinel/nanoCT/20260522.nanoCT/
├── H3K27ac/
│   ├── fragments.tsv.gz
│   ├── possorted_bam.bam
│   ├── peaks.bed
│   └── analysis/
├── H3K27me3/
│   ├── fragments.tsv.gz
│   ├── possorted_bam.bam
│   ├── peaks.bed
│   └── analysis/
├── analysis_05.26/
│   ├── combined_*.h5ad
│   └── ...
├── R_analysis_peaks/
│   ├── output/
│   └── ...
└── SU.analysis.2026.05.22/
    └── Vasso_nanoCT_nanoscope/
        ├── H3K27ac/
        ├── H3K27me3/
        └── pseudobulk/
```

## Typical Workflow

1. **Peak Calling** → `preprocessing/macs3_callpeak.sh`
2. **QC Processing** → `processing/run_to_qc.py`
3. **Core Analysis** → `processing/nanoCT_peak_analysis.R`
4. **SVD Tuning** → `processing/nanoCT_svd_tuning.R`
5. **Clustering** → `processing/nanoCT_multires_clustering.R`
6. **WNN Integration** → `processing/nanoCT_wnn_analysis.R`
7. **scGLUE Integration** → `integration/scglue_integration.py`
8. **BigWig Tracks** → `downstream/create_cluster_bigwigs.py`
9. **Markers** → `downstream/chromatin_marker_differential.py`
10. **Visualization** → `visualization/plot_per_cluster_umap.py`

## Troubleshooting

### scit library not found
```bash
export SCIT_PATH=/data/ebaird/scentinel/nanoCT/20260522.nanoCT/scit_src
```

### Memory issues
```bash
# Increase R memory
options(future.globals.maxSize = 15 * 1024^3)  # 15GB
```

### Missing annotations
```bash
# Download dm6 annotations
wget https://ftp.ensembl.org/pub/release-110/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.46.110.gtf.gz
```
