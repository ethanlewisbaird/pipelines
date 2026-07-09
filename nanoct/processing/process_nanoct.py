#!/usr/bin/env python3
"""
nanoCT Processing Pipeline with scit

Main processing script for single-cell nanoCT data.
Handles: fragments → AnnData → QC → PCA → clustering → UMAP

Usage:
    python process_nanoct.py --fragments <fragments.tsv.gz> --peaks <peaks.bed> --output <output.h5ad>

Environment variables (set via BAIRD):
    NANOCT_FRAGMENTS - path to fragments file
    NANOCT_PEAKS - path to peaks BED file
    NANOCT_OUTPUT - output h5ad path
    NANOCT_MARK - histone mark (H3K27ac or H3K27me3)
    NANOCT_MIN_CELLS - minimum cells per gene (default: 10)
    NANOCT_MIN_GENES - minimum genes per cell (default: 200)
"""

import argparse
import os
import sys
from datetime import datetime

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

# Add scit to path if needed
scit_path = os.environ.get("SCIT_PATH", "")
if scit_path:
    sys.path.insert(0, scit_path)

try:
    import src as scit
except ImportError:
    try:
        import scit
    except ImportError:
        print("ERROR: scit library not found. Set SCIT_PATH environment variable.")
        sys.exit(1)


def parse_args():
    parser = argparse.ArgumentParser(description="Process nanoCT data with scit")
    parser.add_argument("--fragments", required=True, help="Fragments TSV.GZ file")
    parser.add_argument("--peaks", required=True, help="Peaks BED file")
    parser.add_argument("--output", required=True, help="Output h5ad file")
    parser.add_argument("--mark", default="H3K27ac", help="Histone mark name")
    parser.add_argument("--min-cells", type=int, default=10, help="Min cells per gene")
    parser.add_argument("--min-genes", type=int, default=200, help="Min genes per cell")
    parser.add_argument("--n-pcs", type=int, default=50, help="Number of PCs")
    parser.add_argument("--resolution", type=float, default=0.8, help="Clustering resolution")
    parser.add_argument("--output-dir", default="", help="Output directory for plots")
    return parser.parse_args()


def load_fragments(fragments_path, peaks_path):
    """Load fragments and peaks into AnnData."""
    print(f"Loading fragments from: {fragments_path}")
    print(f"Loading peaks from: {peaks_path}")
    
    # Read peaks
    peaks = pd.read_csv(peaks_path, sep="\t", header=None, 
                        names=["chrom", "start", "end"])
    print(f"  Loaded {len(peaks)} peaks")
    
    # Use scit to load fragments
    adata = scit.pp.read_fragments(
        fragments_path,
        peaks,
        file_type="tsv"
    )
    
    print(f"  Created AnnData: {adata.shape[0]} cells x {adata.shape[1]} peaks")
    return adata


def quality_control(adata, min_cells=10, min_genes=200, output_dir=""):
    """Run quality control filtering."""
    print("Running quality control...")
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    
    # Filter
    print(f"  Before filtering: {adata.shape[0]} cells x {adata.shape[1]} genes")
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    print(f"  After filtering: {adata.shape[0]} cells x {adata.shape[1]} genes")
    
    # Plot QC metrics
    if output_dir:
        fig, axes = plt.subplots(1, 3, figsize=(12, 4))
        sc.pl.violin(adata, "n_genes_by_counts", ax=axes[0], show=False)
        sc.pl.violin(adata, "total_counts", ax=axes[1], show=False)
        axes[2].hist(adata.obs["n_genes_by_counts"], bins=50)
        axes[2].set_xlabel("Genes per cell")
        axes[2].set_ylabel("Count")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "qc_metrics.png"), dpi=150)
        plt.close()
        print(f"  Saved QC plot to: {output_dir}/qc_metrics.png")
    
    return adata


def process_adata(adata, n_pcs=50, resolution=0.8, output_dir=""):
    """Process AnnData: normalize → PCA → neighbors → clusters → UMAP."""
    print("Processing AnnData...")
    
    # Normalize
    print("  Normalizing...")
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    
    # Variable genes
    print("  Finding variable genes...")
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    
    # PCA
    print(f"  Running PCA ({n_pcs} components)...")
    sc.tl.pca(adata, n_comps=n_pcs)
    
    # Neighbors
    print("  Computing neighbors...")
    sc.pp.neighbors(adata, n_pcs=n_pcs)
    
    # Clustering
    print(f"  Clustering (resolution={resolution})...")
    sc.tl.leiden(adata, resolution=resolution)
    n_clusters = len(adata.obs["leiden"].unique())
    print(f"  Found {n_clusters} clusters")
    
    # UMAP
    print("  Computing UMAP...")
    sc.tl.umap(adata)
    
    # Plot UMAP
    if output_dir:
        fig, ax = plt.subplots(figsize=(8, 8))
        sc.pl.umap(adata, color="leiden", ax=ax, show=False, 
                   title=f"Clusters (res={resolution})")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "umap_clusters.png"), dpi=150)
        plt.close()
        print(f"  Saved UMAP to: {output_dir}/umap_clusters.png")
    
    return adata


def main():
    args = parse_args()
    
    # Create output directory
    output_dir = args.output_dir or os.path.dirname(args.output)
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 60)
    print("nanoCT Processing Pipeline")
    print("=" * 60)
    print(f"Start time: {datetime.now()}")
    print(f"Fragments: {args.fragments}")
    print(f"Peaks: {args.peaks}")
    print(f"Mark: {args.mark}")
    print(f"Output: {args.output}")
    print()
    
    # Load data
    adata = load_fragments(args.fragments, args.peaks)
    adata.obs["mark"] = args.mark
    
    # Quality control
    adata = quality_control(adata, args.min_cells, args.min_genes, output_dir)
    
    # Process
    adata = process_adata(adata, args.n_pcs, args.resolution, output_dir)
    
    # Save
    print(f"Saving to: {args.output}")
    adata.write(args.output)
    
    print()
    print("=" * 60)
    print("Pipeline complete!")
    print(f"End time: {datetime.now()}")
    print(f"Output: {args.output}")
    print(f"Cells: {adata.shape[0]}")
    print(f"Peaks: {adata.shape[1]}")
    print(f"Clusters: {len(adata.obs['leiden'].unique())}")
    print("=" * 60)


if __name__ == "__main__":
    main()
