#!/usr/bin/env python3
"""
nanoCT Processing Pipeline

Main pipeline for single-cell nanoCT analysis using scit library.
Consolidates: QC → processing → clustering → UMAP → peak analysis

Usage:
    # Full pipeline
    python nanoct_pipeline.py --fragments frags.tsv.gz --peaks peaks.bed --output output.h5ad
    
    # QC only (for inspection)
    python nanoct_pipeline.py --fragments frags.tsv.gz --peaks peaks.bed --output output.h5ad --mode qc
    
    # Finalize (after SVD tuning)
    python nanoct_pipeline.py --fragments frags.tsv.gz --peaks peaks.bed --output output.h5ad --mode finalize --remove-pc 0

Environment variables:
    SCIT_PATH - path to scit library (default: /data/ebaird/scentinel/nanoCT/20260522.nanoCT/scit_src)
"""

import argparse
import os
import sys
from datetime import datetime

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import anndata as ad
import scipy.sparse as sps

# Add scit to path
scit_path = os.environ.get("SCIT_PATH", "/data/ebaird/scentinel/nanoCT/20260522.nanoCT/scit_src")
sys.path.insert(0, scit_path)
import src as scit

# Drosophila dm6 chromosome sizes
CHROM_SIZES = {
    '2L': 23513712, '2R': 25286936, '3L': 28110227, '3R': 32079331,
    '4': 1348131, 'X': 23542271, 'Y': 3667352,
    'mitochondrion_genome': 19524, 'rDNA': 76973,
}


def parse_args():
    parser = argparse.ArgumentParser(description="nanoCT Processing Pipeline")
    parser.add_argument("--fragments", required=True, help="Fragments TSV.GZ file")
    parser.add_argument("--peaks", required=True, help="Peaks BED file")
    parser.add_argument("--output", required=True, help="Output h5ad file")
    parser.add_argument("--mode", default="full", choices=["full", "qc", "process", "finalize"],
                        help="Pipeline mode")
    parser.add_argument("--remove-pc", type=int, default=0, help="PC to remove (for finalize mode)")
    parser.add_argument("--resolution", type=float, default=0.8, help="Clustering resolution")
    parser.add_argument("--output-dir", default="", help="Output directory for plots")
    return parser.parse_args()


def load_peaks(peaks_path):
    """Load peaks BED file."""
    return pl.read_csv(peaks_path, separator='\t', has_header=False,
                       new_columns=['chr', 'start', 'end'])


def count_fragments(fragments_path, peaks_df, batch_size=400_000):
    """Count fragments in peaks."""
    # Build peak index
    peak_idx = {}
    for chrom in peaks_df['chr'].unique().to_list():
        sub = peaks_df.filter(pl.col('chr') == chrom).sort('start')
        peak_idx[chrom] = (sub['start'].to_numpy(), sub['end'].to_numpy())
    
    peak_ids = peaks_df['peak_id'].to_list()
    n_peaks = len(peak_ids)
    
    # Read fragments
    import gzip, io
    with gzip.open(fragments_path, 'rb') as gz:
        buf = io.BytesIO(b''.join(l for l in gz if not l.startswith(b'#')))
    
    df_full = pl.read_csv(buf, separator='\t', has_header=False,
                          new_columns=['chr', 'start', 'end', 'bc', 'readSupport'])
    
    bcs = np.sort(df_full['bc'].unique().to_numpy())
    bc_to_row = {b: i for i, b in enumerate(bcs)}
    
    rows_list, cols_list, data_list = [], [], []
    
    for batch_start in range(0, df_full.height, batch_size):
        batch = df_full.slice(batch_start, batch_size)
        bc_arr = batch['bc'].to_numpy()
        start_arr = batch['start'].to_numpy()
        chr_arr = batch['chr'].to_numpy()
        
        for i in range(len(bc_arr)):
            chrom = chr_arr[i]
            if chrom not in peak_idx:
                continue
            starts, ends = peak_idx[chrom]
            pos = start_arr[i]
            cand = np.searchsorted(starts, pos, side='right') - 1
            if cand >= 0 and pos < ends[cand]:
                rows_list.append(bc_to_row[bc_arr[i]])
                cols_list.append(cand)
                data_list.append(1)
    
    # Build sparse matrix
    from scipy.sparse import coo_matrix
    matrix = coo_matrix((data_list, (rows_list, cols_list)), 
                        shape=(len(bcs), n_peaks)).tocsr()
    
    return matrix, bcs, peak_ids


def run_qc(adata, output_dir):
    """Run QC and generate plots."""
    print("Running QC...")
    
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    
    # Save QC plots
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    sc.pl.violin(adata, "n_genes_by_counts", ax=axes[0], show=False)
    sc.pl.violin(adata, "total_counts", ax=axes[1], show=False)
    axes[2].hist(adata.obs["n_genes_by_counts"], bins=50)
    axes[2].set_xlabel("Genes per cell")
    axes[2].set_ylabel("Count")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "qc_metrics.png"), dpi=150)
    plt.close()
    
    print(f"  QC plot saved to: {output_dir}/qc_metrics.png")
    return adata


def run_processing(adata, resolution=0.8):
    """Run main processing: normalize → PCA → neighbors → clusters → UMAP."""
    print("Processing...")
    
    # Normalize
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    
    # Variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    
    # PCA
    sc.tl.pca(adata, n_comps=50)
    
    # Neighbors
    sc.pp.neighbors(adata, n_pcs=50)
    
    # Clustering
    sc.tl.leiden(adata, resolution=resolution)
    n_clusters = len(adata.obs["leiden"].unique())
    print(f"  Found {n_clusters} clusters")
    
    # UMAP
    sc.tl.umap(adata)
    
    return adata


def run_finalize(adata, remove_pc=0, output_dir=""):
    """Finalize: remove depth-correlated PC, re-cluster."""
    print(f"Finalizing (removing PC{remove_pc})...")
    
    # Remove PC
    if remove_pc > 0 and f"X_pca" in adata.obsm:
        pca = adata.obsm["X_pca"].copy()
        pca[:, remove_pc-1] = 0
        adata.obsm["X_pca"] = pca
    
    # Re-cluster
    sc.pp.neighbors(adata, use_rep="X_pca")
    sc.tl.leiden(adata)
    sc.tl.umap(adata)
    
    # Plot
    if output_dir:
        fig, ax = plt.subplots(figsize=(8, 8))
        sc.pl.umap(adata, color="leiden", ax=ax, show=False, 
                   title=f"Final Clusters (PC{remove_pc} removed)")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "umap_final.png"), dpi=150)
        plt.close()
    
    return adata


def main():
    args = parse_args()
    
    output_dir = args.output_dir or os.path.dirname(args.output)
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 60)
    print("nanoCT Processing Pipeline")
    print("=" * 60)
    print(f"Start time: {datetime.now()}")
    print(f"Mode: {args.mode}")
    print(f"Fragments: {args.fragments}")
    print(f"Peaks: {args.peaks}")
    print(f"Output: {args.output}")
    print()
    
    # Load peaks
    peaks_df = load_peaks(args.peaks)
    peaks_df = peaks_df.with_columns(
        (pl.col('chr') + ':' + pl.col('start').cast(pl.Utf8) + '-' + pl.col('end').cast(pl.Utf8)).alias('peak_id')
    )
    print(f"Loaded {len(peaks_df)} peaks")
    
    # Count fragments
    print("Counting fragments in peaks...")
    matrix, barcodes, peak_ids = count_fragments(args.fragments, peaks_df)
    print(f"  Matrix: {matrix.shape[0]} cells x {matrix.shape[1]} peaks")
    
    # Create AnnData
    adata = ad.AnnData(
        X=matrix,
        obs=pl.DataFrame({"barcode": barcodes}).to_pandas().set_index("barcode"),
        var=pl.DataFrame({"peak_id": peak_ids}).to_pandas().set_index("peak_id")
    )
    
    # Run pipeline based on mode
    if args.mode in ["full", "qc"]:
        adata = run_qc(adata, output_dir)
    
    if args.mode in ["full", "process"]:
        adata = run_processing(adata, args.resolution)
    
    if args.mode == "finalize":
        adata = run_finalize(adata, args.remove_pc, output_dir)
    
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
    print("=" * 60)


if __name__ == "__main__":
    main()
