#!/usr/bin/env python3
"""
Common QC utilities for BAIRD pipelines.

Usage:
    from common.qc import run_qc, plot_qc_metrics
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt


def run_qc(adata, 
           min_genes=200, 
           min_cells=3,
           max_pct_mito=20,
           max_genes=5000):
    """Run standard QC filtering."""
    
    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    n_before = adata.n_obs
    
    # Filter cells
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    # Filter by mito percentage
    adata = adata[adata.obs.pct_counts_mt < max_pct_mito, :]
    
    # Filter by max genes
    adata = adata[adata.obs.n_genes_by_counts < max_genes, :]
    
    n_after = adata.n_obs
    print(f"QC: {n_before} → {n_after} cells ({n_before - n_after} removed)")
    
    return adata


def plot_qc_metrics(adata, output_dir=".", prefix="qc"):
    """Generate QC plots."""
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    sc.pl.violin(adata, 'n_genes_by_counts', ax=axes[0], show=False)
    axes[0].set_title('Genes per cell')
    
    sc.pl.violin(adata, 'total_counts', ax=axes[1], show=False)
    axes[1].set_title('Counts per cell')
    
    sc.pl.violin(adata, 'pct_counts_mt', ax=axes[2], show=False)
    axes[2].set_title('Mito %')
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{prefix}_violin.png", dpi=150, bbox_inches='tight')
    plt.close()
    
    # Scatter plot
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax=ax, show=False)
    plt.savefig(f"{output_dir}/{prefix}_scatter.png", dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"QC plots saved to {output_dir}/")


def normalize_and_log(adata, target_sum=1e4):
    """Normalize and log transform."""
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    return adata
