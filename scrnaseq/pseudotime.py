#!/usr/bin/env python3
"""
Pseudotime analysis pipeline for scRNA-seq data.

Usage:
    python pseudotime.py --input <h5ad_file> [--output <output_dir>] [--root_key <obs_column>]

This is a base pipeline. For specific analyses, copy to your analysis
directory and modify as needed.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(description="Pseudotime analysis")
    parser.add_argument("--input", required=True, help="Input h5ad file")
    parser.add_argument("--output", default=".", help="Output directory")
    parser.add_argument("--root_key", default="leiden", help="Column for root cell selection")
    parser.add_argument("--root_value", default=None, help="Specific value to use as root")
    parser.add_argument("--n_neighbors", type=int, default=15, help="Number of neighbors for diffusion")
    parser.add_argument("--n_comps", type=int, default=50, help="Number of diffusion components")
    parser.add_argument("--use_rep", default="X_pca", help="Representation to use")
    return parser.parse_args()


def run_pseudotime(adata, args):
    """Run pseudotime analysis."""
    print(f"Input: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # Compute diffusion map
    print("Computing diffusion map...")
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, use_rep=args.use_rep)
    sc.tl.diffmap(adata, n_comps=args.n_comps)
    
    # Select root cell
    if args.root_value:
        root_mask = adata.obs[args.root_key] == args.root_value
    else:
        # Use the cluster with highest mean diffusion component 1
        cluster_means = adata.obs.groupby(args.root_key)['diffmap_dc0'].mean()
        root_cluster = cluster_means.idxmax()
        root_mask = adata.obs[args.root_key] == root_cluster
        print(f"Auto-selected root cluster: {root_cluster}")
    
    # Set root cell (cell closest to centroid of root cluster)
    root_idx = np.where(root_mask)[0]
    if len(root_idx) == 0:
        print(f"Warning: No cells found for root selection")
        root_idx = [0]
    
    # Use diffusion pseudotime
    adata.uns['iroot'] = root_idx[0]
    sc.tl.dpt(adata, n_branchings=0, n_dcs=10)
    
    print(f"Pseudotime range: {adata.obs['dpt_pseudotime'].min():.3f} - {adata.obs['dpt_pseudotime'].max():.3f}")
    
    return adata


def plot_results(adata, output_dir):
    """Generate plots."""
    os.makedirs(output_dir, exist_ok=True)
    
    # UMAP with pseudotime
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    sc.pl.umap(adata, color='dpt_pseudotime', ax=axes[0], show=False, 
               title='Pseudotime', colorbar_label='Pseudotime')
    sc.pl.umap(adata, color='leiden', ax=axes[1], show=False, 
               title='Clusters')
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/pseudotime_umap.png", dpi=150, bbox_inches='tight')
    plt.close()
    
    # Diffusion map
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.diffmap(adata, color='dpt_pseudotime', ax=ax, show=False,
                  title='Diffusion Map')
    plt.savefig(f"{output_dir}/diffusion_map.png", dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"Plots saved to {output_dir}/")


def main():
    args = parse_args()
    
    # Load data
    print(f"Loading {args.input}...")
    adata = sc.read_h5ad(args.input)
    
    # Run pseudotime
    adata = run_pseudotime(adata, args)
    
    # Save results
    output_h5ad = f"{args.output}/pseudotime.h5ad"
    adata.write_h5ad(output_h5ad)
    print(f"Saved {output_h5ad}")
    
    # Generate plots
    plot_results(adata, args.output)
    
    print("Done!")


if __name__ == "__main__":
    main()
