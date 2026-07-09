#!/usr/bin/env python3
"""
nanoCT-scRNAseq Integration with scGLUE

Integrates nanoCT chromatin data with scRNA-seq for label transfer.

Usage:
    python scglue_integration.py --chrom chromatin.h5ad --rna rna.h5ad --output output_dir

Environment variables:
    NANOCT_DATA_DIR - base data directory
"""

import argparse
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sps
import scglue
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import torch


def parse_args():
    parser = argparse.ArgumentParser(description="scGLUE integration")
    parser.add_argument("--chrom", required=True, help="Chromatin h5ad file")
    parser.add_argument("--rna", required=True, help="RNA h5ad file")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--gtf", default="", help="GTF file for gene annotation")
    return parser.parse_args()


def prepare_rna(rna):
    """Prepare RNA object for scGLUE."""
    print("Preparing RNA...")
    
    if 'counts' in rna.layers:
        rna.X = rna.layers['counts'].copy()
        rna_pp = rna.copy()
        rna_pp.X = rna.layers['data'].copy() if 'data' in rna.layers else rna.X.copy()
        sc.pp.highly_variable_genes(rna_pp, n_top_genes=4000, flavor='cell_ranger')
        sc.tl.pca(rna_pp, n_comps=50, use_highly_variable=True)
        rna.obsm['X_pca'] = rna_pp.obsm['X_pca']
        rna.var['highly_variable'] = rna_pp.var['highly_variable']
        del rna_pp
        RNA_DIST = 'NB'
    else:
        sc.pp.highly_variable_genes(rna, n_top_genes=4000, flavor='cell_ranger')
        sc.tl.pca(rna, n_comps=50, use_highly_variable=True)
        RNA_DIST = 'Normal'
    
    print(f"  RNA distribution: {RNA_DIST}")
    print(f"  RNA HVGs: {rna.var['highly_variable'].sum()}")
    
    return rna, RNA_DIST


def prepare_chrom(chrom):
    """Prepare chromatin object for scGLUE."""
    print("Preparing chromatin...")
    
    sc.pp.highly_variable_genes(chrom, n_top_genes=2000)
    sc.tl.pca(chrom, n_comps=50)
    
    return chrom


def run_scglue(chrom, rna, output_dir, gtf_path=""):
    """Run scGLUE integration."""
    print("Running scGLUE...")
    
    # Configure scGLUE
    scglue.config.setup_dataset(chrom, rna)
    
    # Build graph
    print("  Building guidance graph...")
    guidance = scglue.genomics.build_guidance_graph(
        rna.var_names, chrom.var_names,
        gene_bed=gtf_path if gtf_path else None
    )
    
    # Configure datasets
    chrom.var["id"] = chrom.var_names
    rna.var["id"] = rna.var_names
    
    scglue.models.configure_dataset(
        chrom, "NB", use_highly_variable=True,
        use_layer="counts", use_rep="X_pca"
    )
    scglue.models.configure_dataset(
        rna, "NB", use_highly_variable=True,
        use_layer="counts", use_rep="X_pca"
    )
    
    # Train model
    print("  Training scGLUE model...")
    model = scglue.models.fit_SCGLUE(
        {"chrom": chrom, "rna": rna},
        guidance,
        fit_kwargs={"directory": os.path.join(output_dir, "scglue_model")}
    )
    
    # Extract embeddings
    print("  Extracting embeddings...")
    chrom.obsm["X_glue"] = model.encode_data("chrom", chrom)
    rna.obsm["X_glue"] = model.encode_data("rna", rna)
    
    # Label transfer
    print("  Performing label transfer...")
    from sklearn.neighbors import KNeighborsClassifier
    
    # Get shared cells or use embedding similarity
    if 'leiden' in rna.obs.columns:
        knn = KNeighborsClassifier(n_neighbors=10)
        knn.fit(rna.obsm["X_glue"], rna.obs['leiden'])
        chrom.obs['transferred_labels'] = knn.predict(chrom.obsm["X_glue"])
    
    return chrom, rna


def main():
    args = parse_args()
    
    os.makedirs(args.output, exist_ok=True)
    
    print("=" * 60)
    print("nanoCT-scRNAseq Integration")
    print("=" * 60)
    
    # Load data
    print("Loading data...")
    chrom = ad.read_h5ad(args.chrom)
    rna = ad.read_h5ad(args.rna)
    
    print(f"  Chromatin: {chrom.shape[0]} cells x {chrom.shape[1]} features")
    print(f"  RNA: {rna.shape[0]} cells x {rna.shape[1]} genes")
    
    # Prepare
    rna, rna_dist = prepare_rna(rna)
    chrom = prepare_chrom(chrom)
    
    # Run scGLUE
    chrom, rna = run_scglue(chrom, rna, args.output, args.gtf)
    
    # Save
    print("Saving results...")
    chrom.write(os.path.join(args.output, "chrom_glue.h5ad"))
    rna.write(os.path.join(args.output, "rna_glue.h5ad"))
    
    # UMAP
    print("Generating UMAP...")
    sc.pp.neighbors(chrom, use_rep="X_glue")
    sc.tl.umap(chrom)
    
    fig, ax = plt.subplots(figsize=(8, 8))
    sc.pl.umap(chrom, color="transferred_labels", ax=ax, show=False,
               title="scGLUE Label Transfer")
    plt.tight_layout()
    plt.savefig(os.path.join(args.output, "umap_transferred.png"), dpi=150)
    plt.close()
    
    print()
    print("=" * 60)
    print("Integration complete!")
    print(f"Output: {args.output}")
    print("=" * 60)


if __name__ == "__main__":
    main()
