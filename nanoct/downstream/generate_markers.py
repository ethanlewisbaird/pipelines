#!/usr/bin/env python3
"""
Generate top 50 markers per cluster for H3K27ac and H3K27me3.
Uses the reclustered chromatin object.
"""
import os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

BASE = '/data/ebaird/scentinel/nanoCT/20260522.nanoCT'
OUT_DIR = f'{BASE}/cluster_markers/reclustered'
os.makedirs(OUT_DIR, exist_ok=True)

print("Loading chromatin object...")
chrom = ad.read_h5ad(f'{BASE}/analysis_05.26/combined_dim_reduced.h5ad')
print(f"Shape: {chrom.shape}")
print(f"Clusters: {chrom.obs['leiden'].value_counts().sort_index().to_dict()}")

def get_markers(adata, layer, group_col, n_top=50):
    """Get top markers for each cluster using Mann-Whitney U test."""
    results = []
    clusters = sorted(adata.obs[group_col].unique(), key=lambda x: int(x))
    
    for cluster in clusters:
        print(f"  Processing cluster {cluster}...")
        mask = adata.obs[group_col] == cluster
        cluster_data = adata[mask].layers[layer] if layer in adata.layers else adata[mask].X
        other_data = adata[~mask].layers[layer] if layer in adata.layers else adata[~mask].X
        
        # Get gene names
        if layer in adata.layers:
            gene_names = adata.var_names
        else:
            gene_names = adata.var_names
        
        # Test each gene
        pvals = []
        log2fc = []
        scores = []
        
        for i, gene in enumerate(gene_names):
            cluster_expr = cluster_data[:, i].toarray().flatten() if hasattr(cluster_data[:, i], 'toarray') else cluster_data[:, i].flatten()
            other_expr = other_data[:, i].toarray().flatten() if hasattr(other_data[:, i], 'toarray') else other_data[:, i].flatten()
            
            # Skip if no expression
            if np.sum(cluster_expr) == 0 and np.sum(other_expr) == 0:
                pvals.append(1.0)
                log2fc.append(0.0)
                scores.append(0.0)
                continue
            
            # Mann-Whitney U test
            try:
                stat, pval = mannwhitneyu(cluster_expr, other_expr, alternative='greater')
                pvals.append(pval)
                
                # Log2FC
                mean_cluster = np.mean(cluster_expr)
                mean_other = np.mean(other_expr)
                if mean_other > 0:
                    l2fc = np.log2(mean_cluster / mean_other)
                else:
                    l2fc = 0.0
                log2fc.append(l2fc)
                scores.append(stat)
            except:
                pvals.append(1.0)
                log2fc.append(0.0)
                scores.append(0.0)
        
        # Multiple testing correction
        _, padj, _, _ = multipletests(pvals, method='fdr_bh')
        
        # Create DataFrame
        df = pd.DataFrame({
            'cluster': cluster,
            'rank': range(1, len(gene_names) + 1),
            'bin': gene_names,
            'nearest_gene': gene_names,
            'score': scores,
            'padj': padj,
            'logFC': log2fc,
            'mark': layer.replace('_binary', '').replace('acet', 'H3K27ac').replace('meth', 'H3K27me3')
        })
        
        # Sort by score and take top n
        df = df.sort_values('score', ascending=False).head(n_top)
        df['rank'] = range(1, len(df) + 1)
        results.append(df)
    
    return pd.concat(results, ignore_index=True)

# Process H3K27ac
print("\nProcessing H3K27ac markers...")
if 'acet' in chrom.layers:
    ac_markers = get_markers(chrom, 'acet', 'leiden', n_top=50)
else:
    print("  No 'acet' layer found, using X...")
    ac_markers = get_markers(chrom, 'X', 'leiden', n_top=50)

# Process H3K27me3
print("\nProcessing H3K27me3 markers...")
if 'meth' in chrom.layers:
    me_markers = get_markers(chrom, 'meth', 'leiden', n_top=50)
else:
    print("  No 'meth' layer found, using X...")
    me_markers = get_markers(chrom, 'X', 'leiden', n_top=50)

# Save
ac_markers.to_csv(f'{OUT_DIR}/top50_markers_H3K27ac.csv', index=False)
me_markers.to_csv(f'{OUT_DIR}/top50_markers_H3K27me3.csv', index=False)

print(f"\nSaved {len(ac_markers)} H3K27ac markers")
print(f"Saved {len(me_markers)} H3K27me3 markers")
print("Done!")
