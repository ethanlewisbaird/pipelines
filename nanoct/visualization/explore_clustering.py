"""Explore leiden clustering at multiple resolutions, save summary."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '/data/ebaird/scRNAseq/20260522.nanoCT/scit_src')

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import src as scit

OUT = os.path.dirname(os.path.abspath(__file__))

# Load
adata = ad.read_h5ad(f'{OUT}/combined_dim_reduced.h5ad')
print(f"Loaded adata: {adata.shape[0]} cells x {adata.shape[1]} features")
print(f"Obs columns: {list(adata.obs.columns)}")
print(f"Obsm keys: {list(adata.obsm.keys())}")

# Rebuild KNN graph on multi_spectral embedding
scit.gr.knn(adata, 'X_multi_spectral')
g = scit.gr.neighbor_graph(adata)

# Test multiple resolutions
resolutions = [0.3, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0]

results = {}
for res in resolutions:
    scit.gr.leiden(adata, g, resolution=res, random_seed=42)
    clusters = adata.obs['leiden']
    n_clusters = clusters.nunique()
    counts = clusters.value_counts().sort_index()
    results[res] = {
        'n_clusters': n_clusters,
        'counts': counts,
        'min_size': counts.min(),
        'max_size': counts.max()
    }
    print(f"  res={res:.1f}: {n_clusters} clusters, sizes {counts.min()}-{counts.max()}")

# Summary table
print("\n=== Clustering Summary ===")
print(f"{'Res':>5} | {'Clusters':>8} | {'Min':>5} | {'Max':>5} | {'Cluster sizes':>80}")
print("-" * 100)
for res in resolutions:
    r = results[res]
    sizes = {str(k): int(v) for k, v in sorted(r['counts'].items())}
    print(f"{res:5.1f} | {r['n_clusters']:>8} | {r['min_size']:>5} | {r['max_size']:>5} | {sizes}")

# Save UMAP plots for each resolution
for res in resolutions:
    scit.gr.leiden(adata, g, resolution=res, random_seed=42)
    fig = scit.pl.embedding2d(adata, 'X_umap', 'leiden', 
                              title=f'Leiden res={res}', show=False)
    fig.savefig(f'{OUT}/cluster_res{res:.1f}.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved cluster_res{res:.1f}.png")

# Restore a reasonable default (res=1.0)
scit.gr.leiden(adata, g, resolution=1.0, random_seed=42)
print("\nRestored resolution=1.0 as default in adata")
