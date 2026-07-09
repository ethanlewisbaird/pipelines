#!/usr/bin/env python3
"""
Generate per-cluster UMAP plots from a completed scGLUE chrom h5ad.
Usage: python3 plot_per_cluster_umap.py <input.h5ad> <output_prefix>
"""
import sys, os, numpy as np, anndata as ad
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

h5ad_path = sys.argv[1]
out_prefix = sys.argv[2]

chrom = ad.read_h5ad(h5ad_path)
umap = chrom.obsm['X_umap']

clusters = sorted(chrom.obs['rna_label_transfer'].unique(), key=lambda x: int(x[1:]))
conf_arr = chrom.obs['transfer_confidence'].values.astype(float)
n = len(clusters)
ncols = min(5, n)
nrows = int(np.ceil(n / ncols))
cmap_label = plt.cm.get_cmap('tab20', n)
label_colours = {cl: cmap_label(i) for i, cl in enumerate(clusters)}

fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 3.8),
                          constrained_layout=True)
fig.suptitle(f'RNA label transfer — per cluster', fontsize=14, y=1.01)

for ax_idx, cl in enumerate(clusters):
    ax = axes[ax_idx // ncols][ax_idx % ncols]
    mask = np.array(chrom.obs['rna_label_transfer'] == cl)
    n_cl = mask.sum()
    ax.scatter(umap[~mask, 0], umap[~mask, 1],
               c='#d4d4d4', s=1, alpha=0.3, rasterized=True, linewidths=0)
    sc = ax.scatter(umap[mask, 0], umap[mask, 1],
                    c=conf_arr[mask], cmap='viridis',
                    vmin=0, vmax=1, s=4, alpha=0.7, rasterized=True, linewidths=0)
    ax.set_title(f'{cl} (n={n_cl:,})', color=label_colours[cl],
                 fontsize=10, fontweight='bold')
    ax.set_xticks([]); ax.set_yticks([])

for ax_idx in range(n, nrows * ncols):
    axes[ax_idx // ncols][ax_idx % ncols].set_visible(False)

cbar_ax = fig.add_axes([1.01, 0.15, 0.015, 0.7])
sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=0, vmax=1))
sm.set_array([])
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label('Transfer confidence', fontsize=11)

fig.savefig(f'{out_prefix}_per_cluster_umap.png', dpi=150, bbox_inches='tight')
print(f"Saved {out_prefix}_per_cluster_umap.png")
