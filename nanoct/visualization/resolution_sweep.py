#!/usr/bin/env python3
"""
Leiden resolution sweep on the chromatin UMAP.
Same UMAP coordinates throughout — only clustering changes.
"""
import numpy as np
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

BASE = '/data/ebaird/scRNAseq/20260522.nanoCT'

print("Loading...")
gl = ad.read_h5ad(f'{BASE}/genelevel_chrom.h5ad')
umap = gl.obsm['X_umap']

# scanpy leiden needs uns['neighbors'] pointing to the connectivity matrix
gl.uns['neighbors'] = {'connectivities_key': 'connectivities',
                       'distances_key': 'distances',
                       'params': {'n_neighbors': 15, 'method': 'umap'}}

RESOLUTIONS = [0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]
NCOLS = 4
NROWS = int(np.ceil(len(RESOLUTIONS) / NCOLS))

fig, axes = plt.subplots(NROWS, NCOLS,
                          figsize=(NCOLS * 4.5, NROWS * 4.2),
                          constrained_layout=True)
axes = axes.flatten()

for ax_i, res in enumerate(RESOLUTIONS):
    sc.tl.leiden(gl, resolution=res, key_added=f'leiden_r{res}', random_state=42)
    labels = gl.obs[f'leiden_r{res}'].values
    n_cl   = len(np.unique(labels))
    label_ints = labels.astype(int)

    cmap = plt.cm.get_cmap('tab20' if n_cl <= 20 else 'hsv', n_cl)
    colors = [cmap(i % n_cl) for i in label_ints]

    ax = axes[ax_i]
    ax.scatter(umap[:, 0], umap[:, 1],
               c=colors, s=1.5, alpha=0.6, rasterized=True, linewidths=0)

    # Label each cluster at its centroid
    for cl in np.unique(label_ints):
        mask = label_ints == cl
        cx, cy = umap[mask, 0].mean(), umap[mask, 1].mean()
        ax.text(cx, cy, str(cl), fontsize=5, ha='center', va='center',
                fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.1', facecolor='white',
                          alpha=0.55, edgecolor='none'))

    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f'resolution = {res}  ({n_cl} clusters)', fontsize=10)

for ax_i in range(len(RESOLUTIONS), len(axes)):
    axes[ax_i].set_visible(False)

fig.suptitle('Leiden resolution sweep — chromatin UMAP', fontsize=14, y=1.01)
out = f'{BASE}/leiden_resolution_sweep.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
print(f"Saved {out}")

# Also save cluster counts
import pandas as pd
counts = pd.DataFrame({
    'resolution': RESOLUTIONS,
    'n_clusters': [len(gl.obs[f'leiden_r{r}'].unique()) for r in RESOLUTIONS]
})
counts.to_csv(f'{BASE}/leiden_resolution_sweep_counts.csv', index=False)
print(counts.to_string(index=False))
