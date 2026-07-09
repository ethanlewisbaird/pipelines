#!/usr/bin/env python3
"""
05a: Post-hoc Harmony alignment on SCGLUE X_glue embeddings.

SCGLUE produces X_glue for both modalities but with a residual modality gap.
Harmony treats 'modality' as a batch variable and iteratively shifts cluster
centroids to remove it — no retraining needed.

Inputs:  scglue_chrom.h5ad, scglue_rna.h5ad
Outputs: harmony_chrom.h5ad (X_harmony + updated label transfer)
         harmony_integration_umap.png
"""

import os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize
from scipy.spatial.distance import cdist
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = '/data/ebaird/scRNAseq/20260522.nanoCT'

print("Loading objects...")
chrom = ad.read_h5ad(f'{BASE}/scglue_chrom.h5ad')
rna   = ad.read_h5ad(f'{BASE}/scglue_rna.h5ad')
print(f"Chromatin: {chrom.n_obs} cells | RNA: {rna.n_obs} cells")

# ---------------------------------------------------------------------------
# 1. Stack X_glue embeddings into a single AnnData for Harmony
# ---------------------------------------------------------------------------
chrom_glue = chrom.obsm['X_glue'].astype(np.float32)
rna_glue   = rna.obsm['X_glue'].astype(np.float32)

combined = ad.AnnData(
    X=np.vstack([chrom_glue, rna_glue]),
    obs=pd.DataFrame({
        'modality': ['chrom'] * chrom.n_obs + ['rna'] * rna.n_obs,
        'source_idx': list(range(chrom.n_obs)) + list(range(rna.n_obs)),
    })
)
combined.obsm['X_pca'] = combined.X.copy()   # Harmony expects X_pca key

print(f"Combined shape: {combined.shape}")

# ---------------------------------------------------------------------------
# 2. Run Harmony — treat modality as the batch variable
# ---------------------------------------------------------------------------
print("Running Harmony...")
sc.external.pp.harmony_integrate(combined, key='modality', basis='X_pca',
                                  adjusted_basis='X_harmony',
                                  max_iter_harmony=30, random_state=42)

chrom_harmony = combined.obsm['X_harmony'][:chrom.n_obs]
rna_harmony   = combined.obsm['X_harmony'][chrom.n_obs:]
print(f"X_harmony shape: {chrom_harmony.shape}")

# ---------------------------------------------------------------------------
# 3. Alignment sanity check
# ---------------------------------------------------------------------------
rna_n   = normalize(rna_harmony[:500].astype(np.float32),   norm='l2')
chrom_n = normalize(chrom_harmony[:500].astype(np.float32), norm='l2')
within_rna   = cdist(rna_n,   rna_n,   metric='cosine').mean()
within_chrom = cdist(chrom_n, chrom_n, metric='cosine').mean()
cross_modal  = cdist(chrom_n, rna_n,   metric='cosine').mean()
print(f"\nAlignment after Harmony:")
print(f"  Within-RNA:    {within_rna:.3f}")
print(f"  Within-chrom:  {within_chrom:.3f}")
print(f"  Cross-modal:   {cross_modal:.3f}")
if cross_modal < (within_rna + within_chrom) / 2:
    print("  -> Good alignment")
else:
    print("  -> WARNING: still poor alignment")

# ---------------------------------------------------------------------------
# 4. Re-run k-NN label transfer in Harmony space
# ---------------------------------------------------------------------------
print("\nRe-running k-NN label transfer in Harmony space...")
cluster_col = 'cluster'
K = 15

rna_emb   = normalize(rna_harmony.astype(np.float32),   norm='l2')
chrom_emb = normalize(chrom_harmony.astype(np.float32), norm='l2')

np.random.seed(42)
rna_labels = rna.obs[cluster_col].values
sample_idx = []
for cl in np.unique(rna_labels):
    idx = np.where(rna_labels == cl)[0]
    sample_idx.extend(np.random.choice(idx, min(500, len(idx)), replace=False).tolist())
sample_idx = np.array(sample_idx)

nn = NearestNeighbors(n_neighbors=K, metric='cosine', algorithm='brute', n_jobs=-1)
nn.fit(rna_emb[sample_idx])
distances, indices = nn.kneighbors(chrom_emb)

rna_labels_sub = rna_labels[sample_idx]
transferred, confidence = [], []
for idx_row in indices:
    votes = rna_labels_sub[idx_row]
    unique, counts = np.unique(votes, return_counts=True)
    winner = unique[np.argmax(counts)]
    transferred.append(winner)
    confidence.append(counts.max() / K)

chrom.obsm['X_harmony'] = chrom_harmony
chrom.obs['harmony_label_transfer']  = transferred
chrom.obs['harmony_confidence']      = confidence
chrom.obs['harmony_knn_distance']    = distances.mean(axis=1)

print("Label distribution (Harmony):")
print(chrom.obs['harmony_label_transfer'].value_counts().head(10))
print(f"\nMean confidence: {np.mean(confidence):.3f}  (vs SCGLUE: {chrom.obs['transfer_confidence'].mean():.3f})")
print(f"Mean kNN dist:   {distances.mean():.3f}  (vs SCGLUE: {chrom.obs['mean_knn_distance'].mean():.3f})")

# Crosstab vs original SCGLUE labels
xtab = pd.crosstab(
    chrom.obs['rna_label_transfer'],
    chrom.obs['harmony_label_transfer'],
    normalize='index'
).round(2)
agreement = np.diag(xtab.values).mean() if xtab.shape[0] == xtab.shape[1] else np.nan
print(f"\nLabel agreement (SCGLUE vs Harmony): {agreement:.3f} where 1=identical")

# ---------------------------------------------------------------------------
# 5. UMAP in Harmony space
# ---------------------------------------------------------------------------
print("\nComputing UMAP on Harmony embedding...")
combined.obsm['X_harmony'] = np.vstack([chrom_harmony, rna_harmony])
sc.pp.neighbors(combined, use_rep='X_harmony', n_neighbors=20, metric='cosine')
sc.tl.umap(combined, min_dist=0.3, random_state=42)

chrom.obsm['X_umap_harmony'] = combined.obsm['X_umap'][:chrom.n_obs]
rna_umap_harmony = combined.obsm['X_umap'][chrom.n_obs:]

# ---------------------------------------------------------------------------
# 6. Plot
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 3, figsize=(18, 11), constrained_layout=True)
umap = chrom.obsm['X_umap_harmony']

# Row 1: chromatin cells coloured by leiden / harmony label / confidence
for cl in sorted(chrom.obs['leiden'].astype(str).unique()):
    m = chrom.obs['leiden'].astype(str) == cl
    axes[0, 0].scatter(umap[m, 0], umap[m, 1], s=1, alpha=0.4, label=cl)
axes[0, 0].set_title('Leiden clusters (chromatin)')
axes[0, 0].legend(markerscale=5, fontsize=6)

for cl in sorted(chrom.obs['harmony_label_transfer'].unique()):
    m = chrom.obs['harmony_label_transfer'] == cl
    axes[0, 1].scatter(umap[m, 0], umap[m, 1], s=1, alpha=0.4, label=cl)
axes[0, 1].set_title('Harmony label transfer')
axes[0, 1].legend(markerscale=5, fontsize=6)

sc2 = axes[0, 2].scatter(umap[:, 0], umap[:, 1], s=1,
                          c=chrom.obs['harmony_confidence'], cmap='viridis', alpha=0.5)
fig.colorbar(sc2, ax=axes[0, 2])
axes[0, 2].set_title('Transfer confidence (Harmony)')

# Row 2: joint UMAP (chrom + RNA)
all_umap = combined.obsm['X_umap']
n_ch = chrom.n_obs
axes[1, 0].scatter(all_umap[n_ch:, 0], all_umap[n_ch:, 1], s=0.5, alpha=0.2, c='lightgrey', label='RNA')
axes[1, 0].scatter(all_umap[:n_ch, 0], all_umap[:n_ch, 1], s=0.5, alpha=0.5, c='steelblue', label='Chromatin')
axes[1, 0].set_title('Joint UMAP (modality)')
axes[1, 0].legend(markerscale=8, fontsize=8)

# RNA cells coloured by cluster on joint UMAP
rna_labs = rna.obs[cluster_col].values
for cl in sorted(np.unique(rna_labs)):
    m = rna_labs == cl
    axes[1, 1].scatter(rna_umap_harmony[m, 0], rna_umap_harmony[m, 1], s=0.5, alpha=0.3, label=cl)
axes[1, 1].set_title('RNA cell types (joint UMAP)')
axes[1, 1].legend(markerscale=5, fontsize=6)

# Confidence comparison
axes[1, 2].scatter(chrom.obs['transfer_confidence'], chrom.obs['harmony_confidence'], s=1, alpha=0.2)
axes[1, 2].set_xlabel('SCGLUE confidence')
axes[1, 2].set_ylabel('Harmony confidence')
axes[1, 2].plot([0, 1], [0, 1], 'r--', lw=1)
axes[1, 2].set_title('Confidence comparison')

plt.savefig(f'{BASE}/harmony_integration_umap.png', dpi=150)
print(f"\nSaved {BASE}/harmony_integration_umap.png")

# ---------------------------------------------------------------------------
# 7. Save
# ---------------------------------------------------------------------------
chrom.write_h5ad(f'{BASE}/harmony_chrom.h5ad')
print(f"Saved harmony_chrom.h5ad")
print("\nDone.")
