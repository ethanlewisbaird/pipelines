#!/usr/bin/env python3
"""
05b: Gene-level re-integration of nanoCT H3K27ac with scRNA-seq.

Instead of integrating at bin level (sparse, 24k features), use the
gene-level H3K27ac activity scores from infer_layer (imputed_acet in obsm).
This collapses feature space to ~8-10k genes directly comparable to RNA.

Strategy:
  1. Build a joint AnnData: chromatin gene-level acet (log-norm) + RNA (log-norm)
  2. Select HVGs shared between both
  3. Joint PCA → Harmony (batch = modality)
  4. k-NN label transfer in Harmony space
  5. Compare to SCGLUE results

Inputs:  inferred_chrom.h5ad, scglue_rna.h5ad
Outputs: genelevel_chrom.h5ad, genelevel_integration_umap.png
"""

import os
import numpy as np
import pandas as pd
import scipy.sparse as sps
import anndata as ad
import scanpy as sc
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize
from scipy.spatial.distance import cdist
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = '/data/ebaird/scRNAseq/20260522.nanoCT'

# ---------------------------------------------------------------------------
# 1. Load
# ---------------------------------------------------------------------------
print("Loading inferred chromatin object...")
chrom = ad.read_h5ad(f'{BASE}/inferred_chrom.h5ad')
print(chrom)
print("obsm keys:", list(chrom.obsm.keys()))

print("\nLoading RNA object...")
rna = ad.read_h5ad(f'{BASE}/scglue_rna.h5ad')
print(rna)

# ---------------------------------------------------------------------------
# 2. Extract gene-level H3K27ac (imputed) from chromatin
# ---------------------------------------------------------------------------
# infer_layer stores results in obsm with the gene_names used during infer_layer
# We need to recover those gene names — they were passed as `names` = rna.var_names
# after filtering to genes connected in the guidance graph.
# The shape of obsm['acet_imp'] is (n_chrom_cells, n_genes).
# Gene order matches the `names` array from the infer_layer run.
# Since we used rna.var_names (after guidance-graph filtering) as gene_names,
# the simplest way to recover them is via the RNA object's var_names.
# The guidance graph filtered RNA to connected genes; after SCGLUE rna has those genes.

acet_imp = chrom.obsm['acet_imp']
if sps.issparse(acet_imp):
    acet_imp = acet_imp.toarray()
acet_imp = acet_imp.astype(np.float32)

n_genes = acet_imp.shape[1]
rna_var_names = rna.var_names   # these were the gene_names used in infer_layer
print(f"\nChromatin gene-level acet_imp: {acet_imp.shape}")
print(f"RNA var_names: {len(rna_var_names)}")

if n_genes != len(rna_var_names):
    print(f"WARNING: shape mismatch ({n_genes} vs {len(rna_var_names)}). "
          "Truncating to min.")
    n_use = min(n_genes, len(rna_var_names))
    acet_imp = acet_imp[:, :n_use]
    rna_var_names = rna_var_names[:n_use]

gene_names = rna_var_names

# ---------------------------------------------------------------------------
# 3. Log-normalise H3K27ac gene scores
# ---------------------------------------------------------------------------
# Library-size normalise each cell then log1p (mirrors RNA preprocessing)
totals = acet_imp.sum(axis=1, keepdims=True)
totals[totals == 0] = 1
acet_lognorm = np.log1p(acet_imp / totals * 1e4)
print(f"Acet log-norm: min={acet_lognorm.min():.3f}, max={acet_lognorm.max():.3f}, "
      f"mean={acet_lognorm.mean():.4f}")
print(f"Fraction non-zero: {(acet_lognorm > 0).mean():.3f}")

# ---------------------------------------------------------------------------
# 4. Build AnnData objects in common gene space
# ---------------------------------------------------------------------------
chrom_ad = ad.AnnData(
    X=acet_lognorm,
    obs=chrom.obs.copy(),
    var=pd.DataFrame(index=gene_names),
)
chrom_ad.obsm['X_umap']  = chrom.obsm['X_umap']
chrom_ad.obsm['X_multi_spectral'] = chrom.obsm['X_multi_spectral']

# RNA: use log-normalised counts
rna_X = rna.layers['counts']
if sps.issparse(rna_X):
    rna_X = rna_X.toarray().astype(np.float32)
else:
    rna_X = rna_X.astype(np.float32)
rna_totals = rna_X.sum(axis=1, keepdims=True)
rna_totals[rna_totals == 0] = 1
rna_lognorm = np.log1p(rna_X / rna_totals * 1e4)

rna_ad = ad.AnnData(
    X=rna_lognorm,
    obs=rna.obs.copy(),
    var=pd.DataFrame(index=rna.var_names),
)
rna_ad.obsm['X_umap'] = rna.obsm.get('X_umap', np.zeros((rna.n_obs, 2)))

# Intersect genes
common = chrom_ad.var_names.intersection(rna_ad.var_names)
print(f"\nCommon genes: {len(common)}")
chrom_ad = chrom_ad[:, common].copy()
rna_ad   = rna_ad[:, common].copy()

# ---------------------------------------------------------------------------
# 5. HVG selection on RNA, apply to both
# ---------------------------------------------------------------------------
print("Selecting HVGs on RNA...")
sc.pp.highly_variable_genes(rna_ad, n_top_genes=min(3000, len(common)),
                             flavor='cell_ranger')
hvg_mask = rna_ad.var['highly_variable'].values
print(f"HVGs: {hvg_mask.sum()}")

chrom_hvg = chrom_ad[:, hvg_mask].copy()
rna_hvg   = rna_ad[:, hvg_mask].copy()

# Also check variance of chromatin — skip genes with zero variance
chrom_var = np.var(chrom_hvg.X, axis=0)
rna_var_v = np.var(rna_hvg.X, axis=0)
both_var  = (chrom_var > 0) & (rna_var_v > 0)
print(f"HVGs with variance in both modalities: {both_var.sum()}")
chrom_hvg = chrom_hvg[:, both_var].copy()
rna_hvg   = rna_hvg[:, both_var].copy()

# ---------------------------------------------------------------------------
# 6. Joint PCA
# ---------------------------------------------------------------------------
print("Running joint PCA (n_comps=50)...")
combined = ad.concat([chrom_hvg, rna_hvg], label='modality',
                      keys=['chrom', 'rna'], merge='same')
sc.tl.pca(combined, n_comps=50, svd_solver='auto')
print(f"PCA done: {combined.obsm['X_pca'].shape}")

# ---------------------------------------------------------------------------
# 7. Harmony batch correction (modality = batch)
# ---------------------------------------------------------------------------
print("Running Harmony...")
sc.external.pp.harmony_integrate(combined, key='modality', basis='X_pca',
                                  adjusted_basis='X_harmony',
                                  max_iter_harmony=30, random_state=42)

chrom_harmony = combined.obsm['X_harmony'][:chrom.n_obs]
rna_harmony   = combined.obsm['X_harmony'][chrom.n_obs:]

# ---------------------------------------------------------------------------
# 8. Alignment sanity check
# ---------------------------------------------------------------------------
rna_n   = normalize(rna_harmony[:500].astype(np.float32),   norm='l2')
chrom_n = normalize(chrom_harmony[:500].astype(np.float32), norm='l2')
within_rna   = cdist(rna_n,   rna_n,   metric='cosine').mean()
within_chrom = cdist(chrom_n, chrom_n, metric='cosine').mean()
cross_modal  = cdist(chrom_n, rna_n,   metric='cosine').mean()
print(f"\nAlignment sanity check (gene-level Harmony):")
print(f"  Within-RNA:    {within_rna:.3f}")
print(f"  Within-chrom:  {within_chrom:.3f}")
print(f"  Cross-modal:   {cross_modal:.3f}")
if cross_modal < (within_rna + within_chrom) / 2:
    print("  -> Good alignment")
else:
    print("  -> WARNING: poor alignment")

# ---------------------------------------------------------------------------
# 9. k-NN label transfer in Harmony space
# ---------------------------------------------------------------------------
print("\nRunning k-NN label transfer (K=15, cosine)...")
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

chrom.obsm['X_genelevel_harmony'] = chrom_harmony
chrom.obs['genelevel_label_transfer'] = transferred
chrom.obs['genelevel_confidence']     = confidence
chrom.obs['genelevel_knn_distance']   = distances.mean(axis=1)

print("Label distribution (gene-level Harmony):")
print(chrom.obs['genelevel_label_transfer'].value_counts().head(10))
print(f"\nMean confidence: {np.mean(confidence):.3f}  (vs SCGLUE: {chrom.obs['transfer_confidence'].mean():.3f})")
print(f"Mean kNN dist:   {distances.mean():.3f}  (vs SCGLUE: {chrom.obs['mean_knn_distance'].mean():.3f})")

print("\nLeiden -> gene-level label fractions:")
xtab = pd.crosstab(
    chrom.obs['leiden'].astype(str),
    chrom.obs['genelevel_label_transfer'],
    normalize='index'
).round(2)
print(xtab)

# ---------------------------------------------------------------------------
# 10. UMAP on joint Harmony embedding
# ---------------------------------------------------------------------------
print("\nComputing joint UMAP...")
sc.pp.neighbors(combined, use_rep='X_harmony', n_neighbors=20, metric='cosine')
sc.tl.umap(combined, min_dist=0.3, random_state=42)

chrom.obsm['X_umap_genelevel'] = combined.obsm['X_umap'][:chrom.n_obs]
rna_umap = combined.obsm['X_umap'][chrom.n_obs:]

# ---------------------------------------------------------------------------
# 11. Plot
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 3, figsize=(18, 11), constrained_layout=True)
umap_ch = chrom.obsm['X_umap_genelevel']
all_umap = combined.obsm['X_umap']
n_ch = chrom.n_obs

# Joint: modality
axes[0, 0].scatter(all_umap[n_ch:, 0], all_umap[n_ch:, 1], s=0.5, alpha=0.15, c='lightgrey', label='RNA')
axes[0, 0].scatter(all_umap[:n_ch, 0], all_umap[:n_ch, 1], s=0.8, alpha=0.6, c='steelblue', label='Chromatin')
axes[0, 0].set_title('Joint UMAP — modality')
axes[0, 0].legend(markerscale=8, fontsize=8)

# RNA cells coloured by cluster
for cl in sorted(np.unique(rna_labels)):
    m = rna_labels == cl
    axes[0, 1].scatter(rna_umap[m, 0], rna_umap[m, 1], s=0.5, alpha=0.3, label=cl)
axes[0, 1].set_title('RNA cell types (joint UMAP)')
axes[0, 1].legend(markerscale=5, fontsize=6)

# Chromatin: gene-level label
for cl in sorted(chrom.obs['genelevel_label_transfer'].unique()):
    m = chrom.obs['genelevel_label_transfer'] == cl
    axes[0, 2].scatter(umap_ch[m, 0], umap_ch[m, 1], s=1, alpha=0.5, label=cl)
axes[0, 2].set_title('Gene-level label transfer')
axes[0, 2].legend(markerscale=5, fontsize=6)

# Chromatin: leiden
for cl in sorted(chrom.obs['leiden'].astype(str).unique()):
    m = chrom.obs['leiden'].astype(str) == cl
    axes[1, 0].scatter(umap_ch[m, 0], umap_ch[m, 1], s=1, alpha=0.5, label=cl)
axes[1, 0].set_title('Leiden clusters (chromatin)')
axes[1, 0].legend(markerscale=5, fontsize=6)

# Confidence
sc2 = axes[1, 1].scatter(umap_ch[:, 0], umap_ch[:, 1], s=1,
                          c=chrom.obs['genelevel_confidence'], cmap='viridis', alpha=0.5)
fig.colorbar(sc2, ax=axes[1, 1])
axes[1, 1].set_title('Transfer confidence (gene-level)')

# Confidence comparison vs SCGLUE
axes[1, 2].scatter(chrom.obs['transfer_confidence'],
                   chrom.obs['genelevel_confidence'], s=1, alpha=0.2)
axes[1, 2].set_xlabel('SCGLUE confidence')
axes[1, 2].set_ylabel('Gene-level Harmony confidence')
axes[1, 2].plot([0, 1], [0, 1], 'r--', lw=1)
axes[1, 2].set_title('Confidence: SCGLUE vs gene-level')

plt.savefig(f'{BASE}/genelevel_integration_umap.png', dpi=150)
print(f"\nSaved {BASE}/genelevel_integration_umap.png")

# ---------------------------------------------------------------------------
# 12. Save
# ---------------------------------------------------------------------------
chrom.write_h5ad(f'{BASE}/genelevel_chrom.h5ad')
print("Saved genelevel_chrom.h5ad")
print("\nDone.")
