#!/usr/bin/env python3
"""
SCGLUE integration for ATAC universe peaks ±3kb chromatin object.

Adapts scglue_integration_v3.py to use the peaks-based chromatin data
(combined_ATACuniverse_binary.h5ad) instead of the 5kb bins.

Key changes:
  - Loads combined_ATACuniverse_binary.h5ad instead of combined_dim_reduced.h5ad
  - Parses chr/start/end from peak_ids (format: chr2L:12717-12917)
  - Uses acet layer for TF-IDF LSI embedding
"""

import os
import re
import numpy as np
import pandas as pd
import scipy.sparse as sps
import anndata as ad
import scanpy as sc
import scglue
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import torch

n_cpus = int(os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count()))
torch.set_num_threads(n_cpus)
print(f"PyTorch using {n_cpus} CPU threads")

BASE      = '/data/ebaird/scRNAseq/20260522.nanoCT'
GTF       = f'{BASE}/colleague_analysis/tutorial/dmel-all-r6.59.snapatac2.gtf'
OUT_DIR   = f'{BASE}/scglue_peaks'
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
print("Loading peaks chromatin object...")
chrom = ad.read_h5ad(f'{BASE}/analysis_05.26/combined_ATACuniverse_binary.h5ad')
print(chrom)

# Parse chr/start/end from peak_ids (format: chr2L:12717-12917)
def parse_peak_id(pid):
    m = re.match(r'(chr\w+):(\d+)-(\d+)', str(pid))
    if m:
        return m.group(1), int(m.group(2)), int(m.group(3))
    return None, None, None

chr_list, start_list, end_list = [], [], []
for pid in chrom.var_names:
    c, s, e = parse_peak_id(pid)
    chr_list.append(c)
    start_list.append(s)
    end_list.append(e)
chrom.var['chr']   = chr_list
chrom.var['start'] = start_list
chrom.var['end']   = end_list

# Drop any peaks that couldn't be parsed
valid = chrom.var['chr'].notna()
if (~valid).sum() > 0:
    print(f"  Dropping {(~valid).sum()} unparseable peaks")
    chrom = chrom[:, valid].copy()
chrom.var['chr']   = chrom.var['chr'].astype(str)
chrom.var['start'] = chrom.var['start'].astype(int)
chrom.var['end']   = chrom.var['end'].astype(int)
print(f"  Parsed {chrom.n_vars} peaks, {chrom.var['chr'].nunique()} chromosomes")

print("\nLoading RNA object...")
rna_full = ad.read_h5ad(f'{BASE}/rna_SCT.h5ad')
print(rna_full)

# ---------------------------------------------------------------------------
# 2. Subsample RNA
# ---------------------------------------------------------------------------
np.random.seed(42)
cluster_col = 'cluster' if 'cluster' in rna_full.obs.columns else 'seurat_clusters'
labels_full = rna_full.obs[cluster_col].values
n_target    = 20000
unique_cls  = np.unique(labels_full)
per_cl      = max(1, n_target // len(unique_cls))

keep_idx = []
for cl in unique_cls:
    idx = np.where(labels_full == cl)[0]
    keep_idx.extend(np.random.choice(idx, min(per_cl, len(idx)), replace=False).tolist())
keep_idx = np.array(sorted(keep_idx))

rna = rna_full[keep_idx].copy()
del rna_full
print(f"\nRNA subsampled: {rna.n_obs} cells ({len(unique_cls)} clusters, ≤{per_cl}/cluster)")

# ---------------------------------------------------------------------------
# 3. Prepare RNA
# ---------------------------------------------------------------------------
print("RNA layers:", list(rna.layers.keys()))
rna.X = rna.layers['counts'].copy()
rna.X = rna.X.astype(np.float32)
if sps.issparse(rna.X):
    sc.pp.normalize_total(rna, target_sum=1e4)
    sc.pp.log1p(rna)
else:
    rna.X = rna.X / rna.X.sum(axis=1, keepdims=True) * 1e4
    rna.X = np.log1p(rna.X)

sc.pp.highly_variable_genes(rna, n_top_genes=4000, flavor='cell_ranger')
sc.tl.pca(rna, n_comps=50, use_highly_variable=True, svd_solver='auto')
print(f"RNA HVGs: {rna.var['highly_variable'].sum()}")
RNA_DIST = 'Normal'
print(f"RNA distribution: {RNA_DIST}")

# ---------------------------------------------------------------------------
# 4. Prepare chromatin — H3K27ac TF-IDF LSI
# ---------------------------------------------------------------------------
chrom.X = chrom.layers['acet'].copy()
if sps.issparse(chrom.X):
    chrom.X = chrom.X.astype(np.float32)

from sklearn.utils.extmath import randomized_svd
from scipy.sparse import diags
from scipy.stats import pearsonr

print("Computing H3K27ac TF-IDF LSI...")
X_acet = chrom.X if sps.issparse(chrom.X) else sps.csr_matrix(chrom.X)

cell_totals = np.array(X_acet.sum(axis=1)).flatten()
cell_totals[cell_totals == 0] = 1
TF    = diags(1.0 / cell_totals) @ X_acet
n_cells_c = X_acet.shape[0]
bin_counts = np.array((X_acet > 0).sum(axis=0)).flatten()
idf   = np.log1p(n_cells_c / (1.0 + bin_counts))
TFIDF = TF @ diags(idf)

U, S, Vt = randomized_svd(TFIDF, n_components=50, random_state=42, n_iter=5)
X_lsi = U * S

depth_corr, _ = pearsonr(np.log1p(cell_totals), X_lsi[:, 0])
print(f"  LSI component 0 vs depth correlation: {depth_corr:.3f}")
if abs(depth_corr) > 0.5:
    print("  -> Dropping component 0 (depth-correlated)")
    X_lsi = X_lsi[:, 1:]

chrom.obsm['X_lsi_acet'] = X_lsi.astype(np.float32)
print(f"  H3K27ac LSI embedding: {chrom.obsm['X_lsi_acet'].shape}")

# ---------------------------------------------------------------------------
# 5. Parse GTF
# ---------------------------------------------------------------------------
def parse_gtf_gene_coords(gtf_path):
    coords = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'gene':
                continue
            chrom_g, start, end, strand, attrs = (
                fields[0], int(fields[3]), int(fields[4]), fields[6], fields[8]
            )
            name = None
            for key in ('gene_symbol', 'gene_name', 'gene_id'):
                for part in attrs.split(';'):
                    part = part.strip()
                    if part.startswith(key + ' "'):
                        name = part.split('"')[1]
                        break
                if name:
                    break
            if name is None:
                continue
            tss    = start if strand == '+' else end
            chrom_p = chrom_g if chrom_g.startswith('chr') else 'chr' + chrom_g
            coords[name] = (chrom_p, start, end, tss, strand)
    return coords

print("Parsing GTF...")
gene_coords = parse_gtf_gene_coords(GTF)
print(f"Genes in GTF: {len(gene_coords)}")

chrs, starts, ends, promoters, found = [], [], [], [], []
for gene in rna.var_names:
    if gene in gene_coords:
        ch, s, e, tss, strand = gene_coords[gene]
        chrs.append(ch); starts.append(s); ends.append(e)
        promoters.append(tss); found.append(True)
    else:
        chrs.append(None); starts.append(None); ends.append(None)
        promoters.append(None); found.append(False)

rna.var['chr']      = chrs
rna.var['start']    = starts
rna.var['end']      = ends
rna.var['promoter'] = promoters
rna.var['in_gtf']   = found
print(f"RNA genes in GTF: {sum(found)}/{len(found)}")

rna = rna[:, rna.var['in_gtf']].copy()
rna.var['chr']      = rna.var['chr'].astype(str)
rna.var['start']    = rna.var['start'].astype(int)
rna.var['end']      = rna.var['end'].astype(int)
rna.var['promoter'] = rna.var['promoter'].astype(int)

# ---------------------------------------------------------------------------
# 6. Build guidance graph — 150 kb window
# ---------------------------------------------------------------------------
WINDOW        = 150_000
PROMOTER_HALF = 500

print(f"\nBuilding guidance graph (window={WINDOW//1000}kb)...")

def build_guidance_graph_python(rna_var, chrom_var, window=WINDOW, promoter_half=PROMOTER_HALF):
    G = nx.DiGraph()
    bin_names  = chrom_var.index.to_numpy()
    bin_chrs   = chrom_var['chr'].to_numpy()
    bin_starts = chrom_var['start'].to_numpy(dtype=int)
    bin_ends   = chrom_var['end'].to_numpy(dtype=int)
    bin_mids   = (bin_starts + bin_ends) // 2

    bin_by_chr = {}
    for i, (ch, bn) in enumerate(zip(bin_chrs, bin_names)):
        bin_by_chr.setdefault(ch, []).append(i)

    gene_names = rna_var.index.to_numpy()
    gene_chrs  = rna_var['chr'].to_numpy()
    gene_tss   = rna_var['promoter'].to_numpy(dtype=int)

    for gene, ch, tss in zip(gene_names, gene_chrs, gene_tss):
        if ch not in bin_by_chr:
            continue
        idxs = bin_by_chr[ch]
        mids = bin_mids[idxs]
        dists = np.abs(mids - tss)
        for ci in np.where(dists <= window)[0]:
            bin_idx = idxs[ci]
            d = int(dists[ci])
            w = float(np.power((d + promoter_half) / promoter_half, -0.75))
            G.add_edge(gene, bin_names[bin_idx],
                       dist=d, weight=w, type='dist', sign=1)
    return G

dist_graph = build_guidance_graph_python(rna.var, chrom.var)
print(f"Guidance graph: {dist_graph.number_of_nodes()} nodes, {dist_graph.number_of_edges()} edges")

guidance = scglue.graph.compose_multigraph(dist_graph, dist_graph.reverse())

all_features = list(rna.var_names) + list(chrom.var_names)
rna.var['connected']   = [guidance.has_node(v) for v in rna.var_names]
chrom.var['connected'] = [guidance.has_node(v) for v in chrom.var_names]
print(f"RNA genes connected: {rna.var['connected'].sum()}/{rna.n_vars}")
print(f"Chrom peaks connected: {chrom.var['connected'].sum()}/{chrom.n_vars}")

rna   = rna[:, rna.var['connected']].copy()
chrom = chrom[:, chrom.var['connected']].copy()

all_ids = list(rna.var_names) + list(chrom.var_names)
guidance = guidance.subgraph(nodes=all_ids).copy()
for item in all_ids:
    guidance.add_edge(item, item, weight=1.0, type='self-loop')
nx.set_edge_attributes(guidance, 1, 'sign')

scglue.graph.check_graph(guidance, [rna, chrom])
print("Graph check passed")

# ---------------------------------------------------------------------------
# 7. Configure and train SCGLUE
# ---------------------------------------------------------------------------
print("\nConfiguring SCGLUE datasets...")
scglue.models.configure_dataset(
    rna, RNA_DIST,
    use_highly_variable=True,
    use_rep='X_pca'
)
scglue.models.configure_dataset(
    chrom, 'NB',
    use_highly_variable=False,
    use_rep='X_lsi_acet'
)

print("Training SCGLUE model...")
MODEL_DIR = f'{OUT_DIR}/model'
os.makedirs(MODEL_DIR, exist_ok=True)
glue = scglue.models.fit_SCGLUE(
    {'rna': rna, 'chrom': chrom},
    guidance,
    fit_kws={'directory': MODEL_DIR, 'data_batch_size': 512, 'max_epochs': 150},
    init_kws={'latent_dim': 20, 'h_depth': 2, 'h_dim': 256}
)
glue.save(f'{MODEL_DIR}/model.dill')
print("Model saved")

# ---------------------------------------------------------------------------
# 8. Extract embeddings
# ---------------------------------------------------------------------------
print("\nExtracting embeddings...")
rna.obsm['X_glue']   = glue.encode_data('rna',   rna)
chrom.obsm['X_glue'] = glue.encode_data('chrom', chrom)

# Alignment sanity check
from scipy.spatial.distance import cdist
from sklearn.preprocessing import normalize as sk_normalize

rna_emb_check   = sk_normalize(rna.obsm['X_glue'][:500].astype(np.float32),   norm='l2')
chrom_emb_check = sk_normalize(chrom.obsm['X_glue'][:500].astype(np.float32), norm='l2')
within_rna   = cdist(rna_emb_check,   rna_emb_check,   metric='cosine').mean()
within_chrom = cdist(chrom_emb_check, chrom_emb_check, metric='cosine').mean()
cross_modal  = cdist(chrom_emb_check, rna_emb_check,   metric='cosine').mean()
print(f"\nAlignment sanity check:")
print(f"  Within-RNA cosine distance:    {within_rna:.3f}")
print(f"  Within-chrom cosine distance:  {within_chrom:.3f}")
print(f"  Cross-modal cosine distance:   {cross_modal:.3f}")
if cross_modal < 0.5 * (within_rna + within_chrom):
    print("  -> Good alignment: cross-modal distance << within-modal")
else:
    print("  -> WARNING: poor alignment — label transfer may be unreliable")

# ---------------------------------------------------------------------------
# 9. k-NN label transfer
# ---------------------------------------------------------------------------
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize

cluster_col = 'cluster' if 'cluster' in rna.obs.columns else 'seurat_clusters'
print(f"Using RNA cluster column: '{cluster_col}'")

K = 15
rna_emb   = normalize(rna.obsm['X_glue'].astype(np.float32),   norm='l2')
chrom_emb = normalize(chrom.obsm['X_glue'].astype(np.float32), norm='l2')

labels = rna.obs[cluster_col].values

nn = NearestNeighbors(n_neighbors=K, metric='cosine', algorithm='brute', n_jobs=-1)
nn.fit(rna_emb)
distances, indices = nn.kneighbors(chrom_emb)

transferred, confidence = [], []
for idx_row in indices:
    votes = labels[idx_row]
    unique, counts = np.unique(votes, return_counts=True)
    winner = unique[np.argmax(counts)]
    transferred.append(winner)
    confidence.append(counts.max() / K)

chrom.obs['rna_label_transfer']  = transferred
chrom.obs['transfer_confidence'] = confidence
chrom.obs['mean_knn_distance']   = distances.mean(axis=1)

print("\nTransfer label distribution:")
print(chrom.obs['rna_label_transfer'].value_counts())
print(f"\nMean transfer confidence: {np.mean(confidence):.3f}")
print(f"Mean KNN cosine distance: {distances.mean():.3f}")

# Leiden -> RNA fractions
frac = pd.crosstab(chrom.obs['leiden'], chrom.obs['rna_label_transfer'],
                   normalize='index').round(2)
print("\nLeiden -> RNA label fractions:")
print(frac.to_string())

# ---------------------------------------------------------------------------
# 10. Summary UMAP
# ---------------------------------------------------------------------------
umap = chrom.obsm['X_umap']
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

sc_kwargs = dict(s=3, alpha=0.6, rasterized=True, linewidths=0)

axes[0].scatter(umap[:, 0], umap[:, 1], c='#cccccc', **sc_kwargs)
axes[0].set_title('All chrom cells'); axes[0].set_xticks([]); axes[0].set_yticks([])

unique_rna = sorted(chrom.obs['rna_label_transfer'].unique(), key=lambda x: int(x[1:]))
cmap_l = plt.cm.get_cmap('tab20', len(unique_rna))
for i, cl in enumerate(unique_rna):
    m = chrom.obs['rna_label_transfer'] == cl
    axes[1].scatter(umap[m, 0], umap[m, 1], c=[cmap_l(i)], label=cl, **sc_kwargs)
axes[1].legend(markerscale=3, fontsize=7, ncol=2, loc='best')
axes[1].set_title('RNA label transfer (peaks)')
axes[1].set_xticks([]); axes[1].set_yticks([])

sc2 = axes[2].scatter(umap[:, 0], umap[:, 1],
                      c=chrom.obs['transfer_confidence'], cmap='viridis',
                      vmin=0, vmax=1, **sc_kwargs)
plt.colorbar(sc2, ax=axes[2], label='Confidence')
axes[2].set_title('Transfer confidence')
axes[2].set_xticks([]); axes[2].set_yticks([])

fig.savefig(f'{BASE}/scglue_peaks_integration_umap.png', dpi=150, bbox_inches='tight')
print(f"\nSaved scglue_peaks_integration_umap.png")

# ---------------------------------------------------------------------------
# 11. Per-cluster UMAP
# ---------------------------------------------------------------------------
clusters = sorted(chrom.obs['rna_label_transfer'].unique(), key=lambda x: int(x[1:]))
conf_arr = chrom.obs['transfer_confidence'].values.astype(float)
n        = len(clusters)
ncols    = 5
nrows    = int(np.ceil(n / ncols))
cmap_label = plt.cm.get_cmap('tab20', n)
label_colours = {cl: cmap_label(i) for i, cl in enumerate(clusters)}

fig2, axes2 = plt.subplots(nrows, ncols,
                            figsize=(ncols * 4, nrows * 3.8),
                            constrained_layout=True)
fig2.suptitle('RNA label transfer (peaks) — per cluster', fontsize=14, y=1.01)

for ax_idx, cl in enumerate(clusters):
    ax = axes2[ax_idx // ncols][ax_idx % ncols]
    mask = np.array(chrom.obs['rna_label_transfer'] == cl)
    n_cl = mask.sum()
    ax.scatter(umap[~mask, 0], umap[~mask, 1],
               c='#d4d4d4', s=1, alpha=0.3, rasterized=True, linewidths=0)
    sc3 = ax.scatter(umap[mask, 0], umap[mask, 1],
                     c=conf_arr[mask], cmap='viridis',
                     vmin=0, vmax=1, s=4, alpha=0.7, rasterized=True, linewidths=0)
    ax.set_title(f'{cl} (n={n_cl:,})', color=label_colours[cl],
                 fontsize=10, fontweight='bold')
    ax.set_xticks([]); ax.set_yticks([])

for ax_idx in range(n, nrows * ncols):
    axes2[ax_idx // ncols][ax_idx % ncols].set_visible(False)

cbar_ax = fig2.add_axes([1.01, 0.15, 0.015, 0.7])
sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=0, vmax=1))
sm.set_array([])
cbar = fig2.colorbar(sm, cax=cbar_ax)
cbar.set_label('Transfer confidence', fontsize=11)

fig2.savefig(f'{BASE}/scglue_peaks_per_cluster_umap.png', dpi=150, bbox_inches='tight')
print(f"Saved scglue_peaks_per_cluster_umap.png")

# ---------------------------------------------------------------------------
# 12. Save
# ---------------------------------------------------------------------------
chrom.write_h5ad(f'{OUT_DIR}/chrom.h5ad')
rna.write_h5ad(f'{OUT_DIR}/rna.h5ad')
print(f"Saved {OUT_DIR}/chrom.h5ad and rna.h5ad")
print("\nDone.")
