#!/usr/bin/env python3
"""
SCGLUE integration: nanoCT H3K27ac chromatin (workshop object) + SCENTINELsep24 scRNA-seq.

Inputs:
  - combined_dim_reduced.h5ad  : 17,704 cells, bins as vars (chr-prefixed), layers acet/meth
  - rna_SCT.h5ad               : SCT assay from Seurat RDS (corrected counts + log-normalized data)
  - dmel-all-r6.59.snapatac2.gtf

Outputs:
  - scglue_model/              : trained SCGLUE model
  - scglue_chrom.h5ad          : chromatin cells with X_glue embedding + label transfer
  - scglue_rna.h5ad            : RNA cells with X_glue embedding
"""

import os
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

# Use all available CPUs
n_cpus = int(os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count()))
torch.set_num_threads(n_cpus)
print(f"PyTorch using {n_cpus} CPU threads")

BASE   = '/data/ebaird/scRNAseq/20260522.nanoCT'
GTF    = f'{BASE}/colleague_analysis/tutorial/dmel-all-r6.59.snapatac2.gtf'
MODEL_DIR = f'{BASE}/scglue_model'
os.makedirs(MODEL_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
print("Loading chromatin object...")
chrom = ad.read_h5ad(f'{BASE}/analysis_05.26/combined_dim_reduced.h5ad')
print(chrom)

print("\nLoading RNA object...")
rna = ad.read_h5ad(f'{BASE}/rna_SCT.h5ad')
print(rna)

# ---------------------------------------------------------------------------
# 2. Prepare RNA — SCT assay: counts slot = corrected counts (NB), data = log-norm (PCA)
# ---------------------------------------------------------------------------
# X from reconstruct_h5ad.py is the 'data' slot (log-normalized) — stored as X
# 'counts' layer = SCT corrected counts (integer, suitable for NB model)
print("RNA layers:", list(rna.layers.keys()))

if 'counts' in rna.layers:
    # Use corrected counts as X for NB model
    rna.X = rna.layers['counts'].copy()
    print(f"Using SCT corrected counts as X  (dtype: {rna.X.dtype})")
    # Compute PCA on log-normalized data for use_rep
    rna_pp = rna.copy()
    rna_pp.X = rna.layers['data'].copy() if 'data' in rna.layers else rna.X.copy()
    sc.pp.highly_variable_genes(rna_pp, n_top_genes=4000, flavor='cell_ranger')
    sc.tl.pca(rna_pp, n_comps=50, use_highly_variable=True, svd_solver='auto')
    rna.obsm['X_pca'] = rna_pp.obsm['X_pca']
    rna.var['highly_variable'] = rna_pp.var['highly_variable']
    del rna_pp
    print(f"RNA HVGs: {rna.var['highly_variable'].sum()}")
    RNA_DIST = 'NB'
else:
    # Fallback: only log-normalized data available, use Normal distribution
    print(f"No counts layer — using X as-is with Normal distribution")
    if 'X_pca' not in rna.obsm:
        sc.pp.highly_variable_genes(rna, n_top_genes=4000, flavor='cell_ranger')
        sc.tl.pca(rna, n_comps=50, use_highly_variable=True, svd_solver='auto')
    else:
        print("Using pre-computed X_pca from Seurat")
    RNA_DIST = 'Normal'

rna.X = rna.X.astype(np.float32) if not sps.issparse(rna.X) else rna.X.astype(np.float32)
print(f"RNA distribution: {RNA_DIST}")

# ---------------------------------------------------------------------------
# 3. Prepare chromatin — use acet layer as the primary ATAC-like signal
# ---------------------------------------------------------------------------
# Set X to acet counts (raw bin counts, already integers from fragment loading)
chrom.X = chrom.layers['acet'].copy()
if sps.issparse(chrom.X):
    chrom.X = chrom.X.astype(np.float32)

# chrom var already has chr, start, end columns from the workshop object
print("Chromatin var columns:", list(chrom.var.columns))

# ---------------------------------------------------------------------------
# 4. Parse GTF to get gene TSS coordinates for RNA var
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
            tss = start if strand == '+' else end
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
        chrs.append(ch)
        starts.append(s)
        ends.append(e)
        promoters.append(tss)
        found.append(True)
    else:
        chrs.append(None); starts.append(None); ends.append(None)
        promoters.append(None); found.append(False)

rna.var['chr']      = chrs
rna.var['start']    = starts
rna.var['end']      = ends
rna.var['promoter'] = promoters
rna.var['in_gtf']   = found
print(f"RNA genes in GTF: {sum(found)}/{len(found)}")

# Filter RNA to genes with coordinates
rna = rna[:, rna.var['in_gtf']].copy()
rna.var['chr']      = rna.var['chr'].astype(str)
rna.var['start']    = rna.var['start'].astype(int)
rna.var['end']      = rna.var['end'].astype(int)
rna.var['promoter'] = rna.var['promoter'].astype(int)

# ---------------------------------------------------------------------------
# 5. Build guidance graph (pure Python — avoids bedtools version requirement)
# ---------------------------------------------------------------------------
print("\nBuilding guidance graph...")

WINDOW = 20000   # bp upstream/downstream of promoter
PROMOTER_HALF = 500  # bp for promoter definition

def build_guidance_graph_python(rna_var, chrom_var, window=WINDOW, promoter_half=PROMOTER_HALF):
    """
    Link each RNA gene to chromatin bins whose midpoint falls within
    `window` bp of the gene's TSS. Edge weight: power-law decay.
    Returns a directed NetworkX graph (gene -> bin).
    """
    G = nx.DiGraph()

    # Index bins by chromosome
    bin_by_chr = {}
    bin_names = chrom_var.index.to_numpy()
    bin_chrs  = chrom_var['chr'].to_numpy()
    bin_starts = chrom_var['start'].to_numpy(dtype=int)
    bin_ends   = chrom_var['end'].to_numpy(dtype=int)
    bin_mids   = (bin_starts + bin_ends) // 2

    for i, (ch, bn) in enumerate(zip(bin_chrs, bin_names)):
        bin_by_chr.setdefault(ch, []).append(i)

    gene_names  = rna_var.index.to_numpy()
    gene_chrs   = rna_var['chr'].to_numpy()
    gene_tss    = rna_var['promoter'].to_numpy(dtype=int)

    edge_count = 0
    for gene, ch, tss in zip(gene_names, gene_chrs, gene_tss):
        if ch not in bin_by_chr:
            continue
        idxs = bin_by_chr[ch]
        mids = bin_mids[idxs]
        dists = np.abs(mids - tss)
        close = np.where(dists <= window)[0]
        for ci in close:
            bin_idx = idxs[ci]
            d = int(dists[ci])
            w = float(np.power((d + promoter_half) / promoter_half, -0.75))
            G.add_edge(gene, bin_names[bin_idx],
                       dist=d, weight=w, type='dist', sign=1)
            edge_count += 1

    return G

dist_graph = build_guidance_graph_python(rna.var, chrom.var)
print(f"Guidance graph: {dist_graph.number_of_nodes()} nodes, {dist_graph.number_of_edges()} edges")

# Bidirectional: compose with reverse
guidance = scglue.graph.compose_multigraph(dist_graph, dist_graph.reverse())

# Keep only nodes present in our datasets, add self-loops
all_features = list(rna.var_names) + list(chrom.var_names)
rna.var['connected']   = [guidance.has_node(v) for v in rna.var_names]
chrom.var['connected'] = [guidance.has_node(v) for v in chrom.var_names]
print(f"RNA genes connected: {rna.var['connected'].sum()}/{rna.n_vars}")
print(f"Chrom bins connected: {chrom.var['connected'].sum()}/{chrom.n_vars}")

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
# 6. Configure and train SCGLUE
# ---------------------------------------------------------------------------
print("\nConfiguring SCGLUE datasets...")
scglue.models.configure_dataset(
    rna, RNA_DIST,
    use_highly_variable=(RNA_DIST == 'NB'),
    use_rep='X_pca'
)
scglue.models.configure_dataset(
    chrom, 'NB',
    use_highly_variable=False,
    use_rep='X_multi_spectral'
)

print("Training SCGLUE model...")
glue = scglue.models.fit_SCGLUE(
    {'rna': rna, 'chrom': chrom},
    guidance,
    fit_kws={'directory': MODEL_DIR, 'data_batch_size': 512, 'max_epochs': 150},
    init_kws={'latent_dim': 20, 'h_depth': 2, 'h_dim': 256}
)
glue.save(f'{MODEL_DIR}/model.dill')
print("Model saved")

# ---------------------------------------------------------------------------
# 7. Extract joint embeddings
# ---------------------------------------------------------------------------
print("\nExtracting embeddings...")
rna.obsm['X_glue']   = glue.encode_data('rna',   rna)
chrom.obsm['X_glue'] = glue.encode_data('chrom', chrom)

# ---------------------------------------------------------------------------
# 7b. Sanity check — cross-modal alignment quality
# ---------------------------------------------------------------------------
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
# 8. k-NN label transfer in SCGLUE latent space
# ---------------------------------------------------------------------------
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize

# Determine RNA cluster column
cluster_col = 'cluster' if 'cluster' in rna.obs.columns else rna.obs.columns[0]
print(f"Using RNA cluster column: '{cluster_col}'")

K = 15
rna_emb   = normalize(rna.obsm['X_glue'].astype(np.float32),   norm='l2')
chrom_emb = normalize(chrom.obsm['X_glue'].astype(np.float32), norm='l2')

# Subsample RNA to ≤500 cells per cluster for speed
np.random.seed(42)
labels = rna.obs[cluster_col].values
sample_idx = []
for cl in np.unique(labels):
    idx = np.where(labels == cl)[0]
    sample_idx.extend(np.random.choice(idx, min(500, len(idx)), replace=False).tolist())
sample_idx = np.array(sample_idx)

nn = NearestNeighbors(n_neighbors=K, metric='cosine', algorithm='brute', n_jobs=-1)
nn.fit(rna_emb[sample_idx])
distances, indices = nn.kneighbors(chrom_emb)

rna_labels_sub = labels[sample_idx]
transferred, confidence = [], []
for idx_row in indices:
    votes = rna_labels_sub[idx_row]
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

# Crosstab
import pandas as pd
xtab = pd.crosstab(
    chrom.obs['leiden'].astype(str),
    chrom.obs['rna_label_transfer'],
    normalize='index'
).round(2)
print("\nLeiden -> RNA label fractions:")
print(xtab)

# ---------------------------------------------------------------------------
# 9. Plot
# ---------------------------------------------------------------------------
umap = chrom.obsm['X_umap']
fig, axes = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=True)

for cl in sorted(chrom.obs['leiden'].astype(str).unique()):
    mask = chrom.obs['leiden'].astype(str) == cl
    axes[0].scatter(umap[mask, 0], umap[mask, 1], s=1, label=cl, alpha=0.5)
axes[0].set_title('Leiden clusters (chromatin)')
axes[0].legend(markerscale=5, fontsize=7)

for cl in sorted(chrom.obs['rna_label_transfer'].unique()):
    mask = chrom.obs['rna_label_transfer'] == cl
    axes[1].scatter(umap[mask, 0], umap[mask, 1], s=1, label=cl, alpha=0.5)
axes[1].set_title('RNA label transfer (SCGLUE)')
axes[1].legend(markerscale=5, fontsize=7)

sc2 = axes[2].scatter(umap[:, 0], umap[:, 1], s=1,
                      c=chrom.obs['transfer_confidence'], cmap='viridis', alpha=0.5)
fig.colorbar(sc2, ax=axes[2])
axes[2].set_title('Transfer confidence')

plt.savefig(f'{BASE}/scglue_integration_umap.png', dpi=150)
print(f"\nSaved plot to {BASE}/scglue_integration_umap.png")

# ---------------------------------------------------------------------------
# 10. Save
# ---------------------------------------------------------------------------
chrom.write_h5ad(f'{BASE}/scglue_chrom.h5ad')
rna.write_h5ad(f'{BASE}/scglue_rna.h5ad')
print(f"Saved scglue_chrom.h5ad and scglue_rna.h5ad")
print("\nDone.")
