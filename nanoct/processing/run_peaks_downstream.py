"""
Full downstream analysis on opt2 low-q MACS3 peaks:
  1. Build peak→gene regulatory links from TSS proximity (dm6)
  2. Run embedding → binarization → gene inference → marker calling
  3. Compare marker specificity with bins analysis

Goal: show whether peaks give a better analysis than 5kb bins.
"""
import sys
sys.path.insert(0, '/data/ebaird/scRNAseq/20260522.nanoCT/scit_src')

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import pandas as pd
import gzip, io, os
import anndata as ad
import scipy.sparse as sps
import src as scit
from src.tools._regulation import RegInference, get_regulatory_matrix

BASE     = "/data/ebaird/scRNAseq/20260522.nanoCT/SU.analysis.2026.05.22/Vasso_nanoCT_nanoscope"
LOWQ_DIR = "/data/ebaird/scRNAseq/20260522.nanoCT/analysis_05.26/macs3_lowq"
OUT      = "/data/ebaird/scRNAseq/20260522.nanoCT/analysis_05.26"
PLOT_DIR = f"{OUT}/pipeline_plots"
BINS_H5  = f"{OUT}/combined_dim_reduced.h5ad"

# ── helpers (fragment counting) ───────────────────────────────────────────────
def load_bed(path):
    return pl.read_csv(path, separator='\t', has_header=False,
                       new_columns=['chr', 'start', 'end', 'name'])

def merge_peaks(dfs):
    combined = pl.concat([d.select(['chr','start','end']) for d in dfs]).sort(['chr','start'])
    rows, cur_chr, cur_start, cur_end = [], None, None, None
    for chrom, start, end in combined.iter_rows():
        if chrom != cur_chr or start > cur_end:
            if cur_chr is not None: rows.append((cur_chr, cur_start, cur_end))
            cur_chr, cur_start, cur_end = chrom, start, end
        else: cur_end = max(cur_end, end)
    if cur_chr is not None: rows.append((cur_chr, cur_start, cur_end))
    merged = pl.DataFrame(rows, schema={'chr':pl.Utf8,'start':pl.Int64,'end':pl.Int64}, orient='row')
    return merged.with_columns(
        (pl.col('chr')+':'+pl.col('start').cast(pl.Utf8)+'-'+pl.col('end').cast(pl.Utf8)).alias('peak_id')
    )

def _build_peak_index(peaks_df):
    idx = {}
    for chrom in peaks_df['chr'].unique().to_list():
        sub = peaks_df.filter(pl.col('chr') == chrom).sort('start')
        idx[chrom] = (sub['start'].to_numpy(), sub['end'].to_numpy())
    return idx

def _pos_to_peak(pos, starts, ends):
    cand = np.searchsorted(starts, pos, side='right') - 1
    valid = (cand >= 0) & (pos < ends[np.clip(cand, 0, len(ends)-1)])
    return np.where(valid, cand, -1)

def count_fragments_in_peaks(fragments_path, peaks_df, batch_size=400_000):
    peak_idx = _build_peak_index(peaks_df)
    peak_ids = peaks_df['peak_id'].to_list()
    n_peaks  = len(peak_ids)
    read_kw  = dict(separator='\t', has_header=False,
                    new_columns=['chr','start','end','bc','readSupport'])
    with gzip.open(fragments_path, 'rb') as gz:
        buf = io.BytesIO(b''.join(l for l in gz if not l.startswith(b'#')))
    df_full   = pl.read_csv(buf, **read_kw)
    n_rows    = df_full.height
    bcs       = np.sort(df_full['bc'].unique().to_numpy())
    bc_to_row = {b: i for i, b in enumerate(bcs)}
    rows_list, cols_list, data_list = [], [], []
    for batch_start in range(0, n_rows, batch_size):
        batch     = df_full.slice(batch_start, batch_size)
        bc_arr    = batch['bc'].to_numpy()
        start_arr = batch['start'].to_numpy()
        end_arr   = batch['end'].to_numpy()
        chr_arr   = batch['chr'].to_numpy()
        for chrom in np.unique(chr_arr):
            if chrom not in peak_idx: continue
            starts_p, ends_p = peak_idx[chrom]
            mask    = chr_arr == chrom
            bc_rows = np.array([bc_to_row[b] for b in bc_arr[mask]], dtype=np.int32)
            pc_s    = _pos_to_peak(start_arr[mask], starts_p, ends_p)
            pc_e    = _pos_to_peak(end_arr[mask],   starts_p, ends_p)
            hit_s   = pc_s >= 0
            rows_list.append(bc_rows[hit_s]);  cols_list.append(pc_s[hit_s])
            data_list.append(np.ones(hit_s.sum(), dtype=np.uint32))
            hit_e   = (pc_e >= 0) & (pc_e != pc_s)
            rows_list.append(bc_rows[hit_e]);  cols_list.append(pc_e[hit_e])
            data_list.append(np.ones(hit_e.sum(), dtype=np.uint32))
        if batch_start % (batch_size * 5) == 0:
            print(f"    {100*batch_start/n_rows:.0f}%...", end='\r')
    print()
    rows_arr = np.concatenate(rows_list).astype(np.int32)
    cols_arr = np.concatenate(cols_list).astype(np.int32)
    data_arr = np.concatenate(data_list).astype(np.uint32)
    X = sps.coo_matrix((data_arr, (rows_arr, cols_arr)), shape=(len(bcs), n_peaks)).tocsr()
    a = ad.AnnData(X); a.obs.index = bcs; a.var.index = peak_ids
    return a


# ── Step 1: build peak → gene regulatory links ────────────────────────────────
print("=== Step 1: Building peak → gene regulatory links ===")
print("  Loading RNA reference...")
rna = ad.read_h5ad(f"{OUT}/rna.h5ad")
print(f"  RNA: {rna.n_vars} genes")

# Load MACS3 low-q peaks
ac_lowq = load_bed(f"{LOWQ_DIR}/H3K27ac_lowq_sorted.bed")
me_lowq = load_bed(f"{LOWQ_DIR}/H3K27me3_lowq_sorted.bed")
union_peaks = merge_peaks([ac_lowq, me_lowq])
print(f"  MACS3 peaks: {len(union_peaks)}")

# Build peak→gene links using TSS proximity
# RNA var uses chr-prefixed names (chr2L); peaks use bare names (2L) — normalise peaks→chr-prefix
# Links format (must match scit.ld.regulatory_links exactly):
#   col1-3: chr, promoter-500, promoter+500  (matches RNA var)
#   col4-6: peak chr, start, end
#   col7-9: score, pval, qval
MAX_DIST = 100_000  # bp from TSS — covers most cis-regulatory elements
PROMOTER_SIZE = 500

print(f"  Building proximity links (max distance: {MAX_DIST//1000} kb)...")
genes_df = rna.var[['chr', 'start', 'end', 'promoter']].copy()
genes_df = genes_df[genes_df['chr'].notna() & genes_df['promoter'].notna()].copy()
genes_df['prom_start'] = np.clip(genes_df['promoter'].astype(int) - PROMOTER_SIZE, 0, None)
genes_df['prom_end']   = genes_df['promoter'].astype(int) + PROMOTER_SIZE

# Build sorted index of peaks per chromosome — add 'chr' prefix to match RNA var
peaks_by_chr = {}
for chrom in union_peaks['chr'].unique().to_list():
    sub = union_peaks.filter(pl.col('chr') == chrom).sort('start')
    mid = ((sub['start'] + sub['end']) / 2).to_numpy()
    chr_key = chrom if chrom.startswith('chr') else f'chr{chrom}'
    peaks_by_chr[chr_key] = (sub['start'].to_numpy(), sub['end'].to_numpy(),
                              mid, [chr_key]*len(sub), sub['peak_id'].to_list())

link_rows = []
for gene_name, row in genes_df.iterrows():
    chrom = row['chr']
    if chrom not in peaks_by_chr: continue
    tss = int(row['promoter'])
    prom_s = int(row['prom_start'])
    prom_e = int(row['prom_end'])
    starts_p, ends_p, mids, chrs_p, ids_p = peaks_by_chr[chrom]

    # find peaks whose midpoint is within MAX_DIST of TSS
    lo = np.searchsorted(mids, tss - MAX_DIST, side='left')
    hi = np.searchsorted(mids, tss + MAX_DIST, side='right')
    for j in range(lo, hi):
        dist = abs(mids[j] - tss)
        score = np.exp(-dist / 50_000)  # distance-weighted score
        link_rows.append((chrom, prom_s, prom_e,
                          chrs_p[j], int(starts_p[j]), int(ends_p[j]),
                          score, 0.001, 0.01))

links_df = pd.DataFrame(link_rows,
    columns=['chr', 'start', 'end', 'chr2', 'start2', 'end2', 'score', 'pval', 'qval'])
print(f"  Total peak→gene links: {len(links_df):,}")
print(f"  Unique genes with links: {links_df.groupby(['chr','start','end']).ngroups:,}")

# Save links file in bin.links format
peak_links_path = f"{OUT}/peak.links"
links_df.to_csv(peak_links_path, sep='\t', header=False, index=False)
print(f"  Saved: {peak_links_path}")


# ── Step 2: rebuild opt2 adata with full pipeline ─────────────────────────────
print("\n=== Step 2: opt2 embedding (proper pipeline) ===")
print("  Counting H3K27ac...")
ac = count_fragments_in_peaks(f"{BASE}/H3K27ac/fragments.tsv.gz", union_peaks)
print("  Counting H3K27me3...")
me = count_fragments_in_peaks(f"{BASE}/H3K27me3/fragments.tsv.gz", union_peaks)

adata = scit.tl.stack_adata([ac, me], ['acet', 'meth'])
scit.tl.add_metadata(adata)
adata = scit.tl.filter(adata, ['acet','meth'],
                       min_obs_counts=[20,15], max_obs_counts=[3000,3000], return_purged=True)
scit.tl.add_metadata(adata)
adata = scit.tl.filter(adata, ['acet','meth'], min_var_counts=[5,5], return_purged=True)
scit.tl.filter(adata, ['acet','meth'], max_var_counts=[5000,5000], return_purged=True)
scit.tl.add_metadata(adata)
print(f"  After QC: {adata.n_obs} cells × {adata.n_vars} peaks")

print("  Spectral embedding...")
scit.em.multiview_spectral(adata, ['acet','meth'])
scit.tl.add_metadata(adata)
scit.tl.remove_pc(adata, 'X_multi_spectral', 0)  # PC0 is depth-correlated (r=-0.98)

scit.em.umap(adata, 'X_multi_spectral', n_neighbors=5)
scit.gr.knn(adata, 'X_multi_spectral', n_neighbors=5)
g = scit.gr.neighbor_graph(adata)
scit.gr.leiden(adata, g, 1)
n_cl = adata.obs['leiden'].nunique()
print(f"  Leiden clusters: {n_cl}")


# ── Step 3: binarization ──────────────────────────────────────────────────────
print("\n=== Step 3: Binarization ===")
# mirror notebook: IDF + log + L2, with a threshold to create binary layers
lc_acet = scit.tl.make_layer_config('acet', idf_transform=True, log_transform=True,
                                     normalize_with_l2_norm=True,
                                     feature_active_threshold=0.2)
lc_meth = scit.tl.make_layer_config('meth', idf_transform=True, log_transform=True,
                                     normalize_with_l2_norm=True,
                                     feature_active_threshold=0.2)
print("  Layer configs created")


# ── Step 4: load regulatory links and build reg matrix ────────────────────────
print("\n=== Step 4: Peak → gene regulatory matrix ===")
reg = scit.ld.regulatory_links(peak_links_path, rna)
print(f"  {reg}")

# get gene names that appear in the links
gene_names = rna.var_names.to_numpy()
# peak IDs in adata use bare chrom names (2L:x-y); links file uses chr prefix (chr2L:x-y)
# build a chr-prefixed version for the regulatory matrix lookup, then strip back
raw_peak_ids = adata.var_names.to_numpy()
def _add_chr(pid):
    chrom, coords = pid.split(':', 1)
    return pid if chrom.startswith('chr') else f'chr{pid}'
peak_ids_chr = np.array([_add_chr(p) for p in raw_peak_ids])
R = get_regulatory_matrix(reg, gene_names, peak_ids_chr)
print(f"  Regulatory matrix: {R.shape[0]} genes × {R.shape[1]} peaks")
print(f"  Non-zero links: {R.nnz:,}")
# subset to genes with at least one linked peak
linked_mask = R.sum(axis=1).A1 > 0
R_linked    = R[linked_mask]
names_linked = gene_names[linked_mask]
print(f"  Genes with ≥1 linked peak: {linked_mask.sum():,}")


# ── Step 5: gene-level inference ─────────────────────────────────────────────
print("\n=== Step 5: Gene-level inference ===")
# infer_layer needs the reg_matrix columns to align with adata.var_names order
# R_linked columns correspond to peak_ids_chr; adata uses raw_peak_ids — same order, OK
scit.tl.infer_layer(adata, R_linked, lc_acet, names_linked, 'acet')
scit.tl.infer_layer(adata, R_linked, lc_meth, names_linked, 'meth')
scit.tl.infer_layer(adata, R_linked, lc_acet, names_linked, 'acet_binary')
scit.tl.infer_layer(adata, R_linked, lc_meth, names_linked, 'meth_binary')
gene_names_acet = np.array(adata.uns['gene_names_acet_binary'])
gene_names_meth = np.array(adata.uns['gene_names_meth_binary'])
print(f"  Gene-level layers created: {len(gene_names_acet)} acet genes, {len(gene_names_meth)} meth genes")
print(f"  Layers: {list(adata.layers.keys())}")


# ── Step 6: marker calling per Leiden cluster ─────────────────────────────────
print("\n=== Step 6: Marker calling ===")
# create layer configs for marker calling WITH gene feature names (must be after infer_layer)
# gene-level inferred values are continuous — need a threshold to binarize for marker calling
# median non-zero value of the inferred acet matrix used as threshold
acet_mat = adata.obsm['acet_binary']
acet_nz = acet_mat.data if sps.issparse(acet_mat) else acet_mat[acet_mat > 0]
gene_threshold = float(np.percentile(acet_nz, 50))
print(f"  Gene activity threshold (p50 of non-zero values): {gene_threshold:.4f}")

lc_acet_bin = scit.tl.make_layer_config('acet_binary', in_obsm=True,
                                         feature_names=gene_names_acet,
                                         feature_active_threshold=gene_threshold)
lc_meth_bin = scit.tl.make_layer_config('meth_binary', in_obsm=True,
                                         feature_names=gene_names_meth,
                                         feature_active_threshold=gene_threshold)

# pick genes active in ≥1% of cells (reduces feature space, avoids MWU overflow)
min_cells = max(10, int(adata.n_obs * 0.01))
scit.tl.layer_pick_features(adata, lc_acet_bin, filter_min=min_cells, after_binarization=True)
scit.tl.layer_pick_features(adata, lc_meth_bin, filter_min=min_cells, after_binarization=True)

markers_acet = scit.tl.enriched_in_group(adata, lc_acet_bin, 'leiden', binarize=True)
markers_meth = scit.tl.enriched_in_group(adata, lc_meth_bin, 'leiden', binarize=True)

scit.tl.filter_top_markers(markers_acet, n_top=20, alpha=0.05, keep_only_significant=True)
scit.tl.filter_top_markers(markers_meth, n_top=20, alpha=0.05, keep_only_significant=True)

# print top 5 acet markers per cluster
print("\n  Top H3K27ac gene markers per cluster:")
# top_markers is a list of (cluster_name, DataFrame) tuples
for cl_name, df in (markers_acet.top_markers or []):
    if df is not None and len(df) > 0:
        top5 = df['name'].to_list()[:5]
        print(f"    Cluster {cl_name}: {', '.join(top5)}")

# save markers to TSV
rows = []
for cl_name, df in (markers_acet.top_markers or []):
    if df is not None:
        for rank, row in enumerate(df.iter_rows(named=True)):
            rows.append({'cluster': cl_name, 'rank': rank+1, 'gene': row['name'],
                         'score': row['score'], 'pval': row['pval'], 'mark': 'H3K27ac'})
for cl_name, df in (markers_meth.top_markers or []):
    if df is not None:
        for rank, row in enumerate(df.iter_rows(named=True)):
            rows.append({'cluster': cl_name, 'rank': rank+1, 'gene': row['name'],
                         'score': row['score'], 'pval': row['pval'], 'mark': 'H3K27me3'})
pd.DataFrame(rows).to_csv(f"{PLOT_DIR}/opt2_markers.tsv", sep='\t', index=False)
print(f"\n  Markers saved: {PLOT_DIR}/opt2_markers.tsv")


# ── Step 7: UMAP coloured by marker expression ────────────────────────────────
print("\n=== Step 7: Saving UMAPs ===")

def _umap_leiden(adata, title, path):
    fig, ax = plt.subplots(figsize=(5, 5))
    coords = adata.obsm['X_umap']
    leiden = adata.obs['leiden'].astype('category')
    for cl in leiden.cat.categories:
        mask = leiden == cl
        ax.scatter(coords[mask, 0], coords[mask, 1], s=2, alpha=0.6, label=cl)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel('UMAP 1'); ax.set_ylabel('UMAP 2')
    ax.legend(markerscale=4, bbox_to_anchor=(1,1), loc='upper left', fontsize=6)
    plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close('all')
    print(f"  saved {os.path.basename(path)}")

_umap_leiden(adata, f"opt2: low-q MACS3 peaks\n{adata.n_vars} peaks  |  {n_cl} Leiden clusters",
             f"{PLOT_DIR}/opt2_final_umap.png")

# gene expression UMAPs for a few known Drosophila markers
known_markers = ['grh', 'wg', 'hh', 'en', 'ptc', 'ci', 'ap', 'ey', 'dac']
gene_arr = gene_names_acet
acet_bin_layer = adata.obsm.get('acet_binary')  # infer_layer stores in obsm, not layers
if acet_bin_layer is not None:
    n_plot = min(9, len(known_markers))
    fig, axes = plt.subplots(3, 3, figsize=(12, 10))
    coords = adata.obsm['X_umap']
    for idx, gene in enumerate(known_markers[:n_plot]):
        ax = axes[idx // 3, idx % 3]
        gene_match = np.where(gene_arr == gene)[0]
        if len(gene_match) == 0:
            ax.set_title(f"{gene} (not found)", fontsize=8)
            ax.axis('off')
            continue
        gi = gene_match[0]
        if sps.issparse(acet_bin_layer):
            expr = np.array(acet_bin_layer[:, gi].todense()).flatten()
        else:
            expr = acet_bin_layer[:, gi]
        sc = ax.scatter(coords[:, 0], coords[:, 1], c=expr, s=1, cmap='Reds', alpha=0.7)
        ax.set_title(f"H3K27ac: {gene}", fontsize=8)
        ax.set_xlabel('UMAP 1', fontsize=7); ax.set_ylabel('UMAP 2', fontsize=7)
        plt.colorbar(sc, ax=ax, shrink=0.7)
    plt.suptitle("Known Drosophila marker genes — H3K27ac activity", fontsize=10)
    plt.tight_layout()
    plt.savefig(f"{PLOT_DIR}/opt2_marker_genes_umap.png", dpi=150, bbox_inches='tight')
    plt.close('all')
    print(f"  saved opt2_marker_genes_umap.png")

print("\n=== DOWNSTREAM DONE ===")
