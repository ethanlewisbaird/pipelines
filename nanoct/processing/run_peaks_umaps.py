import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import anndata as ad
import polars as pl
import scipy.sparse as sps
import gzip, io, os
import sys
sys.path.insert(0, '/data/ebaird/scRNAseq/20260522.nanoCT/scit_src')
import src as scit
import scanpy as sc

OUT = "/data/ebaird/scRNAseq/20260522.nanoCT/analysis_05.26"
BASE = "/data/ebaird/scRNAseq/20260522.nanoCT/SU.analysis.2026.05.22/Vasso_nanoCT_nanoscope"
PLOT_DIR = f"{OUT}/peaks_umaps"
os.makedirs(PLOT_DIR, exist_ok=True)

# ── helpers ──
def load_bed(path, add_chr=False):
    df = pl.read_csv(path, separator='\t', has_header=False,
                     new_columns=['chr', 'start', 'end', 'name'],
                     schema_overrides={'start': pl.Int64, 'end': pl.Int64})
    if add_chr:
        df = df.with_columns(pl.when(~pl.col('chr').str.starts_with('chr')).then(pl.lit('chr') + pl.col('chr')).otherwise(pl.col('chr')).alias('chr'))
    return df

def merge_peaks(dfs):
    combined = pl.concat([d.select(['chr','start','end']) for d in dfs]).sort(['chr','start'])
    rows, cur_chr, cur_start, cur_end = [], None, None, None
    for chrom, start, end in combined.iter_rows():
        if chrom != cur_chr or start > cur_end:
            if cur_chr is not None: rows.append((cur_chr, cur_start, cur_end))
            cur_chr, cur_start, cur_end = chrom, start, end
        else: cur_end = max(cur_end, end)
    if cur_chr is not None: rows.append((cur_chr, cur_start, cur_end))
    merged = pl.DataFrame(rows, schema={'chr': pl.Utf8, 'start': pl.Int64, 'end': pl.Int64}, orient='row')
    return merged.with_columns((pl.col('chr')+':'+pl.col('start').cast(pl.Utf8)+'-'+pl.col('end').cast(pl.Utf8)).alias('peak_id'))

def _build_peak_index(peaks_df):
    peak_ids = peaks_df['peak_id'].to_list()
    id_to_global = {pid: i for i, pid in enumerate(peak_ids)}
    idx = {}
    for chrom in peaks_df['chr'].unique().to_list():
        sub = peaks_df.filter(pl.col('chr') == chrom).sort('start')
        global_idx = np.array([id_to_global[p] for p in sub['peak_id'].to_list()], dtype=np.int32)
        idx[chrom] = (sub['start'].to_numpy(), sub['end'].to_numpy(), global_idx)
    return idx

def _pos_to_peak(pos, starts, ends):
    cand = np.searchsorted(starts, pos, side='right') - 1
    valid = (cand >= 0) & (pos < ends[np.clip(cand, 0, len(ends)-1)])
    return np.where(valid, cand, -1)

def count_fragments_in_peaks(fragments_path, peaks_df, batch_size=400_000):
    peak_idx  = _build_peak_index(peaks_df)
    peak_ids  = peaks_df['peak_id'].to_list()
    n_peaks   = len(peak_ids)
    read_kw   = dict(separator='\t', has_header=False, new_columns=['chr','start','end','bc','readSupport'])
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
            starts_p, ends_p, global_idx = peak_idx[chrom]
            mask    = chr_arr == chrom
            bc_rows = np.array([bc_to_row[b] for b in bc_arr[mask]], dtype=np.int32)
            pc_s    = _pos_to_peak(start_arr[mask], starts_p, ends_p)
            pc_e    = _pos_to_peak(end_arr[mask],   starts_p, ends_p)
            hit_s   = pc_s >= 0
            rows_list.append(bc_rows[hit_s]);  cols_list.append(global_idx[pc_s[hit_s]])
            data_list.append(np.ones(hit_s.sum(), dtype=np.uint32))
            hit_e   = (pc_e >= 0) & (pc_e != pc_s)
            rows_list.append(bc_rows[hit_e]);  cols_list.append(global_idx[pc_e[hit_e]])
            data_list.append(np.ones(hit_e.sum(), dtype=np.uint32))
        print(f"  {100*batch_start/n_rows:.0f}%...", end='\r')
    print()
    rows_arr = np.concatenate(rows_list).astype(np.int32)
    cols_arr = np.concatenate(cols_list).astype(np.int32)
    data_arr = np.concatenate(data_list).astype(np.uint32)
    X = sps.coo_matrix((data_arr, (rows_arr, cols_arr)), shape=(len(bcs), n_peaks)).tocsr()
    a = ad.AnnData(X)
    a.obs.index = bcs
    a.var.index = peak_ids
    return a

def save(path):
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close('all')
    print(f"  saved {os.path.basename(path)}")

# ── load peaks ──
print("Loading MACS3 lowq peaks...")
ac_peaks = load_bed(f"{OUT}/macs3_lowq/H3K27ac_lowq_sorted.bed")
me_peaks = load_bed(f"{OUT}/macs3_lowq/H3K27me3_lowq_sorted.bed")
union_peaks = merge_peaks([ac_peaks, me_peaks])
print(f"  {ac_peaks.height} ac + {me_peaks.height} me → {union_peaks.height} merged")

# ── count fragments ──
print("\nCounting H3K27ac...")
ac = count_fragments_in_peaks(f"{BASE}/H3K27ac/fragments.tsv.gz", union_peaks)
print(f"  {ac.n_obs} x {ac.n_vars}")

print("Counting H3K27me3...")
me = count_fragments_in_peaks(f"{BASE}/H3K27me3/fragments.tsv.gz", union_peaks)
print(f"  {me.n_obs} x {me.n_vars}")

# ── stack ──
print("\nStacking...")
adata = scit.tl.stack_adata([ac, me], ['acet', 'meth'])
print(f"  {adata.n_obs} cells x {adata.n_vars} peaks")

# ── QC histograms ──
scit.tl.add_metadata(adata)
scit.set_defaults(figsize=(8,3))
scit.pl.cell_counts_histogram(adata, xminmax=(0,10), label_exp=True)
save(f"{PLOT_DIR}/01_cell_counts_raw.png")

scit.pl.feature_counts_histogram(adata, xminmax=(0,7), label_exp=True)
save(f"{PLOT_DIR}/02_feature_counts_raw.png")

# ── filter cells ──
print("\nFiltering cells (min_obs=[20,15], max_obs=[3000,3000])...")
adata = scit.tl.filter(adata, ['acet', 'meth'],
                       min_obs_counts=[20, 15],
                       max_obs_counts=[3000, 3000],
                       return_purged=True)
print(f"  {adata.n_obs} cells remaining")

scit.tl.add_metadata(adata)
scit.pl.cell_counts_histogram(adata)
save(f"{PLOT_DIR}/03_cell_counts_filtered.png")

# ── filter features ──
print("\nFiltering features (min_var=[5,5])...")
adata = scit.tl.filter(adata, ['acet', 'meth'],
                       min_var_counts=[5, 5],
                       return_purged=True)
print(f"  {adata.n_vars} features remaining")

scit.tl.add_metadata(adata)
scit.pl.feature_counts_histogram(adata)
save(f"{PLOT_DIR}/04_feature_counts_filtered.png")

# ── spectral embedding ──
print("\nMultiview spectral embedding...")
eigenvalues = scit.em.multiview_spectral(adata, ['acet', 'meth'])
print(f"  eigenvalues: {np.round(eigenvalues[:8], 3)}")

scit.tl.add_metadata(adata)

# ── depth correlation ──
scit.pl.depth_corr(adata, 'X_multi_spectral', 'acet')
save(f"{PLOT_DIR}/05_depth_corr_acet.png")
scit.pl.depth_corr(adata, 'X_multi_spectral', 'meth')
save(f"{PLOT_DIR}/06_depth_corr_meth.png")

# Print correlations
spectral = adata.obsm['X_multi_spectral']
for mark, col in [('acet', 'acet_log_total_counts'), ('meth', 'meth_log_total_counts')]:
    depths = adata.obs[col].values
    corrs = [np.corrcoef(spectral[:, i], depths)[0, 1] for i in range(spectral.shape[1])]
    top = np.argmax(np.abs(corrs))
    print(f"  {mark} depth corr: " + " ".join([f"PC{i}={c:.2f}" for i, c in enumerate(corrs[:8])]))
    print(f"    → highest |corr| at PC{top} ({corrs[top]:.3f})")

# ── remove depth-correlated PC ──
print("\nRemoving PC0 (depth correlation)...")
scit.tl.remove_pc(adata, 'X_multi_spectral', 0)

# ── UMAP + Leiden ──
print("\nClustering...")
scit.gr.knn(adata, 'X_multi_spectral', n_neighbors=25)
g = scit.gr.neighbor_graph(adata)
sc.pp.neighbors(adata, use_rep="X_multi_spectral", n_neighbors=20)
sc.tl.umap(adata, min_dist=0.01, spread=1)
scit.gr.leiden(adata, g, 0.8)

scit.set_defaults(figsize=(6,5))
scit.pl.embedding2d(adata, 'X_umap', 'leiden', title='MACS3 peaks — Leiden clusters')
save(f"{PLOT_DIR}/07_umap_leiden.png")

scit.pl.embedding2d(adata, 'X_umap', 'acet_log_total_counts', title='Log acetylation counts')
save(f"{PLOT_DIR}/08_umap_acet_counts.png")

scit.pl.embedding2d(adata, 'X_umap', 'meth_log_total_counts', title='Log methylation counts')
save(f"{PLOT_DIR}/09_umap_meth_counts.png")

# ── compare with 5kb bins ──
print("\nComparing with 5kb-bin clustering...")
from sklearn.metrics import adjusted_rand_score
adata_bins = ad.read_h5ad(f"{OUT}/combined_dim_reduced.h5ad")
shared_bcs = np.intersect1d(adata_bins.obs_names, adata.obs_names)
bins_sub = adata_bins[shared_bcs]
peaks_sub = adata[shared_bcs]
ari = adjusted_rand_score(bins_sub.obs['leiden'].values, peaks_sub.obs['leiden'].values)
print(f"  Shared cells: {len(shared_bcs)}")
print(f"  Bins clusters: {bins_sub.obs['leiden'].nunique()}, Peaks clusters: {peaks_sub.obs['leiden'].nunique()}")
print(f"  ARI: {ari:.3f}")

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, ad_obj, title in [(axes[0], bins_sub, f"5kb bins ({bins_sub.obs['leiden'].nunique()} cl)"),
                           (axes[1], peaks_sub, f"MACS3 peaks ({peaks_sub.obs['leiden'].nunique()} cl, ARI={ari:.2f})")]:
    umap = ad_obj.obsm['X_umap']
    ax.scatter(umap[:, 0], umap[:, 1], c=ad_obj.obs['leiden'].astype(int), s=3, cmap='tab20', rasterized=True)
    ax.set_title(title, fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])
plt.tight_layout()
save(f"{PLOT_DIR}/10_comparison_bins_vs_peaks.png")

# ── save ──
adata.write_h5ad(f'{OUT}/combined_MACS3peaks_dim_reduced.h5ad')
print(f"\nAll plots saved to {PLOT_DIR}/")
print("DONE")
