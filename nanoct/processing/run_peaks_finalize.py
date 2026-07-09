"""
Finalize peaks pipeline: remove depth-correlated PC0 (confirmed from depth_corr plots),
then UMAP + Leiden, then compare all three peak sets vs 5 kb bins.
"""
import sys
sys.path.insert(0, '/data/ebaird/scRNAseq/20260522.nanoCT/scit_src')

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import gzip, io, os
import anndata as ad
import scipy.sparse as sps
import src as scit
from sklearn.metrics import adjusted_rand_score

BASE     = "/data/ebaird/scRNAseq/20260522.nanoCT/SU.analysis.2026.05.22/Vasso_nanoCT_nanoscope"
LOWQ_DIR = "/data/ebaird/scRNAseq/20260522.nanoCT/analysis_05.26/macs3_lowq"
ORIG_DIR = "/data/ebaird/scRNAseq/20260522.nanoCT/R_analysis_peaks/macs3_bigwig"
OUT      = "/data/ebaird/scRNAseq/20260522.nanoCT/analysis_05.26"
PLOT_DIR = f"{OUT}/pipeline_plots"
BINS_H5  = f"{OUT}/combined_dim_reduced.h5ad"

CHROM_SIZES = {
    '2L': 23513712, '2R': 25286936, '3L': 28110227, '3R': 32079331,
    '4': 1348131, 'X': 23542271, 'Y': 3667352,
    'mitochondrion_genome': 19524, 'rDNA': 76973,
}

# ── helpers (same as pipeline script) ────────────────────────────────────────
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

def extend_and_merge(dfs, pad_bp):
    combined = pl.concat([d.select(['chr','start','end']) for d in dfs])
    combined = combined.with_columns([
        (pl.col('start') - pad_bp).clip(0).alias('start'),
        pl.struct(['chr','end']).map_elements(
            lambda r: min(r['end'] + pad_bp, CHROM_SIZES.get(r['chr'], r['end'] + pad_bp)),
            return_dtype=pl.Int64
        ).alias('end'),
    ]).sort(['chr','start'])
    return merge_peaks([combined.with_columns(pl.lit('x').alias('name'))])

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


def run_full_pipeline(label, union_peaks,
                      min_obs_counts, max_obs_counts,
                      min_var_counts, max_var_counts,
                      remove_pc_index=0):
    """Full pipeline through UMAP + Leiden, returns adata."""
    prefix = f"{PLOT_DIR}/{label}"
    print(f"\n{'='*65}")
    print(f"  {label}  |  {len(union_peaks)} peaks  |  removing PC{remove_pc_index}")
    print(f"{'='*65}")

    print("  Counting H3K27ac...")
    ac = count_fragments_in_peaks(f"{BASE}/H3K27ac/fragments.tsv.gz", union_peaks)
    print("  Counting H3K27me3...")
    me = count_fragments_in_peaks(f"{BASE}/H3K27me3/fragments.tsv.gz", union_peaks)

    adata = scit.tl.stack_adata([ac, me], ['acet', 'meth'])
    print(f"  Stacked: {adata.n_obs} cells × {adata.n_vars} peaks")

    scit.tl.add_metadata(adata)
    adata = scit.tl.filter(adata, ['acet', 'meth'],
                           min_obs_counts=min_obs_counts,
                           max_obs_counts=max_obs_counts,
                           return_purged=True)
    scit.tl.add_metadata(adata)
    adata = scit.tl.filter(adata, ['acet', 'meth'],
                           min_var_counts=min_var_counts,
                           return_purged=True)
    scit.tl.filter(adata, ['acet', 'meth'],
                   max_var_counts=max_var_counts,
                   return_purged=True)
    scit.tl.add_metadata(adata)
    print(f"  After QC: {adata.n_obs} cells × {adata.n_vars} peaks")

    print("  Running multiview spectral...")
    scit.em.multiview_spectral(adata, ['acet', 'meth'])
    scit.tl.add_metadata(adata)

    # cell 23: remove depth-correlated PC
    print(f"  Removing PC{remove_pc_index}...")
    scit.tl.remove_pc(adata, 'X_multi_spectral', remove_pc_index)

    # cell 24
    print("  UMAP + Leiden...")
    scit.em.umap(adata, 'X_multi_spectral', n_neighbors=5)
    scit.gr.knn(adata, 'X_multi_spectral', n_neighbors=5)
    g = scit.gr.neighbor_graph(adata)
    scit.gr.leiden(adata, g, 1)

    n_cl = adata.obs['leiden'].nunique()
    print(f"  Leiden clusters: {n_cl}")

    # save UMAP
    fig, ax = plt.subplots(figsize=(5, 5))
    coords = adata.obsm['X_umap']
    leiden = adata.obs['leiden'].astype('category')
    for cl in leiden.cat.categories:
        mask = leiden == cl
        ax.scatter(coords[mask, 0], coords[mask, 1], s=2, alpha=0.6, label=cl)
    ax.set_title(f"{label}\n{len(union_peaks)} peaks  |  {n_cl} Leiden clusters")
    ax.set_xlabel('UMAP 1'); ax.set_ylabel('UMAP 2')
    ax.legend(markerscale=4, bbox_to_anchor=(1, 1), loc='upper left', fontsize=6)
    plt.tight_layout()
    plt.savefig(f"{prefix}_07_umap.png", dpi=150, bbox_inches='tight')
    plt.close('all')
    print(f"    saved {label}_07_umap.png")

    return adata


# ── define peak sets (same thresholds as pipeline script) ─────────────────────
ac_orig = load_bed(f"{ORIG_DIR}/H3K27ac_macs3_peaks_dm6_sorted.bed")
me_orig = load_bed(f"{ORIG_DIR}/H3K27me3_macs3_peaks_dm6_sorted.bed")
ac_lowq = load_bed(f"{LOWQ_DIR}/H3K27ac_lowq_sorted.bed")
me_lowq = load_bed(f"{LOWQ_DIR}/H3K27me3_lowq_sorted.bed")

configs = {
    "opt1_orig_macs3": (
        merge_peaks([ac_orig, me_orig]),
        [1, 1], [30, 30], [1, 1], [200, 200], 0,
    ),
    "opt2_lowq_macs3": (
        merge_peaks([ac_lowq, me_lowq]),
        [20, 15], [3000, 3000], [5, 5], [5000, 5000], 0,
    ),
    "opt3_pm500bp": (
        extend_and_merge([ac_orig, me_orig], 500),
        [2, 2], [100, 100], [1, 1], [500, 500], 0,
    ),
    "opt3_pm1000bp": (
        extend_and_merge([ac_orig, me_orig], 1000),
        [2, 2], [150, 150], [1, 1], [800, 800], 0,
    ),
}

results = {}
for label, (peaks, min_obs, max_obs, min_var, max_var, pc_idx) in configs.items():
    results[label] = run_full_pipeline(label, peaks, min_obs, max_obs, min_var, max_var, pc_idx)

# ── comparison figure: bins + all three peak sets ─────────────────────────────
print("\n\nBuilding comparison figure...")
adata_bins = ad.read_h5ad(BINS_H5)

# find cells in all four objects
shared = set(adata_bins.obs_names.tolist())
for ad_ in results.values():
    shared &= set(ad_.obs_names.tolist())
shared = sorted(shared)
print(f"Cells shared across all objects: {len(shared)}")

n_panels = 1 + len(results)
fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, 5))

def _scatter(ax, adata_, shared_bcs, title):
    sub    = adata_[shared_bcs]
    coords = sub.obsm['X_umap']
    leiden = sub.obs['leiden'].astype('category')
    for cl in leiden.cat.categories:
        mask = leiden == cl
        ax.scatter(coords[mask, 0], coords[mask, 1], s=2, alpha=0.6, label=cl)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel('UMAP 1'); ax.set_ylabel('UMAP 2')
    ax.legend(markerscale=3, fontsize=5, bbox_to_anchor=(1, 1), loc='upper left')

bins_sub = adata_bins[shared]
n_bins_cl = bins_sub.obs['leiden'].nunique()
_scatter(axes[0], adata_bins, shared, f"5 kb bins\n({n_bins_cl} clusters)")
print(f"\n  {'Feature set':<22}  {'Features':>9}  {'Cells':>7}  {'Clusters':>9}  {'ARI vs bins':>12}")
print(f"  {'-'*22}  {'-'*9}  {'-'*7}  {'-'*9}  {'-'*12}")
print(f"  {'5 kb bins':<22}  {adata_bins.n_vars:>9,}  {adata_bins.n_obs:>7,}  {n_bins_cl:>9}  {'(reference)':>12}")

for i, (label, ad_) in enumerate(results.items()):
    sub = ad_[shared]
    ari = adjusted_rand_score(bins_sub.obs['leiden'].values, sub.obs['leiden'].values)
    n_cl = sub.obs['leiden'].nunique()
    _scatter(axes[i+1], ad_, shared, f"{label}\n({n_cl} cl, ARI={ari:.3f})")
    print(f"  {label:<22}  {ad_.n_vars:>9,}  {ad_.n_obs:>7,}  {n_cl:>9}  {ari:>12.3f}")

plt.suptitle("Peaks vs bins — all options (proper pipeline, PC0 removed)", y=1.02)
plt.tight_layout()
out_path = f"{PLOT_DIR}/comparison_all_options.png"
plt.savefig(out_path, dpi=150, bbox_inches='tight')
plt.close('all')
print(f"\n  Comparison saved: {out_path}")
print("\n=== FINALIZE DONE ===")
