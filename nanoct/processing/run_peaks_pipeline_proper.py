"""
Peaks-vs-bins pipeline following the original nanoCT_workshop.ipynb exactly:
  add_metadata → cell filter → feature filter → soft feature filter
  → multiview_spectral → depth_corr plots → STOP (user picks PC to remove)

Runs all three peak sets: low-q MACS3, ±500 bp extended, ±1000 bp extended.
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

BASE     = "/data/ebaird/scRNAseq/20260522.nanoCT/SU.analysis.2026.05.22/Vasso_nanoCT_nanoscope"
LOWQ_DIR = "/data/ebaird/scRNAseq/20260522.nanoCT/analysis_05.26/macs3_lowq"
ORIG_DIR = "/data/ebaird/scRNAseq/20260522.nanoCT/R_analysis_peaks/macs3_bigwig"
OUT      = "/data/ebaird/scRNAseq/20260522.nanoCT/analysis_05.26"
PLOT_DIR = f"{OUT}/pipeline_plots"
os.makedirs(PLOT_DIR, exist_ok=True)

CHROM_SIZES = {
    '2L': 23513712, '2R': 25286936, '3L': 28110227, '3R': 32079331,
    '4': 1348131, 'X': 23542271, 'Y': 3667352,
    'mitochondrion_genome': 19524, 'rDNA': 76973,
}

# ── helpers ───────────────────────────────────────────────────────────────────
def save(path):
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close('all')
    print(f"    saved {os.path.basename(path)}")

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
        else:
            cur_end = max(cur_end, end)
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
    return merge_peaks([combined.rename({'chr':'chr','start':'start','end':'end'})
                        .with_columns(pl.lit('x').alias('name'))])

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
    peak_idx  = _build_peak_index(peaks_df)
    peak_ids  = peaks_df['peak_id'].to_list()
    n_peaks   = len(peak_ids)
    read_kw   = dict(separator='\t', has_header=False,
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


# ── pipeline function — mirrors notebook exactly ──────────────────────────────
def run_pipeline(label, union_peaks,
                 min_obs_counts, max_obs_counts,
                 min_var_counts, max_var_counts):
    """
    Follows nanoCT_workshop.ipynb cells 13-22 exactly.
    Stops before remove_pc and saves depth_corr plots.
    """
    prefix = f"{PLOT_DIR}/{label}"
    print(f"\n{'='*65}")
    print(f"  {label}  |  {len(union_peaks)} peaks")
    print(f"{'='*65}")

    # count fragments
    print("  Counting H3K27ac...")
    ac = count_fragments_in_peaks(f"{BASE}/H3K27ac/fragments.tsv.gz", union_peaks)
    print("  Counting H3K27me3...")
    me = count_fragments_in_peaks(f"{BASE}/H3K27me3/fragments.tsv.gz", union_peaks)

    # stack
    adata = scit.tl.stack_adata([ac, me], ['acet', 'meth'])
    print(f"  Stacked: {adata.n_obs} cells × {adata.n_vars} peaks")

    # --- cell 13: add_metadata + cell counts histogram ---
    scit.tl.add_metadata(adata)
    scit.set_defaults(figsize=(10, 2))
    scit.pl.cell_counts_histogram(adata, xminmax=(0, 10), label_exp=True)
    save(f"{prefix}_01_cell_counts_raw.png")

    # --- cell 14: feature counts histogram ---
    scit.pl.feature_counts_histogram(adata, xminmax=(0, 7), label_exp=True)
    save(f"{prefix}_02_feature_counts_raw.png")

    # --- cell 15: hard obs filter ---
    print(f"  Filtering cells: min_obs={min_obs_counts}, max_obs={max_obs_counts}")
    adata = scit.tl.filter(adata, ['acet', 'meth'],
                           min_obs_counts=min_obs_counts,
                           max_obs_counts=max_obs_counts,
                           return_purged=True)
    print(f"  After cell filter: {adata.n_obs} cells")

    # --- cell 16: post-filter cell histogram ---
    scit.tl.add_metadata(adata)
    scit.pl.cell_counts_histogram(adata)
    save(f"{prefix}_03_cell_counts_filtered.png")

    # --- cell 17: hard feature filter ---
    print(f"  Filtering features: min_var={min_var_counts}")
    adata = scit.tl.filter(adata, ['acet', 'meth'],
                           min_var_counts=min_var_counts,
                           return_purged=True)
    print(f"  After feature filter: {adata.n_vars} peaks")

    # --- cell 18: soft feature filter (sets exclude column for spectral) ---
    print(f"  Soft feature filter: max_var={max_var_counts}")
    scit.tl.filter(adata, ['acet', 'meth'],
                   max_var_counts=max_var_counts,
                   return_purged=True)
    n_excluded = adata.var['exclude'].sum() if 'exclude' in adata.var.columns else 0
    print(f"  Features soft-excluded from spectral: {n_excluded}")

    # --- cell 19: post-filter feature histogram ---
    scit.tl.add_metadata(adata)
    scit.pl.feature_counts_histogram(adata)
    save(f"{prefix}_04_feature_counts_filtered.png")

    # --- cell 20: multiview spectral ---
    print("  Running multiview spectral...")
    eigenvalues = scit.em.multiview_spectral(adata, ['acet', 'meth'])
    print(f"  Eigenvalues: {np.round(eigenvalues[:8], 3)}")

    # --- cell 21: add_metadata ---
    scit.tl.add_metadata(adata)

    # --- cell 22: depth_corr plots — key diagnostic ---
    print("  Saving depth_corr plots...")
    scit.pl.depth_corr(adata, 'X_multi_spectral', 'acet')
    save(f"{prefix}_05_depth_corr_acet.png")

    scit.pl.depth_corr(adata, 'X_multi_spectral', 'meth')
    save(f"{prefix}_06_depth_corr_meth.png")

    # also print correlation values for each PC
    spectral = adata.obsm['X_multi_spectral']
    for mark, col in [('acet', 'acet_log_total_counts'), ('meth', 'meth_log_total_counts')]:
        depths = adata.obs[col].values
        corrs  = [np.corrcoef(spectral[:, i], depths)[0, 1] for i in range(spectral.shape[1])]
        top    = np.argmax(np.abs(corrs))
        print(f"  Depth corr ({mark}): " +
              "  ".join([f"PC{i}={c:.2f}" for i, c in enumerate(corrs[:8])]))
        print(f"    → highest |corr| at PC{top} ({corrs[top]:.3f})")

    print(f"\n  *** Inspect {label}_05/06 depth_corr plots, then set remove_pc index ***")
    return adata


# ── define peak sets ──────────────────────────────────────────────────────────
ac_orig = load_bed(f"{ORIG_DIR}/H3K27ac_macs3_peaks_dm6_sorted.bed")
me_orig = load_bed(f"{ORIG_DIR}/H3K27me3_macs3_peaks_dm6_sorted.bed")
ac_lowq = load_bed(f"{LOWQ_DIR}/H3K27ac_lowq_sorted.bed")
me_lowq = load_bed(f"{LOWQ_DIR}/H3K27me3_lowq_sorted.bed")

peak_sets = {
    # label: (union_peaks, min_obs, max_obs, min_var, max_var)
    # Thresholds chosen to match bins notebook spirit:
    #   - obs: raw counts (not log); inspect histograms to confirm
    #   - var: scaled relative to feature count range
    "opt2_lowq_macs3": (
        merge_peaks([ac_lowq, me_lowq]),
        [20, 15],           # min_obs_counts — same as bins notebook
        [3000, 3000],       # max_obs_counts — same as bins notebook
        [5, 5],             # min_var_counts — lower than bins (fewer peaks, each more focal)
        [5000, 5000],       # max_var_counts — soft exclude very high-coverage peaks
    ),
    "opt3_pm500bp": (
        extend_and_merge([ac_orig, me_orig], 500),
        [2, 2],             # min_obs_counts — counts are low with narrow peaks
        [100, 100],         # max_obs_counts
        [1, 1],             # min_var_counts
        [500, 500],         # max_var_counts
    ),
    "opt3_pm1000bp": (
        extend_and_merge([ac_orig, me_orig], 1000),
        [2, 2],
        [150, 150],
        [1, 1],
        [800, 800],
    ),
}

results = {}
for label, (peaks, min_obs, max_obs, min_var, max_var) in peak_sets.items():
    results[label] = run_pipeline(label, peaks, min_obs, max_obs, min_var, max_var)

print("\n\n=== ALL DONE ===")
print(f"Plots saved to: {PLOT_DIR}/")
print("\nFor each option, inspect the _05/_06 depth_corr plots to pick the PC index,")
print("then run run_peaks_finalize.py with those indices.")
