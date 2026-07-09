"""
Run nanoCT MACS3 peaks pipeline up to the first QC inspection point.
Saves all plots to analysis_05.26/qc_plots/.
"""
import sys
sys.path.insert(0, '/data/ebaird/scRNAseq/20260522.nanoCT/scit_src')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import anndata as ad
import scipy.sparse as sps
import polars as pl
import gzip, io, os

import src as scit

# ── paths ─────────────────────────────────────────────────────────────────────
BASE  = "/data/ebaird/scRNAseq/20260522.nanoCT/SU.analysis.2026.05.22/Vasso_nanoCT_nanoscope"
PEAKS = "/data/ebaird/scRNAseq/20260522.nanoCT/R_analysis_peaks/macs3_bigwig"
OUT   = "/data/ebaird/scRNAseq/20260522.nanoCT/analysis_05.26"
PLOT_DIR = f"{OUT}/qc_plots"
os.makedirs(PLOT_DIR, exist_ok=True)

def save(name):
    path = f"{PLOT_DIR}/{name}.png"
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close('all')
    print(f"  saved → {path}")

# ── 1. build union peak set ────────────────────────────────────────────────────
print("\n=== Building union MACS3 peak set ===")

def load_bed(path):
    return pl.read_csv(
        path, separator='\t', has_header=False,
        new_columns=['chr', 'start', 'end', 'name']
    )

def merge_peaks(dfs):
    combined = pl.concat([d.select(['chr', 'start', 'end']) for d in dfs]).sort(['chr', 'start'])
    rows = []
    cur_chr, cur_start, cur_end = None, None, None
    for chrom, start, end in combined.iter_rows():
        if chrom != cur_chr or start > cur_end:
            if cur_chr is not None:
                rows.append((cur_chr, cur_start, cur_end))
            cur_chr, cur_start, cur_end = chrom, start, end
        else:
            cur_end = max(cur_end, end)
    if cur_chr is not None:
        rows.append((cur_chr, cur_start, cur_end))
    merged = pl.DataFrame(rows, schema={'chr': pl.Utf8, 'start': pl.Int64, 'end': pl.Int64})
    return merged.with_columns(
        (pl.col('chr') + ':' + pl.col('start').cast(pl.Utf8) + '-' + pl.col('end').cast(pl.Utf8)).alias('peak_id')
    )

ac_peaks = load_bed(f"{PEAKS}/H3K27ac_macs3_peaks_dm6_sorted.bed")
me_peaks = load_bed(f"{PEAKS}/H3K27me3_macs3_peaks_dm6_sorted.bed")
union_peaks = merge_peaks([ac_peaks, me_peaks])
print(f"H3K27ac: {len(ac_peaks)} peaks")
print(f"H3K27me3: {len(me_peaks)} peaks")
print(f"Union after merging: {len(union_peaks)} peaks")

# ── 2. fragment counting ───────────────────────────────────────────────────────
def _build_peak_index(peaks_df):
    idx = {}
    for chrom in peaks_df['chr'].unique().to_list():
        sub = peaks_df.filter(pl.col('chr') == chrom).sort('start')
        idx[chrom] = (
            sub['start'].to_numpy(),
            sub['end'].to_numpy(),
            sub['peak_id'].to_list(),
        )
    return idx

def _pos_to_peak(pos, starts, ends):
    cand = np.searchsorted(starts, pos, side='right') - 1
    valid = (cand >= 0) & (pos < ends[np.clip(cand, 0, len(ends) - 1)])
    return np.where(valid, cand, -1)

def count_fragments_in_peaks(fragments_path, peaks_df, batch_size=400_000):
    peak_idx  = _build_peak_index(peaks_df)
    peak_ids  = peaks_df['peak_id'].to_list()
    n_peaks   = len(peak_ids)

    read_kwargs = dict(
        separator='\t', has_header=False,
        new_columns=['chr', 'start', 'end', 'bc', 'readSupport']
    )
    print(f"  reading {fragments_path} ...")
    with gzip.open(fragments_path, 'rb') as gz:
        buf = io.BytesIO(b''.join(l for l in gz if not l.startswith(b'#')))
    df_full = pl.read_csv(buf, **read_kwargs)
    n_rows  = df_full.height

    bcs      = np.sort(df_full['bc'].unique().to_numpy())
    bc_to_row = {b: i for i, b in enumerate(bcs)}
    n_cells  = len(bcs)
    print(f"  {n_rows:,} fragments  |  {n_cells:,} barcodes  |  {n_peaks} peaks")

    rows_list, cols_list, data_list = [], [], []

    for batch_start in range(0, n_rows, batch_size):
        batch     = df_full.slice(batch_start, batch_size)
        bc_arr    = batch['bc'].to_numpy()
        start_arr = batch['start'].to_numpy()
        end_arr   = batch['end'].to_numpy()
        chr_arr   = batch['chr'].to_numpy()

        for chrom in np.unique(chr_arr):
            if chrom not in peak_idx:
                continue
            starts_p, ends_p, _ = peak_idx[chrom]
            mask   = chr_arr == chrom
            bc_sub = bc_arr[mask]
            bc_rows = np.array([bc_to_row[b] for b in bc_sub], dtype=np.int32)

            pc_start = _pos_to_peak(start_arr[mask], starts_p, ends_p)
            pc_end   = _pos_to_peak(end_arr[mask],   starts_p, ends_p)

            hit_s = pc_start >= 0
            rows_list.append(bc_rows[hit_s]);  cols_list.append(pc_start[hit_s])
            data_list.append(np.ones(hit_s.sum(), dtype=np.uint32))

            hit_e = (pc_end >= 0) & (pc_end != pc_start)
            rows_list.append(bc_rows[hit_e]);  cols_list.append(pc_end[hit_e])
            data_list.append(np.ones(hit_e.sum(), dtype=np.uint32))

        if batch_start % (batch_size * 5) == 0:
            pct = 100 * batch_start / n_rows
            print(f"  {pct:.0f}% ...", end='\r')

    print()
    rows_arr = np.concatenate(rows_list).astype(np.int32)
    cols_arr = np.concatenate(cols_list).astype(np.int32)
    data_arr = np.concatenate(data_list).astype(np.uint32)

    X = sps.coo_matrix((data_arr, (rows_arr, cols_arr)), shape=(n_cells, n_peaks)).tocsr()
    adata = ad.AnnData(X)
    adata.obs.index = bcs
    adata.var.index = peak_ids
    return adata

print("\n=== Counting H3K27ac ===")
ac = count_fragments_in_peaks(f"{BASE}/H3K27ac/fragments.tsv.gz", union_peaks)

print("\n=== Counting H3K27me3 ===")
me = count_fragments_in_peaks(f"{BASE}/H3K27me3/fragments.tsv.gz", union_peaks)

# ── 3. stack ──────────────────────────────────────────────────────────────────
print("\n=== Stacking into multimodal AnnData ===")
adata_peaks = scit.tl.stack_adata([ac, me], ['acet', 'meth'])
print(adata_peaks)
print(f"Cells after barcode intersection: {adata_peaks.n_obs}")

# ── 4. QC — first inspection point ────────────────────────────────────────────
print("\n=== QC: feature counts histogram ===")
scit.tl.add_metadata(adata_peaks)
scit.set_defaults(figsize=(10, 2))

scit.pl.feature_counts_histogram(adata_peaks, xminmax=(0, 7), label_exp=True)
save("01_feature_counts_histogram")

# Also show raw per-mark total count distributions so thresholds are clear
fig, axes = plt.subplots(1, 2, figsize=(12, 3))
for ax, col, mark in [
    (axes[0], 'acet_log_total_counts', 'H3K27ac'),
    (axes[1], 'meth_log_total_counts', 'H3K27me3'),
]:
    vals = adata_peaks.obs[col].values
    ax.hist(vals, bins=60, color='steelblue', edgecolor='none')
    ax.set_xlabel(f'log total counts ({mark})')
    ax.set_ylabel('cells')
    ax.set_title(f'{mark} per-cell count distribution\n'
                 f'median={np.median(vals):.2f}  p5={np.percentile(vals,5):.2f}  p98={np.percentile(vals,98):.2f}')
fig.suptitle('Raw per-cell counts (MACS3 peaks) — use these to set min/max thresholds', y=1.02)
plt.tight_layout()
save("02_per_cell_count_distributions")

# Per-peak count distributions (to guide min_var_counts)
fig, axes = plt.subplots(1, 2, figsize=(12, 3))
for ax, layer, mark in [
    (axes[0], 'acet', 'H3K27ac'),
    (axes[1], 'meth', 'H3K27me3'),
]:
    peak_counts = np.array(adata_peaks.layers[layer].sum(axis=0)).flatten()
    ax.hist(peak_counts, bins=60, color='coral', edgecolor='none')
    ax.set_xlabel(f'total counts per peak ({mark})')
    ax.set_ylabel('peaks')
    n_zero = (peak_counts == 0).sum()
    ax.set_title(f'{mark}: {adata_peaks.n_vars} peaks\n'
                 f'{n_zero} with 0 counts, median={np.median(peak_counts):.1f}')
fig.suptitle('Per-peak count distributions — use these to set min_var_counts', y=1.02)
plt.tight_layout()
save("03_per_peak_count_distributions")

print("\n=== Summary ===")
print(f"adata_peaks: {adata_peaks}")
print(f"\nPer-cell stats (log counts):")
for mark in ['acet', 'meth']:
    col = f'{mark}_log_total_counts'
    v = adata_peaks.obs[col].values
    print(f"  {mark}: min={v.min():.2f}  p5={np.percentile(v,5):.2f}  "
          f"median={np.median(v):.2f}  p98={np.percentile(v,98):.2f}  max={v.max():.2f}")

print(f"\nPer-peak stats (raw counts):")
for layer in ['acet', 'meth']:
    pc = np.array(adata_peaks.layers[layer].sum(axis=0)).flatten()
    print(f"  {layer}: n_peaks={len(pc)}  zeros={(pc==0).sum()}  "
          f"median={np.median(pc):.1f}  p5={np.percentile(pc,5):.1f}  p95={np.percentile(pc,95):.1f}")

print(f"\nPlots saved to: {PLOT_DIR}")
print("DONE — ready for threshold inspection.")
