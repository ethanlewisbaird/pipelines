import sys; sys.path.insert(0, '/data/ebaird/scRNAseq/20260522.nanoCT/scit_src')
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np, polars as pl, gzip, io, anndata as ad, scipy.sparse as sps
import src as scit

PEAKS = "/data/ebaird/scRNAseq/20260522.nanoCT/R_analysis_peaks/macs3_bigwig"
BASE  = "/data/ebaird/scRNAseq/20260522.nanoCT/SU.analysis.2026.05.22/Vasso_nanoCT_nanoscope"
OUT   = "/data/ebaird/scRNAseq/20260522.nanoCT/analysis_05.26"
PLOT_DIR = f"{OUT}/qc_plots"

def load_bed(path):
    return pl.read_csv(path, separator='\t', has_header=False, new_columns=['chr','start','end','name'])

def merge_peaks(dfs):
    combined = pl.concat([d.select(['chr','start','end']) for d in dfs]).sort(['chr','start'])
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
    merged = pl.DataFrame(rows, schema={'chr': pl.Utf8, 'start': pl.Int64, 'end': pl.Int64}, orient='row')
    return merged.with_columns(
        (pl.col('chr') + ':' + pl.col('start').cast(pl.Utf8) + '-' + pl.col('end').cast(pl.Utf8)).alias('peak_id')
    )

def _build_peak_index(peaks_df):
    idx = {}
    for chrom in peaks_df['chr'].unique().to_list():
        sub = peaks_df.filter(pl.col('chr') == chrom).sort('start')
        idx[chrom] = (sub['start'].to_numpy(), sub['end'].to_numpy(), sub['peak_id'].to_list())
    return idx

def _pos_to_peak(pos, starts, ends):
    cand = np.searchsorted(starts, pos, side='right') - 1
    valid = (cand >= 0) & (pos < ends[np.clip(cand, 0, len(ends) - 1)])
    return np.where(valid, cand, -1)

def count_fragments_in_peaks(fragments_path, peaks_df, batch_size=400_000):
    peak_idx = _build_peak_index(peaks_df)
    peak_ids = peaks_df['peak_id'].to_list()
    n_peaks  = len(peak_ids)
    read_kwargs = dict(separator='\t', has_header=False,
                       new_columns=['chr', 'start', 'end', 'bc', 'readSupport'])
    with gzip.open(fragments_path, 'rb') as gz:
        buf = io.BytesIO(b''.join(l for l in gz if not l.startswith(b'#')))
    df_full = pl.read_csv(buf, **read_kwargs)
    n_rows  = df_full.height
    bcs     = np.sort(df_full['bc'].unique().to_numpy())
    bc_to_row = {b: i for i, b in enumerate(bcs)}
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

print("Building union peaks...")
union_peaks = merge_peaks([
    load_bed(f"{PEAKS}/H3K27ac_macs3_peaks_dm6_sorted.bed"),
    load_bed(f"{PEAKS}/H3K27me3_macs3_peaks_dm6_sorted.bed"),
])
print(f"Union: {len(union_peaks)} peaks")

print("Counting H3K27ac...")
ac = count_fragments_in_peaks(f"{BASE}/H3K27ac/fragments.tsv.gz", union_peaks)
print("Counting H3K27me3...")
me = count_fragments_in_peaks(f"{BASE}/H3K27me3/fragments.tsv.gz", union_peaks)

print("Stacking...")
adata_peaks = scit.tl.stack_adata([ac, me], ['acet', 'meth'])
scit.tl.add_metadata(adata_peaks)

ac_total = adata_peaks.obs['acet_total_counts'].values.astype(float)
me_total = adata_peaks.obs['meth_total_counts'].values.astype(float)

# --- knee plots ---
fig, axes = plt.subplots(2, 2, figsize=(14, 8))
fig.suptitle('Cell vs background — MACS3 peaks\n(pick min_obs_counts at the knee)', fontsize=11)

for col_i, (counts, mark) in enumerate([(ac_total, 'H3K27ac'), (me_total, 'H3K27me3')]):
    sorted_desc = np.sort(counts)[::-1]
    counts_nz   = counts[counts > 0]

    ax = axes[0, col_i]
    ax.semilogy(np.arange(1, len(sorted_desc) + 1), sorted_desc + 0.1, lw=0.8, color='steelblue')
    ax.set_xlabel('Barcode rank')
    ax.set_ylabel('Total peak counts (log)')
    ax.set_title(f'{mark} — rank plot')
    for rank, color, label in [(1725, 'red', '~1725 (CellRanger ac)'),
                               (1781, 'orange', '~1781 (CellRanger me)')]:
        if len(sorted_desc) >= rank:
            v = sorted_desc[rank - 1]
            ax.axvline(rank,  color=color, lw=1, ls='--', label=f'{label}: count={v:.0f}')
            ax.axhline(v, color=color, lw=0.6, ls=':')
    ax.legend(fontsize=7)

    ax2 = axes[1, col_i]
    ax2.hist(counts_nz, bins=80, color='steelblue', edgecolor='none')
    ax2.set_xlabel(f'Total peak counts (non-zero barcodes, n={len(counts_nz):,})')
    ax2.set_ylabel('Barcodes')
    ax2.set_title(f'{mark} — count distribution (non-zero only)')
    for pct, ls in [(5, '--'), (50, '-'), (95, '--')]:
        v = np.percentile(counts_nz, pct)
        ax2.axvline(v, ls=ls, color='gray', lw=0.9)
        ax2.text(v + 0.5, ax2.get_ylim()[1] * 0.85, f'p{pct}={v:.0f}', fontsize=7, color='gray')

plt.tight_layout()
path = f'{PLOT_DIR}/04_knee_plots.png'
plt.savefig(path, dpi=150, bbox_inches='tight')
plt.close('all')
print(f"saved {path}")

# --- cell-vs-background overlay using bin-called cells ---
print("Loading bin adata for reference barcodes...")
adata_bins = ad.read_h5ad(f"{OUT}/combined_dim_reduced.h5ad")
bin_bcs    = set(adata_bins.obs_names.tolist())
all_bcs    = adata_peaks.obs_names.tolist()
is_cell    = np.array([b in bin_bcs for b in all_bcs])
print(f"  {is_cell.sum()} barcodes overlap with bin-called cells")

fig, axes = plt.subplots(1, 2, figsize=(14, 4))
fig.suptitle('Peak counts: cells from bin analysis (red) vs remaining barcodes (blue)', fontsize=10)
for ax, counts, mark in [(axes[0], ac_total, 'H3K27ac acet'), (axes[1], me_total, 'H3K27me3 meth')]:
    xlim = np.percentile(counts[is_cell], 99) * 1.3 if is_cell.any() else 200
    ax.hist(counts[~is_cell], bins=60, color='steelblue', alpha=0.5, label='background', density=True, range=(0, xlim))
    ax.hist(counts[is_cell],  bins=60, color='red',       alpha=0.7, label='cell (from bins analysis)', density=True, range=(0, xlim))
    ax.set_xlabel('Total peak counts')
    ax.set_ylabel('Density')
    ax.set_title(mark)
    ax.legend()
plt.tight_layout()
path = f'{PLOT_DIR}/05_cells_vs_background.png'
plt.savefig(path, dpi=150, bbox_inches='tight')
plt.close('all')
print(f"saved {path}")

# --- summary stats ---
print("\n=== Reference counts in bin-called cells ===")
for counts, mark in [(ac_total, 'acet'), (me_total, 'meth')]:
    c = counts[is_cell]
    print(f"  {mark}: n={is_cell.sum()}  min={c.min():.0f}  p5={np.percentile(c,5):.0f}  "
          f"median={np.median(c):.0f}  p95={np.percentile(c,95):.0f}  max={c.max():.0f}")
print("\nDONE")
