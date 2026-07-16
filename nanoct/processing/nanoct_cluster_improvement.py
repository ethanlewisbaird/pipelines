#!/usr/bin/env python3
"""
nanoCT Clustering Improvement Pipeline

Improves clustering via better peak/feature engineering and
dimensionality reduction with proper signal separation.

Improvements over current approach:
  [PEAKS]
  1. ATAC universe + DiffBind consensus peaks (union, weighted)
  2. Dispersion-based variable peak selection (not just min-count)
  3. Multi-resolution: both narrow (2kb) and broad (10kb) windows

  [DIMENSIONALITY REDUCTION]
  4. Per-mark TF-IDF → LSI (standard scATAC-seq practice)
  5. CCA-based joint embedding (shared + mark-specific factors)
  6. Systematic PC quality assessment (remove ALL technical PCs)
  7. Per-mark library-size normalization before combining

Output:
  - Full data: combined_improved.h5ad  (all cells, dispersion-filtered peaks)
  - Comparison UMAPs: Before vs after, per-mark signals, cluster stats
  - QC report: peak selection diagnostics, PC quality, cluster metrics

Usage:
    python nanoct_cluster_improvement.py

Environment variables:
    NANOCT_DATA_DIR  - data directory override
    SCIT_PATH        - scit src path
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import anndata as ad
import polars as pl
import scipy.sparse as sps
from scipy.sparse import diags, vstack
from scipy.stats import pearsonr, poisson, nbinom
import gzip, io, os, sys, warnings, time
warnings.filterwarnings('ignore')

# ── Paths ──────────────────────────────────────────────────────────────────
DATA_DIR = os.environ.get("NANOCT_DATA_DIR",
                          "/data/ebaird/scentinel/nanoCT/20260522.nanoCT")
OUT      = os.path.join(DATA_DIR, "analysis_05.26")
BASE     = os.path.join(DATA_DIR, "SU.analysis.2026.05.22", "Vasso_nanoCT_nanoscope")
DIFFBIND_FILE = "/data/vtheodorou/fastq/EnsemblBDGP6.46_Diff.paper/ATAC.CUT.TAG/DiffBind/DB.T0.Nw.K27.vs.Nw.K4.CUT.TAG.txt"
ATAC_UNIVERSE_FILE = "/data/vtheodorou/2025.05.ATAC.tumor/ATAC.universe.4.homer.txt"

PLOT_DIR = os.path.join(OUT, "cluster_improvement_plots")
os.makedirs(PLOT_DIR, exist_ok=True)

CHROM_SIZES = {
    '2L': 23513712, '2R': 25286936, '3L': 28110227, '3R': 32079331,
    '4': 1348131, 'X': 23542271, 'Y': 3667352, 'mitochondrion_genome': 19524
}
MAJOR_CHROMS = set('chr' + c for c in CHROM_SIZES if c != 'mitochondrion_genome')

# ============================================================================
#  PART 1 — PEAK & FEATURE ENGINEERING
# ============================================================================

def load_diffbind_peaks():
    """Load DiffBind CUT&Tag differential peaks."""
    print("\n[1a] Loading DiffBind peaks...")
    raw = pl.read_csv(
        DIFFBIND_FILE, separator='\t', has_header=True, quote_char='"',
        new_columns=['seqnames','start','end','width','strand',
                     'Conc','Conc_K27me3','Conc_K4me3','Fold','pvalue','FDR'],
        schema_overrides={'start': pl.Int64, 'end': pl.Int64}
    )
    raw = raw.with_columns(pl.col('seqnames').str.replace_all('"', '').alias('seqnames'))
    df = raw.select([
        pl.col('seqnames').alias('chr'),
        pl.col('start').cast(pl.Int64),
        pl.col('end').cast(pl.Int64),
        pl.col('Fold').abs().alias('score'),  # fold change magnitude as weight
        pl.col('FDR'),
    ]).with_columns((pl.lit('chr') + pl.col('chr')).alias('chr'))
    # Keep significant ones
    n_all = df.height
    df = df.filter((pl.col('chr').is_in(MAJOR_CHROMS)) & (pl.col('FDR') < 0.05))
    print(f"  DiffBind significant peaks (FDR<0.05): {df.height}/{n_all}")
    return df

def load_atac_universe_peaks():
    """Load ATAC-seq universe peaks."""
    print("\n[1b] Loading ATAC universe peaks...")
    # ATAC universe format: chr start end ...
    raw = pl.read_csv(
        ATAC_UNIVERSE_FILE, separator='\t', has_header=False,
        new_columns=['chr','start','end','name','score','strand','thickStart','thickEnd','rgb'],
        schema_overrides={'start': pl.Int64, 'end': pl.Int64}
    )
    df = raw.select(['chr','start','end']).with_columns(
        # Ensure chr prefix
        pl.when(pl.col('chr').str.starts_with('chr'))
          .then(pl.col('chr'))
          .otherwise(pl.lit('chr') + pl.col('chr')).alias('chr')
    ).filter(pl.col('chr').is_in(MAJOR_CHROMS))
    print(f"  ATAC universe peaks on major chroms: {df.height}")
    return df

def extend_peaks(df, half_w=3000):
    """Extend peaks ±half_w from center. Returns new df with peak_id."""
    chr_sizes = {('chr' + k): v for k, v in CHROM_SIZES.items()}
    return df.with_columns([
        ((pl.col('start') + pl.col('end')) // 2).alias('center'),
    ]).with_columns([
        (pl.col('center') - half_w).clip(0).alias('start'),
        (pl.col('center') + half_w).alias('end'),
    ]).with_columns([
        pl.struct(['chr','end']).map_elements(
            lambda r: min(r['end'], chr_sizes.get(r['chr'], r['end'])),
            return_dtype=pl.Int64
        ).alias('end'),
    ]).with_columns(
        (pl.col('chr') + ':' + pl.col('start').cast(pl.Utf8) + '-' + pl.col('end').cast(pl.Utf8)).alias('peak_id')
    ).drop('center')

def merge_overlapping(df):
    """Merge overlapping windows."""
    df = df.sort(['chr', 'start'])
    merged = []
    cur_chr, cur_start, cur_end = None, None, None
    for row in df.iter_rows(named=True):
        c, s, e = row['chr'], row['start'], row['end']
        if c != cur_chr or s > cur_end:
            if cur_chr is not None:
                merged.append((cur_chr, cur_start, cur_end))
            cur_chr, cur_start, cur_end = c, s, e
        else:
            cur_end = max(cur_end, e)
    if cur_chr is not None:
        merged.append((cur_chr, cur_start, cur_end))
    return pl.DataFrame(merged, schema={'chr': pl.Utf8, 'start': pl.Int64, 'end': pl.Int64}, orient='row')

def build_consensus_peaks():
    """
    Build three peak sets for comparison:
    1. DiffBind-only extended ±3kb, merged (current best)
    2. ATAC universe extended ±3kb, merged
    3. Consensus: union of (1) and (2), filtered to ATAC-open regions only
    """
    # Load sources
    db = load_diffbind_peaks()
    atac = load_atac_universe_peaks()

    # --- Set A: DiffBind only ---
    print("\n[1c] Building peak sets...")
    set_a = extend_peaks(db, half_w=3000)
    set_a = merge_overlapping(set_a)
    set_a = set_a.with_columns(
        (pl.col('chr') + ':' + pl.col('start').cast(pl.Utf8) + '-' +
         pl.col('end').cast(pl.Utf8)).alias('peak_id')
    )
    print(f"  Set A (DiffBind, ±3kb merged): {set_a.height} peaks")

    # --- Set B: ATAC universe ---
    set_b = extend_peaks(atac, half_w=3000)
    set_b = merge_overlapping(set_b)
    set_b = set_b.with_columns(
        (pl.col('chr') + ':' + pl.col('start').cast(pl.Utf8) + '-' +
         pl.col('end').cast(pl.Utf8)).alias('peak_id')
    )
    print(f"  Set B (ATAC universe, ±3kb merged): {set_b.height} peaks")

    # --- Set C: Consensus = ATAC + DiffBind union ---
    # Union: take all ATAC peaks, add DiffBind peaks not overlapping ATAC
    print("  Building consensus (ATAC + DiffBind union)...")

    # Merge all together
    union_all = pl.concat([
        set_a.select(['chr','start','end']),
        set_b.select(['chr','start','end'])
    ]).unique().sort(['chr','start'])

    # Re-merge overlapping
    set_c = merge_overlapping(union_all)
    set_c = set_c.with_columns(
        (pl.col('chr') + ':' + pl.col('start').cast(pl.Utf8) + '-' +
         pl.col('end').cast(pl.Utf8)).alias('peak_id')
    )
    print(f"  Set C (consensus ATAC+DiffBind union, merged): {set_c.height} peaks")

    widths = set_c['end'] - set_c['start']
    print(f"  Consensus width: min={widths.min()}, median={widths.median():.0f}, max={widths.max()}")

    return set_a, set_b, set_c


# ============================================================================
#  PART 2 — FRAGMENT COUNTING (shared utility)
# ============================================================================

def _build_peak_index(peaks_df):
    """Build chromosome-indexed peak lookup."""
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

def count_fragments(fragments_path, peaks_df, batch_size=400_000):
    """Count fragment endpoints in peaks. Returns AnnData with cells x peaks."""
    peak_idx = _build_peak_index(peaks_df)
    peak_ids = peaks_df['peak_id'].to_list()
    n_peaks = len(peak_ids)

    with gzip.open(fragments_path, 'rb') as gz:
        buf = io.BytesIO(b''.join(l for l in gz if not l.startswith(b'#')))
    df_full = pl.read_csv(buf, separator='\t', has_header=False,
                          new_columns=['chr','start','end','bc','readSupport'])

    bcs = np.sort(df_full['bc'].unique().to_numpy())
    bc_to_row = {b: i for i, b in enumerate(bcs)}
    rows_list, cols_list, data_list = [], [], []
    n_rows = df_full.height

    for batch_start in range(0, n_rows, batch_size):
        batch = df_full.slice(batch_start, batch_size)
        bc_arr = batch['bc'].to_numpy()
        start_arr = batch['start'].to_numpy()
        end_arr = batch['end'].to_numpy()
        chr_arr = np.array([
            'chr' + c if not c.startswith('chr') else c
            for c in batch['chr'].to_numpy()
        ])
        for chrom in np.unique(chr_arr):
            if chrom not in peak_idx:
                continue
            starts_p, ends_p, global_idx = peak_idx[chrom]
            mask = chr_arr == chrom
            bc_rows = np.array([bc_to_row[b] for b in bc_arr[mask]], dtype=np.int32)
            pc_s = _pos_to_peak(start_arr[mask], starts_p, ends_p)
            pc_e = _pos_to_peak(end_arr[mask], starts_p, ends_p)
            hit_s = pc_s >= 0
            rows_list.append(bc_rows[hit_s])
            cols_list.append(global_idx[pc_s[hit_s]])
            data_list.append(np.ones(hit_s.sum(), dtype=np.uint32))
            hit_e = (pc_e >= 0) & (pc_e != pc_s)
            rows_list.append(bc_rows[hit_e])
            cols_list.append(global_idx[pc_e[hit_e]])
            data_list.append(np.ones(hit_e.sum(), dtype=np.uint32))
        progress = 100 * batch_start / n_rows
        if progress % 10 < 1:
            print(f"    counting: {progress:.0f}%", end='\r')
    print()

    rows_arr = np.concatenate(rows_list).astype(np.int32)
    cols_arr = np.concatenate(cols_list).astype(np.int32)
    data_arr = np.concatenate(data_list).astype(np.uint32)
    X = sps.coo_matrix((data_arr, (rows_arr, cols_arr)),
                        shape=(len(bcs), n_peaks)).tocsr()
    adata = ad.AnnData(X)
    adata.obs.index = bcs
    adata.var.index = peak_ids
    return adata


# ============================================================================
#  PART 3 — VARIABLE PEAK SELECTION (dispersion-based)
# ============================================================================

def select_variable_peaks(ac_adata, me_adata, n_top=10000):
    """
    Select variable peaks using normalized dispersion (similar to
    Seurat's FindVariableFeatures but for sparse count data).

    For each mark separately, compute:
        mean = per-peak mean count
        disp = per-peak variance / mean (overdispersion)
    Then rank peaks by mean-binned dispersion residual, select top N.
    Take the union of variable peaks across marks.
    """
    print("\n[2] Variable peak selection (dispersion-based)...")

    def _var_peaks(adata, n_top=8000, n_bins=20):
        """Find top N variable peaks by dispersion."""
        counts = adata.X
        if sps.issparse(counts):
            means = np.array(counts.mean(axis=0)).flatten()
            # Variance for sparse: E[X^2] - E[X]^2
            sq = counts.multiply(counts)
            sq_mean = np.array(sq.mean(axis=0)).flatten()
            variances = sq_mean - means ** 2
        else:
            means = counts.mean(axis=0)
            variances = counts.var(axis=0, ddof=0)

        # Avoid division by zero
        means = np.maximum(means, 1e-9)
        dispersion = variances / means

        # Bin by mean expression
        valid = means > 0
        mean_log = np.log10(means[valid])
        disp_log = np.log10(dispersion[valid])

        # Bin
        bin_edges = np.percentile(mean_log, np.linspace(0, 100, n_bins + 1))
        bin_idx = np.digitize(mean_log, bin_edges[:-1])

        # Compute bin-specific z-score of dispersion
        disp_z = np.full(len(means), -np.inf)
        for b in range(1, n_bins + 1):
            mask = bin_idx == b
            if mask.sum() < 2:
                continue
            d_bin = disp_log[mask]
            d_mean = np.mean(d_bin)
            d_std = np.std(d_bin, ddof=0)
            if d_std > 0:
                disp_z_full = (np.log10(dispersion) - d_mean) / d_std
                disp_z[valid] = np.where(bin_idx == b,
                                         disp_z_full[valid], disp_z[valid])

        # Select top N
        top_idx = np.argsort(-disp_z)[:n_top]
        return top_idx, means, dispersion, disp_z

    print("  Computing variable peaks for H3K27ac...")
    ac_idx, ac_means, ac_disp, ac_z = _var_peaks(ac_adata, n_top=8000)
    print(f"    H3K27ac variable peaks: {len(ac_idx)}")

    print("  Computing variable peaks for H3K27me3...")
    me_idx, me_means, me_disp, me_z = _var_peaks(me_adata, n_top=8000)
    print(f"    H3K27me3 variable peaks: {len(me_idx)}")

    # Union of variable peaks
    ac_set = set(me_adata.var_names[ac_idx])
    me_set = set(me_adata.var_names[me_idx])
    union = ac_set | me_set
    intersect = ac_set & me_set
    print(f"    Union: {len(union)}, Intersection: {len(intersect)}")

    # Select top N from union by max dispersion z-score
    all_var_names = me_adata.var_names
    union_idx = [i for i, v in enumerate(all_var_names) if v in union]
    union_z = np.maximum(
        np.full(len(all_var_names), -np.inf),
        ac_z if hasattr(ac_z, '__len__') and len(ac_z) == len(all_var_names) else -np.inf
    )
    # Use me_z if available
    if hasattr(me_z, '__len__') and len(me_z) == len(all_var_names):
        for i in range(len(all_var_names)):
            union_z[i] = max(union_z[i], me_z[i] if i < len(me_z) else -np.inf)

    # Top N
    top_n = min(n_top, len(union_idx))
    sorted_union = sorted(union_idx, key=lambda i: -union_z[i])[:top_n]
    var_names_selected = all_var_names[sorted_union]

    print(f"    Selected top {len(var_names_selected)} variable peaks")

    # Plot dispersion
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, means_arr, disp_arr, z_arr, title in [
        (axes[0], ac_means, ac_disp, ac_z, 'H3K27ac'),
        (axes[1], me_means, me_disp, me_z, 'H3K27me3'),
    ]:
        valid = means_arr > 0
        ax.scatter(np.log10(means_arr[valid]), np.log10(disp_arr[valid]),
                   s=1, alpha=0.3, c='gray', rasterized=True)
        # Highlight selected
        sel_mask = np.zeros(len(means_arr), dtype=bool)
        sel_mask[ac_idx if 'ac' in title.lower() else me_idx] = True
        ax.scatter(np.log10(means_arr[sel_mask & valid]),
                   np.log10(disp_arr[sel_mask & valid]),
                   s=2, c='red', alpha=0.5, rasterized=True)
        ax.set_xlabel('log10(mean)')
        ax.set_ylabel('log10(dispersion)')
        ax.set_title(f'{title} — variable peaks (red)')
    plt.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, 'variable_peaks_dispersion.png'), dpi=150)
    plt.close()
    print(f"  Saved variable_peaks_dispersion.png")

    return var_names_selected


# ============================================================================
#  PART 4 — DIMENSIONALITY REDUCTION
# ============================================================================

def per_mark_tfidf_lsi(adata, layer_name, n_components=50):
    """
    Standard scATAC-seq TF-IDF → LSI for a single mark.
    Returns the LSI embedding and depth correlation diagnostics.
    """
    print(f"\n[3a] TF-IDF → LSI for {layer_name}...")
    X = adata.layers[layer_name] if layer_name in adata.layers else adata.X
    if sps.issparse(X):
        X = X.astype(np.float64)
    else:
        X = sps.csr_matrix(X)

    # TF-IDF
    cell_totals = np.array(X.sum(axis=1)).flatten()
    cell_totals[cell_totals == 0] = 1
    tf = diags(1.0 / cell_totals) @ X

    n_cells = X.shape[0]
    bin_counts = np.array((X > 0).sum(axis=0)).flatten()
    idf = np.log1p(n_cells / (1.0 + bin_counts))
    tfidf = tf @ diags(idf)

    # SVD on TF-IDF
    from sklearn.utils.extmath import randomized_svd
    U, S, Vt = randomized_svd(tfidf, n_components=n_components,
                               random_state=42, n_iter=8)
    lsi = U * S

    # Depth correlation for each component
    depth = np.log1p(cell_totals)
    depth_corrs = np.array([pearsonr(depth, lsi[:, i])[0] for i in range(lsi.shape[1])])

    print(f"  LSI shape: {lsi.shape}")
    print(f"  Depth correlations: "
          f"PC0={depth_corrs[0]:.3f}, PC1={depth_corrs[1]:.3f}, "
          f"PC2={depth_corrs[2]:.3f}, PC3={depth_corrs[3]:.3f}")

    return lsi, depth_corrs


def cca_joint_embedding(ac_lsi, me_lsi, n_shared=20, n_unique=10):
    """
    CCA-based joint embedding that separates shared and mark-specific signals.

    Shared: CCA of the two LSIs (canonical correlation)
    Unique: Residual variance after removing shared components

    Returns concatenated [shared, ac_unique, me_unique] embedding.
    """
    print("\n[3b] CCA joint embedding (shared + mark-specific)...")

    from sklearn.cross_decomposition import CCA

    n_min = min(ac_lsi.shape[1], me_lsi.shape[1], 50)
    n_shared = min(n_shared, n_min - 1)

    # CCA
    cca = CCA(n_components=n_shared, max_iter=500, tol=1e-6)
    ac_c, me_c = cca.fit_transform(ac_lsi[:, :n_min], me_lsi[:, :n_min])

    # Canonical correlations
    can_corrs = []
    for i in range(n_shared):
        corr = pearsonr(ac_c[:, i], me_c[:, i])[0]
        can_corrs.append(corr)
    print(f"  Canonical correlations: "
          f"{[f'{c:.3f}' for c in can_corrs[:6]]}...")

    # Signal-specific residuals: regress out shared from each
    from sklearn.linear_model import LinearRegression

    def residualize(X, shared):
        """Remove shared signal from X."""
        lr = LinearRegression().fit(shared, X)
        pred = lr.predict(shared)
        return X - pred

    n_unique = min(n_unique, ac_lsi.shape[1] - n_shared, me_lsi.shape[1] - n_shared)

    ac_resid = residualize(ac_lsi[:, :n_min], ac_c)
    me_resid = residualize(me_lsi[:, :n_min], me_c)

    # PCA on residuals for mark-specific signals
    from sklearn.decomposition import PCA
    ac_pca = PCA(n_components=n_unique).fit_transform(ac_resid)
    me_pca = PCA(n_components=n_unique).fit_transform(me_resid)

    # Concatenate
    joint = np.concatenate([ac_c, ac_pca, me_pca], axis=1)
    print(f"  Joint embedding: {joint.shape} "
          f"(shared={n_shared}, ac_unique={n_unique}, me_unique={n_unique})")

    return joint, np.array(can_corrs)


def spectral_embedding(adata, layers, n_components=50):
    """Multiview spectral embedding (current method, for comparison)."""
    print("\n[3c] Multiview spectral embedding (baseline)...")
    sys.path.insert(0, os.environ.get("SCIT_PATH",
        "/data/ebaird/scentinel/nanoCT/20260522.nanoCT/scit_src"))
    import src as scit

    scit.em.multiview_spectral(adata, layers)
    scit.tl.add_metadata(adata)
    scit.tl.remove_pc(adata, 'X_multi_spectral', 0)

    embed = adata.obsm['X_multi_spectral']
    print(f"  Spectral embedding: {embed.shape}")
    return embed


def assess_pc_quality(adata, lsi_embedding, layer_name, plot_name):
    """
    Comprehensive PC quality assessment:
    1. Depth correlation
    2. Variance explained
    3. Cluster separation index (silhouette on each PC)
    Returns a quality score for each component.
    """
    print(f"\n[3d] PC quality assessment for {layer_name}...")

    n_pcs = lsi_embedding.shape[1]
    X = adata.layers[layer_name] if layer_name in adata.layers else adata.X
    cell_totals = np.array(X.sum(axis=1)).flatten()
    depth = np.log1p(cell_totals)

    # Metrics
    depth_corrs = np.array([
        pearsonr(depth, lsi_embedding[:, i])[0] for i in range(n_pcs)
    ])
    var_explained = np.var(lsi_embedding, axis=0)
    var_frac = var_explained / var_explained.sum()

    # Identify technical PCs: |depth correlation| > 0.3
    tech_pcs = np.where(np.abs(depth_corrs) > 0.3)[0].tolist()
    print(f"  Technical PCs (|depth corr| > 0.3): {tech_pcs}")
    print(f"  Top depth corrs: {[f'PC{i}:{depth_corrs[i]:.2f}' for i in range(min(10,n_pcs))]}")

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    axes[0].bar(range(n_pcs), np.abs(depth_corrs))
    axes[0].axhline(0.3, color='red', linestyle='--', alpha=0.5)
    axes[0].set_xlabel('PC')
    axes[0].set_ylabel('|Depth correlation|')
    axes[0].set_title(f'{layer_name} PC depth correlation')

    axes[1].plot(range(n_pcs), var_frac, 'o-')
    axes[1].set_xlabel('PC')
    axes[1].set_ylabel('Variance fraction')
    axes[1].set_title(f'{layer_name} variance explained')

    axes[2].scatter(lsi_embedding[:, 0], lsi_embedding[:, 1],
                    c=depth, s=2, cmap='viridis', rasterized=True)
    axes[2].set_xlabel('PC1')
    axes[2].set_ylabel('PC2')
    axes[2].set_title(f'{layer_name} PC1/2 colored by depth')
    plt.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, f'pc_quality_{plot_name}.png'), dpi=150)
    plt.close()
    print(f"  Saved pc_quality_{plot_name}.png")

    return {'depth_corrs': depth_corrs, 'var_frac': var_frac, 'tech_pcs': tech_pcs}


# ============================================================================
#  PART 5 — CLUSTERING & EVALUATION
# ============================================================================

def run_clustering(adata, embedding_key, resolution=0.8, n_neighbors=20,
                   min_dist=0.1, spread=1.0):
    """
    Run UMAP + Leiden clustering on a given embedding.
    Returns the adata with obsm X_umap_{suffix} and obs leiden_{suffix}.
    """
    import scanpy as sc

    emb = adata.obsm[embedding_key]
    suffix = embedding_key.replace('X_', '')

    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=embedding_key)
    sc.tl.umap(adata, min_dist=min_dist, spread=spread)

    sc.tl.leiden(adata, resolution=resolution, key_added=f'leiden_{suffix}')

    n_clusters = adata.obs[f'leiden_{suffix}'].nunique()
    print(f"  {suffix}: {n_clusters} clusters (res={resolution}, k={n_neighbors})")

    return adata


def evaluate_clustering(adata, embedding_key, suffix, plot_dir):
    """
    Evaluate cluster quality: silhouette score, stability, biological separation.
    """
    from sklearn.metrics import silhouette_score, silhouette_samples

    emb = adata.obsm[embedding_key]
    labels = adata.obs[f'leiden_{suffix}'].values
    n_clusters = len(np.unique(labels))

    metrics = {}
    metrics['n_clusters'] = n_clusters

    # Silhouette score on the embedding itself
    if n_clusters > 1 and n_clusters < len(labels):
        sil = silhouette_score(emb, labels, metric='euclidean')
        metrics['silhouette'] = sil
        print(f"  Silhouette score ({suffix}): {sil:.3f}")
    else:
        metrics['silhouette'] = np.nan

    # Per-mark separation: do clusters separate by H3K27ac / H3K27me3 depth?
    for layer in ['acet', 'meth']:
        if layer in adata.layers:
            counts = np.array(adata.layers[layer].sum(axis=1)).flatten()
            # ANOVA-like: between-cluster variance / total variance
            grand_mean = np.mean(counts)
            between = 0
            total = 0
            for lbl in np.unique(labels):
                mask = labels == lbl
                grp = counts[mask]
                between += len(grp) * (np.mean(grp) - grand_mean) ** 2
                total += np.sum((grp - grand_mean) ** 2)
            sep = between / max(total, 1e-9)
            metrics[f'{layer}_separation'] = sep
            print(f"  {layer} separation: {sep:.3f}")

    # UMAP plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 11))
    umap_key = f'X_umap'  # after sc.tl.umap was run

    if umap_key in adata.obsm:
        umap = adata.obsm[umap_key]

        # Clusters
        from matplotlib.patches import Patch
        colors = plt.cm.tab20(np.linspace(0, 1, n_clusters))
        for i in range(n_clusters):
            mask = labels == str(i)
            axes[0, 0].scatter(umap[mask, 0], umap[mask, 1],
                               c=[colors[i % 20]], s=3, rasterized=True)
        axes[0, 0].set_title(f'{suffix} — {n_clusters} clusters')
        axes[0, 0].set_xticks([]); axes[0, 0].set_yticks([])

        # H3K27ac
        if 'acet' in adata.layers:
            ac = np.array(adata.layers['acet'].sum(axis=1)).flatten()
            sc1 = axes[0, 1].scatter(umap[:, 0], umap[:, 1], c=np.log1p(ac),
                                      s=2, cmap='viridis', rasterized=True)
            axes[0, 1].set_title('H3K27ac signal')
            plt.colorbar(sc1, ax=axes[0, 1])

        # H3K27me3
        if 'meth' in adata.layers:
            me = np.array(adata.layers['meth'].sum(axis=1)).flatten()
            sc2 = axes[1, 0].scatter(umap[:, 0], umap[:, 1], c=np.log1p(me),
                                      s=2, cmap='magma', rasterized=True)
            axes[1, 0].set_title('H3K27me3 signal')
            plt.colorbar(sc2, ax=axes[1, 0])

        # Acet/Meth ratio
        if 'acet' in adata.layers and 'meth' in adata.layers:
            ratio = ac / (ac + me + 1)
            sc3 = axes[1, 1].scatter(umap[:, 0], umap[:, 1], c=ratio,
                                      s=2, cmap='RdBu_r', vmin=0, vmax=1,
                                      rasterized=True)
            axes[1, 1].set_title('Acet/(Acet+Meth) ratio')
            plt.colorbar(sc3, ax=axes[1, 1])

    plt.tight_layout()
    fig.savefig(os.path.join(plot_dir, f'umap_{suffix}.png'), dpi=150)
    plt.close()

    return metrics


def resolution_sweep(adata, embedding_key, resolutions, suffix=''):
    """Sweep resolution and report cluster counts + stability."""
    import scanpy as sc

    results = []
    for res in resolutions:
        suffix_key = f'leiden_{suffix}_r{res}' if suffix else f'leiden_r{res}'
        sc.tl.leiden(adata, resolution=res, key_added=suffix_key)
        n = adata.obs[suffix_key].nunique()
        results.append({'resolution': res, 'n_clusters': n})
        print(f"  res={res:.1f}: {n} clusters")

    return pd.DataFrame(results)


# ============================================================================
#  MAIN
# ============================================================================

def main():
    print("=" * 70)
    print("  nanoCT Clustering Improvement Pipeline")
    print("=" * 70)
    t0 = time.time()

    # ------------------------------------------------------------------
    # 1. Peak engineering
    # ------------------------------------------------------------------
    set_a, set_b, set_c = build_consensus_peaks()

    # Use consensus set (C) for primary analysis, set A as baseline comparison
    primary_peaks = set_c
    baseline_peaks = set_a

    print(f"\nUsing consensus peaks: {primary_peaks.height}")
    print(f"Baseline (DiffBind only): {baseline_peaks.height}")

    # ------------------------------------------------------------------
    # 2. Count fragments in consensus peaks
    # ------------------------------------------------------------------
    print("\n--- Counting fragments (consensus peaks) ---")
    print("H3K27ac...")
    ac = count_fragments(f"{BASE}/H3K27ac/fragments.tsv.gz", primary_peaks)
    print(f"  H3K27ac: {ac.n_obs} cells x {ac.n_vars} peaks")

    print("H3K27me3...")
    me = count_fragments(f"{BASE}/H3K27me3/fragments.tsv.gz", primary_peaks)
    print(f"  H3K27me3: {me.n_obs} cells x {me.n_vars} peaks")

    # Filter cells present in both
    common_cells = sorted(set(ac.obs_names) & set(me.obs_names))
    ac = ac[common_cells].copy()
    me = me[common_cells].copy()
    print(f"\nCommon cells: {len(common_cells)}")

    # ------------------------------------------------------------------
    # 3. Cell QC filtering
    # ------------------------------------------------------------------
    print("\n--- Cell QC ---")
    for name, ad in [('H3K27ac', ac), ('H3K27me3', me)]:
        totals = np.array(ad.X.sum(axis=1)).flatten()
        print(f"  {name}: median={np.median(totals):.0f}, "
              f"p5={np.percentile(totals,5):.0f}, "
              f"p95={np.percentile(totals,95):.0f}")

    # Filter: remove cells with very low or very high counts
    ac_totals = np.array(ac.X.sum(axis=1)).flatten()
    me_totals = np.array(me.X.sum(axis=1)).flatten()
    keep = (
        (ac_totals >= 20) & (ac_totals <= 3000) &
        (me_totals >= 15) & (me_totals <= 3000)
    )
    print(f"  Cells passing QC: {keep.sum()}/{len(keep)}")
    ac = ac[keep].copy()
    me = me[keep].copy()

    # Also filter peaks: keep peaks with >=5 counts in at least one mark
    ac_peak_counts = np.array((ac.X > 0).sum(axis=0)).flatten()
    me_peak_counts = np.array((me.X > 0).sum(axis=0)).flatten()
    peak_keep = (ac_peak_counts >= 3) | (me_peak_counts >= 3)
    print(f"  Peaks passing min-count filter: {peak_keep.sum()}/{len(peak_keep)}")
    ac = ac[:, peak_keep].copy()
    me = me[:, peak_keep].copy()

    # ------------------------------------------------------------------
    # 4. Variable peak selection
    # ------------------------------------------------------------------
    var_peaks = select_variable_peaks(ac, me, n_top=10000)
    var_peak_idx = [i for i, v in enumerate(me.var_names) if v in set(var_peaks)]
    print(f"  Variable peaks retained: {len(var_peak_idx)}")

    # Subset
    ac = ac[:, var_peak_idx].copy()
    me = me[:, var_peak_idx].copy()

    # ------------------------------------------------------------------
    # 5. Create multimodal AnnData
    # ------------------------------------------------------------------
    print("\n--- Building multimodal AnnData ---")
    # Stack with proper per-mark normalization
    # Normalize each mark to counts per million (CPM) before combining
    from scipy.sparse import csr_matrix

    def normalize_cpm(X):
        """Normalize to counts per million."""
        X = X.astype(np.float64) if sps.issparse(X) else csr_matrix(X, dtype=np.float64)
        totals = np.array(X.sum(axis=1)).flatten()
        totals[totals == 0] = 1
        return diags(1e6 / totals) @ X

    # Keep raw counts in layers, normalized in X
    adata = ad.AnnData(
        X=normalize_cpm(ac.X) + normalize_cpm(me.X),  # combined normalized
        obs=ac.obs.copy(),
        var=ac.var.copy(),
        layers={
            'acet': ac.X.copy(),
            'meth': me.X.copy(),
            'acet_cpm': normalize_cpm(ac.X),
            'meth_cpm': normalize_cpm(me.X),
        }
    )
    print(f"  Combined: {adata.n_obs} cells x {adata.n_vars} peaks")

    # ------------------------------------------------------------------
    # 6. Dimensionality reduction — multiple strategies
    # ------------------------------------------------------------------

    # Strategy A: Multiview spectral embedding (baseline)
    print("\n" + "=" * 60)
    print("  Strategy A: Multiview Spectral (baseline)")
    print("=" * 60)
    sys.path.insert(0, os.environ.get("SCIT_PATH",
        "/data/ebaird/scentinel/nanoCT/20260522.nanoCT/scit_src"))
    import src as scit
    scit.em.multiview_spectral(adata, ['acet', 'meth'])
    scit.tl.add_metadata(adata)
    scit.tl.remove_pc(adata, 'X_multi_spectral', 0)
    adata.obsm['X_spectral'] = adata.obsm['X_multi_spectral'].copy()
    print(f"  Spectral embedding: {adata.obsm['X_spectral'].shape}")

    # Strategy B: Per-mark LSI → CCA joint embedding
    print("\n" + "=" * 60)
    print("  Strategy B: LSI + CCA (shared + specific factors)")
    print("=" * 60)
    ac_lsi, ac_dc = per_mark_tfidf_lsi(adata, 'acet', n_components=50)
    me_lsi, me_dc = per_mark_tfidf_lsi(adata, 'meth', n_components=50)

    # PC quality assessment
    ac_qc = assess_pc_quality(adata, ac_lsi, 'acet', 'acet')
    me_qc = assess_pc_quality(adata, me_lsi, 'meth', 'meth')

    # Keep only non-technical PCs (|depth corr| <= 0.3)
    ac_keep = np.abs(ac_dc) <= 0.3
    me_keep = np.abs(me_dc) <= 0.3
    n_ac_keep = ac_keep.sum()
    n_me_keep = me_keep.sum()
    print(f"  Keeping {n_ac_keep}/{len(ac_keep)} H3K27ac PCs, "
          f"{n_me_keep}/{len(me_keep)} H3K27me3 PCs")

    if n_ac_keep >= 5 and n_me_keep >= 5:
        ac_lsi_clean = ac_lsi[:, ac_keep]
        me_lsi_clean = me_lsi[:, me_keep]

        # CCA joint embedding
        n_shared = min(15, n_ac_keep - 1, n_me_keep - 1)
        n_unique = min(8, n_ac_keep - n_shared, n_me_keep - n_shared, 8)
        if n_shared > 0:
            joint_emb, can_corrs = cca_joint_embedding(
                ac_lsi_clean, me_lsi_clean,
                n_shared=n_shared,
                n_unique=max(1, n_unique)
            )
            adata.obsm['X_cca_joint'] = joint_emb
        else:
            print("  WARNING: Not enough non-technical PCs for CCA, falling back to concat")
            adata.obsm['X_cca_joint'] = np.concatenate([ac_lsi_clean, me_lsi_clean], axis=1)
    else:
        print("  WARNING: Too many technical PCs, using concatenated LSI")
        adata.obsm['X_cca_joint'] = np.concatenate([ac_lsi, me_lsi], axis=1)

    # Strategy C: Concatenated clean LSI (simple alternative)
    print("\n" + "=" * 60)
    print("  Strategy C: Concatenated clean LSI")
    print("=" * 60)
    ac_clean = ac_lsi[:, ac_keep] if n_ac_keep > 0 else ac_lsi
    me_clean = me_lsi[:, me_keep] if n_me_keep > 0 else me_lsi
    adata.obsm['X_lsi_concat'] = np.concatenate([ac_clean, me_clean], axis=1)
    print(f"  Concat LSI: {adata.obsm['X_lsi_concat'].shape}")

    # ------------------------------------------------------------------
    # 7. Clustering with each strategy
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("  Clustering comparison")
    print("=" * 60)

    strategies = {
        'spectral': 'X_spectral',
        'cca_joint': 'X_cca_joint',
        'lsi_concat': 'X_lsi_concat',
    }

    all_metrics = []
    for name, emb_key in strategies.items():
        print(f"\n--- {name} ---")
        adata = run_clustering(adata, emb_key, resolution=0.8, n_neighbors=20)
        metrics = evaluate_clustering(adata, emb_key, name, PLOT_DIR)
        all_metrics.append({'strategy': name, **metrics})

        # Resolution sweep
        print(f"  Resolution sweep:")
        res_df = resolution_sweep(
            adata, emb_key,
            resolutions=[0.3, 0.5, 0.8, 1.2, 1.6, 2.0],
            suffix=name
        )
        res_df.to_csv(os.path.join(PLOT_DIR, f'resolution_sweep_{name}.csv'), index=False)

    # Summary comparison
    metrics_df = pd.DataFrame(all_metrics)
    metrics_df.to_csv(os.path.join(PLOT_DIR, 'strategy_comparison.csv'), index=False)
    print("\n" + "=" * 60)
    print("  Strategy Comparison Summary")
    print("=" * 60)
    print(metrics_df.to_string(index=False))

    # ------------------------------------------------------------------
    # 8. Save
    # ------------------------------------------------------------------
    out_path = os.path.join(OUT, "combined_improved.h5ad")
    print(f"\nSaving to: {out_path}")
    adata.write_h5ad(out_path)
    print(f"  {adata.n_obs} cells x {adata.n_vars} peaks")

    # ------------------------------------------------------------------
    # 9. Combined comparison UMAP
    # ------------------------------------------------------------------
    print("\n--- Generating comparison figure ---")
    fig, axes = plt.subplots(2, 3, figsize=(18, 11))
    for idx, (name, emb_key) in enumerate(strategies.items()):
        row, col = divmod(idx, 3)
        umap_key = 'X_umap'
        if umap_key in adata.obsm:
            umap = adata.obsm[umap_key]
            label_key = f'leiden_{name}'
            if label_key in adata.obs:
                labels = adata.obs[label_key].values
                n_cl = len(np.unique(labels))
                for i in range(n_cl):
                    mask = labels == str(i)
                    axes[row, col].scatter(
                        umap[mask, 0], umap[mask, 1],
                        c=[plt.cm.tab20(i % 20)], s=3, rasterized=True
                    )
                axes[row, col].set_title(f'{name} ({n_cl} clusters)')
                axes[row, col].set_xticks([]); axes[row, col].set_yticks([])

    # Bottom row: signal overlays on best strategy
    best = list(strategies.keys())[1]  # cca_joint is our hypothesis
    if f'leiden_{best}' in adata.obs and 'X_umap' in adata.obsm:
        umap = adata.obsm['X_umap']
        ac_t = np.array(adata.layers['acet'].sum(axis=1)).flatten()
        me_t = np.array(adata.layers['meth'].sum(axis=1)).flatten()
        sc1 = axes[1, 0].scatter(umap[:, 0], umap[:, 1], c=np.log1p(ac_t),
                                  s=2, cmap='viridis', rasterized=True)
        axes[1, 0].set_title(f'{best}: H3K27ac')
        plt.colorbar(sc1, ax=axes[1, 0])
        sc2 = axes[1, 1].scatter(umap[:, 0], umap[:, 1], c=np.log1p(me_t),
                                  s=2, cmap='magma', rasterized=True)
        axes[1, 1].set_title(f'{best}: H3K27me3')
        plt.colorbar(sc2, ax=axes[1, 1])
        sc3 = axes[1, 2].scatter(umap[:, 0], umap[:, 1],
                                  c=ac_t / (ac_t + me_t + 1), s=2,
                                  cmap='RdBu_r', vmin=0, vmax=1, rasterized=True)
        axes[1, 2].set_title(f'{best}: Acet/Meth ratio')
        plt.colorbar(sc3, ax=axes[1, 2])

    plt.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, 'comparison_overview.png'), dpi=150)
    plt.close()
    print(f"Saved comparison_overview.png")

    elapsed = time.time() - t0
    print(f"\n{'=' * 60}")
    print(f"  Done! Elapsed: {elapsed / 60:.1f} minutes")
    print(f"  Output: {out_path}")
    print(f"  Plots:  {PLOT_DIR}/")
    print(f"{'=' * 60}")


if __name__ == '__main__':
    main()
