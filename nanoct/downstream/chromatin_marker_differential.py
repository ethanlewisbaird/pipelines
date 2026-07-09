"""
Marker bin analysis for scGLUE chromatin clusters 0-5.
1. Top 10 marker bins per cluster (scanpy rank_genes_groups) for acet & meth
2. Pairwise differential: cluster 4 vs 5, cluster 5 vs 1
3. Map significant bins to nearest genes via dm6 GTF
"""
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.sparse import issparse
from scipy.stats import mannwhitneyu, ranksums
import gzip
import os
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

BASE = '/data/ebaird/scRNAseq/20260522.nanoCT'
GTF = f'{BASE}/dm6.refGene.gtf.gz'
OUT = f'{BASE}/cluster_markers'
os.makedirs(OUT, exist_ok=True)

# ── 1. Load chromatin object ──
print("Loading chromatin object...", flush=True)
adata = sc.read_h5ad(f'{BASE}/scglue_chrom.h5ad')
clusters_of_interest = [0, 1, 2, 3, 4, 5]
adata = adata[adata.obs['leiden'].isin(clusters_of_interest)].copy()
print(f"Cells in clusters 0-5: {adata.n_obs}", flush=True)

# ── 2. Load GTF and build bin→gene mapping ──
print("Loading GTF and mapping bins to genes...", flush=True)
gene_tss = {}
with gzip.open(GTF, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9 or parts[2] != 'transcript':
            continue
        chrom = parts[0]
        if not chrom.startswith('chr'):
            chrom = 'chr' + chrom
        if chrom not in ['chr2L','chr2R','chr3L','chr3R','chrX','chrY','chr4']:
            continue
        attr = parts[8]
        gn = None
        for a in attr.split(';'):
            a = a.strip()
            if a.startswith('gene_name "'):
                gn = a.split('"')[1]
                break
        if gn is None:
            continue
        tss = int(parts[3]) if parts[6] == '+' else int(parts[4])
        if chrom not in gene_tss:
            gene_tss[chrom] = []
        gene_tss[chrom].append((tss, gn))

for chrom in gene_tss:
    gene_tss[chrom].sort(key=lambda x: x[0])

def nearest_gene(chrom, pos):
    if chrom not in gene_tss or not gene_tss[chrom]:
        return 'intergenic'
    genes = gene_tss[chrom]
    lo, hi = 0, len(genes) - 1
    while lo < hi:
        mid = (lo + hi) // 2
        if genes[mid][0] < pos:
            lo = mid + 1
        else:
            hi = mid
    best_dist = float('inf')
    best_gene = 'intergenic'
    for i in range(max(0, lo - 2), min(len(genes), lo + 3)):
        dist = abs(genes[i][0] - pos)
        if dist < best_dist:
            best_dist = dist
            best_gene = genes[i][1]
    return best_gene

# Map all bins at once
bin_chroms = adata.var['chr'].values
bin_centers = ((adata.var['start'] + adata.var['end']) // 2).values
bin_to_gene = []
for i in range(len(bin_chroms)):
    c = bin_chroms[i]
    if not c.startswith('chr'):
        c = 'chr' + c
    bin_to_gene.append(nearest_gene(c, bin_centers[i]))
adata.var['nearest_gene'] = bin_to_gene
n_mapped = sum(1 for g in bin_to_gene if g != 'intergenic')
print(f"  {n_mapped}/{len(bin_to_gene)} bins mapped to genes", flush=True)

# ── 3. Marker analysis using scanpy rank_genes_groups ──
def run_markers(layer_name, mark_label):
    print(f"\n=== Markers for {mark_label} ===", flush=True)
    saved_X = adata.X.copy()
    adata.X = adata.layers[layer_name].copy()
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', n_genes=20, use_raw=False)
    
    results = []
    for cluster in clusters_of_interest:
        cluster_str = str(cluster)
        names = adata.uns['rank_genes_groups']['names'][cluster_str]
        scores = adata.uns['rank_genes_groups']['scores'][cluster_str]
        pvals = adata.uns['rank_genes_groups']['pvals'][cluster_str]
        pvals_adj = adata.uns['rank_genes_groups']['pvals_adj'][cluster_str]
        logfoldchanges = adata.uns['rank_genes_groups']['logfoldchanges'][cluster_str]
        
        for i in range(len(names)):
            bin_name = names[i]
            bin_idx = list(adata.var_names).index(bin_name)
            results.append({
                'mark': mark_label,
                'cluster': cluster,
                'rank': i + 1,
                'bin_name': bin_name,
                'chr': adata.var['chr'].iloc[bin_idx],
                'start': adata.var['start'].iloc[bin_idx],
                'end': adata.var['end'].iloc[bin_idx],
                'nearest_gene': bin_to_gene[bin_idx],
                'score': scores[i],
                'log2FC': logfoldchanges[i],
                'pval': pvals[i],
                'padj': pvals_adj[i],
            })
    
    adata.X = saved_X
    return pd.DataFrame(results)

acet_markers_df = run_markers('acet', 'H3K27ac')
meth_markers_df = run_markers('meth', 'H3K27me3')
markers_df = pd.concat([acet_markers_df, meth_markers_df], ignore_index=True)

# Top 10 per cluster per mark
top10 = markers_df[markers_df['rank'] <= 10].copy()
top10.to_csv(f'{OUT}/top10_markers_per_cluster.csv', index=False)
print(f"\nSaved top markers to {OUT}/top10_markers_per_cluster.csv", flush=True)

# ── 4. Pairwise differential for specific comparisons ──
def pairwise_diff(ca, cb, layer):
    X = adata.layers[layer].toarray() if issparse(adata.layers[layer]) else np.array(adata.layers[layer])
    leiden = adata.obs['leiden'].values
    ma, mb = leiden == ca, leiden == cb
    
    pvals = np.full(X.shape[1], np.nan)
    log2fc = np.full(X.shape[1], np.nan)
    mu_a = np.full(X.shape[1], np.nan)
    mu_b = np.full(X.shape[1], np.nan)
    pct_a = np.full(X.shape[1], np.nan)
    pct_b = np.full(X.shape[1], np.nan)
    
    for i in range(X.shape[1]):
        a, b = X[ma, i], X[mb, i]
        if a.sum() == 0 and b.sum() == 0:
            continue
        try:
            _, p = mannwhitneyu(a, b, alternative='two-sided')
        except ValueError:
            continue
        pvals[i] = p
        mu_a[i], mu_b[i] = a.mean(), b.mean()
        pct_a[i] = (a > 0).mean() * 100
        pct_b[i] = (b > 0).mean() * 100
        log2fc[i] = np.log2((a.mean() + 1) / (b.mean() + 1))
    
    df = pd.DataFrame({
        'cluster_a': ca, 'cluster_b': cb,
        'bin_name': adata.var_names,
        'chr': bin_chroms,
        'start': adata.var['start'].values,
        'end': adata.var['end'].values,
        'nearest_gene': bin_to_gene,
        'mean_a': mu_a, 'mean_b': mu_b,
        'pct_a': pct_a, 'pct_b': pct_b,
        'log2FC': log2fc, 'pval': pvals,
    }).dropna(subset=['pval'])
    df['padj'] = (df['pval'] * len(df)).clip(upper=1.0)
    return df.sort_values('pval')

comparisons = [
    (5, 1, 'acet', 'H3K27ac'),
    (4, 5, 'acet', 'H3K27ac'),
    (5, 1, 'meth', 'H3K27me3'),
    (4, 5, 'meth', 'H3K27me3'),
]

all_diff = []
for ca, cb, layer, mark in comparisons:
    print(f"\n=== Differential: cluster {ca} vs {cb} ({mark}) ===", flush=True)
    df = pairwise_diff(ca, cb, layer)
    df['mark'] = mark
    all_diff.append(df)
    n_sig = (df['padj'] < 0.05).sum()
    print(f"  Significant bins (padj<0.05): {n_sig}", flush=True)

diff_df = pd.concat(all_diff, ignore_index=True)
diff_df.to_csv(f'{OUT}/pairwise_differential.csv', index=False)
print(f"\nSaved pairwise differential to {OUT}/pairwise_differential.csv", flush=True)

# ── 5. Pretty print summary ──
print("\n" + "="*80, flush=True)
print("TOP 10 MARKER BINS PER CLUSTER", flush=True)
print("="*80, flush=True)

for mark_label in ['H3K27ac', 'H3K27me3']:
    print(f"\n--- {mark_label} ---", flush=True)
    for cluster in clusters_of_interest:
        top = top10[(top10['mark'] == mark_label) & (top10['cluster'] == cluster)]
        print(f"  Cluster {cluster}:", flush=True)
        for _, row in top.iterrows():
            print(f"    rank {int(row['rank']):2d}: {row['nearest_gene']:20s}  ({row['bin_name']})  "
                  f"log2FC={row['log2FC']:.2f}  padj={row['padj']:.2e}", flush=True)

print("\n" + "="*80, flush=True)
print("PAIRWISE DIFFERENTIAL (top 20 by p-value per comparison)", flush=True)
print("="*80, flush=True)

for ca, cb, layer, mark in comparisons:
    print(f"\n--- Cluster {ca} vs {cb} ({mark}) ---", flush=True)
    df = diff_df[(diff_df['cluster_a'] == ca) & (diff_df['cluster_b'] == cb) & (diff_df['mark'] == mark)]
    n_sig = (df['padj'] < 0.05).sum()
    print(f"  Significant: {n_sig} bins (padj<0.05)", flush=True)
    top20 = df.head(20)
    for _, row in top20.iterrows():
        direction = "↑" if row['log2FC'] > 0 else "↓"
        print(f"  {direction} {row['nearest_gene']:20s}  ({row['bin_name']})  "
              f"log2FC={row['log2FC']:+.2f}  pct_a={row['pct_a']:.0f}%  pct_b={row['pct_b']:.0f}%  "
              f"padj={row['padj']:.2e}", flush=True)

print(f"\nAll results saved to {OUT}/", flush=True)
