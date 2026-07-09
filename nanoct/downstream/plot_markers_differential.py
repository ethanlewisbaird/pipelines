"""
Generate plots for marker and differential analysis.
1. Heatmap of top 10 marker bins per cluster (acet + meth)
2. Volcano plots for pairwise: 4 vs 5, 5 vs 1 (acet + meth)
3. Dotplot of marker genes per cluster
"""
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.sparse import issparse
import gzip, os

plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.bbox'] = 'tight'

BASE = '/data/ebaird/scRNAseq/20260522.nanoCT'
OUT = f'{BASE}/cluster_markers/figures'
MARKERS_CSV = f'{BASE}/cluster_markers/top10_markers_per_cluster.csv'
DIFF_CSV = f'{BASE}/cluster_markers/pairwise_differential.csv'
os.makedirs(OUT, exist_ok=True)

# ── 1. Heatmap of top 10 marker bins ──
print("Generating heatmaps...", flush=True)
adata = sc.read_h5ad(f'{BASE}/scglue_chrom.h5ad')
clusters_oi = [0, 1, 2, 3, 4, 5]
adata = adata[adata.obs['leiden'].isin(clusters_oi)].copy()

markers_df = pd.read_csv(MARKERS_CSV)

for mark_label, layer_name in [('H3K27ac', 'acet'), ('H3K27me3', 'meth')]:
    top_bins = []
    for c in clusters_oi:
        top = markers_df[(markers_df['mark'] == mark_label) & (markers_df['cluster'] == c)]
        top = top.sort_values('rank').head(10)
        top_bins.extend(top['bin_name'].tolist())
    
    gene_labels = []
    for b in top_bins:
        row = markers_df[markers_df['bin_name'] == b].iloc[0]
        label = f"{row['nearest_gene']} ({b.split(':')[0]}:{b.split(':')[1].split('-')[0]})"
        gene_labels.append(label)
    
    # Get expression data for these bins
    X = adata.layers[layer_name]
    if issparse(X):
        X = X.toarray()
    
    # Find bin indices
    bin_indices = [list(adata.var_names).index(b) for b in top_bins]
    subX = X[:, bin_indices]
    
    # Normalize per cell (CPM-like)
    row_sums = subX.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    subX = subX / row_sums * 10000
    
    # Sort cells by cluster
    order = adata.obs['leiden'].argsort().values
    subX = subX[order]
    clusters_sorted = adata.obs['leiden'].values[order]
    
    # Plot
    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.imshow(subX.T, aspect='auto', cmap='Reds' if mark_label == 'H3K27ac' else 'Blues',
                   interpolation='nearest', vmin=0)
    
    # Cluster boundaries
    unique_clusters = sorted(adata.obs['leiden'].unique())
    boundaries = []
    for c in unique_clusters:
        boundaries.append(np.sum(clusters_sorted == c))
    boundaries = np.cumsum(boundaries)
    
    for b in boundaries[:-1]:
        ax.axvline(b - 0.5, color='gray', linewidth=0.5, linestyle='--')
    
    # Labels
    ax.set_yticks(range(len(top_bins)))
    ax.set_yticklabels(gene_labels, fontsize=7)
    ax.set_xlabel('Cells (sorted by cluster)')
    ax.set_title(f'{mark_label} — Top 10 marker bins per cluster')
    
    # Cluster labels on top
    tick_positions = []
    start = 0
    for c, end in enumerate(boundaries):
        tick_positions.append((start + end - 1) / 2)
        start = end
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(tick_positions)
    ax2.set_xticklabels([f'L{c}' for c in unique_clusters], fontsize=8)
    ax2.set_xlabel('Leiden cluster')
    
    plt.colorbar(im, ax=ax, label='Norm. counts', shrink=0.6)
    plt.tight_layout()
    plt.savefig(f'{OUT}/{mark_label}_marker_heatmap.png')
    plt.close()
    print(f"  Saved {OUT}/{mark_label}_marker_heatmap.png", flush=True)

# ── 2. Volcano plots ──
print("\nGenerating volcano plots...", flush=True)
diff_df = pd.read_csv(DIFF_CSV)

comparisons = [
    ('5', '1', 'H3K27ac'),
    ('4', '5', 'H3K27ac'),
    ('5', '1', 'H3K27me3'),
    ('4', '5', 'H3K27me3'),
]

colors = {'H3K27ac': ('#d62728', '#ffb3b3'), 'H3K27me3': ('#1f77b4', '#b3d4ff')}

for ca, cb, mark in comparisons:
    df = diff_df[(diff_df['cluster_a'] == ca) & (diff_df['cluster_b'] == cb) & (diff_df['mark'] == mark)]
    if len(df) == 0:
        continue
    
    # Top genes for labeling
    top_n = df.dropna(subset=['padj']).nsmallest(15, 'padj')
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Non-significant
    ns = df[df['padj'] >= 0.05]
    ax.scatter(ns['log2FC'], -np.log10(ns['pval']), c='gray', alpha=0.3, s=3, label='NS')
    
    # Significant (padj < 0.05)
    sig = df[df['padj'] < 0.05]
    up = sig[sig['log2FC'] > 0]
    down = sig[sig['log2FC'] < 0]
    ax.scatter(up['log2FC'], -np.log10(up['pval']), c=colors[mark][0], alpha=0.5, s=5, label=f'↑ in {ca} (n={len(up)})')
    ax.scatter(down['log2FC'], -np.log10(down['pval']), c=colors[mark][1], alpha=0.5, s=5, label=f'↓ in {ca} (n={len(down)})')
    
    # Label top genes
    for _, row in top_n.iterrows():
        gene = row['nearest_gene']
        if gene == 'intergenic':
            continue
        ax.annotate(gene, (row['log2FC'], -np.log10(row['pval'])),
                    fontsize=6, alpha=0.8, ha='center', va='bottom')
    
    ax.axhline(-np.log10(0.05 / len(df)), color='red', linestyle='--', linewidth=0.5, alpha=0.5, label='Bonferroni')
    ax.axvline(0, color='black', linewidth=0.5)
    ax.set_xlabel('log2(FC)')
    ax.set_ylabel('-log10(p-value)')
    ax.set_title(f'{mark} — Cluster {ca} vs {cb}')
    ax.legend(fontsize=7, loc='upper right')
    
    plt.tight_layout()
    plt.savefig(f'{OUT}/{mark}_volcano_{ca}_vs_{cb}.png')
    plt.close()
    print(f"  Saved {OUT}/{mark}_volcano_{ca}_vs_{cb}.png", flush=True)

# ── 3. Dotplot of selected marker genes ──
print("\nGenerating dotplot...", flush=True)
# Pick the top marker gene per cluster per mark (the most differentially expressed)
for mark_label, layer_name in [('H3K27ac', 'acet'), ('H3K27me3', 'meth')]:
    top_genes = []
    for c in clusters_oi:
        top = markers_df[(markers_df['mark'] == mark_label) & (markers_df['cluster'] == c)]
        if len(top) > 0:
            gene = top.iloc[0]['nearest_gene']
            if gene != 'intergenic':
                top_genes.append(gene)
    
    # Also add top pairwise differential genes
    for ca, cb, m in comparisons:
        if m != mark_label: continue
        sub = diff_df[(diff_df['cluster_a'] == ca) & (diff_df['cluster_b'] == cb) & (diff_df['mark'] == m)]
        if len(sub) == 0: continue
        for _, row in sub.head(5).iterrows():
            g = row['nearest_gene']
            if g != 'intergenic' and g not in top_genes:
                top_genes.append(g)
    
    if len(top_genes) == 0:
        continue
    
    # For dotplot, find bins near these genes and extract expression
    # Simple approach: get mean expression per cluster for each gene's nearest bin
    gene_bins = {}
    for g in top_genes:
        rows = markers_df[markers_df['nearest_gene'] == g]
        if len(rows) > 0:
            gene_bins[g] = rows.iloc[0]['bin_name']
    
    bin_indices = [list(adata.var_names).index(b) for b in gene_bins.values()]
    bin_names = list(gene_bins.values())
    gene_names = list(gene_bins.keys())
    
    X = adata.layers[layer_name]
    if issparse(X):
        X = X.toarray()
    subX = X[:, bin_indices]
    
    # Per cluster mean and pct
    leiden = adata.obs['leiden'].values
    means = []
    pcts = []
    for c in clusters_oi:
        mask = leiden == c
        means.append(np.array(subX[mask].mean(axis=0)).flatten())
        pcts.append(np.array((subX[mask] > 0).mean(axis=0)).flatten() * 100)
    means = np.array(means)
    pcts = np.array(pcts)
    
    # Normalize means per gene (0-1)
    for j in range(means.shape[1]):
        col = means[:, j]
        if col.max() > col.min():
            means[:, j] = (col - col.min()) / (col.max() - col.min())
    
    fig, ax = plt.subplots(figsize=(max(8, len(gene_names) * 0.5), 4))
    for gi, gn in enumerate(gene_names):
        for ci, c in enumerate(clusters_oi):
            size = pcts[ci, gi]
            color = means[ci, gi]
            ax.scatter(gi, ci, s=size * 10, c=[color], cmap='Reds' if mark_label == 'H3K27ac' else 'Blues',
                       vmin=0, vmax=1, edgecolors='gray', linewidth=0.5)
    
    ax.set_yticks(range(len(clusters_oi)))
    ax.set_yticklabels([f'L{c}' for c in clusters_oi])
    ax.set_xticks(range(len(gene_names)))
    ax.set_xticklabels(gene_names, rotation=45, ha='right', fontsize=8)
    ax.set_title(f'{mark_label} — Selected marker genes')
    ax.set_xlabel('Gene')
    ax.set_ylabel('Leiden cluster')
    
    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=5, label='10%'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=10, label='50%'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=15, label='100%'),
    ]
    ax.legend(handles=legend_elements, title='% cells', fontsize=6, title_fontsize=7, loc='upper right')
    
    plt.tight_layout()
    plt.savefig(f'{OUT}/{mark_label}_dotplot.png')
    plt.close()
    print(f"  Saved {OUT}/{mark_label}_dotplot.png", flush=True)

print("\nAll plots done!", flush=True)
