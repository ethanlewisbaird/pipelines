#!/usr/bin/env python3
"""
Generate deeptools heatmaps for unmerged clusters (0-5) using reclustered bigWigs.
Same layout as merged_pairwise/deeptools but with 6 clusters instead of 4.
"""
import os
import subprocess
import numpy as np
import pandas as pd
import pyBigWig
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import tempfile

BASE = '/data/ebaird/scentinel/nanoCT/20260522.nanoCT'
BW_DIR = f'{BASE}/cluster_bigwigs/reclustered_leiden'
MARKER_DIR = f'{BASE}/cluster_markers/reclustered'
OUT_DIR = f'{BASE}/cluster_markers/reclustered/deeptools'
os.makedirs(OUT_DIR, exist_ok=True)

GTF = f'{BASE}/dm6.refGene.gtf.gz'
CHROM_SIZES = f'{BASE}/dm6.chrom.sizes'

CLUSTERS = ['c0', 'c1', 'c2', 'c4', 'c5']
CLUSTER_NAMES = ['Cluster 0', 'Cluster 1', 'Cluster 2', 'Cluster 4', 'Cluster 5']
CLUSTER_COLOURS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#9467bd', '#8c564b']

RNA_GROUPS = {
    'metabolic': ['c11', 'c2'],
    'progenitor_diff': ['c6', 'c8', 'c7', 'c13', 'c14'],
    'stem_like': ['c5', 'c0', 'c1', 'c3', 'c4', 'c10']
}



def get_tss_from_bin(bin_str, gene_name):
    chrom, rest = bin_str.split(':')
    start, end = rest.split('-')
    center = (int(start) + int(end)) // 2
    return chrom, center

def create_bed_from_markers(marker_df, bed_path, flank=3000):
    regions = []
    for _, row in marker_df.iterrows():
        gene = row['nearest_gene']
        chrom, tss = get_tss_from_bin(row['bin'], gene)
        start = max(0, tss - flank)
        end = tss + flank
        regions.append(f"{chrom}\t{start}\t{end}\t{gene}\t{row.get('score', 0)}\t+")
    with open(bed_path, 'w') as f:
        f.write('\n'.join(regions) + '\n')
    return len(regions)

def run_compute_matrix(bed_path, bigwig_list, labels, out_prefix, region_size=5000):
    mat_path = f'{out_prefix}.gz'
    compute_matrix = '/data/ebaird/miniconda3/envs/nanoCT/bin/computeMatrix'
    cmd = [
        compute_matrix, 'reference-point',
        '--referencePoint', 'TSS',
        '-S'] + bigwig_list + [
        '-R', bed_path,
        '-a', str(region_size),
        '-b', str(region_size),
        '--binSize', '100',
        '--missingDataAsZero',
        '-o', mat_path,
        '--samplesLabel'] + labels + [
        '--numberOfProcessors', '8'
    ]
    subprocess.run(cmd, check=True)
    return mat_path

def plot_custom_heatmap(mat_path, out_path, title, sort_by_idx=None, sort_ascending=False):
    from deeptools.heatmapper import heatmapper
    hm = heatmapper()
    hm.read_matrix_file(mat_path)

    n_samples = len(hm.matrix.sample_labels)
    n_genes = hm.matrix.matrix.shape[0]
    n_bins = hm.matrix.matrix.shape[1] // n_samples

    all_data = []
    for s in range(n_samples):
        start = s * n_bins
        end = (s + 1) * n_bins
        all_data.append(hm.matrix.matrix[:, start:end])

    gene_labels = [r[2] for r in hm.matrix.regions]

    if sort_by_idx is not None and sort_by_idx < n_samples:
        sort_scores = all_data[sort_by_idx].mean(axis=1)
        order = np.argsort(sort_scores)[::-1] if not sort_ascending else np.argsort(sort_scores)
        all_data = [d[order] for d in all_data]
        gene_labels = [gene_labels[i] for i in order]

    n_clusters = len(CLUSTERS)
    fig = plt.figure(figsize=(28, 28))
    gs = GridSpec(n_clusters, 2, height_ratios=[1]*n_clusters, width_ratios=[1, 1],
                  hspace=0.3, wspace=0.3)

    ac_cmap = 'YlOrRd'
    me_cmap = 'Blues'

    for c_idx in range(n_clusters):
        ax_ac = fig.add_subplot(gs[c_idx, 0])
        ax_me = fig.add_subplot(gs[c_idx, 1])

        ac_data = all_data[c_idx * 2]
        me_data = all_data[c_idx * 2 + 1]

        vmax_ac = min(5, np.percentile(ac_data[ac_data > 0], 95)) if np.any(ac_data > 0) else 5
        vmax_me = min(5, np.percentile(me_data[me_data > 0], 95)) if np.any(me_data > 0) else 5

        ax_ac.imshow(ac_data, aspect='auto', cmap=ac_cmap, vmin=0, vmax=vmax_ac,
                     interpolation='nearest')
        ax_me.imshow(me_data, aspect='auto', cmap=me_cmap, vmin=0, vmax=vmax_me,
                     interpolation='nearest')

        ax_ac.set_title(f'{CLUSTER_NAMES[c_idx]} — H3K27ac', fontsize=10)
        ax_me.set_title(f'{CLUSTER_NAMES[c_idx]} — H3K27me3', fontsize=10)

        if c_idx == 0:
            ax_ac.set_ylabel('Genes', fontsize=9)
            ax_me.set_ylabel('Genes', fontsize=9)

        for ax in [ax_ac, ax_me]:
            ax.set_xticks([0, n_bins//2, n_bins-1])
            ax.set_xticklabels(['-5kb', 'TSS', '+5kb'], fontsize=8)
            ax.tick_params(axis='y', labelsize=6)

        if c_idx == 0:
            ax_ac.set_yticks(range(len(gene_labels)))
            ax_ac.set_yticklabels(gene_labels, fontsize=6)

    fig.suptitle(title, fontsize=14, y=1.01)
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {out_path}")

def png_to_jpg(png_path, jpg_path, quality=95):
    from PIL import Image
    img = Image.open(png_path)
    img = img.convert('RGB')
    img.save(jpg_path, 'JPEG', quality=quality)
    os.remove(png_path)
    print(f"Converted to {jpg_path}")

def plot_with_profiles(mat_path, out_path, title, sort_by_idx=None, sort_ascending=False):
    from deeptools.heatmapper import heatmapper
    hm = heatmapper()
    hm.read_matrix_file(mat_path)

    n_samples = len(hm.matrix.sample_labels)
    n_genes = hm.matrix.matrix.shape[0]
    n_bins = hm.matrix.matrix.shape[1] // n_samples

    all_data = []
    for s in range(n_samples):
        start = s * n_bins
        end = (s + 1) * n_bins
        all_data.append(hm.matrix.matrix[:, start:end])

    gene_labels = [r[2] for r in hm.matrix.regions]

    if sort_by_idx is not None and sort_by_idx < n_samples:
        sort_scores = all_data[sort_by_idx].mean(axis=1)
        order = np.argsort(sort_scores)[::-1] if not sort_ascending else np.argsort(sort_scores)
        all_data = [d[order] for d in all_data]
        gene_labels = [gene_labels[i] for i in order]

    n_clusters = len(CLUSTERS)
    fig = plt.figure(figsize=(28, 32))
    gs = GridSpec(n_clusters + 1, 2, height_ratios=[1]*(n_clusters+1), width_ratios=[1, 1],
                  hspace=0.3, wspace=0.3)

    ac_cmap = 'YlOrRd'
    me_cmap = 'Blues'
    x = np.linspace(-5000, 5000, n_bins)

    ax_ac_profile = fig.add_subplot(gs[0, 0])
    ax_me_profile = fig.add_subplot(gs[0, 1])

    for c_idx in range(n_clusters):
        ac_mean = np.clip(all_data[c_idx * 2], 0, 5).mean(axis=0)
        me_mean = np.clip(all_data[c_idx * 2 + 1], 0, 5).mean(axis=0)
        ax_ac_profile.plot(x, ac_mean, color=CLUSTER_COLOURS[c_idx], linewidth=2,
                          label=CLUSTER_NAMES[c_idx])
        ax_me_profile.plot(x, me_mean, color=CLUSTER_COLOURS[c_idx], linewidth=2,
                          label=CLUSTER_NAMES[c_idx])

    for ax, mark in [(ax_ac_profile, 'H3K27ac'), (ax_me_profile, 'H3K27me3')]:
        ax.axvline(0, color='black', linestyle='--', linewidth=0.5)
        ax.set_title(f'{mark} — Average signal', fontsize=12)
        ax.set_xlabel('Distance from TSS', fontsize=10)
        ax.set_ylabel('Mean signal', fontsize=10)
        ax.legend(fontsize=8, loc='upper left')
        ax.set_xlim(-5000, 5000)

    for c_idx in range(n_clusters):
        ax_ac = fig.add_subplot(gs[c_idx + 1, 0])
        ax_me = fig.add_subplot(gs[c_idx + 1, 1])

        ac_data = all_data[c_idx * 2]
        me_data = all_data[c_idx * 2 + 1]

        vmax_ac = min(5, np.percentile(ac_data[ac_data > 0], 95)) if np.any(ac_data > 0) else 5
        vmax_me = min(5, np.percentile(me_data[me_data > 0], 95)) if np.any(me_data > 0) else 5

        ax_ac.imshow(ac_data, aspect='auto', cmap=ac_cmap, vmin=0, vmax=vmax_ac,
                     interpolation='nearest')
        ax_me.imshow(me_data, aspect='auto', cmap=me_cmap, vmin=0, vmax=vmax_me,
                     interpolation='nearest')

        ax_ac.set_title(f'{CLUSTER_NAMES[c_idx]} — H3K27ac', fontsize=10)
        ax_me.set_title(f'{CLUSTER_NAMES[c_idx]} — H3K27me3', fontsize=10)

        if c_idx == 0:
            ax_ac.set_ylabel('Genes', fontsize=9)
            ax_me.set_ylabel('Genes', fontsize=9)

        for ax in [ax_ac, ax_me]:
            ax.set_xticks([0, n_bins//2, n_bins-1])
            ax.set_xticklabels(['-5kb', 'TSS', '+5kb'], fontsize=8)
            ax.tick_params(axis='y', labelsize=6)

        if c_idx == 0:
            ax_ac.set_yticks(range(len(gene_labels)))
            ax_ac.set_yticklabels(gene_labels, fontsize=6)

    fig.suptitle(title, fontsize=14, y=1.01)
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {out_path}")

def main():
    print("=== Generating deeptools heatmaps for unmerged clusters (excluding cluster 3) ===")

    ac_markers = pd.read_csv(f'{MARKER_DIR}/top10_markers_H3K27ac.csv')
    me_markers = pd.read_csv(f'{MARKER_DIR}/top10_markers_H3K27me3.csv')

    for cluster_idx in [0, 1, 2, 4, 5]:  # Skip cluster 3
        cluster_name = f'c{cluster_idx}'
        print(f"\n--- Processing {cluster_name} ---")

        ac_cl = ac_markers[ac_markers['cluster'] == cluster_idx]
        me_cl = me_markers[me_markers['cluster'] == cluster_idx]

        if len(ac_cl) == 0:
            print(f"  No H3K27ac markers for {cluster_name}, skipping")
            continue

        ac_bed = f'{OUT_DIR}/{cluster_name}_H3K27ac_markers_TSS.bed'
        me_bed = f'{OUT_DIR}/{cluster_name}_H3K27me3_markers_TSS.bed'

        n_ac = create_bed_from_markers(ac_cl, ac_bed)
        n_me = create_bed_from_markers(me_cl, me_bed)
        print(f"  Created BED files: {n_ac} acetyl markers, {n_me} methyl markers")

        ac_bws = [f'{BW_DIR}/H3K27ac/{c}_H3K27ac.bw' for c in CLUSTERS]
        me_bws = [f'{BW_DIR}/H3K27me3/{c}_H3K27me3.bw' for c in CLUSTERS]
        all_bws = []
        for c in CLUSTERS:
            all_bws.append(f'{BW_DIR}/H3K27ac/{c}_H3K27ac.bw')
            all_bws.append(f'{BW_DIR}/H3K27me3/{c}_H3K27me3.bw')

        all_labels = []
        for c in CLUSTERS:
            all_labels.append(f'{c} H3K27ac')
            all_labels.append(f'{c} H3K27me3')

        combined_bed = f'{OUT_DIR}/{cluster_name}_combined_markers_TSS.bed'
        with open(combined_bed, 'w') as f:
            with open(ac_bed) as fa:
                f.write(fa.read())
            with open(me_bed) as fm:
                f.write(fm.read())

        mat_path = f'{OUT_DIR}/matrix_combined_{cluster_name}_5kbpTSS.gz'
        if not os.path.exists(mat_path):
            print(f"  Computing matrix for {cluster_name}...")
            run_compute_matrix(combined_bed, all_bws, all_labels,
                             f'{OUT_DIR}/matrix_combined_{cluster_name}_5kbpTSS')
        else:
            print(f"  Matrix already exists for {cluster_name}")

        # Acetyl marker heatmap
        jpg_path = f'{OUT_DIR}/heatmap_{cluster_name}_acetyl_marker_genes_sorted.jpg'
        if not os.path.exists(jpg_path):
            print(f"  Plotting {cluster_name} acetyl heatmap...")
            png_path = jpg_path.replace('.jpg', '.png')
            plot_with_profiles(mat_path, png_path,
                             f'{CLUSTER_NAMES[CLUSTERS.index(cluster_name)]} acetyl marker genes (n={n_ac})',
                             sort_by_idx=CLUSTERS.index(cluster_name) * 2)
            png_to_jpg(png_path, jpg_path)

        # Methyl marker heatmap
        me_mat_path = f'{OUT_DIR}/matrix_methyl_{cluster_name}_5kbpTSS.gz'
        if not os.path.exists(me_mat_path):
            print(f"  Computing methyl matrix for {cluster_name}...")
            run_compute_matrix(me_bed, all_bws, all_labels,
                             f'{OUT_DIR}/matrix_methyl_{cluster_name}_5kbpTSS')
        
        jpg_path = f'{OUT_DIR}/heatmap_{cluster_name}_methyl_marker_genes_sorted.jpg'
        if not os.path.exists(jpg_path):
            print(f"  Plotting {cluster_name} methyl heatmap...")
            png_path = jpg_path.replace('.jpg', '.png')
            plot_with_profiles(me_mat_path, png_path,
                             f'{CLUSTER_NAMES[CLUSTERS.index(cluster_name)]} methyl marker genes (n={n_me})',
                             sort_by_idx=CLUSTERS.index(cluster_name) * 2 + 1)
            png_to_jpg(png_path, jpg_path)

    print("\n=== Processing RNA group markers ===")
    merged_deeptools = f'{BASE}/cluster_markers/merged_pairwise/deeptools'
    for group_name in RNA_GROUPS.keys():
        print(f"\n--- RNA group: {group_name} ---")
        bed_path = f'{merged_deeptools}/rna_{group_name}_markers_TSS.bed'

        if not os.path.exists(bed_path):
            print(f"  BED file not found: {bed_path}")
            continue

        all_bws = []
        all_labels = []
        for c in CLUSTERS:
            all_bws.append(f'{BW_DIR}/H3K27ac/{c}_H3K27ac.bw')
            all_bws.append(f'{BW_DIR}/H3K27me3/{c}_H3K27me3.bw')
            all_labels.append(f'{c} H3K27ac')
            all_labels.append(f'{c} H3K27me3')

        mat_path = f'{OUT_DIR}/matrix_combined_rna_{group_name}_5kbpTSS.gz'
        if not os.path.exists(mat_path):
            print(f"  Computing matrix for {group_name}...")
            run_compute_matrix(bed_path, all_bws, all_labels,
                             f'{OUT_DIR}/matrix_combined_rna_{group_name}_5kbpTSS')
        else:
            print(f"  Matrix already exists for {group_name}")

        jpg_path = f'{OUT_DIR}/heatmap_rna_{group_name}_marker_genes_sorted.jpg'
        if not os.path.exists(jpg_path):
            print(f"  Plotting {group_name} heatmap...")
            png_path = jpg_path.replace('.jpg', '.png')
            plot_with_profiles(mat_path, png_path,
                             f'RNA {group_name} marker genes')
            png_to_jpg(png_path, jpg_path)

    print("\n=== Processing ATAC differential regions ===")
    atac_files = {
        'increase_2_3': '/data/vtheodorou/2025.05.ATAC.tumor/clusters.2.3. increase during differentiation.bed',
        'decrease_5_6': '/data/vtheodorou/2025.05.ATAC.tumor/clusters.5.6 decrease during differentiation.bed'
    }

    for name, bed_file in atac_files.items():
        print(f"\n--- ATAC: {name} ---")

        all_bws = []
        all_labels = []
        for c in CLUSTERS:
            all_bws.append(f'{BW_DIR}/H3K27ac/{c}_H3K27ac.bw')
            all_bws.append(f'{BW_DIR}/H3K27me3/{c}_H3K27me3.bw')
            all_labels.append(f'{c} H3K27ac')
            all_labels.append(f'{c} H3K27me3')

        mat_path = f'{OUT_DIR}/matrix_atac_{name}_5kbpTSS.gz'
        if not os.path.exists(mat_path):
            print(f"  Computing matrix for {name}...")
            run_compute_matrix(bed_file, all_bws, all_labels,
                             f'{OUT_DIR}/matrix_atac_{name}_5kbpTSS')
        else:
            print(f"  Matrix already exists for {name}")

        jpg_path = f'{OUT_DIR}/heatmap_atac_{name}.jpg'
        if not os.path.exists(jpg_path):
            print(f"  Plotting {name} heatmap...")
            png_path = jpg_path.replace('.jpg', '.png')
            plot_with_profiles(mat_path, png_path,
                             f'ATAC {name} regions')
            png_to_jpg(png_path, jpg_path)

    print("\n=== All done! ===")

if __name__ == '__main__':
    main()
