"""
Create RPKM-normalized per-cluster bigWig coverage tracks.
Uses merged clusters (1+5 combined) from combined_dim_reduced_merged.h5ad.

RPKM = fragments per kilobase of bin per million fragments in cluster.
"""
import scanpy as sc
import pyBigWig
import gzip
import os
from collections import defaultdict

BASE = '/data/ebaird/scRNAseq/20260522.nanoCT'
FRAG_DIR = f'{BASE}/SU.analysis.2026.05.22/Vasso_nanoCT_nanoscope'
OUT = f'{BASE}/cluster_bigwigs/reclustered_merged_rpkm'
CHROM_SIZES = f'{BASE}/dm6.chrom.sizes'

CHROM_MAP = {
    '2L': 'chr2L', '2R': 'chr2R', '3L': 'chr3L', '3R': 'chr3R',
    'X': 'chrX', 'Y': 'chrY', '4': 'chr4'
}
UCSC_CHROMS = list(CHROM_MAP.values())

chrom_sizes = {}
with open(CHROM_SIZES) as f:
    for line in f:
        chrom, size = line.strip().split()
        chrom_sizes[chrom] = int(size)

print("Loading merged chromatin object...", flush=True)
adata = sc.read_h5ad(f'{BASE}/analysis_05.26/combined_dim_reduced_merged.h5ad')
clusters_of_interest = sorted(adata.obs['leiden_merged'].unique())
print(f"Clusters: {clusters_of_interest}", flush=True)

bc_to_cluster = {}
for cluster in clusters_of_interest:
    bcs = adata.obs_names[adata.obs['leiden_merged'] == cluster].tolist()
    for bc in bcs:
        bc_to_cluster[bc] = cluster
total_bcs = len(bc_to_cluster)
cluster_counts = {c: sum(1 for v in bc_to_cluster.values() if v == c) for c in clusters_of_interest}
print(f"Loaded {total_bcs} barcodes: {cluster_counts}", flush=True)

os.makedirs(f'{OUT}/H3K27ac', exist_ok=True)
os.makedirs(f'{OUT}/H3K27me3', exist_ok=True)

BIN_SIZE = 5000  # 5kb bins for RPKM normalization

def process_fragments(frag_path, mark_label):
    print(f"\n=== Processing {mark_label} ===", flush=True)
    
    # Count fragments per bin per cluster
    bin_counts = {c: defaultdict(int) for c in clusters_of_interest}
    cluster_total = defaultdict(int)
    
    n_total = 0
    n_matched = 0
    
    with gzip.open(frag_path, 'rt') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 5:
                continue
            frag_chrom = parts[0]
            ucsc_chrom = CHROM_MAP.get(frag_chrom)
            if ucsc_chrom is None:
                continue
            start = int(parts[1])
            end = int(parts[2])
            barcode = parts[3]
            
            n_total += 1
            cluster = bc_to_cluster.get(barcode)
            if cluster is not None:
                n_matched += 1
                # Assign fragment to bins it overlaps
                bin_start = (start // BIN_SIZE) * BIN_SIZE
                bin_end = ((end // BIN_SIZE) + 1) * BIN_SIZE
                for b in range(bin_start, bin_end, BIN_SIZE):
                    bin_counts[cluster][(ucsc_chrom, b)] += 1
                cluster_total[cluster] += 1
    
    print(f"  Total fragments: {n_total}", flush=True)
    print(f"  Matched: {n_matched}", flush=True)
    print(f"  Per cluster: {dict(cluster_total)}", flush=True)
    
    for c in clusters_of_interest:
        out_path = f'{OUT}/{mark_label}/c{c}_{mark_label}_rpkm.bw'
        print(f"  Cluster {c} (n={cluster_total[c]}): building RPKM bigWig ...", flush=True)
        
        bw = pyBigWig.open(out_path, 'w')
        bw.addHeader([(chrom, chrom_sizes[chrom]) for chrom in UCSC_CHROMS])
        
        # Group bins by chrom
        chrom_bins = defaultdict(list)
        for (chrom, pos), count in bin_counts[c].items():
            chrom_bins[chrom].append((pos, count))
        
        for chrom in UCSC_CHROMS:
            bins = sorted(chrom_bins.get(chrom, []))
            if not bins:
                continue
            
            starts = []
            ends = []
            values = []
            
            for pos, count in bins:
                # RPKM = count / (bin_length_kb * total_fragments_M)
                bin_len_kb = BIN_SIZE / 1000.0
                total_M = cluster_total[c] / 1e6
                rpkm = count / (bin_len_kb * total_M) if total_M > 0 else 0
                
                starts.append(pos)
                ends.append(min(pos + BIN_SIZE, chrom_sizes[chrom]))
                values.append(float(rpkm))
            
            if starts:
                bw.addEntries([chrom] * len(starts), starts, ends=ends, values=values)
        
        bw.close()
        size_mb = os.path.getsize(out_path) / 1e6
        print(f"    Done ({size_mb:.1f} MB)", flush=True)
    
    del bin_counts

process_fragments(f'{FRAG_DIR}/H3K27ac/fragments.tsv.gz', 'H3K27ac')
process_fragments(f'{FRAG_DIR}/H3K27me3/fragments.tsv.gz', 'H3K27me3')

print("\nDone! All RPKM bigWigs created.")
