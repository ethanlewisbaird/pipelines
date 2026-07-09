"""
Create per-cluster bigWig coverage tracks for scGLUE leiden clusters 0-5.
Dual-mark: H3K27ac + H3K27me3 from Cell Ranger ATAC fragment files.

Interval-based approach: build coverage as sorted intervals, merge overlapping,
then write bigWig. Avoids full-genome per-base arrays.
"""
import scanpy as sc
import pyBigWig
import gzip
import os
from collections import defaultdict

BASE = '/data/ebaird/scRNAseq/20260522.nanoCT'
FRAG_DIR = f'{BASE}/SU.analysis.2026.05.22/Vasso_nanoCT_nanoscope'
OUT = f'{BASE}/cluster_bigwigs/scglue_leiden'
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

print("Loading chromatin object...", flush=True)
adata = sc.read_h5ad(f'{BASE}/scglue_chrom.h5ad')
clusters_of_interest = [0, 1, 2, 3, 4, 5]

bc_to_cluster = {}
for cluster in clusters_of_interest:
    bcs = adata.obs_names[adata.obs['leiden'] == cluster].tolist()
    for bc in bcs:
        bc_to_cluster[bc] = cluster
total_bcs = len(bc_to_cluster)
cluster_counts = {c: sum(1 for v in bc_to_cluster.values() if v == c) for c in clusters_of_interest}
print(f"Loaded {total_bcs} barcodes from clusters 0-5: {cluster_counts}", flush=True)

os.makedirs(f'{OUT}/H3K27ac', exist_ok=True)
os.makedirs(f'{OUT}/H3K27me3', exist_ok=True)

def merge_intervals(intervals):
    """Merge overlapping intervals with same depth into a bedGraph line."""
    if not intervals:
        return []
    intervals.sort(key=lambda x: (x[0], x[1]))
    merged = []
    cur_start, cur_end, cur_val = intervals[0]
    for start, end, val in intervals[1:]:
        if start <= cur_end and val == cur_val:
            cur_end = max(cur_end, end)
        else:
            merged.append((cur_start, cur_end, cur_val))
            cur_start, cur_end, cur_val = start, end, val
    merged.append((cur_start, cur_end, cur_val))
    return merged

def process_fragments(frag_path, mark_label):
    """Single-pass fragment → per-cluster interval accumulation → bigWigs."""
    print(f"\n=== Processing {mark_label} ===", flush=True)
    
    # Store intervals as list of (start, end) per cluster per chrom
    intervals = {}
    for c in clusters_of_interest:
        intervals[c] = {chrom: [] for chrom in UCSC_CHROMS}
    
    n_total = 0
    n_matched = 0
    n_matched_cluster = defaultdict(int)
    
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
                n_matched_cluster[cluster] += 1
                intervals[cluster][ucsc_chrom].append((start, end))
    
    print(f"  Total fragments: {n_total}", flush=True)
    print(f"  Matched: {n_matched}", flush=True)
    print(f"  Per cluster: {dict(n_matched_cluster)}", flush=True)
    
    for c in clusters_of_interest:
        out_path = f'{OUT}/{mark_label}/c{c}_{mark_label}.bw'
        print(f"  Cluster {c}: building bigWig ...", flush=True)
        
        bw = pyBigWig.open(out_path, 'w')
        bw.addHeader([(chrom, chrom_sizes[chrom]) for chrom in UCSC_CHROMS])
        
        for chrom in UCSC_CHROMS:
            ivs = intervals[c][chrom]
            if not ivs:
                continue
            # Sort by start position
            ivs.sort(key=lambda x: x[0])
            
            # Sweep-line to compute depth at each position
            # Use event-based approach: +1 at starts, -1 at ends
            events = []
            for start, end in ivs:
                events.append((start, 1))
                events.append((end, -1))
            events.sort(key=lambda x: (x[0], x[1]))
            
            # Build bedGraph: iterate sweep, record depth segments
            depth = 0
            prev_pos = None
            starts, ends, values = [], [], []
            for pos, delta in events:
                if depth > 0 and prev_pos is not None and pos > prev_pos:
                    starts.append(prev_pos)
                    ends.append(pos)
                    values.append(float(depth))
                depth += delta
                prev_pos = pos
            
            # Write in chunks per chromosome
            if starts:
                bw.addEntries([chrom] * len(starts), starts, ends=ends, values=values)
        
        bw.close()
        print(f"    Done ({os.path.getsize(out_path)} B)", flush=True)
    
    del intervals

process_fragments(f'{FRAG_DIR}/H3K27ac/fragments.tsv.gz', 'H3K27ac')
process_fragments(f'{FRAG_DIR}/H3K27me3/fragments.tsv.gz', 'H3K27me3')

print("\nDone! All bigWigs created.")
