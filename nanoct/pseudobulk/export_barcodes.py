import pandas as pd
import os

def export_barcodes(csv_path, output_dir, prefix):
    df = pd.read_csv(csv_path)
    for cluster in df['Cluster'].unique():
        cluster_barcodes = df[df['Cluster'] == cluster]['Barcode']
        output_path = os.path.join(output_dir, f"{prefix}_cluster_{cluster}.txt")
        cluster_barcodes.to_csv(output_path, index=False, header=False)
        print(f"Exported {len(cluster_barcodes)} barcodes to {output_path}")

# H3K27ac
export_barcodes(
    "/data/ebaird/scRNAseq/20260522.nanoCT/SU.analysis.2026.05.22/Vasso_nanoCT_nanoscope/H3K27ac/analysis/clustering/graphclust/clusters.csv",
    "/data/ebaird/scRNAseq/20260522.nanoCT/SU.analysis.2026.05.22/Vasso_nanoCT_nanoscope/pseudobulk/H3K27ac",
    "H3K27ac"
)

# H3K27me3
export_barcodes(
    "/data/ebaird/scRNAseq/20260522.nanoCT/SU.analysis.2026.05.22/Vasso_nanoCT_nanoscope/H3K27me3/analysis/clustering/graphclust/clusters.csv",
    "/data/ebaird/scRNAseq/20260522.nanoCT/SU.analysis.2026.05.22/Vasso_nanoCT_nanoscope/pseudobulk/H3K27me3",
    "H3K27me3"
)
