#!/usr/bin/env Rscript
#===================================================================
# Single-cell nanoCUT&Tag Analysis — Peak-based (Signac/Seurat)
# SCENTINEL Project — Delidakis Lab, IMBB-FORTH
# Drosophila melanogaster (dm6) — H3K27ac + H3K27me3
#
# This script uses called peaks (not 5kb bins) as features.
# Integrates both histone marks for joint analysis.
#===================================================================

# ---- 0. Setup ----
library(Signac)
library(Seurat)
library(GenomicRanges)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(future)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(Matrix)
set.seed(42)
options(future.globals.maxSize = 15 * 1024^3)
plan("multicore", workers = 8)

BASE_DIR <- "/data/ebaird/scRNAseq/20260522.nanoCT/SU.analysis.2026.05.22/Vasso_nanoCT_nanoscope"
OUT_DIR  <- "/data/ebaird/scRNAseq/20260522.nanoCT/R_analysis_peaks/output"
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Gene Annotations ----
message("=== 1. Loading dm6 gene annotations from GTF ===")

gtf_file <- "/data/ebaird/scRNAseq/20260522.nanoCT/R_analysis_peaks/dm6_genes.gtf.gz"
annotations_raw <- rtracklayer::import(gtf_file)

# UCSC GTF has transcript-level features. Take transcript ranges and
# collapse to gene-level ranges (min start, max end per gene).
tx <- annotations_raw[annotations_raw$type == "transcript"]

# Split by gene_id FIRST, then compute ranges
tx_by_gene <- split(tx, tx$gene_id)
gene_ranges <- unlist(reduce(range(tx_by_gene)))
gene_ranges$gene_id <- names(gene_ranges)

# Get gene names
gn <- unique(mcols(tx)[, c("gene_id", "gene_name")])
m <- match(gene_ranges$gene_id, gn$gene_id)
gene_ranges$gene_name <- gn$gene_name[m]

annotations <- gene_ranges

# Add required metadata columns for Signac
annotations$gene_biotype <- "protein_coding"

# Strip "chr" prefix to match fragment file naming
seqlevels(annotations) <- gsub("^chr", "", seqlevels(annotations))

# Keep only standard chromosomes
seqlevels_to_keep <- c("2L", "2R", "3L", "3R", "4", "X", "Y")
annotations <- keepSeqlevels(annotations,
                              intersect(seqlevels(annotations), seqlevels_to_keep),
                              pruning.mode = "coarse")
# Ensure gene_name column exists (UCSC uses gene_id for Ensembl ID and gene_name separately)
if (is.null(annotations$gene_name) && !is.null(annotations$gene_id)) {
  annotations$gene_name <- annotations$gene_id
}
# Get TSS positions (1bp at transcription start)
tss_positions <- resize(annotations, width = 1, fix = "start")
message("  Annotations: ", length(annotations), " genes")
message("  TSS positions: ", length(tss_positions))

# ---- 2. Load data per mark ----
message("=== 2. Loading marks: H3K27ac + H3K27me3 ===")

marks <- c("H3K27ac", "H3K27me3")
seurat_list <- list()

for (mark in marks) {
  message("  Loading ", mark, "...")
  path <- file.path(BASE_DIR, mark)
  
  frag_file  <- file.path(path, "fragments.tsv.gz")
  peak_file  <- file.path(path, "peaks.bed")
  meta_file  <- file.path(path, "singlecell.csv")
  
  # Load peaks (skip CellRanger headers)
  peaks <- read.table(peak_file, col.names = c("chr", "start", "end"), comment.char = "#")
  peaks_gr <- makeGRangesFromDataFrame(peaks)
  
  # Load CellRanger metadata
  metadata <- read.csv(meta_file, header = TRUE, row.names = 1)
  # Exclude NO_BARCODE row
  metadata <- metadata[rownames(metadata) != "NO_BARCODE", ]
  
  # Filter: use is__cell_barcode for clean cells, or at least passed_filters > 500
  # CellRanger calls ~1725 cells per mark as true barcodes
  # passed_filters > 500 gives ~3000-4000 cells with meaningful signal
  cells_pass <- rownames(metadata)[metadata$passed_filters > 500]
  message("    Cells with >500 passing fragments: ", length(cells_pass))
  
  # Fragment object
  frags <- CreateFragmentObject(path = frag_file, cells = cells_pass)
  
  # Count matrix (peaks x cells)
  counts <- FeatureMatrix(fragments = frags, features = peaks_gr, 
                          cells = cells_pass, verbose = FALSE)
  message("    Count matrix: ", nrow(counts), " peaks x ", ncol(counts), " cells")
  
  # ChromatinAssay — no aggressive built-in filtering, we'll QC manually
  assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    fragments = frags,
    annotation = annotations
  )
  
  # Seurat object
  obj <- CreateSeuratObject(counts = assay, assay = "peaks")
  obj$mark <- mark
  
  # Add CellRanger QC metadata (select useful columns)
  qc_cols <- intersect(c("total", "duplicate", "chimeric", "unmapped", "lowmapq",
                          "mitochondrial", "passed_filters", "is__cell_barcode",
                          "peak_region_fragments", "peak_region_cutsites",
                          "TSS_fragments", "promoter_region_fragments",
                          "blacklist_region_fragments", "on_target_fragments"),
                       colnames(metadata))
  for (col in qc_cols) {
    obj <- AddMetaData(obj, metadata[Cells(obj), col], col.name = col)
  }
  
  # Compute QC metrics
  obj <- TSSEnrichment(obj, tss.positions = tss_positions, fast = TRUE)
  obj <- NucleosomeSignal(obj)
  
  # Compute FRiP (Fraction of Reads in Peaks)
  obj$FRiP <- obj$peak_region_fragments / obj$passed_filters
  
  seurat_list[[mark]] <- obj
}

# ---- 3. QC Diagnostics ----
message("=== 3. QC diagnostics ===")

# Per-mark QC plots
for (mark in marks) {
  obj <- seurat_list[[mark]]
  
  p1 <- VlnPlot(obj, features = "nCount_peaks", pt.size = 0.01) + 
    ggtitle(paste0(mark, " — Fragments per cell")) + NoLegend()
  p2 <- VlnPlot(obj, features = "TSS.enrichment", pt.size = 0.01) + 
    ggtitle(paste0(mark, " — TSS Enrichment")) + NoLegend()
  p3 <- VlnPlot(obj, features = "nucleosome_signal", pt.size = 0.01) + 
    ggtitle(paste0(mark, " — Nucleosome Signal")) + NoLegend()
  p4 <- VlnPlot(obj, features = "FRiP", pt.size = 0.01) + 
    ggtitle(paste0(mark, " — FRiP")) + NoLegend()
  
  ggsave(file.path(FIG_DIR, paste0("QC_", mark, ".png")), 
         wrap_plots(p1, p2, p3, p4, ncol = 2), width = 12, height = 10)
  
  # Print QC summary
  message("  ", mark, " QC summary:")
  message("    Median fragments/cell: ", median(obj$nCount_peaks))
  message("    Median TSS enrichment: ", median(obj$TSS.enrichment, na.rm = TRUE))
  message("    Median nucleosome signal: ", median(obj$nucleosome_signal, na.rm = TRUE))
  message("    Median FRiP: ", median(obj$FRiP, na.rm = TRUE))
}

# ---- 4. Cell Filtering ----
message("=== 4. Cell filtering ===")

# CUT&Tag is inherently sparser and noisier than ATAC-seq.
# Use very lenient thresholds — only remove extreme outliers.
# Heavy filtering will be done downstream by clustering.

filtered_list <- list()
for (mark in marks) {
  obj <- seurat_list[[mark]]
  
  # Get percentiles for guidance
  ncounts <- obj$nCount_peaks
  tss <- obj$TSS.enrichment[!is.na(obj$TSS.enrichment)]
  
  lo_counts <- quantile(ncounts, 0.02)
  hi_counts <- quantile(ncounts, 0.98)
  lo_tss <- ifelse(length(tss) > 0, quantile(tss, 0.02, na.rm = TRUE), 0)
  
  message("  ", mark, " count thresholds: ", round(lo_counts), " - ", round(hi_counts))
  message("  ", mark, " TSS p2: ", round(lo_tss, 2))
  
  # Keep: exclude extreme outliers (bottom/top 2% by counts)
  # Don't filter on TSS yet — let's just remove obvious noise
  cells_keep <- colnames(obj)[
    obj$nCount_peaks >= lo_counts &
    obj$nCount_peaks <= hi_counts
  ]
  
  filtered_list[[mark]] <- subset(obj, cells = cells_keep)
  message("  ", mark, " cells after filtering: ", length(cells_keep), 
          " (", round(100 * length(cells_keep) / ncol(obj), 1), "% retained)")
}

# ---- 5. Multi-modal Integration ----
message("=== 5. Multi-modal integration ===")

# Find common cells across marks
common_cells <- Reduce(intersect, lapply(filtered_list, Cells))
message("  Common cells: ", length(common_cells))

# Build combined object: first mark as base, second as additional assay
combined <- filtered_list[["H3K27ac"]][, common_cells]
combined <- RenameAssays(combined, peaks = "H3K27ac")

# Add H3K27me3 as second assay — subset the object first, then extract assay
me_subset <- subset(filtered_list[["H3K27me3"]], cells = common_cells)
combined[["H3K27me3"]] <- me_subset[["peaks"]]

# Add gene annotations for gene activity
Annotation(combined[["H3K27ac"]]) <- annotations
Annotation(combined[["H3K27me3"]]) <- annotations

message("  Combined object: ", ncol(combined), " cells")
message("    H3K27ac: ", nrow(combined[["H3K27ac"]]), " peaks")
message("    H3K27me3: ", nrow(combined[["H3K27me3"]]), " peaks")

# ---- 6. Normalization + Dimension Reduction ----
message("=== 6. TF-IDF normalization and SVD ===")

# H3K27ac
DefaultAssay(combined) <- "H3K27ac"
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = "q5")  # Keep top 95% features
combined <- RunSVD(combined, n = 50, reduction.name = "svd_H3K27ac")

# H3K27me3
DefaultAssay(combined) <- "H3K27me3"
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = "q5")
combined <- RunSVD(combined, n = 30, reduction.name = "svd_H3K27me3")

# Depth correlation diagnostics
p1 <- DepthCor(combined, reduction = "svd_H3K27ac") + ggtitle("H3K27ac SVD depth correlation")
p2 <- DepthCor(combined, reduction = "svd_H3K27me3") + ggtitle("H3K27me3 SVD depth correlation")
ggsave(file.path(FIG_DIR, "SVD_depth_correlation.png"), wrap_plots(p1, p2), width = 12, height = 5)

# Determine how many SVD components to remove (usually the first)
# Signac's DepthCor shows correlation of each component with depth
# Component 1 is almost always depth-correlated for ATAC/CUT&Tag

# ---- 7. Joint Embedding ----
message("=== 7. Joint embedding ===")

# Concatenate SVD reductions (removing depth-correlated component 1)
joint_embedding <- cbind(
  Embeddings(combined, "svd_H3K27ac")[, 2:30],   # Remove PC1
  Embeddings(combined, "svd_H3K27me3")[, 2:15]   # Remove PC1, fewer dims for smaller feature set
)

colnames(joint_embedding) <- paste0("joint_", seq_len(ncol(joint_embedding)))
combined[["joint"]] <- CreateDimReducObject(
  embeddings = joint_embedding,
  key = "joint_",
  assay = DefaultAssay(combined)
)

# UMAP on joint embedding
combined <- RunUMAP(combined, reduction = "joint", dims = 1:ncol(joint_embedding),
                    reduction.name = "umap_joint", n.neighbors = 30, min.dist = 0.3)

# Clustering on joint graph
combined <- FindNeighbors(combined, reduction = "joint", dims = 1:ncol(joint_embedding),
                          graph.name = "joint_snn")
combined <- FindClusters(combined, graph.name = "joint_snn", resolution = 0.5,
                         algorithm = 1)  # 1 = Louvain

p <- DimPlot(combined, reduction = "umap_joint", group.by = "seurat_clusters", 
             label = TRUE, pt.size = 0.3) + 
  ggtitle("Joint H3K27ac + H3K27me3 (peaks)") + NoLegend()
ggsave(file.path(FIG_DIR, "umap_joint_clusters.png"), p, width = 8, height = 7)

message("  Joint clusters found: ", length(unique(combined$seurat_clusters)))
print(table(combined$seurat_clusters))

# ---- 8. Per-mark UMAPs ----
message("=== 8. Per-mark visualizations ===")

# H3K27ac only UMAP
combined <- RunUMAP(combined, reduction = "svd_H3K27ac", dims = 2:30,
                    reduction.name = "umap_H3K27ac", n.neighbors = 30)
p_ac <- DimPlot(combined, reduction = "umap_H3K27ac", group.by = "seurat_clusters",
                label = TRUE, pt.size = 0.3) + ggtitle("H3K27ac only")

# H3K27me3 only UMAP
combined <- RunUMAP(combined, reduction = "svd_H3K27me3", dims = 2:20,
                    reduction.name = "umap_H3K27me3", n.neighbors = 30)
p_me <- DimPlot(combined, reduction = "umap_H3K27me3", group.by = "seurat_clusters",
                label = TRUE, pt.size = 0.3) + ggtitle("H3K27me3 only")

ggsave(file.path(FIG_DIR, "umap_per_mark.png"), wrap_plots(p_ac, p_me, ncol = 2), 
       width = 16, height = 7)

# Mark contribution to joint embedding
p_counts <- FeaturePlot(combined, reduction = "umap_joint",
                        features = c("nCount_peaks", "TSS.enrichment", "FRiP"),
                        ncol = 3, pt.size = 0.2)
ggsave(file.path(FIG_DIR, "umap_joint_QCfeatures.png"), p_counts, width = 18, height = 5)

# ---- 9. Differential Peak Accessibility ----
message("=== 9. Differential peak accessibility ===")

# Find peaks that differ between clusters for each mark
Idents(combined) <- "seurat_clusters"

for (mark in marks) {
  message("  ", mark, " differential peaks...")
  DefaultAssay(combined) <- mark
  
  # Find all markers (each cluster vs rest)
  markers <- FindAllMarkers(
    combined,
    only.pos = TRUE,
    min.pct = 0.05,
    logfc.threshold = 0.25,
    test.use = "LR",  # Logistic regression — accounts for latent vars
    latent.vars = "nCount_peaks",
    verbose = FALSE
  )
  
  if (nrow(markers) > 0) {
    write.csv(markers, file.path(OUT_DIR, paste0(mark, "_differential_peaks.csv")), 
              row.names = FALSE)
    
    # Top markers per cluster
    top_markers <- markers %>%
      group_by(cluster) %>%
      slice_max(n = 20, order_by = avg_log2FC)
    
    write.csv(top_markers, file.path(OUT_DIR, paste0(mark, "_top20_peaks.csv")), 
              row.names = FALSE)
    
    message("    Total DE peaks: ", nrow(markers))
    print(table(markers$cluster))
  } else {
    message("    No significantly DE peaks found with current thresholds")
  }
}

# ---- 10. Gene Activity ----
message("=== 10. Gene activity scores ===")

# Use annotation stored in the assay (no need to pass features explicitly)
DefaultAssay(combined) <- "H3K27ac"
gene_activities <- GeneActivity(combined)
combined[["ACTIVITY_H3K27ac"]] <- CreateAssayObject(counts = gene_activities)

DefaultAssay(combined) <- "H3K27me3"
gene_activities_me <- GeneActivity(combined)
combined[["ACTIVITY_H3K27me3"]] <- CreateAssayObject(counts = gene_activities_me)

# Normalize gene activities
DefaultAssay(combined) <- "ACTIVITY_H3K27ac"
combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)

DefaultAssay(combined) <- "ACTIVITY_H3K27me3"
combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)

# Feature plots for known Drosophila neural markers
neural_markers <- c("ase", "D", "pros", "elav", "repo", "gem", "wg", "hh", "ptc")
neural_markers <- intersect(neural_markers, rownames(combined[["ACTIVITY_H3K27ac"]]))
if (length(neural_markers) > 0) {
  DefaultAssay(combined) <- "ACTIVITY_H3K27ac"
  p <- FeaturePlot(combined, reduction = "umap_joint", features = neural_markers,
                   ncol = 3, pt.size = 0.2, max.cutoff = "q95")
  ggsave(file.path(FIG_DIR, "gene_activity_neural_markers.png"), p, 
         width = 18, height = 4 * ceiling(length(neural_markers) / 3))
}

# ---- 11. Save ----
message("=== 11. Saving results ===")

saveRDS(combined, file.path(OUT_DIR, "combined_peak_analysis.rds"))
message("  Object saved to: ", file.path(OUT_DIR, "combined_peak_analysis.rds"))
message("  Plots saved to: ", FIG_DIR)
message("=== Done ===")