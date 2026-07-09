#!/usr/bin/env Rscript
#===================================================================
# WNN (Weighted Nearest Neighbor) Integration for nanoCT
# Uses Seurat v5 WNN to jointly embed H3K27ac and H3K27me3
# instead of concatenated SVD. WNN learns per-cell, per-modality
# weights that can better preserve modality-specific signal.
#===================================================================

library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(future)

set.seed(42)
options(future.globals.maxSize = 15 * 1024^3)
plan("multicore", workers = 8)

OUT_DIR  <- "/data/ebaird/scRNAseq/20260522.nanoCT/R_analysis_peaks/output"
FIG_DIR  <- file.path(OUT_DIR, "figures_wnn")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

message("=== Loading data for WNN analysis ===")

BASE_DIR <- "/data/ebaird/scRNAseq/20260522.nanoCT/SU.analysis.2026.05.22/Vasso_nanoCT_nanoscope"
gtf_file <- "/data/ebaird/scRNAseq/20260522.nanoCT/R_analysis_peaks/dm6_genes.gtf.gz"

annotations_raw <- rtracklayer::import(gtf_file)
tx <- annotations_raw[annotations_raw$type == "transcript"]
tx_by_gene <- split(tx, tx$gene_id)
gene_ranges <- unlist(reduce(range(tx_by_gene)))
gene_ranges$gene_id <- names(gene_ranges)
gn <- unique(mcols(tx)[, c("gene_id", "gene_name")])
m <- match(gene_ranges$gene_id, gn$gene_id)
gene_ranges$gene_name <- gn$gene_name[m]
annotations <- gene_ranges
annotations$gene_biotype <- "protein_coding"
seqlevels(annotations) <- gsub("^chr", "", seqlevels(annotations))
seqlevels_to_keep <- c("2L", "2R", "3L", "3R", "4", "X", "Y")
annotations <- keepSeqlevels(annotations,
                              intersect(seqlevels(annotations), seqlevels_to_keep),
                              pruning.mode = "coarse")
tss_positions <- resize(annotations, width = 1, fix = "start")

marks <- c("H3K27ac", "H3K27me3")
seurat_list <- list()

for (mark in marks) {
  message("  Loading ", mark, "...")
  path <- file.path(BASE_DIR, mark)
  frag_file  <- file.path(path, "fragments.tsv.gz")
  peak_file  <- file.path(path, "peaks.bed")
  meta_file  <- file.path(path, "singlecell.csv")
  
  peaks <- read.table(peak_file, col.names = c("chr", "start", "end"), comment.char = "#")
  peaks_gr <- makeGRangesFromDataFrame(peaks)
  metadata <- read.csv(meta_file, header = TRUE, row.names = 1)
  metadata <- metadata[rownames(metadata) != "NO_BARCODE", ]
  cells_pass <- rownames(metadata)[metadata$passed_filters > 500]
  
  frags <- CreateFragmentObject(path = frag_file, cells = cells_pass)
  counts <- FeatureMatrix(fragments = frags, features = peaks_gr, cells = cells_pass, verbose = FALSE)
  
  assay <- CreateChromatinAssay(counts = counts, sep = c(":", "-"), fragments = frags, annotation = annotations)
  obj <- CreateSeuratObject(counts = assay, assay = "peaks")
  obj$mark <- mark
  
  qc_cols <- intersect(c("total", "duplicate", "passed_filters", "is__cell_barcode",
                          "peak_region_fragments", "TSS_fragments"),
                       colnames(metadata))
  for (col in qc_cols) {
    obj <- AddMetaData(obj, metadata[Cells(obj), col], col.name = col)
  }
  obj <- TSSEnrichment(obj, tss.positions = tss_positions, fast = TRUE)
  obj <- NucleosomeSignal(obj)
  obj$FRiP <- obj$peak_region_fragments / obj$passed_filters
  
  ncounts <- obj$nCount_peaks
  lo_counts <- quantile(ncounts, 0.02)
  hi_counts <- quantile(ncounts, 0.98)
  cells_keep <- colnames(obj)[obj$nCount_peaks >= lo_counts & obj$nCount_peaks <= hi_counts]
  obj <- subset(obj, cells = cells_keep)
  
  seurat_list[[mark]] <- obj
  message("    Cells: ", ncol(obj))
}

common_cells <- Reduce(intersect, lapply(seurat_list, Cells))
message("\n  Common cells: ", length(common_cells))

# Build multi-assay object
combined <- seurat_list[["H3K27ac"]][, common_cells]
combined <- RenameAssays(combined, peaks = "H3K27ac")
me_subset <- subset(seurat_list[["H3K27me3"]], cells = common_cells)
combined[["H3K27me3"]] <- me_subset[["peaks"]]
Annotation(combined[["H3K27ac"]]) <- annotations
Annotation(combined[["H3K27me3"]]) <- annotations

# Normalize and find variable features for each assay independently
DefaultAssay(combined) <- "H3K27ac"
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = "q5")

DefaultAssay(combined) <- "H3K27me3"
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = "q5")

# Run SVD on each assay independently
DefaultAssay(combined) <- "H3K27ac"
combined <- RunSVD(combined, n = 50)
combined <- RunUMAP(combined, dims = 2:30, reduction.name = "umap_H3K27ac_wnn",
                    reduction.key = "umapH3K27ac_")

DefaultAssay(combined) <- "H3K27me3"
combined <- RunSVD(combined, n = 30)
combined <- RunUMAP(combined, dims = 2:15, reduction.name = "umap_H3K27me3_wnn",
                    reduction.key = "umapH3K27me3_")

# WNN: compute multimodal neighbors
# k = number of neighbors from each modality to use
combined <- FindMultiModalNeighbors(
  combined,
  reduction.list = list("svd_H3K27ac", "svd_H3K27me3"),
  dims.list = list(2:30, 2:15),
  k.nn = 20,
  knn.graph.name = "wknn",
  weighted.graph.name = "wknnw",
  verbose = TRUE
)

# WNN UMAP
combined <- RunUMAP(combined, nn.name = "wknn", reduction.name = "umap_wnn",
                    reduction.key = "umapWnn_")

# Per-modality weighted UMAP
combined <- RunUMAP(combined, nn.name = "wknn", reduction.name = "umap_wnn_weighted",
                    reduction.key = "umapWnnW_", weight.by.modality = TRUE)

# Clustering using WNN graph at multiple resolutions
resolutions <- c(0.5, 0.8, 1.0, 1.2, 1.5, 2.0)

message("\n=== WNN Multi-resolution clustering ===")

wnn_results <- data.frame()

for (res in resolutions) {
  message("  Resolution: ", res)
  
  combined <- FindClusters(combined, graph.name = "wknnw", resolution = res,
                           algorithm = 1, verbose = FALSE)
  
  n_clust <- length(unique(combined$seurat_clusters))
  clust_sizes <- table(combined$seurat_clusters)
  
  message("    Clusters: ", n_clust)
  message("    Sizes: ", paste(clust_sizes, collapse = ", "))
  
  col_name <- paste0("wnn_res_", gsub("\\.", "p", res))
  combined[[col_name]] <- combined$seurat_clusters
  
  wnn_results <- rbind(wnn_results, data.frame(
    resolution = res, n_clusters = n_clust,
    largest_pct = max(clust_sizes) / sum(clust_sizes) * 100,
    smallest_size = min(clust_sizes)
  ))
  
  # UMAP colored by clusters
  p <- DimPlot(combined, reduction = "umap_wnn", group.by = col_name,
               label = TRUE, pt.size = 0.3) +
    ggtitle(paste0("WNN res=", res, " (", n_clust, " clusters)")) +
    NoLegend()
  ggsave(file.path(FIG_DIR, paste0("umap_wnn_res", gsub("\\.", "p", res), ".png")),
         p, width = 8, height = 7)
}

# Also plot weighted UMAP at best resolution
p_weighted <- DimPlot(combined, reduction = "umap_wnn_weighted",
                      group.by = "seurat_clusters", label = TRUE, pt.size = 0.3) +
  ggtitle(paste0("WNN Weighted UMAP (res=", tail(resolutions, 1), ")"))
ggsave(file.path(FIG_DIR, "umap_wnn_weighted.png"), p_weighted, width = 8, height = 7)

# Per-mark UMAPS colored by WNN clusters
p_ac <- DimPlot(combined, reduction = "umap_H3K27ac_wnn", group.by = "seurat_clusters",
                label = TRUE, pt.size = 0.3) + ggtitle("H3K27ac UMAP (WNN clusters)")
p_me <- DimPlot(combined, reduction = "umap_H3K27me3_wnn", group.by = "seurat_clusters",
                label = TRUE, pt.size = 0.3) + ggtitle("H3K27me3 UMAP (WNN clusters)")
ggsave(file.path(FIG_DIR, "umap_per_mark_wnn.png"), wrap_plots(p_ac, p_me, ncol = 2),
       width = 16, height = 7)

# WNN weight visualization
p_w <- FeaturePlot(combined, reduction = "umap_wnn",
                   features = c("weighted.nn_H3K27ac", "weighted.nn_H3K27me3"),
                   ncol = 2, pt.size = 0.2)
ggsave(file.path(FIG_DIR, "wnn_weights.png"), p_w, width = 14, height = 6)

message("\n=== WNN Resolution Summary ===")
print(wnn_results)
write.csv(wnn_results, file.path(OUT_DIR, "wnn_resolution_summary.csv"), row.names = FALSE)

# Save final WNN object
saveRDS(combined, file.path(OUT_DIR, "combined_wnn_analysis.rds"))

message("\n=== Done ===")
