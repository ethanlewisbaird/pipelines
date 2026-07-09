#!/usr/bin/env Rscript
#===================================================================
# nanoCT Signac Analysis with MACS3 Peaks
# Replaces CellRanger peaks with MACS3-called peaks, then runs
# the full pipeline: TF-IDF, SVD, joint embedding, multi-res clustering.
#===================================================================

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
PEAK_DIR <- "/data/ebaird/scRNAseq/20260522.nanoCT/R_analysis_peaks/macs3_peaks"
OUT_DIR  <- "/data/ebaird/scRNAseq/20260522.nanoCT/R_analysis_peaks/output_macs3"
FIG_DIR  <- file.path(OUT_DIR, "figures")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Gene Annotations (same as before) ----
message("=== 1. Loading dm6 gene annotations ===")
gtf_file <- "/data/ebaird/scRNAseq/20260522.nanoCT/R_analysis_peaks/dm6_genes.gtf.gz"
annotations_raw <- rtracklayer::import(gtf_file)
tx <- annotations_raw[annotations_raw$type == "transcript"]
tx_by_gene <- split(tx, tx$gene_id)
gene_ranges <- unlist(reduce(range(tx_by_gene)))
gene_ranges$gene_id <- names(gene_ranges)
gn <- unique(mcols(tx)[, c("gene_id", "gene_name")])
m_idx <- match(gene_ranges$gene_id, gn$gene_id)
gene_ranges$gene_name <- gn$gene_name[m_idx]
annotations <- gene_ranges
annotations$gene_biotype <- "protein_coding"
seqlevels(annotations) <- gsub("^chr", "", seqlevels(annotations))
annotations <- keepSeqlevels(annotations, c("2L","2R","3L","3R","4","X","Y"),
                              pruning.mode = "coarse")
tss_positions <- resize(annotations, width = 1, fix = "start")
message("  Annotations: ", length(annotations), " genes")

# ---- 2. Load MACS3 peaks ----
message("=== 2. Loading MACS3 peaks ===")

macs3_peaks_ac <- read.table(file.path(PEAK_DIR, "H3K27ac_macs3_peaks.bed"),
                              col.names = c("chr", "start", "end"))
macs3_peaks_me <- read.table(file.path(PEAK_DIR, "H3K27me3_macs3_peaks.bed"),
                              col.names = c("chr", "start", "end"))

message("  MACS3 H3K27ac peaks: ", nrow(macs3_peaks_ac))
message("  MACS3 H3K27me3 peaks: ", nrow(macs3_peaks_me))

# Compare with CellRanger
cr_peaks_ac <- read.table(file.path(BASE_DIR, "H3K27ac/peaks.bed"),
                           col.names = c("chr", "start", "end"), comment.char = "#")
cr_peaks_me <- read.table(file.path(BASE_DIR, "H3K27me3/peaks.bed"),
                           col.names = c("chr", "start", "end"), comment.char = "#")
message("  CellRanger H3K27ac peaks: ", nrow(cr_peaks_ac))
message("  CellRanger H3K27me3 peaks: ", nrow(cr_peaks_me))

# Peak width distributions
macs3_width_ac <- macs3_peaks_ac$end - macs3_peaks_ac$start
macs3_width_me <- macs3_peaks_me$end - macs3_peaks_me$start
cr_width_ac <- cr_peaks_ac$end - cr_peaks_ac$start
cr_width_me <- cr_peaks_me$end - cr_peaks_me$start

message("  MACS3 H3K27ac median width: ", median(macs3_width_ac))
message("  CellRanger H3K27ac median width: ", median(cr_width_ac))
message("  MACS3 H3K27me3 median width: ", median(macs3_width_me))
message("  CellRanger H3K27me3 median width: ", median(cr_width_me))

# Peak overlap
macs3_gr_ac <- makeGRangesFromDataFrame(macs3_peaks_ac)
macs3_gr_me <- makeGRangesFromDataFrame(macs3_peaks_me)
cr_gr_ac <- makeGRangesFromDataFrame(cr_peaks_ac)
cr_gr_me <- makeGRangesFromDataFrame(cr_peaks_me)

overlap_ac <- sum(countOverlaps(cr_gr_ac, macs3_gr_ac, maxgap = 100) > 0)
overlap_me <- sum(countOverlaps(cr_gr_me, macs3_gr_me, maxgap = 100) > 0)
message("  H3K27ac overlap (CR peaks within 100bp of MACS3): ", overlap_ac, "/", nrow(cr_peaks_ac))
message("  H3K27me3 overlap (CR peaks within 100bp of MACS3): ", overlap_me, "/", nrow(cr_peaks_me))

# Save peak comparison
peak_comparison <- data.frame(
    mark = c("H3K27ac", "H3K27me3"),
    cellranger_peaks = c(nrow(cr_peaks_ac), nrow(cr_peaks_me)),
    macs3_peaks = c(nrow(macs3_peaks_ac), nrow(macs3_peaks_me)),
    cr_median_width = c(median(cr_width_ac), median(cr_width_me)),
    macs3_median_width = c(median(macs3_width_ac), median(macs3_width_me)),
    overlap_100bp = c(overlap_ac, overlap_me)
)
write.csv(peak_comparison, file.path(OUT_DIR, "peak_comparison.csv"), row.names = FALSE)

# ---- 3. Load data with MACS3 peaks ----
message("=== 3. Building count matrices with MACS3 peaks ===")

marks <- c("H3K27ac", "H3K27me3")
peak_files <- list(H3K27ac = macs3_peaks_ac, H3K27me3 = macs3_peaks_me)
seurat_list <- list()

for (mark in marks) {
    message("  Loading ", mark, "...")
    path <- file.path(BASE_DIR, mark)
    frag_file <- file.path(path, "fragments.tsv.gz")
    meta_file <- file.path(path, "singlecell.csv")
    
    peaks_df <- peak_files[[mark]]
    peaks_gr <- makeGRangesFromDataFrame(peaks_df)
    
    metadata <- read.csv(meta_file, header = TRUE, row.names = 1)
    metadata <- metadata[rownames(metadata) != "NO_BARCODE", ]
    cells_pass <- rownames(metadata)[metadata$passed_filters > 500]
    message("    Cells with >500 passing fragments: ", length(cells_pass))
    
    frags <- CreateFragmentObject(path = frag_file, cells = cells_pass)
    counts <- FeatureMatrix(fragments = frags, features = peaks_gr,
                            cells = cells_pass, verbose = FALSE)
    message("    Count matrix: ", nrow(counts), " peaks x ", ncol(counts), " cells")
    
    assay_obj <- CreateChromatinAssay(
        counts = counts, sep = c(":", "-"),
        fragments = frags, annotation = annotations
    )
    obj <- CreateSeuratObject(counts = assay_obj, assay = "peaks")
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
    
    # Lenient QC filter
    ncounts <- obj$nCount_peaks
    lo_counts <- quantile(ncounts, 0.02)
    hi_counts <- quantile(ncounts, 0.98)
    cells_keep <- colnames(obj)[obj$nCount_peaks >= lo_counts & obj$nCount_peaks <= hi_counts]
    obj <- subset(obj, cells = cells_keep)
    message("    Cells after filtering: ", ncol(obj))
    
    seurat_list[[mark]] <- obj
}

# ---- 4. Multi-modal Integration ----
message("=== 4. Multi-modal integration ===")
common_cells <- Reduce(intersect, lapply(seurat_list, Cells))
message("  Common cells: ", length(common_cells))

combined <- seurat_list$H3K27ac[, common_cells]
combined <- RenameAssays(combined, peaks = "H3K27ac")
me_subset <- subset(seurat_list$H3K27me3, cells = common_cells)
combined$H3K27me3 <- me_subset$peaks
Annotation(combined$H3K27ac) <- annotations
Annotation(combined$H3K27me3) <- annotations

message("  Combined: ", ncol(combined), " cells")
message("    H3K27ac: ", nrow(combined$H3K27ac), " peaks")
message("    H3K27me3: ", nrow(combined$H3K27me3), " peaks")

# ---- 5. TF-IDF + SVD ----
message("=== 5. TF-IDF normalization and SVD ===")
DefaultAssay(combined) <- "H3K27ac"
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = "q5")
combined <- RunSVD(combined, n = 50, reduction.name = "svd_H3K27ac")

DefaultAssay(combined) <- "H3K27me3"
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = "q5")
combined <- RunSVD(combined, n = 30, reduction.name = "svd_H3K27me3")

# DepthCor
p1 <- DepthCor(combined, reduction = "svd_H3K27ac") + ggtitle("H3K27ac SVD depth corr (MACS3)")
p2 <- DepthCor(combined, reduction = "svd_H3K27me3") + ggtitle("H3K27me3 SVD depth corr (MACS3)")
ggsave(file.path(FIG_DIR, "SVD_depth_correlation.png"), wrap_plots(p1, p2), width = 12, height = 5)

# ---- 6. Joint Embedding ----
message("=== 6. Joint embedding ===")
joint_embedding <- cbind(
    Embeddings(combined, "svd_H3K27ac")[, 2:30],
    Embeddings(combined, "svd_H3K27me3")[, 2:15]
)
colnames(joint_embedding) <- paste0("joint_", seq_len(ncol(joint_embedding)))
combined[["joint"]] <- CreateDimReducObject(
    embeddings = joint_embedding, key = "joint_", assay = DefaultAssay(combined)
)
combined <- RunUMAP(combined, reduction = "joint", dims = 1:ncol(joint_embedding),
                    reduction.name = "umap_joint", n.neighbors = 30, min.dist = 0.3)
combined <- FindNeighbors(combined, reduction = "joint", dims = 1:ncol(joint_embedding),
                          graph.name = "joint_snn")

# ---- 7. Multi-resolution clustering ----
message("=== 7. Multi-resolution clustering ===")
resolutions <- c(0.5, 0.8, 1.0, 1.2, 1.5, 2.0)

results <- data.frame()
for (res in resolutions) {
    message("  Resolution: ", res)
    combined <- FindClusters(combined, graph.name = "joint_snn", resolution = res,
                             algorithm = 1, verbose = FALSE)
    n_clust <- length(unique(combined$seurat_clusters))
    clust_sizes <- table(combined$seurat_clusters)
    message("    Clusters: ", n_clust, " sizes: ", paste(clust_sizes, collapse = ", "))
    
    results <- rbind(results, data.frame(
        resolution = res, n_clusters = n_clust,
        largest_pct = max(clust_sizes) / sum(clust_sizes) * 100,
        smallest_size = min(clust_sizes)
    ))
    
    col_name <- paste0("res_", gsub("\\.", "p", res))
    combined[[col_name]] <- combined$seurat_clusters
    
    p <- DimPlot(combined, reduction = "umap_joint", group.by = col_name,
                 label = TRUE, pt.size = 0.3) +
        ggtitle(paste0("MACS3 res=", res, " (", n_clust, " clusters)")) + NoLegend()
    ggsave(file.path(FIG_DIR, paste0("umap_res", gsub("\\.", "p", res), ".png")),
           p, width = 8, height = 7)
}

# ---- 8. Per-mark UMAPs ----
message("=== 8. Per-mark visualizations ===")
combined <- RunUMAP(combined, reduction = "svd_H3K27ac", dims = 2:30,
                    reduction.name = "umap_H3K27ac", n.neighbors = 30)
combined <- RunUMAP(combined, reduction = "svd_H3K27me3", dims = 2:15,
                    reduction.name = "umap_H3K27me3", n.neighbors = 30)

p_ac <- DimPlot(combined, reduction = "umap_H3K27ac", group.by = "seurat_clusters",
                label = TRUE, pt.size = 0.3) + ggtitle("H3K27ac only (MACS3)")
p_me <- DimPlot(combined, reduction = "umap_H3K27me3", group.by = "seurat_clusters",
                label = TRUE, pt.size = 0.3) + ggtitle("H3K27me3 only (MACS3)")
ggsave(file.path(FIG_DIR, "umap_per_mark.png"), wrap_plots(p_ac, p_me, ncol = 2),
       width = 16, height = 7)

# ---- 9. Save ----
message("=== 9. Saving results ===")
saveRDS(combined, file.path(OUT_DIR, "combined_macs3_analysis.rds"))
write.csv(results, file.path(OUT_DIR, "resolution_summary_macs3.csv"), row.names = FALSE)

message("\n=== Resolution Summary (MACS3 peaks) ===")
print(results)
message("\n=== Done ===")
