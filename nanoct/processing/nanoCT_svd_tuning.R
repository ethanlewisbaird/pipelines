#!/usr/bin/env Rscript
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
FIG_DIR  <- file.path(OUT_DIR, "figures_svdtune")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
BASE_DIR <- "/data/ebaird/scRNAseq/20260522.nanoCT/SU.analysis.2026.05.22/Vasso_nanoCT_nanoscope"
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
annotations <- keepSeqlevels(annotations, c("2L","2R","3L","3R","4","X","Y"), pruning.mode = "coarse")
tss_positions <- resize(annotations, width = 1, fix = "start")
marks <- c("H3K27ac", "H3K27me3")
seurat_list <- list()
for (mark in marks) {
  path <- file.path(BASE_DIR, mark)
  peaks <- read.table(file.path(path, "peaks.bed"), col.names = c("chr","start","end"), comment.char = "#")
  peaks_gr <- makeGRangesFromDataFrame(peaks)
  metadata <- read.csv(file.path(path, "singlecell.csv"), header = TRUE, row.names = 1)
  metadata <- metadata[rownames(metadata) != "NO_BARCODE", ]
  cells_pass <- rownames(metadata)[metadata$passed_filters > 500]
  frags <- CreateFragmentObject(path = file.path(path, "fragments.tsv.gz"), cells = cells_pass)
  counts <- FeatureMatrix(fragments = frags, features = peaks_gr, cells = cells_pass, verbose = FALSE)
  assay_obj <- CreateChromatinAssay(counts, sep = c(":","-"), fragments = frags, annotation = annotations)
  obj <- CreateSeuratObject(counts = assay_obj, assay = "peaks"); obj$mark <- mark
  obj <- TSSEnrichment(obj, tss.positions = tss_positions, fast = TRUE)
  obj <- NucleosomeSignal(obj)
  lo_c <- quantile(obj$nCount_peaks, 0.02); hi_c <- quantile(obj$nCount_peaks, 0.98)
  obj <- subset(obj, cells = colnames(obj)[obj$nCount_peaks >= lo_c & obj$nCount_peaks <= hi_c])
  seurat_list[[mark]] <- obj
}
common_cells <- Reduce(intersect, lapply(seurat_list, Cells))
combined <- seurat_list$H3K27ac[, common_cells]
combined <- RenameAssays(combined, peaks = "H3K27ac")
me_subset <- subset(seurat_list$H3K27me3, cells = common_cells)
combined$H3K27me3 <- me_subset$peaks
Annotation(combined$H3K27ac) <- annotations; Annotation(combined$H3K27me3) <- annotations
DefaultAssay(combined) <- "H3K27ac"; combined <- RunTFIDF(); combined <- FindTopFeatures(combined, min.cutoff = "q5"); combined <- RunSVD(combined, n = 50)
DefaultAssay(combined) <- "H3K27me3"; combined <- RunTFIDF(); combined <- FindTopFeatures(combined, min.cutoff = "q5"); combined <- RunSVD(combined, n = 30)
p1 <- DepthCor(combined, reduction = "svd_H3K27ac") + ggtitle("H3K27ac depth corr")
p2 <- DepthCor(combined, reduction = "svd_H3K27me3") + ggtitle("H3K27me3 depth corr")
ggsave(file.path(FIG_DIR, "depth_cor.png"), wrap_plots(p1, p2), width = 12, height = 5)
configs <- list(
  list(name = "orig",    ac_dims = 2:30, me_dims = 2:15, n_ac = 29, n_me = 14),
  list(name = "wide",    ac_dims = 2:40, me_dims = 2:20, n_ac = 39, n_me = 19),
  list(name = "narrow",  ac_dims = 2:15, me_dims = 2:10, n_ac = 14, n_me = 9),
  list(name = "ac_only", ac_dims = 2:30, me_dims = NULL, n_ac = 29, n_me = 0)
)
config_results <- data.frame()
for (cfg in configs) {
  je <- cbind(Embeddings(combined, "svd_H3K27ac")[, cfg$ac_dims, drop = FALSE],
              if (!is.null(cfg$me_dims)) Embeddings(combined, "svd_H3K27me3")[, cfg$me_dims, drop = FALSE] else NULL)
  colnames(je) <- paste0("joint_", seq_len(ncol(je)))
  cr <- combined
  cr[["joint"]] <- CreateDimReducObject(embeddings = je, key = "joint_", assay = DefaultAssay(combined))
  cr <- RunUMAP(cr, reduction = "joint", dims = 1:ncol(je), reduction.name = paste0("umap_", cfg$name), n.neighbors = 30, min.dist = 0.3)
  cr <- FindNeighbors(cr, reduction = "joint", dims = 1:ncol(je), graph.name = paste0("snn_", cfg$name))
  cr <- FindClusters(cr, graph.name = paste0("snn_", cfg$name), resolution = 1.0, algorithm = 1, verbose = FALSE)
  nc <- length(unique(cr$seurat_clusters))
  message("Config ", cfg$name, ": ", nc, " clusters, ", ncol(je), " dims")
  p <- DimPlot(cr, reduction = paste0("umap_", cfg$name), group.by = "seurat_clusters", label = TRUE, pt.size = 0.3) +
    ggtitle(paste0(cfg$name, " (", nc, " cl, ", ncol(je), " dims)")) + NoLegend()
  ggsave(file.path(FIG_DIR, paste0("umap_svdtune_", cfg$name, ".png")), p, width = 8, height = 7)
  config_results <- rbind(config_results, data.frame(config = cfg$name, n_dims = ncol(je), n_clusters = nc))
}
write.csv(config_results, file.path(OUT_DIR, "svd_config_summary.csv"), row.names = FALSE)
saveRDS(combined, file.path(OUT_DIR, "combined_for_svdtune.rds"))
message("=== SVD tuning done ===")
