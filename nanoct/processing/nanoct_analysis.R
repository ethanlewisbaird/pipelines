#!/usr/bin/env Rscript
# nanoCT R Analysis Pipeline
# Unified script for Signac/Seurat analysis
#
# Usage:
#   Rscript nanoct_analysis.R --mode peak_analysis
#   Rscript nanoct_analysis.R --mode wnn
#   Rscript nanoct_analysis.R --mode svd_tune
#
# Environment variables:
#   NANOCT_DATA_DIR - data directory
#   NANOCT_MARK - histone mark (H3K27ac or H3K27me3)
#   NANOCT_OUTPUT_DIR - output directory

library(Signac)
library(Seurat)
library(GenomicRanges)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(future)
library(ggplot2)
library(patchwork)
library(dplyr)

set.seed(42)
options(future.globals.maxSize = 15 * 1024^3)
plan("multicore", workers = 8)

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
mode <- "peak_analysis"
for (i in seq_along(args)) {
  if (args[i] == "--mode" && i < length(args)) mode <- args[i + 1]
}

# Configuration
data_dir <- Sys.getenv("NANOCT_DATA_DIR", "/data/ebaird/scentinel/nanoCT/20260522.nanoCT")
output_dir <- Sys.getenv("NANOCT_OUTPUT_DIR", file.path(data_dir, "analysis_R_output"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load dm6 annotations
message("Loading dm6 annotations...")
gtf_file <- file.path(data_dir, "R_analysis_peaks/dm6_genes.gtf.gz")
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

# Load data
marks <- c("H3K27ac", "H3K27me3")
base_dir <- file.path(data_dir, "SU.analysis.2026.05.22/Vasso_nanoCT_nanoscope")

load_mark <- function(mark) {
  message("  Loading ", mark, "...")
  path <- file.path(base_dir, mark)
  frag_file <- file.path(path, "fragments.tsv.gz")
  peak_file <- file.path(path, "peaks.bed")
  
  peaks <- read.table(peak_file, col.names = c("chr", "start", "end"))
  peaks_gr <- makeGRangesFromDataFrame(peaks)
  
  frags <- CreateFragmentObject(path = frag_file, cells = Cells(readRDS(file.path(path, "seurat.rds"))))
  counts <- FeatureMatrix(fragments = frags, features = peaks_gr, verbose = FALSE)
  
  assay <- CreateChromatinAssay(counts = counts, sep = c(":", "-"), 
                                 fragments = frags, annotation = annotations)
  obj <- CreateSeuratObject(counts = assay, assay = "peaks")
  obj$mark <- mark
  
  # QC
  obj <- TSSEnrichment(obj, tss.positions = tss_positions, fast = TRUE)
  obj <- NucleosomeSignal(obj)
  
  return(obj)
}

# Main analysis based on mode
if (mode == "peak_analysis") {
  message("=== Peak-based Analysis ===")
  
  # Load both marks
  seu_list <- lapply(marks, load_mark)
  
  # Merge
  seu <- merge(seu_list[[1]], seu_list[[2]])
  
  # Process
  message("Processing...")
  seu <- FindTopFeatures(seu, min.cutoff = "q0")
  seu <- RunTFIDF(seu)
  seu <- RunSVD(seu)
  seu <- RunUMAP(seu, reduction = "lsi", dims = 2:30)
  seu <- FindNeighbors(seu, reduction = "lsi", dims = 2:30)
  seu <- FindClusters(seu, resolution = 0.8)
  
  # Save
  saveRDS(seu, file.path(output_dir, "seurat_peak_analysis.rds"))
  
  # Plots
  pdf(file.path(output_dir, "umap_peak_analysis.pdf"), width = 10, height = 8)
  DimPlot(seu, label = TRUE) + ggtitle("Peak-based Clustering")
  dev.off()
  
} else if (mode == "wnn") {
  message("=== WNN Multi-modal Integration ===")
  
  # Load both marks separately
  seu_list <- lapply(marks, load_mark)
  
  # Process each
  for (i in seq_along(seu_list)) {
    seu_list[[i]] <- FindTopFeatures(seu_list[[i]], min.cutoff = "q0")
    seu_list[[i]] <- RunTFIDF(seu_list[[i]])
    seu_list[[i]] <- RunSVD(seu_list[[i]])
  }
  
  # WNN integration
  seu <- FindMultiModalNeighbors(seu_list[[1]], seu_list[[2]], 
                                  dims.list = list(2:30, 2:30))
  seu <- RunUMAP(seu, nn.name = "weighted.nn")
  seu <- FindClusters(seu, resolution = 0.8)
  
  # Save
  saveRDS(seu, file.path(output_dir, "seurat_wnn.rds"))
  
  # Plots
  pdf(file.path(output_dir, "umap_wnn.pdf"), width = 10, height = 8)
  DimPlot(seu, label = TRUE) + ggtitle("WNN Clustering")
  dev.off()
  
} else if (mode == "svd_tune") {
  message("=== SVD Tuning ===")
  
  # Load and process
  seu <- load_mark("H3K27ac")
  seu <- FindTopFeatures(seu, min.cutoff = "q0")
  seu <- RunTFIDF(seu)
  seu <- RunSVD(seu)
  
  # Depth correlation
  depth_corr <- DepthCor(seu)
  
  # Save plot
  pdf(file.path(output_dir, "depth_correlation.pdf"), width = 8, height = 6)
  print(depth_corr)
  dev.off()
  
  message("Check depth_correlation.pdf to decide which PC to remove")
}

message("Done! Results in: ", output_dir)
