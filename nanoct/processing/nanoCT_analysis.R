#!/usr/bin/env Rscript
# nanoCT Analysis Pipeline
# Full analysis script for single-cell nanoCT (Drosophila dm6)
# Uses environment variables for configuration
#
# Environment variables:
#   NANOCT_DATA_DIR - base data directory
#   NANOCT_MARK - histone mark (H3K27ac or H3K27me3)
#   NANOCT_H5AD - path to processed h5ad file
#   NANOCT_OUTPUT_DIR - output directory
#   NANOCT_RESOLUTION - clustering resolution (default: 0.8)

library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(stringr)
library(ggplot2)
library(gghalves)
library(ggpubr)
library(ComplexUpset)
library(regioneR)
library(scales)
library(ggVennDiagram)
library(BSgenome.Dmelanogaster.UCSC.dm6)

# Increase memory for future
options(future.globals.maxSize = 10 * 1024^3) # 10GB
plan("multicore", workers = 4)

# Source helper functions
script_dir <- Sys.getenv("NANOCT_SCRIPT_DIR", ".")
source(file.path(script_dir, "functions_scCT.R"))

# Read configuration from environment variables
data_dir <- Sys.getenv("NANOCT_DATA_DIR", "/data/ebaird/scentinel/nanoCT/20260522.nanoCT")
mark <- Sys.getenv("NANOCT_MARK", "H3K27ac")
h5ad_path <- Sys.getenv("NANOCT_H5AD", "")
output_dir <- Sys.getenv("NANOCT_OUTPUT_DIR", file.path(data_dir, "analysis_R_output"))
resolution <- as.numeric(Sys.getenv("NANOCT_RESOLUTION", "0.8"))

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Configuration:\n")
cat("  DATA_DIR:", data_dir, "\n")
cat("  MARK:", mark, "\n")
cat("  H5AD:", h5ad_path, "\n")
cat("  OUTPUT_DIR:", output_dir, "\n")
cat("  RESOLUTION:", resolution, "\n")

# 1. Get Annotations for dm6
message("Fetching dm6 annotations...")
options(timeout = 300)
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "dmelanogaster_gene_ensembl", 
                       host = "https://www.ensembl.org")
genes <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position", "strand", 
                 "external_gene_name", "gene_biotype", "ensembl_gene_id"), 
  mart = ensembl
)

# Convert to GRanges
annotations <- GRanges(
  seqnames = genes$chromosome_name,
  ranges = IRanges(start = genes$start_position, end = genes$end_position),
  strand = ifelse(genes$strand == 1, "+", "-"),
  gene_name = genes$external_gene_name,
  gene_biotype = genes$gene_biotype,
  gene_id = genes$ensembl_gene_id
)

cat("Loaded", length(annotations), "gene annotations\n")

# 2. Load Data
message("Loading data...")
base_dir <- file.path(data_dir, "SU.analysis.2026.05.22/Vasso_nanoCT_nanoscope")
marks <- c("H3K27ac", "H3K27me3")

# Load fragments for specified mark
fragments_file <- file.path(base_dir, mark, "fragments.tsv.gz")
if (!file.exists(fragments_file)) {
  # Try alternative path
  fragments_file <- file.path(base_dir, mark, "possorted_fragments.tsv.gz")
}

cat("Loading fragments from:", fragments_file, "\n")

# Create ChromatinAssay
counts <- Read10X(file.path(base_dir, mark, "analysis"))
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = "dm6",
  fragments = fragments_file,
  min.cells = 10,
  annotation = annotations
)

# Create Seurat object
seu <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)

cat("Created Seurat object:", ncol(seu), "cells,", nrow(seu), "peaks\n")

# 3. Quality Control
message("Running quality control...")

# Compute nucleosome signal
seu <- NucleosomeSignal(object = seu)

# Compute TSS enrichment
seu <- TSSEnrichment(object = seu, fast = TRUE)

# Add blacklist ratio
seu$blacklist_ratio <- CountsInRegion(
  object = seu, 
  assay = "peaks", 
  regions = blacklist_hg38
)

# Plot QC
pdf(file.path(output_dir, paste0("qc_", mark, ".pdf")), width = 12, height = 8)
VlnPlot(
  object = seu,
  features = c("nCount_peaks", "nFeature_peaks", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4
)
dev.off()

# Filter
seu <- subset(
  x = seu,
  subset = nCount_peaks > 1000 &
    nFeature_peaks > 500 &
    TSS.enrichment > 1 &
    nucleosome_signal < 2
)

cat("After QC filtering:", ncol(seu), "cells\n")

# 4. Dimensionality Reduction
message("Running dimensionality reduction...")

# LSI
seu <- FindTopFeatures(seu, min.cutoff = "q0")
seu <- RunTFIDF(seu)
seu <- RunSVD(seu)

# UMAP
seu <- RunUMAP(object = seu, reduction = "lsi", dims = 2:30)
seu <- FindNeighbors(object = seu, reduction = "lsi", dims = 2:30)
seu <- FindClusters(object = seu, resolution = resolution)

cat("Found", length(levels(Idents(seu))), "clusters\n")

# 5. Generate Plots
message("Generating plots...")

# UMAP
pdf(file.path(output_dir, paste0("umap_", mark, ".pdf")), width = 10, height = 8)
DimPlot(seu, label = TRUE) + 
  ggtitle(paste(mark, "- Clusters (res =", resolution, ")"))
dev.off()

# TSS enrichment
pdf(file.path(output_dir, paste0("tss_enrichment_", mark, ".pdf")), width = 10, height = 8)
TSSPlot(seu) + ggtitle(paste(mark, "TSS Enrichment"))
dev.off()

# Fragment length
pdf(file.path(output_dir, paste0("fragment_length_", mark, ".pdf")), width = 10, height = 8)
FragmentHistogram(object = seu) + ggtitle(paste(mark, "Fragment Length"))
dev.off()

# 6. Save
output_rds <- file.path(output_dir, paste0("seurat_", mark, ".rds"))
cat("Saving Seurat object to:", output_rds, "\n")
saveRDS(seu, file = output_rds)

cat("\nDone! Results saved to:", output_dir, "\n")
cat("Files:\n")
list.files(output_dir, pattern = mark)
