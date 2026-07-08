# Subset Reanalysis Pipeline
# Uses environment variables for configuration
# Set these before running:
#   DATA_DIR - main data directory
#   SEU_FILE - path to Seurat object (relative to DATA_DIR)
#   SUBSET_NAME - name for this subset analysis
#   CLUSTERS - comma-separated cluster IDs to subset
#   OUT_DIR - output directory (optional, defaults to DATA_DIR/subset_reanalysis/TIMESTAMP)

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(harmony)

# Read configuration from environment variables
mainDir <- Sys.getenv("DATA_DIR", "/data/ebaird/scentinel/scRNAseq/analysis/SCENTINELsep24")
seu_file <- Sys.getenv("SEU_FILE", "QC_clustering/20250624/merged_clusters.rds")
subset_name <- Sys.getenv("SUBSET_NAME", "subset_analysis")
clusters_str <- Sys.getenv("CLUSTERS", "12,11,1,3")
out_dir <- Sys.getenv("OUT_DIR", "")

# Parse clusters
clusters <- as.numeric(strsplit(clusters_str, ",")[[1]])

# Set up output directory
if (out_dir == "") {
  repDir <- paste0(mainDir, "/subset_reanalysis/", format(Sys.time(), "%Y%m%d_%H%M%S"), "/")
} else {
  repDir <- out_dir
}
dir.create(repDir, recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(repDir, "figs/"), showWarnings = FALSE)
dir.create(paste0(repDir, "tables/"), showWarnings = FALSE)

cat("Configuration:\n")
cat("  DATA_DIR:", mainDir, "\n")
cat("  SEU_FILE:", seu_file, "\n")
cat("  SUBSET_NAME:", subset_name, "\n")
cat("  CLUSTERS:", clusters, "\n")
cat("  OUT_DIR:", repDir, "\n")

# Load Seurat object
cat("\nLoading Seurat object...\n")
seu <- readRDS(file = paste0(mainDir, "/", seu_file))
cat("Loaded:", ncol(seu), "cells,", nrow(seu), "genes\n")

# Subset clusters
cat("Subsetting clusters:", clusters, "\n")
seu_sub <- subset(seu, merged_clusters %in% clusters)
cat("Subset:", ncol(seu_sub), "cells\n")

# SCTransform
cat("Running SCTransform...\n")
seu_sub <- SCTransform(seu_sub, vars.to.regress = "percent.mt", verbose = FALSE)

# PCA
cat("Running PCA...\n")
seu_sub <- RunPCA(seu_sub, npcs = 30)

# Elbow plot
pdf(paste0(repDir, "figs/elbow_", subset_name, "subset.pdf"))
ElbowPlot(seu_sub, ndims = 30)
dev.off()

# Dims sweep
cat("Running dims sweep...\n")
dims_sweep <- lapply(seq(5, 25, by = 5), function(d) {
  seu_sub <- FindNeighbors(seu_sub, dims = 1:d)
  seu_sub <- FindClusters(seu_sub, resolution = 1.0)
  seu_sub <- RunUMAP(seu_sub, dims = 1:d)
  DimPlot(seu_sub, label = TRUE) + ggtitle(paste("dims =", d))
})
pdf(paste0(repDir, "figs/UMAP_dims_sweep.pdf"), width = 15, height = 10)
wrap_plots(dims_sweep, ncol = 3)
dev.off()

# Resolution sweep
cat("Running resolution sweep...\n")
res_sweep <- lapply(seq(0.4, 1.6, by = 0.2), function(r) {
  seu_sub <- FindClusters(seu_sub, resolution = r)
  DimPlot(seu_sub, label = TRUE) + ggtitle(paste("resolution =", r))
})
pdf(paste0(repDir, "figs/UMAP_resolution_sweep.pdf"), width = 15, height = 10)
wrap_plots(res_sweep, ncol = 3)
dev.off()

# Final clustering (default: dims=16, res=0.6)
cat("Final clustering (dims=16, res=0.6)...\n")
seu_sub <- FindNeighbors(seu_sub, dims = 1:16)
seu_sub <- FindClusters(seu_sub, resolution = 0.6)
seu_sub <- RunUMAP(seu_sub, dims = 1:16)

# UMAP plots
pdf(paste0(repDir, "figs/UMAP.", subset_name, ".pdf"), width = 10, height = 8)
DimPlot(seu_sub, label = TRUE) + ggtitle(paste(subset_name, "- New Clusters"))
dev.off()

pdf(paste0(repDir, "figs/UMAP.", subset_name, "_split.pdf"), width = 15, height = 8)
DimPlot(seu_sub, split.by = "merged_clusters", label = TRUE) + ggtitle(paste(subset_name, "- Split by Original"))
dev.off()

# Save subset object
cat("Saving subset object...\n")
saveRDS(seu_sub, file = paste0(repDir, subset_name, "_subset.rds"))

# Find markers
cat("Finding markers...\n")
Idents(seu_sub) <- "seurat_clusters"
all_markers <- FindAllMarkers(seu_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save markers
write.csv(all_markers, paste0(repDir, "tables/allMarkers_", subset_name, ".csv"), row.names = FALSE)

# Top 10 markers
top10 <- all_markers %>% group_by(cluster) %>% slice_max(n = 10, order_by = avg_log2FC)
write.csv(top10, paste0(repDir, "tables/top10Markers_", subset_name, ".csv"), row.names = FALSE)

# Dot plot
pdf(paste0(repDir, "figs/top10markers.dotplot_", subset_name, ".pdf"), width = 15, height = 10)
DotPlot(seu_sub, features = unique(top10$gene)) + RotatedAxis()
dev.off()

# Add new clusters back to full object
cat("Adding new clusters to full object...\n")
seu$new_clusters <- seu$merged_clusters
seu$new_clusters[Cells(seu_sub)] <- paste0("sub_", seu_sub$seurat_clusters)

# Full UMAP with new clusters
pdf(paste0(repDir, "figs/", subset_name, "_recluster_full_umap.pdf"), width = 12, height = 10)
DimPlot(seu, group.by = "new_clusters", label = TRUE) + ggtitle("Full Object with Subset Reclustering")
dev.off()

cat("\nDone! Results saved to:", repDir, "\n")
cat("Files:\n")
list.files(repDir, recursive = TRUE)
