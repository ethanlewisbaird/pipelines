#!/usr/bin/env Rscript
### Marker gene calculation and plots for SCENTINELsep24
### Saves to date-stamped directory to avoid overwriting existing results

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(future)

options(future.globals.maxSize = 214748364800)
plan("multicore", workers = 10)

### Set directories - use today's date stamp
mainDir <- "/data/ebaird/scentinel/scRNAseq/analysis/SCENTINELsep24/"
qcDir <- paste0(mainDir, "QC_clustering/")
dateStamp <- format(Sys.Date(), "%Y%m%d")
repDir <- paste0(qcDir, dateStamp, "/")
figDir <- paste0(repDir, "figs/")
tabDir <- paste0(repDir, "tables/")

dir.create(repDir, recursive = TRUE, showWarnings = FALSE)
dir.create(figDir, recursive = TRUE, showWarnings = FALSE)
dir.create(tabDir, recursive = TRUE, showWarnings = FALSE)

### Set colours
mycols <- c(1, '#ffffe5','#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506')

cat("=== Loading Seurat object ===\n")
seu <- readRDS(paste0(qcDir, "20250624/merged_clusters.rds"))
cat("Loaded:", paste0(qcDir, "20250624/merged_clusters.rds"), "\n")
cat("Clusters:", length(levels(Idents(seu))), "\n")
cat("Cells:", ncol(seu), "\n")

### Check if marker files already exist in today's directory
markersFile <- paste0(tabDir, "allMarkers_merged_clusters.csv")
top5File <- paste0(tabDir, "top5Markers_merged_clusters.csv")
top10File <- paste0(tabDir, "top10Markers_merged_clusters.csv")

if (file.exists(markersFile)) {
  cat("WARNING: Marker file already exists at:", markersFile, "\n")
  cat("Skipping FindAllMarkers to avoid overwriting.\n")
  all.markers <- read.csv(markersFile, row.names = 1)
} else {
  cat("=== Running FindAllMarkers ===\n")
  DefaultAssay(seu) <- 'SCT'
  Idents(seu) <- 'seurat_clusters'
  
  all.markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)
  write.csv(all.markers, file = markersFile)
  cat("Saved markers to:", markersFile, "\n")
  cat("Total markers found:", nrow(all.markers), "\n")
}

### Get top markers
all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) -> top5

all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

if (!file.exists(top5File)) {
  write.csv(top5, file = top5File)
  cat("Saved top5 markers to:", top5File, "\n")
}

if (!file.exists(top10File)) {
  write.csv(top10, file = top10File)
  cat("Saved top10 markers to:", top10File, "\n")
}

### Generate plots
cat("=== Generating plots ===\n")

DefaultAssay(seu) <- 'SCT'
Idents(seu) <- 'seurat_clusters'

# Top 5 markers heatmap
heatmap5File <- paste0(figDir, "top5markers_heatmap.jpeg")
if (!file.exists(heatmap5File)) {
  jpeg(heatmap5File, quality = 100, width = 2500, height = 1500, res = 150)
  print(DoHeatmap(seu, features = unique(top5$gene)) + NoLegend() + 
    theme(axis.text.y = element_text(size = 8)))
  dev.off()
  cat("Saved:", heatmap5File, "\n")
} else {
  cat("EXISTS:", heatmap5File, "\n")
}

# Top 5 markers dotplot
dotplot5File <- paste0(figDir, "top5markers_dotplot.jpeg")
if (!file.exists(dotplot5File)) {
  jpeg(dotplot5File, quality = 100, width = 2500, height = 1000, res = 150)
  print(
    DotPlot(seu, features = unique(top5$gene), dot.scale = 6) + 
    RotatedAxis() +
    theme(axis.text.x = element_text(size = 5)) +
    scale_color_gradientn(
      colours = c("white", "forestgreen"),
      limits = c(0, 1.5),
      oob = scales::squish
    ) +
    ggtitle("Top 5 Markers per cluster")
  )
  dev.off()
  cat("Saved:", dotplot5File, "\n")
} else {
  cat("EXISTS:", dotplot5File, "\n")
}

# Top 10 markers heatmap
heatmap10File <- paste0(figDir, "top10markers_heatmap.jpeg")
if (!file.exists(heatmap10File)) {
  jpeg(heatmap10File, quality = 100, width = 3000, height = 1800, res = 150)
  print(DoHeatmap(seu, features = unique(top10$gene)) + NoLegend() + 
    theme(axis.text.y = element_text(size = 7)))
  dev.off()
  cat("Saved:", heatmap10File, "\n")
} else {
  cat("EXISTS:", heatmap10File, "\n")
}

# Top 10 markers dotplot
dotplot10File <- paste0(figDir, "top10markers_dotplot.jpeg")
if (!file.exists(dotplot10File)) {
  jpeg(dotplot10File, quality = 100, width = 3800, height = 1000, res = 150)
  print(
    DotPlot(seu, features = unique(top10$gene), dot.scale = 6) + 
    RotatedAxis() +
    theme(axis.text.x = element_text(size = 7)) +
    scale_color_gradientn(
      colours = c("white", "forestgreen"),
      limits = c(0, 1.5),
      oob = scales::squish
    ) +
    ggtitle("Top 10 Markers per cluster")
  )
  dev.off()
  cat("Saved:", dotplot10File, "\n")
} else {
  cat("EXISTS:", dotplot10File, "\n")
}

### Per-genotype dotplots
cat("=== Generating per-genotype dotplots ===\n")
genoDotFile <- paste0(figDir, "top5markers_dotplot_per_genotype.jpeg")
if (!file.exists(genoDotFile)) {
  if ("genotype" %in% colnames(seu@meta.data)) {
    unique_genotypes <- unique(seu$genotype)
    genes_to_plot <- unique(top5$gene)
    
    plots <- list()
    for (gen in unique_genotypes) {
      seu_subset <- subset(seu, subset = genotype == gen)
      p <- DotPlot(seu_subset, features = genes_to_plot, dot.scale = 6) +
        scale_color_gradientn(colours = c("white", "blue"), limits = c(0, 1.5), oob = scales::squish) +
        ggtitle(gen) +
        RotatedAxis() +
        theme(
          axis.text.x = element_text(size = 5),
          plot.title = element_text(hjust = 0.5)
        ) +
        guides(color = guide_colorbar(title = "Expr"),
               size = guide_legend(title = "Pct. Expr"))
      plots[[gen]] <- p
    }
    
    combined <- wrap_plots(plots, ncol = 1)
    
    jpeg(genoDotFile, quality = 100, width = 1800, height = 2000, res = 150)
    print(combined)
    dev.off()
    cat("Saved:", genoDotFile, "\n")
  } else {
    cat("No 'genotype' column found in metadata, skipping per-genotype dotplots\n")
  }
} else {
  cat("EXISTS:", genoDotFile, "\n")
}

cat("=== Marker gene analysis complete ===\n")
cat("Results saved to:", repDir, "\n")
cat("Tables:", tabDir, "\n")
cat("Figures:", figDir, "\n")
