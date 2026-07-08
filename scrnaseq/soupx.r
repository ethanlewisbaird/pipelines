library(SoupX)
library(ggplot2)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
sample_names <- strsplit(args[1], ",")[[1]]  
output_dir <- args[2]
base_sample_dir <- args[3]

scs <- list()
print(sample_names)
for (sample_name in sample_names) {
  
  sample_dir <- file.path(base_sample_dir, sample_name, "outs")
  sample_output_dir <- file.path(output_dir, sample_name)
  
  if (!dir.exists(sample_output_dir)) {
    dir.create(sample_output_dir, recursive = TRUE)
  }
  
  pdf(file.path(sample_output_dir, "Rplots.pdf"))
  
  # Load the data and estimate the soup profile
  sc <- load10X(file.path(sample_dir))
  sc <- autoEstCont(sc, forceAccept = TRUE)  # Estimate the soup contamination
  out <- adjustCounts(sc)  # Adjust the counts
  print(slotNames(sc))
  
  scs[[sample_name]] <- sc  # Store soup profile for later analysis
  
  # Create Seurat objects from the original and adjusted counts
  seurat_obj_original <- CreateSeuratObject(counts = sc$toc)
  seurat_obj_adjusted <- CreateSeuratObject(counts = out)
  
  # Process original data
  seurat_obj_original <- NormalizeData(seurat_obj_original, verbose = FALSE)
  seurat_obj_original <- FindVariableFeatures(seurat_obj_original, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  seurat_obj_original <- ScaleData(seurat_obj_original, verbose = FALSE)
  seurat_obj_original <- RunPCA(seurat_obj_original, npcs = 20, verbose = FALSE)
  
  # Process adjusted data
  seurat_obj_adjusted <- NormalizeData(seurat_obj_adjusted, verbose = FALSE)
  seurat_obj_adjusted <- FindVariableFeatures(seurat_obj_adjusted, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  seurat_obj_adjusted <- ScaleData(seurat_obj_adjusted, verbose = FALSE)
  seurat_obj_adjusted <- RunPCA(seurat_obj_adjusted, npcs = 20, verbose = FALSE)
  
  pca_coords_original <- Embeddings(seurat_obj_original, "pca")[, 1:2]
  pca_coords_adjusted <- Embeddings(seurat_obj_adjusted, "pca")[, 1:2]
  
  data_original <- data.frame(
    PC1 = pca_coords_original[, 1],
    PC2 = pca_coords_original[, 2]
  )
  
  data_adjusted <- data.frame(
    PC1 = pca_coords_adjusted[, 1],
    PC2 = pca_coords_adjusted[, 2]
  )
  
  # Create scatter plots
  p_original <- ggplot(data_original, aes(x = PC1, y = PC2)) +
    geom_point() +
    labs(title = paste("PCA Plot of Original Data for", sample_name), x = "PC1", y = "PC2") +
    theme_minimal()
  
  p_adjusted <- ggplot(data_adjusted, aes(x = PC1, y = PC2)) +
    geom_point() +
    labs(title = paste("PCA Plot of Adjusted Data for", sample_name), x = "PC1", y = "PC2") +
    theme_minimal()
  
  # Save the plots as PDFs
  ggsave(file.path(sample_output_dir, paste0("PCA_Plot_of_Original_Data_", sample_name, ".pdf")), plot = p_original, device = "pdf")
  ggsave(file.path(sample_output_dir, paste0("PCA_Plot_of_Adjusted_Data_", sample_name, ".pdf")), plot = p_adjusted, device = "pdf")
  print(p_original)
  print(p_adjusted)
  
  dev.off()
}

# Integrate soupProfile calculations
res <- data.frame()
allg <- unique(unlist(lapply(names(scs), function(i) rownames(scs[[i]]$soupProfile))))
res <- data.frame(lapply(names(scs), function(i) scs[[i]]$soupProfile[allg, 'est']), row.names = allg)
colnames(res) <- paste0('SoupX_est_', names(scs))

print("Intermediate res dataframe:")
print(head(res))

res$SoupX_est_mean <- rowMeans(res)
res$SoupX_est_variance <- apply(res[, 1:length(scs)], 1, var)

print("Final res dataframe after sorting:")
print(head(res))

res <- res[order(res$SoupX_est_mean, decreasing = TRUE), ]
res$rank <- seq(1, length(allg))

# Save the top 100 ranked soup contamination genes
write.csv(head(res, 100), file = file.path(output_dir, "Top_100_SoupX_Genes.csv"), row.names = TRUE)
