# miRNA Correlation and GO Enrichment Analysis
# For miRNAs: mir-317, mir-9948 (formerly mir-4942), mir-927 (formerly mir-2279), mir-4974

# Load required libraries
library(Seurat)
library(dplyr)
library(clusterProfiler)
library(org.Dm.eg.db)
library(ggplot2)

so <- readRDS("/data/ebaird/scRNAseq/SCENTINELsep24/composition_DEG_signatures/signatures.rds")

# Function to perform correlation and GO enrichment for a miRNA
analyze_mirna_correlation_go <- function(seurat_obj, mirna_name, n_top_genes = 100) {
  cat("Analyzing", mirna_name, "...\n")
  
  # Check if miRNA exists in dataset
  if (!mirna_name %in% rownames(seurat_obj)) {
    cat("Warning:", mirna_name, "not found in dataset. Skipping.\n")
    return(NULL)
  }
  
  # Get normalized expression matrix
  expression_matrix <- LayerData(seurat_obj, assay = "RNA", layer = "data")
  mirna_expression <- expression_matrix[mirna_name, ]
  
  # Calculate correlations with all genes
  cat("Calculating correlations for", mirna_name, "...\n")
  
  # Use only variable genes for efficiency if needed
  variable_genes <- VariableFeatures(seurat_obj)
  if (length(variable_genes) > 0) {
    genes_to_test <- variable_genes
  } else {
    genes_to_test <- rownames(expression_matrix)
  }
  
  # Calculate Spearman correlations
  correlation_results <- sapply(genes_to_test, function(gene) {
    if (gene %in% rownames(expression_matrix)) {
      gene_expression <- expression_matrix[gene, ]
      cor_result <- cor(mirna_expression, gene_expression, method = "spearman", use = "complete.obs")
      return(cor_result)
    } else {
      return(NA)
    }
  })
  
  # Remove NAs and create results data frame
  correlation_results <- correlation_results[!is.na(correlation_results)]
  cor_df <- data.frame(
    gene = names(correlation_results),
    correlation = correlation_results,
    row.names = NULL
  ) %>%
    arrange(correlation)
  
  # Get top anti-correlated genes
  top_anti_cor <- head(cor_df, n_top_genes)
  
  # Perform GO enrichment on anti-correlated genes
  cat("Performing GO enrichment for", mirna_name, "...\n")
  
  # Convert gene symbols to ENTREZ IDs
  entrez_ids <- tryCatch({
    bitr(
      top_anti_cor$gene,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Dm.eg.db
    )
  }, error = function(e) {
    cat("Gene ID conversion failed for", mirna_name, ":", e$message, "\n")
    return(data.frame(SYMBOL = character(), ENTREZID = character()))
  })
  
  if (nrow(entrez_ids) > 5) {
    go_enrichment <- enrichGO(
      gene = entrez_ids$ENTREZID,
      OrgDb = org.Dm.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.1,
      readable = TRUE
    )
  } else {
    go_enrichment <- NULL
    cat("Not enough genes for pathway enrichment for", mirna_name, "\n")
  }
  
  # Return results
  results <- list(
    mirna_name = mirna_name,
    correlation_results = cor_df,
    top_anti_cor = top_anti_cor,
    go_enrichment = go_enrichment
  )
  
  return(results)
}

# Function to plot correlation results
plot_correlation_results <- function(results, seurat_obj) {
  mirna_name <- results$mirna_name
  
  # Plot miRNA expression on UMAP
  p1 <- FeaturePlot(seurat_obj, features = mirna_name) +
    ggtitle(paste(mirna_name, "expression"))
  
  # Plot distribution of correlations
  p2 <- ggplot(results$correlation_results, aes(x = correlation)) +
    geom_histogram(bins = 50, fill = "lightblue", color = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    ggtitle(paste("Correlation distribution for", mirna_name)) +
    xlab("Spearman Correlation") +
    ylab("Count")
  
  # Plot top 3 anti-correlated genes if available
  plot_list <- list(p1, p2)
  
  top_genes <- head(results$top_anti_cor$gene, 3)
  for (gene in top_genes) {
    if (gene %in% rownames(seurat_obj)) {
      p <- FeaturePlot(seurat_obj, features = gene) +
        ggtitle(paste(gene, "(anti-correlated with", mirna_name, ")"))
      plot_list <- c(plot_list, list(p))
    }
  }
  
  # Plot GO enrichment results if available
  if (!is.null(results$go_enrichment) && nrow(results$go_enrichment) > 0) {
    p_go <- dotplot(results$go_enrichment, showCategory = 10) +
      ggtitle(paste("GO Enrichment for genes anti-correlated with", mirna_name))
    plot_list <- c(plot_list, list(p_go))
  }
  
  # Arrange plots
  combined_plot <- patchwork::wrap_plots(plot_list, ncol = 2)
  print(combined_plot)
  
  return(combined_plot)
}

# Main analysis
mirnas_to_analyze <- c("mir-317", "mir-4942", "mir-2279", "mir-4974")
all_results <- list()

# Create output directory
if (!dir.exists("mirna_analysis_results")) {
  dir.create("mirna_analysis_results")
}

for (mirna in mirnas_to_analyze) {
  results <- analyze_mirna_correlation_go(so, mirna, n_top_genes = 100)
  if (!is.null(results)) {
    all_results[[mirna]] <- results
    
    # Save results to files
    write.csv(results$correlation_results, 
              file = paste0("mirna_analysis_results/", mirna, "_all_correlations.csv"),
              row.names = FALSE)
    
    write.csv(results$top_anti_cor, 
              file = paste0("mirna_analysis_results/", mirna, "_top_anti_correlated.csv"),
              row.names = FALSE)
    
    if (!is.null(results$go_enrichment)) {
      write.csv(results$go_enrichment@result, 
                file = paste0("mirna_analysis_results/", mirna, "_go_enrichment.csv"),
                row.names = FALSE)
    }
    
    # Generate and save plots
    plot_file <- paste0("mirna_analysis_results/", mirna, "_analysis_plot.png")
    png(plot_file, width = 12, height = 10, units = "in", res = 300)
    plot_correlation_results(results, so)
    dev.off()
  }
}

# Print summary
cat("\n=== ANALYSIS SUMMARY ===\n")
for (mirna in names(all_results)) {
  cat("\nFor", mirna, ":\n")
  cat("- Number of genes analyzed:", nrow(all_results[[mirna]]$correlation_results), "\n")
  cat("- Range of correlations:", 
      round(min(all_results[[mirna]]$correlation_results$correlation), 3), "to",
      round(max(all_results[[mirna]]$correlation_results$correlation), 3), "\n")
  cat("- Top anti-correlated gene:", 
      all_results[[mirna]]$top_anti_cor$gene[1], 
      "(r =", round(all_results[[mirna]]$top_anti_cor$correlation[1], 3), ")\n")
  
  if (!is.null(all_results[[mirna]]$go_enrichment)) {
    top_go <- head(all_results[[mirna]]$go_enrichment$Description, 3)
    cat("- Top GO terms:", paste(top_go, collapse = ", "), "\n")
  }
  cat("- Results saved to: mirna_analysis_results/", mirna, "_*.csv\n", sep = "")
}

# Create a summary report of top anti-correlated genes
summary_df <- data.frame()
for (mirna in names(all_results)) {
  top_genes <- head(all_results[[mirna]]$top_anti_cor, 10)
  for (i in 1:nrow(top_genes)) {
    summary_df <- rbind(summary_df, data.frame(
      miRNA = mirna,
      Gene = top_genes$gene[i],
      Correlation = round(top_genes$correlation[i], 4),
      Rank = i
    ))
  }
}

write.csv(summary_df, "mirna_analysis_results/summary_top_anti_correlated_genes.csv", row.names = FALSE)
cat("\nSummary of top anti-correlated genes saved to: mirna_analysis_results/summary_top_anti_correlated_genes.csv\n")