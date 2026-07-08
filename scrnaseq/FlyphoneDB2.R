library(FlyPhone)
library(Seurat)

mainDir <- "/data/ebaird/scRNAseq/2025_2026_int/"
objectDir <- paste0(mainDir, "QC_clustering/integrated.rds")
base_output_dir <- paste0(mainDir, "FlyPhone/", format(Sys.time(), "%Y%m%d_%H%M%S"))
dir.create(base_output_dir, recursive = TRUE, showWarnings = FALSE)

seu <- readRDS(objectDir)

seu$cluster <- seu$seurat_clusters

seu$cluster <- paste0("c", seu$cluster)

head(seu$cluster)

saveRDS(seu, objectDir)

RunFlyPhone(
    knowledgebase_version = "Version 1",
    seuratObject = objectDir,
    base_output_dir = base_output_dir#,
    # control_name = "gal",
    # mutant_name = "flp"
)