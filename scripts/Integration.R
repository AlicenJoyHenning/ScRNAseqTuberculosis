# Library
library(Seurat)

cat("Reading in objects...")
# Read in object
PBMC_TB_1 <- readRDS("/home/alicen/Projects/PBMC/R_objects/PBMC_TB_1.rds")
PBMC_TB_2 <- readRDS("/home/alicen/Projects/PBMC/R_objects/PBMC_TB_2.rds")
PBMC_TB_3 <- readRDS("/home/alicen/Projects/PBMC/R_objects/PBMC_TB_3.rds")
PBMC_LTBI_1 <- readRDS("/home/alicen/Projects/PBMC/R_objects/PBMC_LTBI_1.rds")
PBMC_LTBI_2 <- readRDS("/home/alicen/Projects/PBMC/R_objects/PBMC_LTBI_2.rds")
PBMC_HC_1 <- readRDS("/home/alicen/Projects/PBMC/R_objects/PBMC_HC_1.rds")
PBMC_HC_2 <- readRDS("/home/alicen/Projects/PBMC/R_objects/PBMC_HC_2.rds")

# Based on marker expression from UMAP_marker output select samples
sample_list <- c(

  PBMC_TB_1$seuratObject,
  PBMC_TB_2$seuratObject,
  PBMC_TB_3$seuratObject,
  PBMC_LTBI_1$seuratObject,
  PBMC_LTBI_2$seuratObject,
  PBMC_HC_1$seuratObject,
  PBMC_HC_2$seuratObject

)

# Subset of genes used to integrate
cat("Selecting integration features...")
PBMC_features <- SelectIntegrationFeatures(sample_list, nfeatures = 2000)

# Sets of genes used to integrate
cat("Selecting integration anchors...")
PBMC_anchors <- FindIntegrationAnchors(
  object.list = sample_list,
  anchor.features = cellranger_features,
  normalization.method = "LogNormalize" ,
  dims = 1:30,
  k.anchor = 20)

cat("Integrating data...")
# Actual integration
PBMC <- IntegrateData(
  anchorset = PBMC_anchors,
  new.assay.name = "integrated",
  normalization.method = "LogNormalize"
)

cat("Saving data...")
# Save output
saveRDS(PBMC, "/home/alicen/Projects/PBMC/R_objects/PBMC.rds")

cat("Complete!")
