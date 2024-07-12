# Library
library(Seurat)
library(magrittr)

# 1. Read in individual objects for integration ####
cat("Reading in objects...")
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

# 2. Run integration functions #####
# Subset of genes used to integrate
cat("Selecting integration features...")
PBMC_features <- SelectIntegrationFeatures(sample_list, nfeatures = 2000)

# Sets of genes used to integrate
cat("Selecting integration anchors...")
PBMC_anchors <- FindIntegrationAnchors(
  object.list = sample_list,
  anchor.features = PBMC_features,
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

# 3. Reduction with UMAP ####
# Read in output if necessary
#cat("Reading in the integrated object...")
#PBMC <- readRDS("/home/alicen/Projects/PBMC/R_objects/PBMC.rds")

cat("Reducing dimensions...")
PBMC <- ScaleData(PBMC) %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.6, group.singletons = TRUE) %>%
  RunUMAP(dims = 1:30)


cat("Saving UMAP reduced object...")
saveRDS(PBMC, "/home/alicen/Projects/PBMC/R_objects/PBMC_integrated.rds")



# 4. Reduction with TSNE ####
# cat("Reading in the integrated object...\n")
# PBMC <- readRDS("/home/alicen/Projects/PBMC/R_objects/PBMC.rds")
#
#
# cat("Reducing dimensions...\n")
# PBMC <- ScaleData(PBMC) %>%
#   RunPCA() %>%
#   FindNeighbors(dims = 1:30) %>%
#   FindClusters(resolution = 0.6, group.singletons = TRUE) %>%
#   RunTSNE(dims = 1:30)
#
#
# cat("Saving data...\n")
# saveRDS(PBMC, "/home/alicen/Projects/PBMC/R_objects/PBMC_TSNE.rds")


cat("Complete!")
