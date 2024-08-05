# [1] Packages ####
library(Seurat)
library(tidyr)
library(magrittr)

# [2] Read in data ####
cat("Reading in the samples...\n")
PBMC_TB_1 <- readRDS("home/alicen/Projects/PBMC/")
PBMC_TB_2 <- readRDS("home/alicen/Projects/PBMC/")
PBMC_TB_3 <- readRDS("home/alicen/Projects/PBMC/")
PBMC_LTBI_1 <- readRDS("home/alicen/Projects/PBMC/")
PBMC_LTBI_2 <- readRDS("home/alicen/Projects/PBMC/")
PBMC_HC_1 <- readRDS("home/alicen/Projects/PBMC/")
PBMC_HC_2 <- readRDS("home/alicen/Projects/PBMC/")

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


# [3] Merging samples ####
cat("Merging samples...\n")
PBMC_merged <- merge(
  x = samples_list[[1]],
  y = c(
    samples_list[[2]],
    samples_list[[3]],
    samples_list[[4]],
    samples_list[[5]],
    samples_list[[6]],
    samples_list[[7]]
  ),
  add.cell.ids = c( 
    
    "PBMC_TB_1",
    "PBMC_TB_2",
    "PBMC_TB_3",
    "PBMC_LTBI_1",
    "PBMC_LTBI_2",
    "PBMC_HC_1",
    "PBMC_HC_2"

  )
)

# [4] Run reduction algorithms ####
cat("creating UMAP reduction...\n")
PBMC_merged <- ScaleData(PBMC_merged) %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.6, group.singletons = TRUE) %>%
  cat("creating TSNE reduction...")
  RunUMAP(dims = 1:30) %>%
  RunUMAP(dims = 1:30)


# [5] Save #####
cat("Saving object...\n")
saveRDS(PBMC_merged, "/PBMC_merged.rds")
cat("Complete!")
