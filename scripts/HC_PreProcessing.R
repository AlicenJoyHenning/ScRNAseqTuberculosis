 
# 1. PreProcess ####
# Function for healthy controls manual selection of damaged cells

PreProcess <- function(
    filtered_path,
    raw_path,
    nf_csv,
    project_name
) {

  # Timing function
  start_time <- Sys.time()

  # 1. Ambient RNA correction with SoupX ####
  cat("\nBegin pre-processing for", project_name, "...\n")

  # Create soup channel (zipped input files)
  tod <- suppressWarnings(Read10X(raw_path))
  toc <- suppressWarnings(Read10X(filtered_path))
  check <- suppressWarnings(CreateSeuratObject(
    counts = toc,
    project = project_name,
    min.cells = 1))

  # Only run Soup X if more than 1500 cells

  seurat <- NULL

  if (length(Cells(check)) >= 1500) {

    # Run Soup X
    sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)

    # Estimate the soup
    sc <- estimateSoup(sc)

    # clustering (required) of filtered matrix
    seurat_soup <- suppressWarnings(CreateSeuratObject(toc, min.cells = 0))
    seurat_soup <- suppressWarnings(SCTransform(seurat_soup, verbose = FALSE) %>%
                                      RunPCA(verbose = FALSE) %>%
                                      RunUMAP(dims = 1:30, verbose = FALSE) %>%
                                      FindNeighbors(dims = 1:30, verbose = FALSE) %>%
                                      FindClusters(verbose = FALSE))

    # add clusters to the channel
    meta.data <- seurat_soup@meta.data
    umap.embedding <- seurat_soup@reductions$umap@cell.embeddings
    sc <- suppressWarnings(setClusters(sc, setNames(meta.data$seurat_clusters, rownames(meta.data))))
    sc <- suppressWarnings(setDR(sc, umap.embedding, c("UMAP_1", "UMAP_2")))

    # with defined clusters, run the main Soup X function to calculate ambient RNA profile
    sc <- autoEstCont(sc, verbose = FALSE, doPlot = FALSE)

    # Output integer matrix of soup-corrected reads (unzipped output)
    adj.matrix <- suppressWarnings(adjustCounts(sc, roundToInt = T))

    # Output results in seurat object
    seurat <- suppressWarnings(CreateSeuratObject(counts = adj.matrix,
                                                  min.cells = 1,
                                                  project = project_name))

    cat("\u2714  Soup X ambient RNA correction\n")

  } else {

    # Create seurat object straight from filtered output
    seurat <- Read10X(filtered_path)
    seurat <- suppressWarnings(CreateSeuratObject(
      counts = seurat,
      project = project_name,
      min.cells = 1))

  }

  # Print number of cells
  Initial_cell_number <- length(Cells(seurat))
  cat("\tInitial cell number:", Initial_cell_number, "\n")

  # 2. Add nuclear fraction data ####

  # Read in velocyto nuclear fraction csv (NB before filtering cells to merge correctly)
  nf_csv <- read.csv2(nf_csv,
                      sep = ",",
                      col.names = c("barcodes", "nf"))

  # add to seurat object (ensure row order correct)
  seurat@meta.data[['nf']] <- nf_csv$nf[match(rownames(seurat@meta.data), nf_csv$barcodes)]

  cat("\u2714  Velocyto nuclear fraction output added\n")

  # 3. Isolating immune cells ####
  hemo_gene <- c("HBA1", "HBA2", "HBB")
  cd45_gene <- c("PTPRC")
  house_keeping_gene_1 <- c("MT-ND2")
  house_keeping_gene_2 <- c("GAPDH")

  seurat[['ptprc.percent']] <- PercentageFeatureSet(
    object = seurat,
    features = intersect(cd45_gene, rownames(seurat@assays$RNA)),
    assay = "RNA"
  )

  seurat[['house.keeping.1']] <- PercentageFeatureSet(
    object = seurat,
    features = intersect(house_keeping_gene_1, rownames(seurat@assays$RNA)),
    assay = "RNA"
  )

  seurat[['house.keeping.2']] <- PercentageFeatureSet(
    object = seurat,
    features = intersect(house_keeping_gene_2, rownames(seurat@assays$RNA)),
    assay = "RNA"
  )

  # Add metadata column
  seurat[['hbb.percent']] <- PercentageFeatureSet(
    object = seurat,
    features = intersect(hemo_gene, rownames(seurat@assays$RNA)),
    assay = "RNA"
  )

  seurat[['IMC']] <- ifelse(seurat$hbb.percent >=  1 |
                              seurat$ptprc.percent == 0,
                            "non_IMC",
                            "IMC")

  # Print number of cells
  Initial_cell_number <- length(Cells(seurat))
  cat("\tCell number:", Initial_cell_number, "\n")

  # Plot the hemoglobin percent
  df <- as.data.frame(seurat@meta.data)

  # Plot hemo vs CD45
  hemo_feature_QC <- ggplot(df, aes(x = hbb.percent, y = ptprc.percent, color = IMC)) +
    geom_point(size = 0.25) +
    scale_color_manual(values = c("grey","darkred")) +
    xlab("Hemoglobin expression") +
    ylab("CD45 expression") +
    labs(caption = paste(project_name)) +
    theme_classic() +
    theme(
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.4), 'cm'),
      plot.caption = element_text(hjust = 0.5, face = "bold", size = 8),
      axis.title.x = element_text(hjust = 0.5, vjust = -1, face = "bold", size = 8),
      axis.title.y = element_text(hjust = 0.5, vjust = 3, face = "bold", size = 8),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = -1.2, size = 8),
      axis.text.y = element_text(angle = 0, hjust = 1, vjust = -1.2, size = 8),
      legend.title = element_text(face = "bold"),
      legend.position = "below",
      legend.box.background = element_rect(colour = "black"),
      strip.background = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(colour = "black")
    )

  # Plot housekeeping genes
  mt_feature_QC <- ggplot(df, aes(x = house.keeping.1, y = house.keeping.2, color = IMC)) +
    geom_point(size = 0.25) +
    scale_color_manual(values = c("grey","#324776")) +
    xlab("MT-ND1 expression") +
    ylab("GAPDH expression") +
    labs(caption = paste(project_name)) +
    theme_classic() +
    theme(
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.4), 'cm'),
      plot.caption = element_text(hjust = 0.5, face = "bold", size = 8),
      axis.title.x = element_text(hjust = 0.5, vjust = -1, face = "bold", size = 8),
      axis.title.y = element_text(hjust = 0.5, vjust = 3, face = "bold", size = 8),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = -1.2, size = 8),
      axis.text.y = element_text(angle = 0, hjust = 1, vjust = -1.2, size = 8),
      legend.title = element_text(face = "bold"),
      legend.position = "below",
      legend.box.background = element_rect(colour = "black"),
      strip.background = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(colour = "black")
    )

  feature_QC <- ( hemo_feature_QC | mt_feature_QC )

  # Calculate the number of red blood cells & immune cells
  # Account for if no hemoglobin expression
  RBC_number <- NULL

  if (all(seurat$hbb.percent == 0)) {

    # Say there were no RBCs
    RBC_number = 0
    cat("\u2714  No red blood cells found\n")

  } else {

    # Count RBCs if they are present
    RBC <- subset(seurat, hbb.percent > 0)
    RBC_number <- length(Cells(RBC))
    cat("\u2714 ", RBC_number, "Red blood cells removed\n")

  }

  # Filtering step after visualization, this is the seurat object we want
  seurat <- subset(seurat, IMC == "IMC")
  IMC_number <- length(Cells(seurat))
  cat("\u2714 ", IMC_number, "Immune cells isolated\n")

  # 4. Identify damaged cells with DropletQC #####

  # Identify empty droplets
    ed_df <- data.frame(nf = as.numeric(seurat$nf),
                        umi = seurat$nCount_RNA)

    ed_results <- suppressMessages(identify_empty_drops(ed_df))

    # Identify damaged cells
    # Create vector of length same as cell number (runs default with no adjustment on groups of cells, need to fill this column requirement while getting results for that sample as a whole)
    n_cells <- length(Cells(seurat))
    cell_type <- rep(1, n_cells)

    dc_df <- data.frame(
      nf = as.numeric(seurat$nf),
      umi = seurat$nCount_RNA,
      cell_status = ed_results$cell_status,
      cell_type = cell_type
    )

    # Open a connection to a temporary file for writing
    tmp_conn <- file(tempfile(), open = "wt")

    # Redirect standard output and messages to the temporary file
    sink(tmp_conn)
    sink(tmp_conn, type = "message")

    # Your function call
    dc_results <- identify_damaged_cells(dc_df)

    # Reset output redirection
    sink(NULL)
    sink(NULL, type = "message")

    # Close the connection
    close(tmp_conn)

    # Add output to seurat object metadata
    seurat$DQC_cell_status <- dc_results$df$cell_status
    seurat$DQC_cell_status <- ifelse(seurat$DQC_cell_status == "empty_droplet" & seurat$nf >= 0.02, "cell", seurat$DQC_cell_status)
    seurat$DQC_cell_status <- ifelse(seurat$DQC_cell_status == "empty_droplet" & seurat$nCount_RNA >= 2000, "cell", seurat$DQC_cell_status)

    cat("\u2714  Droplet QC predictions\n")

    # 5. Identify damaged cells with Mt-Rb ####
    # Get gene annotations for MT & RB genes
    mt_gene_annotations <- annotations[grep("MT-", annotations$gene_name, perl=TRUE),]
    mt_gene_annotations <- mt_gene_annotations[grepl("protein_coding", mt_gene_annotations$gene_biotype, perl=TRUE),]
    mt_genes <- mt_gene_annotations %>% pull(gene_name)

    # isolate ribosomal genes
    rps_genes <- annotations[grep("RPS", annotations$gene_name, perl=TRUE),]
    rps_genes <- rps_genes[grepl("protein_coding", rps_genes$gene_biotype, perl=TRUE),]
    rps_genes <- rps_genes %>% pull(gene_name)
    rpl_genes <- annotations[grep("RPL", annotations$gene_name, perl=TRUE),]
    rpl_genes <- rpl_genes[grepl("protein_coding", rpl_genes$gene_biotype, perl=TRUE),]
    rpl_genes <- rpl_genes %>% pull(gene_name)
    rb_genes <- c(rps_genes, rpl_genes)

    # combine mt and rb genes
    mt_rb_genes <- c(mt_genes, rb_genes)
    mt_rb_genes <- unique(mt_rb_genes)

    # reduce based on mt & rb genes only
    seurat_mtrb <- subset(seurat, features = intersect(mt_rb_genes, rownames(seurat@assays$RNA)))
    seurat_mtrb <- NormalizeData(seurat_mtrb, verbose = FALSE) %>%
      FindVariableFeatures(verbose = FALSE) %>%
      ScaleData(verbose = FALSE) %>%
      RunPCA(verbose = FALSE) %>%
      FindNeighbors(dims = 1:10, verbose = FALSE) %>%
      FindClusters(resolution = 0.4, verbose = FALSE) %>%
      RunTSNE(dims = 1:10, verbose = FALSE)

    # Annotations to the reduced object based on the unreduced seurat object (all genes)
    DefaultAssay(seurat_mtrb) <- "RNA"
    DefaultAssay(seurat) <- "RNA"

    # Define mitochondrial expression
    seurat_mtrb$mt.percent <- PercentageFeatureSet(
      object =  seurat,
      features = intersect(mt_genes, rownames(seurat@assays$RNA)),
      assay = "RNA"
    )

    # Define mitochondrial expression
    seurat$mt.percent <- PercentageFeatureSet(
      object =  seurat,
      features = intersect(mt_genes, rownames(seurat@assays$RNA)),
      assay = "RNA"
    )

    # Define ribosomal expression
    seurat_mtrb$rb.percent <- PercentageFeatureSet(
      object = seurat,
      features = intersect(rb_genes, rownames(seurat@assays$RNA)),
      assay = "RNA"
    )

    seurat$rb.percent <- PercentageFeatureSet(
      object = seurat,
      features = intersect(rb_genes, rownames(seurat@assays$RNA)),
      assay = "RNA"
    )

    # Automatically find the damaged cell population
    # Calculate the average mt.percent and rb.percent for each cluster
    mt_rb_avg_per_cluster <- seurat_mtrb@meta.data %>%
      group_by(seurat_clusters) %>%
      summarise(
        avg_mt_percent = mean(mt.percent, na.rm = TRUE),
        avg_rb_percent = mean(rb.percent, na.rm = TRUE)
      )

    # Add ranks for each criterion
    mt_rb_avg_per_cluster <- mt_rb_avg_per_cluster %>%
      mutate(rank_mt = rank(-avg_mt_percent),  # Negative for descending order (highest mt to be #1)
             rank_rb = rank(avg_rb_percent),
             combined_rank = rank_mt + rank_rb)

    # Identify the cluster with the best combined rank
    best_cluster <- mt_rb_avg_per_cluster %>%
      arrange(combined_rank) %>%
      slice(1) %>%
      pull(seurat_clusters)

    best_cluster <- as.character(best_cluster)

    # Label all cells belonging to this cluster as "damaged
    seurat_mtrb$seurat_clusters <- ifelse(seurat_mtrb$seurat_clusters == best_cluster,
                                          'damaged_cell',
                                          seurat_mtrb$seurat_clusters)

    # Add Cell QC meta data to object
    # MtRb results
    seurat_mtrb$MtRb <- ifelse(seurat_mtrb$seurat_clusters == "damaged_cell", "damaged_cell", "cell")

    # Add cluster #s to og seurat object for subsetting in next function
    seurat$MtRb <- seurat_mtrb$MtRb
    seurat$MtRb_clusters <- seurat_mtrb$seurat_clusters

    cat("\u2714  MtRb QC predictions\n")

    # View labels in reduced space (MtRb metric)
    clusters <- DimPlot(
      seurat_mtrb, pt.size = 1, group.by = "MtRb", cols = c("cell" = "grey", "damaged_cell" = "#64263E")) +
      labs(caption = "MtRb metric") + NoAxes() +
      theme(
        plot.title = element_blank(),
        plot.caption = element_text(hjust = 0.5, face = "bold", size = 16))

    clusters_labelled <- DimPlot(
      seurat_mtrb, pt.size = 1, label = TRUE) +
      labs(caption = "Labels") + NoAxes() +
      theme(
        plot.title = element_blank(),
        plot.caption = element_text(hjust = 0.5, face = "bold", size = 16))

    # tSNE of mt percent
    mt_plot <- FeaturePlot(seurat_mtrb, features = c("mt.percent"), cols = c("#BCBCBC", "#64263E"), pt.size = 1) +
      NoAxes() + labs(caption = "Mt") +
      theme(
        plot.title = element_blank(),
        plot.caption = element_text(hjust = 0.5, face = "bold", size = 16))

    # tSNE of rb percent
    rb_plot <- FeaturePlot(seurat_mtrb, features = c("rb.percent"), cols = c("#BCBCBC", "#64263E"), pt.size = 1) +
      NoAxes() + labs(caption = "Rb") +
      theme(
        plot.title = element_blank(),
        plot.caption = element_text(hjust = 0.5, face = "bold", size = 16))

    cell_QC <- (clusters | clusters_labelled | mt_plot | rb_plot)

    return(list(
      seuratObject = seurat,
      ImmuneCellFiltering = feature_QC,
      CellQualityFiltering = cell_QC
    ))

  }

# 2. Running control samples ####

# Error with gzipped file redoing separately
PBMC_HC_1 <- PreProcess(
  "PBMC_HC_1/filtered/",
  "PBMC_HC_1/raw/",
  "velocyto/PBMC_HC_1.csv",
  "PBMC_HC_1")

# HC: 2
PBMC_HC_2 <- PreProcess(
  "PBMC_HC_2/filtered/",
  "PBMC_HC_2/raw/",
  "velocyto/PBMC_HC_2.csv",
  "PBMC_HC_2")


PBMC_HC_1$ImmuneCellFiltering / PBMC_HC_2$ImmuneCellFiltering
PBMC_HC_1$CellQualityFiltering / PBMC_HC_2$CellQualityFiltering

# Label all cells belonging to this cluster as "damaged
PBMC_HC_1$seuratObject$MtRb <- ifelse(
  # Already in the damaged cell pop but satisfy criteria
  PBMC_HC_1$seuratObject$MtRb_clusters  == "damaged_cell" &
    PBMC_HC_1$seuratObject$mt.percent >= 10,
  'damaged_cell', 'cell')


# 3. Rest of pre-processing ####

cont_PreProcess <- function(
    seurat,
    project_name
) {
  # 1. Combine cell QC output ####
  seurat$QC <- "cell"
  seurat$QC <- seurat@meta.data %>%
    dplyr::mutate(QC = case_when(
      DQC_cell_status == "cell" & MtRb == "cell" ~ "cell",
      DQC_cell_status == "damaged_cell" & MtRb == "damaged_cell" ~ "DQ_MtRb",
      DQC_cell_status == "damaged_cell" & MtRb == "cell" ~ "DQ",
      DQC_cell_status == "cell" & MtRb == "damaged_cell" ~ "MtRb",
      DQC_cell_status == "empty_droplet" ~ "empty_droplet"
    ))


  # View labels in scatter plot nf vs UMI_count (DropletQC metric)
  df <- seurat@meta.data
  df$nf <- as.numeric(df$nf)
  df$nCount_RNA <- as.numeric(df$nCount_RNA)

  scatter <- ggplot(df, aes(x = nf, y = nCount_RNA, color = QC)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = c("cell" = "grey", "DQ_MtRb" = "red", "DQ" = "#5674B7", "MtRb" = "blue", "empty_droplet" = "#64263E")) +
    coord_trans(y = "log10") + xlab("Nuclear fraction") + ylab("UMI count") + labs(caption = "DropletQC metric") +
    ggtitle(paste(project_name, " cell quality")) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.caption = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title.x = element_text(hjust = 0.5),
      axis.title.y = element_text(hjust = 0.5),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "below",
      legend.box.background = element_rect(colour = "black"),
      strip.background = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black")
    )


  # After visualising, remove damaged cells & empty drops
  Unfiltered_cells <- length(Cells(seurat))
  seurat <- subset(seurat, QC == "cell" | QC == "DQ" | QC == "mt") # only removing empty drops (empty) & those marked as damaged by both metrics
  Final_cells <-  length(Cells(seurat))
  Damaged_cells <- Unfiltered_cells - Final_cells
  cat("\tFiltered cell number:", Final_cells, "\n")


  # 2. DoubletFinder ####

  # Prepare seurat object (required)
  seurat_DF <- NormalizeData(seurat, verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    FindNeighbors(dims = 1:10, verbose = FALSE) %>%
    FindClusters(verbose = FALSE) %>%
    RunUMAP(dims = 1:10, verbose = FALSE)


  cat("\tDoublet Finder running...\n")
  # Open a connection to a temporary file for writing
  tmp_conn <- file(tempfile(), open = "wt")

  # Redirect standard output and messages to the temporary file
  sink(tmp_conn)
  sink(tmp_conn, type = "message")

  # pK Identification (no ground-truth)
  sweep.res.list <- paramSweep(seurat_DF, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)

  pK <- bcmvn %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  pK <- as.numeric(as.character(pK[[1]]))

  # Homotypic Doublet Proportion Estimate
  annotations <- seurat_DF@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075 * nrow(seurat_DF@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  nExp_poi.adj <- nExp_poi.adj - (nExp_poi.adj * (1 / 7))

  # run
  seurat_DF <- doubletFinder(
    seurat_DF,
    PCs = 1:10,
    pN = 0.25,
    pK = pK,
    nExp = nExp_poi.adj,
    reuse.pANN = FALSE,
  )

  # Reset output redirection
  sink(NULL)
  sink(NULL, type = "message")

  # Close the connection
  close(tmp_conn)

  cat("\u2714  Doublet Finder complete\n")

  # Meta data for OG seurat object
  df_col_name <- grep("^DF\\.classifications", colnames(seurat_DF@meta.data), value = TRUE)
  doublet_classifications <- seurat_DF@meta.data[[df_col_name]]
  seurat[['DF']] <- doublet_classifications

  # Print number of damaged cells
  Cells_before_doublets <- length(Cells(seurat))
  seurat <- subset(seurat, DF == "Singlet")
  Remaining_cells <- length(Cells(seurat))
  Doublet_cells <- Cells_before_doublets - Remaining_cells
  cat("\tCell number after doublet removal:", Remaining_cells, "\n")


  # 3. UMAP check to visualize markers #####
  # Prepare seurat object (required)
  seurat <- NormalizeData(seurat, verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE)

  cat(paste("\u2714  Seurat object prepared for integration\n"))

  # No longer altering the same seurat object, this is just for visualisation
  seurat_check <- RunPCA(seurat, verbose = FALSE) %>%
    FindNeighbors(dims = 1:10, verbose = FALSE) %>%
    FindClusters(verbose = FALSE) %>%
    RunUMAP(dims = 1:10, verbose = FALSE)

  markers_plot <- FeaturePlot(seurat_check, features = c("PTPRC", "MS4A1", "CD3E", "NKG7", "CD14", "MARCO"), cols = c("lightgrey", "#324776")) +
    labs(caption = project_name) +
    theme(plot.caption = element_text(hjust = 0.5, vjust = -1.2, size = rel(1), face = "bold", margin = margin(b = 0.5, unit = "cm")))


  plot_clusters <- DimPlot(seurat_check,
                           pt.size = 0.5,
                           label = TRUE,
                           label.size = 4) +
    labs(caption = project_name) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    theme(
      plot.title = element_blank(),
      plot.caption = element_text(hjust = 0.5, size = rel(1), face = "bold"))

  UMAP_markers <- (plot_clusters | markers_plot)

  # 4. Clean & save output #####

  # Filter unnecessary columns
  columns <- c("nf", "ptprc.percent", "hbb.percent", "house.keeping.1", "house.keeping.2", "IMC", "DQC_cell_status", "MtRb", "QC", "DF")

  for (column in columns){
    if (column %in% colnames(seurat@meta.data)) {
      seurat@meta.data$column <- NULL }
    else {
      cat(column, " not found\n") }
  }


  # Timing function
  end_time <- Sys.time()
  time <- (end_time - start_time)
  time_output <- as.character(round(time, digits = 2))
  cat(paste("TIME:", time_output, "\n"))

  PreProcess_output[nrow(PreProcess_output) + 1, ] = list(

    Sample_name   = project_name,                         # Sample_name
    Doublets      = Doublet_cells,                        # Doublets
    Final_cells   = Remaining_cells,                      # final cells
    Time          = time_output                           # Time

  )

  assign("PreProcess_output", PreProcess_output, envir = .GlobalEnv)

  return(list(
    seuratObject = seurat,
    CellQualityFiltering = scatter,
    MarkerVisualisation =  UMAP_markers
  ))

}

PBMC_HC_1 <- PBMC_HC_1$seuratObject
PBMC_HC_2 <- PBMC_HC_2$seuratObject

cont_PreProcess(PBMC_HC_1, "PBMC_HC_1")
cont_PreProcess(PBMC_HC_2, "PBMC_HC_2")

PBMC_HC_1$CellQualityFiltering
PBMC_HC_1$UMAP_markers
PBMC_HC_2$CellQualityFiltering
PBMC_HC_2$UMAP_markers
