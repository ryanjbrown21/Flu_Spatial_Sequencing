### NEED TO USE SPOTlight v 0.1.7 for this to work
### Default version of SPOTlight has removed several functions
# install_github("MarcElosua/SPOTlight@spotlight-0.1.7")

setwd("/scratch/u/mkasmani/SPOTlight_Flu_4timepoints")

### Rearrange libpaths for package installs on the RCC (only needed for SPOTlight because it's installed via install_github)
# .libPaths(.libPaths()[c(2, 1)])
# install_github("MarcElosua/SPOTlight")
### Failed initially b/c of failed install of "imager" package, with error:
# FFTW library not found, please install fftw3 for better FFT support
### Contacted RCC admin to help install FFTW3 in main terminal (don't have permission to do so myself)

library(tidyverse)
library(Seurat)
library(patchwork)
library(ggplot2) 
library(SPOTlight)
library(RColorBrewer)
library(gsubfn)
library(ggcorrplot)


set.seed(123)
scRNAseq_datasets_Visium_final <- readRDS("Flu_pure_RJB.RDS")

DefaultAssay(scRNAseq_datasets_Visium_final) <- "RNA"
# scRNAseq_datasets_Visium_final <- NormalizeData(scRNAseq_datasets_Visium_final)
# scRNAseq_datasets_Visium_final <- ScaleData(scRNAseq_datasets_Visium_final, features = rownames(scRNAseq_datasets_Visium_final))

scRNAseq_datasets_Visium_final <- StashIdent(scRNAseq_datasets_Visium_final, save.name = "cell_type")

cluster_markers_all <- FindAllMarkers(object = scRNAseq_datasets_Visium_final, 
                                      assay = "RNA",
                                      slot = "data",
                                      verbose = TRUE, 
                                      only.pos = TRUE)

write.csv(cluster_markers_all, "Flu_pure_markers_RJB.csv")
cluster_markers_all <- read.csv("Flu_pure_markers_RJB.csv", row.names = 1)

se_sc_down <- downsample_se_obj(se_obj = scRNAseq_datasets_Visium_final,
                                clust_vr = "cell_type",
                                cluster_markers = cluster_markers_all,
                                cl_n = 100,
                                hvg = 2000)

# Load and process spatial data
# Need "Spatial" folder and "filtered_feature_bc_matrix.h5" in the same folder entered into Load10X_Spatial

Day0_Young <- Load10X_Spatial("Visium_data/A1_Day0_Young/", slice = "Day 0 Young")
Day0_Old <- Load10X_Spatial("Visium_data/B1_Day0_Old/", slice = "Day 0 Aged")
Day3_Young <- Load10X_Spatial("Visium_data/A1_Day3_Young/", slice = "Day 3 Young")
Day3_Old <- Load10X_Spatial("Visium_data/B1_Day3_Old/", slice = "Day 3 Aged")
Day7_Young <- Load10X_Spatial("Visium_data/C1_Day7_Young/", slice = "Day 7 Young")
Day7_Old <- Load10X_Spatial("Visium_data/D1_Day7_Old/", slice = "Day 7 Aged")
Day9_Young <- Load10X_Spatial("Visium_data/C1_Day9_Young/", slice = "Day 9 Young")
Day9_Old <- Load10X_Spatial("Visium_data/D1_Day9_Old/", slice = "Day 9 Aged")

# Add sample names to barcodes
Day0_Young <- RenameCells(Day0_Young, add.cell.id = "Day0_Young")
Day0_Old <- RenameCells(Day0_Old, add.cell.id = "Day0_Aged")
Day3_Young <- RenameCells(Day3_Young, add.cell.id = "Day3_Young")
Day3_Old <- RenameCells(Day3_Old, add.cell.id = "Day3_Aged")
Day7_Young <- RenameCells(Day7_Young, add.cell.id = "Day7_Young")
Day7_Old <- RenameCells(Day7_Old, add.cell.id = "Day7_Aged")
Day9_Young <- RenameCells(Day9_Young, add.cell.id = "Day9_Young")
Day9_Old <- RenameCells(Day9_Old, add.cell.id = "Day9_Aged")

# Add columns specifying orig.ident, timepoint, and condition
Day0_Young$orig.ident <- "Day 0 Young"
Day0_Old$orig.ident <- "Day 0 Aged"
Day3_Young$orig.ident <- "Day 3 Young"
Day3_Old$orig.ident <- "Day 3 Aged"
Day7_Young$orig.ident <- "Day 7 Young"
Day7_Old$orig.ident <- "Day 7 Aged"
Day9_Young$orig.ident <- "Day 9 Young"
Day9_Old$orig.ident <- "Day 9 Aged"

Day0_Young$Timepoint <- "Day 0"
Day0_Old$Timepoint <- "Day 0"
Day3_Young$Timepoint <- "Day 3"
Day3_Old$Timepoint <- "Day 3"
Day7_Young$Timepoint <- "Day 7"
Day7_Old$Timepoint <- "Day 7"
Day9_Young$Timepoint <- "Day 9"
Day9_Old$Timepoint <- "Day 9"

Day0_Young$Age <- "Young"
Day0_Old$Age <- "Aged"
Day3_Young$Age <- "Young"
Day3_Old$Age <- "Aged"
Day7_Young$Age <- "Young"
Day7_Old$Age <- "Aged"
Day9_Young$Age <- "Young"
Day9_Old$Age <- "Aged"

# Add pct mito
Day0_Young[["percent.mt"]] <- PercentageFeatureSet(Day0_Young, pattern = "^mt-")
Day0_Old[["percent.mt"]] <- PercentageFeatureSet(Day0_Old, pattern = "^mt-")
Day3_Young[["percent.mt"]] <- PercentageFeatureSet(Day3_Young, pattern = "^mt-")
Day3_Old[["percent.mt"]] <- PercentageFeatureSet(Day3_Old, pattern = "^mt-")
Day7_Young[["percent.mt"]] <- PercentageFeatureSet(Day7_Young, pattern = "^mt-")
Day7_Old[["percent.mt"]] <- PercentageFeatureSet(Day7_Old, pattern = "^mt-")
Day9_Young[["percent.mt"]] <- PercentageFeatureSet(Day9_Young, pattern = "^mt-")
Day9_Old[["percent.mt"]] <- PercentageFeatureSet(Day9_Old, pattern = "^mt-")


DefaultAssay(Day0_Young) <- "Spatial"
DefaultAssay(Day0_Old) <- "Spatial"
DefaultAssay(Day3_Young) <- "Spatial"
DefaultAssay(Day3_Old) <- "Spatial"
DefaultAssay(Day7_Young) <- "Spatial"
DefaultAssay(Day7_Old) <- "Spatial"
DefaultAssay(Day9_Young) <- "Spatial"
DefaultAssay(Day9_Old) <- "Spatial"

library(future)
plan(strategy = "multiprocess", workers = 24)
options(future.globals.maxSize = 4000 * 1024^2)
#this will run all samples in parallel.
#Change workers = 1 if you run out of memory (RAM).

spotlight_ls_d0_young %<-% spotlight_deconvolution(
  se_sc = se_sc_down,
  counts_spatial = Day0_Young@assays$Spatial@counts,
  clust_vr = "cell_type", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = cluster_markers_all, # Dataframe with the marker genes
  cl_n = 100, # number of cells per cell type to use
  hvg = 2000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorization and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold 
)

spotlight_ls_d0_old %<-% spotlight_deconvolution(
  se_sc = se_sc_down,
  counts_spatial = Day0_Old@assays$Spatial@counts,
  clust_vr = "cell_type", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = cluster_markers_all, # Dataframe with the marker genes
  cl_n = 100, # number of cells per cell type to use
  hvg = 2000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorization and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold 
)

spotlight_ls_d3_young %<-% spotlight_deconvolution(
  se_sc = se_sc_down,
  counts_spatial = Day3_Young@assays$Spatial@counts,
  clust_vr = "cell_type", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = cluster_markers_all, # Dataframe with the marker genes
  cl_n = 100, # number of cells per cell type to use
  hvg = 2000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorization and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold 
)

spotlight_ls_d3_old %<-% spotlight_deconvolution(
  se_sc = se_sc_down,
  counts_spatial = Day3_Old@assays$Spatial@counts,
  clust_vr = "cell_type", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = cluster_markers_all, # Dataframe with the marker genes
  cl_n = 100, # number of cells per cell type to use
  hvg = 2000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorization and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold 
)

spotlight_ls_d7_young %<-% spotlight_deconvolution(
  se_sc = se_sc_down,
  counts_spatial = Day7_Young@assays$Spatial@counts,
  clust_vr = "cell_type", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = cluster_markers_all, # Dataframe with the marker genes
  cl_n = 100, # number of cells per cell type to use
  hvg = 2000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorization and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold 
)

spotlight_ls_d7_old %<-% spotlight_deconvolution(
  se_sc = se_sc_down,
  counts_spatial = Day7_Old@assays$Spatial@counts,
  clust_vr = "cell_type", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = cluster_markers_all, # Dataframe with the marker genes
  cl_n = 100, # number of cells per cell type to use
  hvg = 2000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorization and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold 
)

spotlight_ls_d9_young %<-% spotlight_deconvolution(
  se_sc = se_sc_down,
  counts_spatial = Day9_Young@assays$Spatial@counts,
  clust_vr = "cell_type", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = cluster_markers_all, # Dataframe with the marker genes
  cl_n = 100, # number of cells per cell type to use
  hvg = 2000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorization and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold 
)

spotlight_ls_d9_old %<-% spotlight_deconvolution(
  se_sc = se_sc_down,
  counts_spatial = Day9_Old@assays$Spatial@counts,
  clust_vr = "cell_type", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = cluster_markers_all, # Dataframe with the marker genes
  cl_n = 100, # number of cells per cell type to use
  hvg = 2000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorization and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold 
)


#I've wrapped most of the figure creation and data extraction code from:
#https://marcelosua.github.io/SPOTlight/
#as a function here

figs_and_mats <- function(spotlight_obj, visium_obj, visium_sample_image_path_low_res){
  
  output_list <- list()
  
  nmf_mod <- spotlight_obj[[1]]
  decon_mtrx <- spotlight_obj[[2]]
  
  h <- NMF::coef(nmf_mod[[1]])
  rownames(h) <- paste("Topic", 1:nrow(h), sep = "_")
  topic_profile_plts <- SPOTlight::dot_plot_profiles_fun(
    h = h,
    train_cell_clust = nmf_mod[[2]])
  
  output_list[["topic_by_cell_type"]] <- topic_profile_plts[[2]] + ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90), 
    axis.text = ggplot2::element_text(size = 12))
  
  output_list[["topic_within_cell_type"]] <- topic_profile_plts[[1]] + theme(axis.text.x = element_text(angle = 90), 
                                                                             axis.text = element_text(size = 12))
  
  basis_spotlight <- data.frame(NMF::basis(nmf_mod[[1]]))
  
  colnames(basis_spotlight) <- unique(stringr::str_wrap(nmf_mod[[2]], width = 30))
  
  decon_mtrx_sub <- decon_mtrx[, colnames(decon_mtrx) != "res_ss"]
  decon_mtrx_sub[decon_mtrx_sub < 0.08] <- 0
  decon_mtrx <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx[, "res_ss"])
  rownames(decon_mtrx) <- colnames(visium_obj)
  
  decon_df <- decon_mtrx %>%
    data.frame() %>%
    tibble::rownames_to_column("barcodes")
  
  visium_obj@meta.data <- visium_obj@meta.data %>%
    tibble::rownames_to_column("barcodes") %>%
    dplyr::left_join(decon_df, by = "barcodes") %>%
    tibble::column_to_rownames("barcodes")
  
  cell_types_all <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "res_ss")]
  
  output_list[["scatterpie"]] <- SPOTlight::spatial_scatterpie(se_obj = visium_obj,
                                                               cell_types_all = cell_types_all,
                                                               img_path = visium_sample_image_path_low_res,
                                                               pie_scale = 0.4)
  
  #Interaction Graph
  graph_ntw <- SPOTlight::get_spatial_interaction_graph(decon_mtrx = decon_mtrx[, cell_types_all])
  
  deg <- degree(graph_ntw, mode="all")
  
  # Get color palette for difusion
  edge_importance <- E(graph_ntw)$importance
  
  # Select a continuous palette
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'seq',]
  
  # Create a color palette
  getPalette <- colorRampPalette(brewer.pal(9, "YlOrRd"))
  
  # Get how many values we need
  grad_edge <- seq(0, max(edge_importance), 0.1)
  
  # Generate extended gradient palette dataframe
  graph_col_df <- data.frame(value = as.character(grad_edge),
                             color = getPalette(length(grad_edge)),
                             stringsAsFactors = FALSE)
  # Assign color to each edge
  color_edge <- data.frame(value = as.character(round(edge_importance, 1)), stringsAsFactors = FALSE) %>%
    dplyr::left_join(graph_col_df, by = "value") %>%
    dplyr::pull(color)
  
  # Open a pdf file
  output_list[["graph_ntw"]] <- plot(graph_ntw,
                                     # Size of the edge
                                     edge.width = edge_importance,
                                     edge.color = color_edge,
                                     # Size of the buble
                                     vertex.size = deg,
                                     vertex.color = "#cde394",
                                     vertex.frame.color = "white",
                                     vertex.label.color = "black",
                                     vertex.label.family = "Times", # Font family of the label (e.g.“Times”, “Helvetica”)
                                     layout = layout.circle)
  
  #Cell cell correlation matrix
  # Remove cell types not predicted to be on the tissue
  decon_mtrx_sub <- decon_mtrx[, cell_types_all]
  decon_mtrx_sub <- decon_mtrx_sub[, colSums(decon_mtrx_sub) > 0]
  
  # Compute correlation
  decon_cor <- cor(decon_mtrx_sub)
  
  # Compute correlation P-value
  p.mat <- corrplot::cor.mtest(mat = decon_mtrx_sub, conf.level = 0.95)
  
  output_list[["decon_cor"]] <- decon_cor
  
  # Visualize
  output_list[["corrplot"]] <- ggcorrplot::ggcorrplot(
    corr = decon_cor,
    p.mat = p.mat[[1]],
    hc.order = TRUE,
    type = "full",
    insig = "blank",
    lab = TRUE,
    outline.col = "lightgrey",
    method = "square",
    # colors = c("#4477AA", "white", "#BB4444"))
    colors = c("#6D9EC1", "white", "#E46726"),
    title = "Predicted cell-cell proportion correlation",
    legend.title = "Correlation\n(Pearson)") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 22, hjust = 0.5, face = "bold"),
      legend.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(size = 15),
      axis.text.x = ggplot2::element_text(angle = 90),
      axis.text = ggplot2::element_text(size = 18, vjust = 0.5))
  
  return(output_list)
}

d0_young_out <- figs_and_mats(spotlight_obj = spotlight_ls_d0_young, visium_obj = Day0_Young, visium_sample_image_path_low_res = "Visium_data/A1_Day0_Young/spatial/tissue_lowres_image.png")

d0_old_out <- figs_and_mats(spotlight_obj = spotlight_ls_d0_old, visium_obj = Day0_Old, visium_sample_image_path_low_res = "Visium_data/B1_Day0_Old/spatial/tissue_lowres_image.png")

d3_young_out <- figs_and_mats(spotlight_obj = spotlight_ls_d3_young, visium_obj = Day3_Young, visium_sample_image_path_low_res = "Visium_data/A1_Day3_Young/spatial/tissue_lowres_image.png")

d3_old_out <- figs_and_mats(spotlight_obj = spotlight_ls_d3_old, visium_obj = Day3_Old, visium_sample_image_path_low_res = "Visium_data/B1_Day3_Old/spatial/tissue_lowres_image.png")


d7_young_out <- figs_and_mats(spotlight_obj = spotlight_ls_d7_young, visium_obj = Day7_Young, visium_sample_image_path_low_res = "Visium_data/C1_Day7_Young/spatial/tissue_lowres_image.png")

d7_old_out <- figs_and_mats(spotlight_obj = spotlight_ls_d7_old, visium_obj = Day7_Old, visium_sample_image_path_low_res = "Visium_data/D1_Day7_Old/spatial/tissue_lowres_image.png")

d9_young_out <- figs_and_mats(spotlight_obj = spotlight_ls_d9_young, visium_obj = Day9_Young, visium_sample_image_path_low_res = "Visium_data/C1_Day9_Young/spatial/tissue_lowres_image.png")

d9_old_out <- figs_and_mats(spotlight_obj = spotlight_ls_d9_old, visium_obj = Day9_Old, visium_sample_image_path_low_res = "Visium_data/D1_Day9_Old/spatial/tissue_lowres_image.png")


dir.create("SPOTlight_plots_RJB")

saveRDS(spotlight_ls_d0_young, "SPOTlight_plots_RJB/d0_young_spotlight.RDS")
saveRDS(spotlight_ls_d0_old, "SPOTlight_plots_RJB/d0_old_spotlight.RDS")
saveRDS(spotlight_ls_d3_young, "SPOTlight_plots_RJB/d3_young_spotlight.RDS")
saveRDS(spotlight_ls_d3_old, "SPOTlight_plots_RJB/d3_old_spotlight.RDS")

saveRDS(d0_young_out, "SPOTlight_plots_RJB/d0_young_out.RDS")
saveRDS(d0_old_out, "SPOTlight_plots_RJB/d0_old_out.RDS")
saveRDS(d3_young_out, "SPOTlight_plots_RJB/d3_young_out.RDS")
saveRDS(d3_old_out, "SPOTlight_plots_RJB/d3_old_out.RDS")

saveRDS(spotlight_ls_d7_young, "SPOTlight_plots_RJB/d7_young_spotlight.RDS")
saveRDS(spotlight_ls_d7_old, "SPOTlight_plots_RJB/d7_old_spotlight.RDS")
saveRDS(spotlight_ls_d9_young, "SPOTlight_plots_RJB/d9_young_spotlight.RDS")
saveRDS(spotlight_ls_d9_old, "SPOTlight_plots_RJB/d9_old_spotlight.RDS")

saveRDS(d7_young_out, "SPOTlight_plots_RJB/d7_young_out.RDS")
saveRDS(d7_old_out, "SPOTlight_plots_RJB/d7_old_out.RDS")
saveRDS(d9_young_out, "SPOTlight_plots_RJB/d9_young_out.RDS")
saveRDS(d9_old_out, "SPOTlight_plots_RJB/d9_old_out.RDS")



ggsave(filename = "SPOTlight_plots_RJB/d0_young_scatterpie.pdf", plot = d0_young_out$scatterpie, device = "pdf", height = 20, width = 20, units = "in")
ggsave(filename = "SPOTlight_plots_RJB/d0_old_scatterpie.pdf", plot = d0_old_out$scatterpie, device = "pdf", height = 20, width = 20, units = "in")
ggsave(filename = "SPOTlight_plots_RJB/d3_young_scatterpie.pdf", plot = d3_young_out$scatterpie, device = "pdf", height = 20, width = 20, units = "in")
ggsave(filename = "SPOTlight_plots_RJB/d3_old_scatterpie.pdf", plot = d3_old_out$scatterpie, device = "pdf", height = 20, width = 20, units = "in")
ggsave(filename = "SPOTlight_plots_RJB/d0_young_corrplot.pdf", plot = d0_young_out$corrplot, device = "pdf", height = 15, width = 15, units = "in")
ggsave(filename = "SPOTlight_plots_RJB/d0_old_corrplot.pdf", plot = d0_old_out$corrplot, device = "pdf", height = 15, width = 15, units = "in")
ggsave(filename = "SPOTlight_plots_RJB/d3_young_corrplot.pdf", plot = d3_young_out$corrplot, device = "pdf", height = 15, width = 15, units = "in")
ggsave(filename = "SPOTlight_plots_RJB/d3_old_corrplot.pdf", plot = d3_old_out$corrplot, device = "pdf", height = 15, width = 15, units = "in")

ggsave(filename = "SPOTlight_plots_RJB/d7_young_scatterpie.pdf", plot = d7_young_out$scatterpie, device = "pdf", height = 20, width = 20, units = "in")
ggsave(filename = "SPOTlight_plots_RJB/d7_old_scatterpie.pdf", plot = d7_old_out$scatterpie, device = "pdf", height = 20, width = 20, units = "in")
ggsave(filename = "SPOTlight_plots_RJB/d9_young_scatterpie.pdf", plot = d9_young_out$scatterpie, device = "pdf", height = 20, width = 20, units = "in")
ggsave(filename = "SPOTlight_plots_RJB/d9_old_scatterpie.pdf", plot = d9_old_out$scatterpie, device = "pdf", height = 20, width = 20, units = "in")
ggsave(filename = "SPOTlight_plots_RJB/d7_young_corrplot.pdf", plot = d7_young_out$corrplot, device = "pdf", height = 15, width = 15, units = "in")
ggsave(filename = "SPOTlight_plots_RJB/d7_old_corrplot.pdf", plot = d7_old_out$corrplot, device = "pdf", height = 15, width = 15, units = "in")
ggsave(filename = "SPOTlight_plots_RJB/d9_young_corrplot.pdf", plot = d9_young_out$corrplot, device = "pdf", height = 15, width = 15, units = "in")
ggsave(filename = "SPOTlight_plots_RJB/d9_old_corrplot.pdf", plot = d9_old_out$corrplot, device = "pdf", height = 15, width = 15, units = "in")


