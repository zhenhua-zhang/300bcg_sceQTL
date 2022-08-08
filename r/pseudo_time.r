#!/usr/bin/env Rscript
library(Seurat)
library(monocle3)
library(magrittr)
library(tidyverse)
library(SeuratWrappers)

proj_dir <- "~/Documents/projects/wp_bcg_eqtl"

cr_dir <- file.path(proj_dir, "inputs/sc_rnaseq/cellranger")
save_dir <- file.path(proj_dir, "outputs/pseudo_time")

# FIXME:
# 1. How to remove MT-genes?
# 2. How to remove cells with highly expressed MT-genes?

# Load the dataset (10X)
# Q: If combine the CDS, how to distinguish them from each other?
# A: Using a named list of cell_data_set object.

from_scratch <- FALSE
if (from_scratch) {
  tar_pools <- c()
  cds <- list.dirs(cr_dir, recursive = FALSE) %>%
    setNames(nm = str_remove_all(basename(.), "atch|ool|_|count$")) %>%
    lapply(FUN = function(e) {
      load_cellranger_data(e)
    }) %>%
    combine_cds()

  # Preprocess

  # Here we use the first 50 PCs to normalize the data.
  cds <- preprocess_cds(cds, num_dim = 50)

  pdf(file.path(save_dir, "pc_variance_explained.pdf")) # Var. explained by PCs
  plot_pc_variance_explained(cds)
  dev.off()

  # Reduce dimensionality
  cds <- reduce_dimension(cds,
    reduction_method = "UMAP", preprocess_method = "PCA", cores = 4
  )

  pdf(file.path(save_dir, "cells_reduced_raw.pdf")) # Show the cell groups
  plot_cells(cds,
    reduction_method = "UMAP", color_cells_by = "sample", label_cell_groups = F
  )
  dev.off()

  # Remove batches
  # TODO: cite the Batchelor paper.
  cds <- align_cds(cds, num_dim = 50, alignment_group = "sample")
  cds <- reduce_dimension(cds,
    reduction_method = "UMAP", preprocess_method = "Aligned"
  )

  pdf(file.path(save_dir, "cells_reduced_aligned.pdf"))
  plot_cells(cds,
    reduction_method = "UMAP", color_cells_by = "sample", label_cell_groups = F
  )
  dev.off()

  # Group cells into clusters
  cds <- cluster_cells(cds, resolution = 1e-5)

  pdf(file.path(save_dir, "cells_reduced_aligned_clustered.pdf"))
  plot_cells(cds, reduction_method = "UMAP")
  dev.off()

  # Find marker genes
  marker_test_res <- top_markers(cds,
    group_cells_by = "partition", reference_cells = 1000, cores = 4
  )

  top_specific_markers <- marker_test_res %>%
    dplyr::filter(fraction_expressing >= 0.10) %>%
    dplyr::group_by(cell_group) %>%
    dplyr::slice_min(order_by = pseudo_R2, n = 10)

  pdf(file.path(save_dir, "top_markers_expression.pdf"))
  plot_genes_by_group(cds, unique(top_specific_markers$gene_short_name),
    group_cells_by = "partition", ordering_type = "maximal_on_diag"
  )
  dev.off()

  # Assign cell type to the partitions

  ## Manually
  colData(cds)$assigned_cell_type <- as.character(partitions(cds))

  colData(cds)$assigned_cell_type <- dplyr::recode(
    colData(cds)$assigned_cell_type,
    "1" = "a", "2" = "b", "3" = "c", "4" = "d", "5" = "e", "6" = "f", "7" = "g"
  )

  pdf(file.path(save_dir, "cells_reduced_aligned_clustered_assigned.pdf"))
  plot_cells(cds,
    group_cells_by = "partition", color_cells_by = "assigned_cell_type"
  )
  dev.off()

  ## Automatically by Garnett (https://cole-trapnell-lab.github.io/garnett)

  # Learn the trajectory graph
  cds <- learn_graph(cds)

  pdf(file.path(
    save_dir, "cells_reduced_aligned_clustered_assigned_graphed.pdf"
  ))
  plot_cells(cds,
    color_cells_by = "assigned_cell_type", label_groups_by_cluster = F,
    label_leaves = T, label_branch_points = T, label_principal_points = T,
    graph_label_size = 1.5
  )
  dev.off()

  # Pseudotime by given root cells.
  cds <- order_cells(cds, reduction_method = "UMAP", root_pr_nodes = "Y_1")
  pdf(file.path(
    save_dir,
    "cells_reduced_aligned_clustered_assigned_graphed_pseudotime.pdf"
  ))
  plot_cells(cds,
    color_cells_by = "pseudotime",
    label_cell_groups = F,
    label_leaves = F,
    label_branch_points = F,
    graph_label_size = 1.5
  )
  dev.off()

  # Save the monocle3 object
  # FIXME: Not possible to save the required rdd_umap_transform_model_umap.idx
  save_monocle_objects(
    cds = cds, directory_path = file.path(save_dir, "monocle_object"),
    comment = "Test"
  )

  cds <- load_monocle_objects(
    directory_path = file.path(save_dir, "monocle_object")
  )
} else {
  pbmc_path <- file.path(proj_dir, "inputs/sc_rnaseq/bcg4-0712.rds")
  pbmc <- readRDS(pbmc_path)
  DefaultAssay(pbmc) <- "RNA"

  group_by <- "time"
  groups <- switch(group_by == "time",
    c("T0", "T3m"),
    c("LPS", "RPMI")
  )

  color_by <- switch(group_by == "time",
    "stim",
    "time"
  )
  colors <- switch(color_by == "time",
    c("T0", "T3m"),
    c("LPS", "RPMI")
  )

  celltype <- "Monocytes"
  for (per_group in groups) {
    pref <- paste0(celltype, "_", per_group, "_")
    tar_cells <- pbmc@meta.data[, group_by] == per_group &
      pbmc@meta.data$clusters1 == celltype

    sub_pbmc <- pbmc[, tar_cells]
    sub_pbmc_cds <- as.cell_data_set(sub_pbmc)

    pdf(file.path(save_dir, paste0(pref, "reduced_aligned.pdf")))
    plot_cells(sub_pbmc_cds,
      reduction_method = "UMAP", color_cells_by = color_by,
      label_cell_groups = F
    )
    dev.off()
  }

  # Subset the cell types.
  tar_cells <- pbmc@meta.data$clusters1 == celltype
  if (condition %in% c("T0", "T3m")) {
    tar_cells <- tar_cells & pbmc@meta.data$time == condition
  } else if (condition %in% c("LPS", "RPMI")) {
    tar_cells <- tar_cells & pbmc@meta.data$stim == condition
  }
  sub_pbmc <- pbmc[, tar_cells]

  # Transform the Seurat object into a cell_data_set (CDS) object.
  sub_pbmc_cds <- as.cell_data_set(sub_pbmc)
  pdf(file.path(save_dir, paste0(pref, "reduced_aligned.pdf")))
  plot_cells(sub_pbmc_cds,
    reduction_method = "UMAP", color_cells_by = "time", label_cell_groups = F
  )
  dev.off()

  # Group cells into clusters
  sub_pbmc_cds <- cluster_cells(cds = sub_pbmc_cds, reduction_method = "UMAP")
  pdf(file.path(save_dir, paste0(pref, "reduced_aligned_clustered.pdf")))
  plot_cells(sub_pbmc_cds, reduction_method = "UMAP")
  dev.off()

  # Learn the trajectory graph
  sub_pbmc_cds <- learn_graph(sub_pbmc_cds, use_partition = TRUE)
  pdf(file.path(
    save_dir, paste0(pref, "reduced_aligned_clustered_graphed.pdf")
  ))
  plot_cells(sub_pbmc_cds,
    color_cells_by = "time", label_groups_by_cluster = F,
    label_leaves = F, label_branch_points = F, label_principal_points = T,
    graph_label_size = 1.5
  )
  dev.off()

  # Pseudotime by given root cells.
  sub_pbmc_cds <- order_cells(
    sub_pbmc_cds,
    reduction_method = "UMAP", root_pr_nodes = "Y_1"
  )
  pdf(file.path(
    save_dir, paste0(pref, "reduced_aligned_clustered_graphed_pseudotime.pdf")
  ))
  plot_cells(
    cds = sub_pbmc_cds,
    color_cells_by = "pseudotime",
    show_trajectory_graph = FALSE
  )
  dev.off()

  pt_mat <- sub_pbmc_cds@principal_graph_aux@listData$UMAP$pseudotime %>%
    purrr::map_dbl(.f = ~ if_else(is.finite(.x), .x, 15))

  sub_pbmc <- AddMetaData(sub_pbmc, metadata = pt_mat, col.name = "pseudotime")

  g_pstm <- FeaturePlot(sub_pbmc, features = "pseudotime", raster = F) +
    theme_classic()

  g_stim <- FeaturePlot(sub_pbmc, features = "stim", raster = F) +
    theme_classic()
}
