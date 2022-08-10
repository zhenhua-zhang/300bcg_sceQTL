#!/usr/bin/env Rscript
library(magrittr)
library(data.table)
library(dynamicTreeCut)
library(fastcluster)
library(WGCNA)
library(gplots)

#' Load expression and covariates from the disk
load_exp_cov_data <- function(path, subset_by = NULL) {
  exp_tab <- file.path(path, "phenotypes.tsv") %>% data.table::fread(verbose = FALSE)
  proc_exp_tab <- exp_tab %>%
    dplyr::select(feature_id, dplyr::contains(subset_by)) %>%
    dplyr::filter(rowSums(dplyr::select(., -feature_id)) > 0) %>%
    tidyr::pivot_longer(cols = -feature_id, names_to = "sample_id", values_to = "expression") %>%
    tidyr::pivot_wider(id_cols = "sample_id", names_from = "feature_id", values_from = "expression") %>%
    (function(dat) {
      dat <- as.data.frame(dat) %>% dplyr::arrange(sample_id)
      rownames(dat) <- dat$sample_id
      dplyr::select(dat, -sample_id)
    })

  all_ok <- goodSamplesGenes(proc_exp_tab)$allOK
  if (!all_ok) cat("[W]: Not all samples and genes passed the check.\n")

  cov_tab <- file.path(path, "covariates.tsv") %>% data.table::fread(verbose = FALSE)
  proc_cov_tab <- cov_tab %>%
    dplyr::filter(stringr::str_detect(sample_id, subset_by)) %>%
    (function(dat) {
      dat <- as.data.frame(dat) %>% dplyr::arrange(sample_id)
      rownames(dat) <- dat$sample_id
      dplyr::select(dat, -sample_id)
    })

  return(list(exp = proc_exp_tab, cov = proc_cov_tab))
}


#
## Main steps
#
proj_dir <- "/home/zzhang/Documents/projects/wp_bcg_eqtl"
in_dir <- file.path(proj_dir, "outputs/pseudo_bulk")

# Main cell types.
cell_type_vec <- c("Monocytes", "CD4T", "CD8T", "NK", "B") # , "pDC", "mDC")

# All comparisons used in the analysis.
condition_vec <- c("T0_LPS", "T0_RPMI", "T3m_LPS", "T3m_RPMI")

# Network cosntruction
enableWGCNAThreads(8)

# setseed
set.seed(1024)

# Color
myheatcol <- colorpanel(250, "red", "orange", "lemonchiffon")

for (cell_type in cell_type_vec) {
  celltype_dir <- file.path(proj_dir, "outputs/pseudo_bulk/outcomes/coexpression", cell_type)

  for (condition in condition_vec) {
    out_dir <- file.path(celltype_dir, condition)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

    in_dir <- file.path(proj_dir, "inputs/pseudo_bulk", cell_type)
    exp_cov_dat <- load_exp_cov_data(in_dir, condition)

    gene_sets <- list()
    gene_sets$egene <- data.table::fread(file.path(celltype_dir, "scegene.txt"))$gene_name
    gene_sets$random <- base::sample(1:ncol(exp_cov_dat$exp), length(gene_sets$egene))

    # Show the samples tree
    sample_tree <- hclust(dist(exp_cov_dat$exp), method = "average")

    saveto <- file.path(out_dir, "sample_clustering_tree.pdf")
    pdf(saveto, width = 16, height = 7)
    par(cex = 0.6, mar = c(0, 4, 2, 0))
    plot(sample_tree, main = "Sample clustering to detect outliers", sub = NULL, xlab = NULL, cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
    abline(h = 400, col = "red")
    dev.off()

    # Dendrogram and trait
    trait_colors <- numbers2colors(exp_cov_dat$cov, signed = FALSE)
    saveto <- file.path(out_dir, "sample_dendrogram_and_trait_heatmap.pdf")
    pdf(saveto, width = 20, height = 8)
    par(cex = 1.5, mar = c(0, 4, 2, 0))
    plotDendroAndColors(sample_tree, trait_colors, setLayout = TRUE, autoColorHeight = TRUE, groupLabels = names(exp_cov_dat$cov), main = "Sample dendrogram and trait heatmap")
    dev.off()

    ## Determine soft threshold
    power_vec <- c(c(1:10, seq(12, 20, 2)))
    sft <- pickSoftThreshold(exp_cov_dat$exp, verbose = 5, networkType = "unsigned")
    saveto <- file.path(out_dir, "soft_threhold_power.pdf")
    pdf(saveto, width = 7, height = 7)
    plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], xlab = "Soft threshold (power)", ylab = "Scale free topology model fit, singed R^2", type = "n", main = "Scale independence")
    text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = power_vec, cex = 0.9, col = "red")
    abline(h = 0.9, col = "red")
    plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft threshold (power)", ylab = "Mean connectivity", type = "n", main = "Mean connectivity")
    text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = power_vec, cex = 0.9, col = "red")
    dev.off()

    ## Chose the power for the module defination analysis.
    power <- sft$fitIndices %>%
      dplyr::filter(SFT.R.sq < 0.9) %>%
      dplyr::slice_max(SFT.R.sq, n = 1) %>%
      dplyr::slice_max(Power, n = 1) %>%
      dplyr::select(Power)

    ## Create a network
    net <- blockwiseModules(exp_cov_dat$exp,
      power = power, TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0,
      mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMS = TRUE,
      saveTOMFileBase = out_dir, maxBlockSize = 10000, verbose = 3
    )

    gene_tree <- net$dendrogram[[1]]
    module_colors <- labels2colors(net$colors)

    saveto <- file.path(out_dir, "hierarchical_clustering_dendrogram_tree.pdf")
    pdf(saveto, width = 12, height = 6)
    plotDendroAndColors(gene_tree, module_colors[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
    dev.off()

    saveto <- file.path(out_dir, "net_WGCNA.rdata")
    saveRDS(net, file = saveto)

    # Visulizing the network
    diss_TOM <- 1 - TOMsimilarityFromExpr(exp_cov_dat$exp, power = power)
    rownames(diss_TOM) <- colnames(exp_cov_dat$exp)
    colnames(diss_TOM) <- colnames(exp_cov_dat$exp)

    # saveto <- file.path(out_dir, "disstance_TOM.rdata")
    # save(diss_TOM, file = saveto)

    for (per_set in names(gene_sets)) {
      if (per_set == "egene") {
        select <- which(gene_sets[[per_set]] %in% rownames(diss_TOM))
      } else {
        select <- gene_sets[[per_set]]
      }

      select_tom <- diss_TOM[select, select]
      select_tree <- hclust(as.dist(select_tom), method = "average")
      select_colors <- module_colors[select]
      plot_diss <- select_tom^7
      diag(plot_diss) <- NA

      saveto <- file.path(out_dir, paste0("network_heatmap_", per_set, ".pdf"))
      pdf(saveto, width = 7, height = 7)
      TOMplot(plot_diss, select_tree, select_colors,
        col = myheatcol,
        main = paste0("Network heatmap plot ", "(", per_set, ")")
      )
      dev.off()
    }

    gc()
  }
}
