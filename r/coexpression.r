#!/usr/bin/env Rscript
library(magrittr)
library(data.table)
library(WGCNA)

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


#' Assign eQTL genes to co-expression modules
assign_egenes <- function(ec_dat, saveto = ".", power_rsq = 0.9) {
  sample_tree <- hclust(dist(ec_dat$exp), method = "average")
}


#
## IO
#
proj_dir <- "~/Documents/projects/wp_bcg_eqtl"

mode <- "pseudo_bulk"
model <- "normal"
cell_type <- "Monocytes"
comp_pair <- ""

in_dir <- file.path(proj_dir, "inputs", mode, cell_type)
out_dir <- file.path(proj_dir, "outputs", mode, "outcomes", model, "coexpression", cell_type, comp_pair)


#
## Main steps
#
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
exp_cov_dat <- load_exp_cov_data(in_dir, "T0_LPS")

# Show the samples tree
sample_tree <- hclust(dist(exp_cov_dat$exp), method = "average")
saveto <- file.path(out_dir, "sample_clustering_tree.pdf")
pdf(saveto, width = 16, height = 7)
par(cex = 0.6, mar = c(0, 4, 2, 0))
plot(sample_tree, main = "Sample clustering to detect outliers", sub = NULL, xlab = NULL, cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = 400, col = "red")
dev.off()

trait_colors <- numbers2colors(exp_cov_dat$cov, signed = FALSE)
saveto <- file.path(out_dir, "sample_dendrogram_and_trait_heatmap.pdf")
pdf(saveto, width = 20, height = 8)
par(cex = 1.5, mar = c(0, 4, 2, 0))
plotDendroAndColors(sample_tree, trait_colors, setLayout = TRUE, autoColorHeight = TRUE, groupLabels = names(exp_cov_dat$cov), main = "Sample dendrogram and trait heatmap")
dev.off()

# Network cosntruction
enableWGCNAThreads(8)

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
net <- blockwiseModules(exp_cov_dat$exp, power = power, TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0,
                        mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMS = TRUE,
                        saveTOMFileBase = out_dir, maxBlockSize = 10000, verbose = 3)

gene_tree <- net$dendrogram[[1]]
module_colors <- labels2colors(net$colors)

saveto <- file.path(out_dir, "hierarchical_clustering_dendrogram_tree.pdf")
pdf(saveto, width = 12, height = 6)
plotDendroAndColors(gene_tree, module_colors[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()


# Visulizing the network
tar_genes <- data.table::fread(file.path(out_dir, "eQTL_genes.csv")) %>% unlist()

diss_TOM <- 1 - TOMsimilarityFromExpr(exp_cov_dat$exp, power = power)

select <- colnames(exp_cov_dat$exp) %in% tar_genes
select_tom <- diss_TOM[select, select]
select_tree <- hclust(as.dist(select_tom), method = "average")
select_colors <- module_colors[select]
plot_tom <- select_tom ^ 7
diag(plot_tom) <- NA

saveto <- file.path(out_dir, "network_heatmap.pdf")
pdf(saveto, width = 7, height = 7)
TOMplot(plot_tom, select_tree, select_colors, main = "Network heatmap plot, selected genes")
dev.off()
