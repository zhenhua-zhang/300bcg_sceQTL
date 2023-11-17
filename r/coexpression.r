options(stringsAsFactors = FALSE, data.table.verbose = FALSE, datatable.showProgress = FALSE, digits = 2)
suppressPackageStartupMessages({
  # Single-cell RNA-seq
  library(SeuratObject)
  library(Seurat)
  library(Matrix)

  # Correlation analysis
  library(lme4)
  library(lmerTest)

  # Co-expression analysis
  library(WGCNA)
  library(hdWGCNA)

  # Functional analysis
  library(enrichR)

  # Data manipulation
  library(magrittr)
  library(tidyverse)
  library(data.table)

  # Plotting
  library(colorRamp2)
  library(cowplot)
  library(patchwork)
  library(RColorBrewer)
  library(ggrepel)
  library(ggpubr)
  library(ggsci)
})

enableWGCNAThreads(nThreads = 4)


# Working directories
proj_dir <- "/home/zzhang/Documents/projects/wp_bcg_eqtl"
celltype_vec <- c("Monocytes", "CD4T", "CD8T", "NK", "B")
condition_vec <- c("T0_LPS.vs.T0_RPMI", "T3m_LPS.vs.T3m_RPMI", "T3m_RPMI.vs.T0_RPMI") #, "T3m_LPS.vs.T0_RPMI")

target_features <- c("CD55", "SPI1", "STAT1", "ESR1", "FOXK2", "EP300", "CEBPA", "KLF5", "TBX3", "HLF", "FOXP1", "RARA")

per_condition <- "RPMI"
per_cell_type <- "Monocytes"
save_token <- paste(per_condition, per_cell_type, sep = ".")

out_dir <- file.path(proj_dir, "outputs/pseudo_bulk/coexpression/300BCG")
setwd(out_dir)
set.seed(31415926)

#
## Load single-cell RNA-seq results
#
wgcna_300bcg <- file.path(out_dir, paste0(save_token, ".rds"))
if (!file.exists(wgcna_300bcg)) {
  pbmc_300bcg <- readRDS(file.path(proj_dir, "inputs/sc_rnaseq/bcg4-0712.rds"))
  pbmc_300bcg <- pbmc_300bcg[, (!pbmc_300bcg@meta.data$clusters1 %in% c("HSP(T)"))] # Removing HSP(T) cells
  new_cell_lvl <- c("Monocytes", "CD4+ T", "CD8+ T", "NK", "B", "Platelet", "mDC", "pDC", "Undefined")
  pbmc_300bcg@meta.data$clusters1 <- factor(pbmc_300bcg@meta.data$clusters1, levels = new_cell_lvl)
  Idents(pbmc_300bcg) <- "clusters1"
  DefaultAssay(pbmc_300bcg) <- "RNA"


  # Setup for WGCNA
  selected_cells <- as.data.frame(pbmc_300bcg@meta.data) %>% dplyr::filter(stim == per_condition, clusters1 == per_cell_type) %>% rownames()
  subobj <- SetupForWGCNA(pbmc_300bcg[, selected_cells], gene_select = "fraction", fraction = 0.05, wgcna_name = "wgcna_300bcg")

  # Construct metacells
  subobj <- MetacellsByGroups(subobj, group.by = c("clusters1", "ids"), min_cells = 50, max_shared = 15, ident.group = "clusters1")
  subobj <- NormalizeMetacells(subobj)

  # Setup expression matrix
  subobj <- SetDatExpr(subobj, group_name = per_cell_type, group.by = "clusters1", assay = "RNA", slot = "data")

  subobj <- TestSoftPowers(subobj, networkType = "signed")
  plot_list <- PlotSoftPowers(subobj)
  p <- wrap_plots(plot_list, ncol = 2)
  save_to <- file.path(out_dir, paste0(save_token, ".soft_power.pdf"))
  ggsave(save_to, p, width = 8, height = 5)

  # Construct co-expression network
  subobj <- ConstructNetwork(
    subobj, soft_power = 14, tom_outdir = paste0(save_token, ".TOM"), overwrite_tom = TRUE, setDatExpr = FALSE, tom_name = "INH",
    TOMType = "unsigned", TOMDenom = "min"
  )

  save_to <- file.path(out_dir, paste0(save_token, ".dengrogram.pdf"))
  pdf(save_to, width = 8, height = 5)
  PlotDendrogram(subobj, main = "300bcg")
  dev.off()

  modules <- GetModules(subobj)
  modules %>% dplyr::filter(gene_name %in% target_features)

  # Module eigengenes and connectivity
  subobj <- ScaleData(subobj, features = VariableFeatures(subobj))

  # compute all MEs in the full single-cell dataset
  subobj <- ModuleEigengenes(subobj, group.by.vars = c("ids"), verbose = FALSE)

  # Compute module connectivity
  subobj <- ModuleConnectivity(subobj, group.by = c("clusters1"), group_name = "Monocytes")

  # Visualize module connectivity
  p <- PlotKMEs(subobj, ncol = 5, n_hubs = 20)
  save_to <- file.path(out_dir, paste0(save_token, ".model_eigengenes.pdf"))
  ggsave(save_to, p, width = 10, height = 8)

  # module network plot
  ModuleNetworkPlot(subobj, outdir = paste0(save_token, ".module_network"))

  modules <- GetModules(subobj)
  hubgenes <- GetHubGenes(subobj, n_hubs = 200)

  # Functional annotation and enrichment
  dbs <- c("GO_Biological_Process_2021","GO_Cellular_Component_2021","GO_Molecular_Function_2021")
  subobj <- RunEnrichr(subobj, max_genes = Inf)
  enrich_df <- GetEnrichrTable(subobj)

  EnrichrBarPlot(subobj, outdir = paste0(save_token, ".enrich_bar_plot"), n_terms = 20, plot_size = c(9, 7), logscale = TRUE)

  saveRDS(subobj, file.path(out_dir, paste0(save_token, ".rds")))
} else {
  subobj <- readRDS(wgcna_300bcg)
}

modules <- GetModules(subobj) %>%
  dplyr::filter(module %in% "blue") %>%
  dplyr::arrange(kME_blue) %>%
  dplyr::mutate(index = row_number())

hightlight_genes <- modules %>% dplyr::filter(gene_name %in% target_features)

p <- ggplot(modules) +
  geom_point(aes(y = index, x = kME_blue), size = 1) +
  geom_point(aes(y = index, x = kME_blue), data = hightlight_genes, size = 2, color = "red") +
  geom_text_repel(aes(y = index, x = kME_blue, label = gene_name), dplyr::filter(modules, gene_name %in% target_features)) +
  theme_classic()
save_to <- file.path(out_dir, paste0(save_token, ".blue.model_eigengenes.pdf"))
ggsave(save_to, plot = p, width = 4, height = 5)



# Check those papers
# https://www.jci.org/articles/view/42869, BCL6 repress EP300 and CD55

# ------- TO BE REMOVED -------
co_exp_genes <- list.files(file.path(proj_dir, "inputs/ADASTRA/TF")) %>% stringr::str_remove("_HUMAN.tsv") %>% unlist() %>% unique() %>% intersect(rownames(pbmc)) %>% c("CD55", .)
co_exp_list <- list()
for (per_ct in c("Monocytes", "B")) {
  for (per_gt in c("TT", "TC", "CC")) {
    run_id <- paste0(per_ct, "_", per_gt)

    # Obtain subset of cells to estimate co-expression
    per_gt_ids <- gt_per_ind %>% dplyr::filter(rs2564978 %in% per_gt) %>% pull(ids)
    sub_pbmc <- pbmc[, pbmc$clusters1 %in% per_ct & pbmc$ids %in% per_gt_ids]

    # Estimate co-expression
    coexp <- CSCORE(sub_pbmc, genes = co_exp_genes)

    # Adjust p-values
    coexp_est <- coexp$est
    coexp_p_value <- coexp$p_value

    p_matrix_bh <- matrix(0, length(co_exp_genes), length(co_exp_genes))
    rownames(p_matrix_bh) <- colnames(p_matrix_bh) <- colnames(coexp_est)
    p_matrix_bh[upper.tri(p_matrix_bh)] <- p.adjust(coexp_p_value[upper.tri(p_matrix_bh)], method = "BH")
    p_matrix_bh <- p_matrix_bh + t(p_matrix_bh)
    coexp_est[p_matrix_bh > 0.05] <- 0

    coexp$p_value_bh_adj <- p_matrix_bh
    coexp$est_bh_masked <- coexp_est

    co_exp_list[[run_id]] <- coexp
  }
}

plot_tab <- lapply(names(co_exp_list), function(pid) {
  est_masked <- co_exp_list[[pid]]$est_bh_masked["CD55", ]
  p_val_adj <- co_exp_list[[pid]]$p_value_bh_adj["CD55", ]
  p_val <- co_exp_list[[pid]]$p_value["CD55", ]

  data.frame(name = pid, TF = names(est_masked), est_masked = est_masked, p_val = p_val, p_val_adj = p_val_adj)
}) %>%
  Reduce(rbind, .) %>%
  tidyr::separate(name, into = c("celltype", "genotype")) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::arrange(est_masked) %>%
  dplyr::mutate(p_val_adj = -log10(p_val_adj), genotype = factor(genotype, levels = c("TT", "TC", "CC"))) %>%
  dplyr::filter(!is.infinite(p_val_adj) & !is.na(p_val_adj)) %>%
  as.data.table()

top_tfs <- plot_tab %>% dplyr::group_by(celltype, genotype) %>% dplyr::slice_max(p_val_adj, n = 10)
p <- ggplot() +
  geom_point(aes(x = est_masked, y = p_val_adj), plot_tab, size = 0.5) +
  geom_point(aes(x = est_masked, y = p_val_adj), top_tfs, color = "red", size = 2, alpha = 0.75) +
  geom_text_repel(aes(x = est_masked, y = p_val_adj, label = TF), top_tfs, min.segment.length = 0.0) +
  facet_nested_wrap(celltype ~ genotype, ncol = 3, scales = "free") +
  labs(x = "Co-expression", y = "-log10(p-value)") +
  theme_classic() +
  theme(legend.position = "top", strip.background = element_blank(), ggh4x.facet.nestline = element_line(colour = "blue"))
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.coexpression.rank_dot_plot.300BCG.pdf"))
ggsave(save_to, plot = p, width = 10, height = 6)


selected_genes <- (co_exp_list$Monocytes_CC$p_value_bh_adj["CD55", ] < 0.05) %>% purrr::keep(~!is.na(.x) && .x) %>% names()
selected_est <- co_exp_list$Monocytes_CC$est_bh_masked[selected_genes, selected_genes]
selected_stat <- co_exp_list$Monocytes_CC$test_stat[selected_genes, selected_genes]

tom <- WGCNA::adjacency.fromSimilarity(selected_est) %>% WGCNA::TOMsimilarity(TOMType = "unsigned", TOMDenom = "mean", indent = 2)
diss_tom <- 1.0 - tom
rownames(diss_tom) <- colnames(diss_tom) <- rownames(selected_est)

gene_tree <- hclust(as.dist(diss_tom), method = "average")
dynamic_mods <- dynamicTreeCut::cutreeDynamic(
  dendro = gene_tree, distM = diss_tom, deepSplit = 3, pamRespectsDendro = FALSE, minClusterSize = 10
)

names(dynamic_mods) <- rownames(diss_tom)
module_list <- lapply(sort(unique(dynamic_mods)), function(i_k) names(which(dynamic_mods == i_k)))

dynamic_colors <- labels2colors(dynamic_mods)
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.coexpression.wgcna_dendro_tree.300BCG.pdf"))
pdf(save_to, width = 10, height = 5)
plotDendroAndColors(gene_tree, dynamic_colors)
dev.off()

plot_tom <- diss_tom ^ 13
diag(plot_tom) <- NA
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.coexpression.wgcna_tom_heatmap.300BCG.pdf"))
pdf(save_to, width = 5, height = 5)
TOMplot(plot_tom, gene_tree, dynamic_colors, main = "Network heatmap plot, all genes")
dev.off()

cd55_module <- which(dynamic_mods == dynamic_mods["CD55"])
cd55_plot_tom <- plot_tom[cd55_module, cd55_module]
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.coexpression.wgcna_tom_heatmap.CD55_cluster.300BCG.pdf"))
pdf(save_to, width = 5, height = 5)
TOMplot(cd55_plot_tom, gene_tree[cd55_module], dynamic_colors[cd55_module], main = "Network heatmap plot, all genes")
dev.off()
