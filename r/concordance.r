#!/usr/bin/env Rscript
# File: concordance.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Aug 22, 2023
# Updated: Aug 22, 2023


suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
})

proj_dir <- "~/Documents/projects/wp_bcg_eqtl"

celltype_vec <- c("Monocytes", "CD4T", "CD8T", "NK", "B")
condition_vec <- c("T0_LPS", "T0_RPMI", "T3m_LPS", "T3m_RPMI")
comparison_vec <- c("T0_LPS.vs.T0_RPMI", "T3m_LPS.vs.T3m_RPMI", "T3m_RPMI.vs.T0_RPMI")#, "T3m_LPS.vs.T0_RPMI")

evar_tab <- fread(file.path(proj_dir, "outputs/pseudo_bulk/overview/filtered.eGene_eVariants.FDR0.05.csv"))

# QTL concordance, normal
tar_qtl <- evar_tab %>% dplyr::filter(condition == "Common") %>% dplyr::pull(QTL) %>% unique()
cmp_tab <- combn(celltype_vec, 2) %>% t() %>% as.data.frame() %>% dplyr::rename(x = V1, y = V2)
asso_tab <- lapply(celltype_vec, function(pct) {
  cat("[I]: Loading summary statistic for", pct, "\n")
  file.path(proj_dir, "outputs/pseudo_bulk/summary_statistic/normal", pct, "qtl_results_all.txt") %>%
    data.table::fread(showProgress = FALSE) %>%
    dplyr::filter(p_value < 0.5) %>%
    dplyr::mutate(QTL = paste(feature_id, snp_id, sep = "-")) %>%
    dplyr::filter(QTL %in% tar_qtl) %>%
    dplyr::mutate(celltype = pct, condition = "Common")
  }) %>%
  Reduce(rbind, .)

plot_tab <- apply(cmp_tab, 1, function(vec) {
  dplyr::mutate(asso_tab, z_score = beta / beta_se) %>%
    dplyr::select(feature_id, snp_id, celltype, z_score) %>%
    dplyr::filter(celltype %in% vec) %>%
    tidyr::pivot_wider(names_from = celltype, values_from = z_score) %>%
    dplyr::mutate(Celltype_x = vec[1], Celltype_y = vec[2]) %>%
    dplyr::rename(vec)
}) %>%
  Reduce(rbind, .) %>%
  dplyr::filter(!is.na(x) & !is.na(y)) %>%
  dplyr::mutate(Celltype_x = factor(Celltype_x, levels = celltype_vec), Celltype_y = factor(Celltype_y, levels = celltype_vec))

p <- ggplot(data = plot_tab) +
  geom_point(aes(x = x, y = y), size = 0.1, alpha = 0.5) +
  geom_hline(yintercept = c(1.96, -1.96), linetype = "dotted", color = "red") +
  geom_vline(xintercept = c(1.96, -1.96), linetype = "dotted", color = "red") +
  facet_grid(Celltype_x ~ Celltype_y) +
  labs(x = "Z-score", y = "Z-score") +
  theme_bw()

file.path(proj_dir, "outputs/pseudo_bulk/concordance/normal/qtl_concordance.pdf") %>% ggsave(plot = p, width= 9, height = 9)


# QTL concordance, across interactions per cell type
tar_qtl <- dplyr::filter(evar_tab, condition %in% comparison_vec) %>% dplyr::pull(QTL) %>% unique()
cmp_tab <- combn(comparison_vec, 2) %>% t() %>% as.data.frame() %>% dplyr::rename(x = V1, y = V2)
for (pct in celltype_vec) {
  asso_tab <- lapply(comparison_vec, function(pcd) {
    cat("[I]: Loading summary statistic for", pcd, "in", pct, "\n")
    file.path(proj_dir, "outputs/pseudo_bulk/summary_statistic/interaction", pct, pcd, "qtl_results_all.txt") %>%
      data.table::fread(showProgress = FALSE) %>%
      dplyr::filter(p_value < 0.5) %>%
      dplyr::mutate(QTL = paste(feature_id, snp_id, sep = "-")) %>%
      dplyr::filter(QTL %in% tar_qtl) %>%
      dplyr::mutate(celltype = pct, condition = pcd)
  }) %>%
    Reduce(rbind, .)

  plot_tab <- apply(cmp_tab, 1, function(vec) {
    dplyr::mutate(asso_tab, z_score = beta / beta_se) %>%
      dplyr::select(feature_id, snp_id, condition, z_score) %>%
      dplyr::filter(condition %in% vec) %>%
      tidyr::pivot_wider(names_from = condition, values_from = z_score) %>%
      dplyr::mutate(Condition_x = vec[1], Condition_y = vec[2]) %>%
      dplyr::rename(vec)
  }) %>%
    Reduce(rbind, .) %>%
    dplyr::filter(!is.na(x) & !is.na(y)) %>%
    dplyr::mutate(Condition_x = factor(Condition_x, levels = comparison_vec), Condition_y = factor(Condition_y, levels = comparison_vec))

  p <- ggplot(data = plot_tab) +
    geom_point(aes(x = x, y = y), size = 0.1, alpha = 0.5) +
    geom_hline(yintercept = c(1.96, -1.96), linetype = "dotted", color = "red") +
    geom_vline(xintercept = c(1.96, -1.96), linetype = "dotted", color = "red") +
    facet_grid(Condition_x ~ Condition_y) +
    labs(x = "Z-score", y = "Z-score") +
    theme_bw()
  file.path(proj_dir, "outputs/pseudo_bulk/concordance/interaction", paste(pct, "qtl_concordance.pdf", sep = ".")) %>% ggsave(plot = p, width= 4, height = 4)
}


# QTL concordance, across cell-types per interaction
tar_qtl <- dplyr::filter(evar_tab, condition %in% comparison_vec) %>% dplyr::pull(QTL) %>% unique()
cmp_tab <- combn(celltype_vec, 2) %>% t() %>% as.data.frame() %>% dplyr::rename(x = V1, y = V2)
for (pcd in comparison_vec) {
  asso_tab <- lapply(celltype_vec, function(pct) {
    cat("[I]: Loading summary statistic for", pcd, "in", pct, "\n")
    file.path(proj_dir, "outputs/pseudo_bulk/summary_statistic/interaction", pct, pcd, "qtl_results_all.txt") %>%
      data.table::fread(showProgress = FALSE) %>%
      dplyr::filter(p_value < 0.5) %>%
      dplyr::mutate(QTL = paste(feature_id, snp_id, sep = "-")) %>%
      dplyr::filter(QTL %in% tar_qtl) %>%
      dplyr::mutate(celltype = pct, condition = pcd)
  }) %>%
    Reduce(rbind, .)

  plot_tab <- apply(cmp_tab, 1, function(vec) {
    dplyr::mutate(asso_tab, z_score = beta / beta_se) %>%
      dplyr::select(feature_id, snp_id, celltype, z_score) %>%
      dplyr::filter(celltype %in% vec) %>%
      tidyr::pivot_wider(names_from = celltype, values_from = z_score) %>%
      dplyr::mutate(Celltype_x = vec[1], Celltype_y = vec[2]) %>%
      dplyr::rename(vec)
  }) %>%
    Reduce(rbind, .) %>%
    dplyr::filter(!is.na(x) & !is.na(y)) %>%
    dplyr::mutate(Celltype_x = factor(Celltype_x, levels = celltype_vec), Celltype_y = factor(Celltype_y, levels = celltype_vec))

  p <- ggplot(data = plot_tab) +
    geom_point(aes(x = x, y = y), size = 0.1, alpha = 0.5) +
    geom_hline(yintercept = c(1.96, -1.96), linetype = "dotted", color = "red") +
    geom_vline(xintercept = c(1.96, -1.96), linetype = "dotted", color = "red") +
    facet_grid(Celltype_x ~ Celltype_y) +
    labs(x = "Z-score", y = "Z-score") +
    theme_bw()
  file.path(proj_dir, "outputs/pseudo_bulk/concordance/interaction", paste(pcd, "qtl_concordance.pdf", sep = ".")) %>% ggsave(plot = p, width= 6, height = 6)
}
