#!/usr/bin/env Rscript
# File: concordance.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Aug 22, 2023
# Updated: Aug 22, 2023

# Estimate the concordance

options(stringsAsFactors = FALSE, datatable.verbose = FALSE, datatable.showProgress = FALSE, error = NULL)
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggrepel)
  library(ggsci)
})

proj_dir <- "~/Documents/projects/wp_bcg_eqtl"

time_vec <- c("T0", "T3m")
stim_vec <- c("RPMI", "LPS")
mode_vec <- c("normal", "interaction")
celltype_vec <- c("Monocytes", "CD4T", "CD8T", "NK", "B")
condition_vec <- c("T0_LPS", "T0_RPMI", "T3m_LPS", "T3m_RPMI")
comparison_vec <- c("T0_LPS.vs.T0_RPMI", "T3m_LPS.vs.T3m_RPMI", "T3m_RPMI.vs.T0_RPMI") #, "T3m_LPS.vs.T0_RPMI")
col_pattern_vec <- c("normal" = "_Common", "interaction" = ".vs.")
effect_vec <- c("LPS eff. (T0)", "LPS eff. (T3m)", "BCG eff.")
names(effect_vec) <- comparison_vec

#
## Replication in eqtlGen, 1M_bloodNL
#
overwrite <- FALSE
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/replication/normal.eqtlgen.csv")
if (!file.exists(save_to) || overwrite) {
  eqtl_gen_cols <- c(
    "gene_id" = "Gene", "snp_chrom" = "SNPChr", "snp_pos" = "SNPPos", "snp_id" = "SNP", "assessed_allele" = "AssessedAllele",
    "other_allele" = "OtherAllele", "p_value" = "Pvalue", "beta" = "Zscore" # The real meaning of beta is Zscore.
  )
  eqtl_gen_tab <- fread("/vol/projects/BIIM/resources/eqtlGen/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz") %>%
    dplyr::select(dplyr::all_of(eqtl_gen_cols)) %>%
    dplyr::filter(p_value < 0.1)

  eqtl_cols <- c(
    "snp_id", "snp_chrom" = "snp_chromosome", "snp_pos" = "snp_position", "assessed_allele", "maf",
    "gene_symbol" = "feature_id", "gene_id" = "ensembl_gene_id", "p_value", "beta", "beta_se",
    "empirical_feature_p_value", "q_value", "celltype", "condition"
  )
  eqtl_tab <- fread(file.path(proj_dir, "outputs/pseudo_bulk/overview/filtered.eGene_eVariants.FDR0.05.csv")) %>%
    dplyr::filter(condition == "common") %>%
    dplyr::select(dplyr::all_of(eqtl_cols))

  comb_tab <- dplyr::left_join(eqtl_tab, eqtl_gen_tab, by = c("gene_id", "snp_id"), suffix = c(".eqtl", ".eqtlgen"))

  comb_tab %>% data.table::fwrite(save_to)
} else {
  comb_tab <- fread(save_to)
}


# Number of replicated eVariants
plot_tab <- comb_tab %>%
  dplyr::mutate(eqtlgen.rep = !is.na(p_value.eqtlgen) & (beta.eqtl * beta.eqtlgen > 0 & assessed_allele.eqtlgen == assessed_allele.eqtl) | (beta.eqtl * beta.eqtlgen < 0 & assessed_allele.eqtlgen != assessed_allele.eqtl)) %>%
  dplyr::filter(!is.na(p_value.eqtlgen)) %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(`eVariants` = n(), `Replicated` = sum(eqtlgen.rep), `Non-replicated` = `eVariants` - `Replicated`) %>%
  tidyr::pivot_longer(cols = c(`Non-replicated`, `Replicated`), names_to = "Category", values_to = "Counts") %>%
  dplyr::mutate(celltype = factor(celltype, celltype_vec)) %>%
  dplyr::mutate(
    Percent = Counts / `eVariants` * 100,
    label = paste0(formatC(Counts, big.mark = ",", format = "d", width = 6), " (", round(Percent, 1), "%)"),
    label = dplyr::if_else(Category == "Non-replicated", "", label)
  ) %>%
  dplyr::mutate(Category = factor(Category, c("Replicated", "Non-replicated")))

p <- plot_tab %>%
  ggplot() +
  geom_bar(aes(y = celltype, x = Percent, fill = Category), stat = "identity", position = "stack") +
  geom_text(aes(y = celltype, x = 25, label = label), position = "stack") +
  scale_fill_npg() +
  labs(y = "Cell type", x = "Percentage") +
  theme_classic() +
  theme(legend.position = "top")

plot_save_to <- file.path(proj_dir, "outputs/pseudo_bulk/replication/normal.replication.by_existing_evariant.bar_plot.pdf")
ggsave(plot_save_to, plot = p, height = 3, width = 5)


# Concordance of eVariants
plot_tab <- comb_tab %>%
  dplyr::filter(!is.na(p_value.eqtlgen)) %>%
  dplyr::mutate(beta.eqtlgen = dplyr::if_else(assessed_allele.eqtl == assessed_allele.eqtlgen, beta.eqtlgen, -beta.eqtlgen)) %>%
  dplyr::mutate(beta.eqtlgen = sign(beta.eqtlgen) * log10(abs(beta.eqtlgen))) %>%
  dplyr::mutate(beta.eqtl = beta.eqtl / beta_se, celltype = factor(celltype, celltype_vec))

p <- plot_tab %>%
  ggplot(aes(y = beta.eqtlgen, x = beta.eqtl)) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_hline(yintercept = c(log10(1.96), -log10(1.96)), color = "red", linetype = "dotted") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = c(1.96, -1.96), color = "red", linetype = "dotted") +
  geom_vline(xintercept = 0) +
  labs(y = "Sign(Z-score) * Log10(abs(Z-score)) in eQTLGen", x = "Z-score in sc-eQTL") +
  facet_wrap(~ celltype, ncol = 1) +
  theme_bw()

plot_save_to <- file.path(proj_dir, "outputs/pseudo_bulk/replication/normal.concordance.by_existing_evariant.point_plot.pdf")
ggsave(plot_save_to, plot = p, height = 12, width = 3)


# Replication rate
plot_tab <- comb_tab %>%
  dplyr::mutate(eqtlgen.rep = !is.na(p_value.eqtlgen) & (beta.eqtl * beta.eqtlgen > 0 & assessed_allele.eqtlgen == assessed_allele.eqtl) | (beta.eqtl * beta.eqtlgen < 0 & assessed_allele.eqtlgen != assessed_allele.eqtl)) %>%
  dplyr::filter(!is.na(p_value.eqtlgen)) %>%
  dplyr::group_by(celltype, gene_id) %>%
  dplyr::summarise(eVariants = n(), Category = dplyr::if_else(sum(eqtlgen.rep) > 0, "Replicated", "Non-replicated")) %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(n_total = n(), Yes = sum(Category == "Replicated"), No = sum(Category == "Non-replicated")) %>%
  dplyr::mutate(celltype = factor(celltype, celltype_vec)) %>%
  tidyr::pivot_longer(c(Yes, No), names_to = "Category", values_to = "Counts") %>%
  dplyr::mutate(
    Percentage = as.double(Counts) / n_total * 100,
    Label = paste0(formatC(Counts, big.mark = ",", format = "d", width = 6), " (", round(Percentage, 1), "%)"),
    Category = factor(Category, c("Yes", "No"))
  )

p <- ggplot(data = plot_tab) +
  geom_col(aes(x = celltype, y = Counts, fill = Category), position = position_dodge2(width = 1)) +
  geom_text(aes(x = celltype, y = 0, label = Label), hjust = 0, position = position_dodge2(width = 1)) +
  labs(x = "Cell type", y = "Number of eGenes", fill = "Replicated in eQTLGen") +
  scale_fill_npg() +
  theme_classic() +
  theme(legend.position = "top") +
  coord_flip()

plot_save_to <- file.path(proj_dir, "outputs/pseudo_bulk/replication/normal.replication.by_existing_egene.bar_plot.pdf")
ggsave(plot_save_to, plot = p, height = 4, width = 6)



#
## Concordance
#
evar_tab <- fread(file.path(proj_dir, "outputs/pseudo_bulk/overview/filtered.eGene_eVariants.FDR0.05.csv"))

# QTL concordance, normal
tar_qtl <- evar_tab %>% dplyr::filter(condition == "common") %>% dplyr::pull(QTL) %>% unique()
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
