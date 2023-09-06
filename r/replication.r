#!/usr/bin/env Rscript
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: May 21, 2023
# Updated: Jul 03, 2023

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
comparison_vec <- c("T0_LPS.vs.T0_RPMI", "T3m_LPS.vs.T3m_RPMI", "T3m_RPMI.vs.T0_RPMI") #, "T3m_LPS.vs.T0_RPMI")
col_pattern_vec <- c("normal" = "_Common", "interaction" = ".vs.")
effect_vec <- c("LPS eff. (T0)", "LPS eff. (T3m)", "BCG eff.")
names(effect_vec) <- comparison_vec

overwrite <- FALSE
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/replication/normal.eqtlgen_gtex.csv")
final_cols <- c("snp_chromosome", "snp_position", "snp_id", "gene_name", "gene_ensembl_id", "celltype", "BCG.assessed_allele", "BCG.beta", "BCG.p_value", "BCG.q_value", "GTEx.assessed_allele", "GTEx.other_allele", "GTEx.p_value", "GTEx.beta", "eqtlGen.assessed_allele", "eqtlGen.other_allele", "eqtlGen.p_value", "eqtlGen.beta")
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

  gtex_whbl_cols <- c("gene_id", "snp_chrom", "snp_pos", "snp_id", "assessed_allele", "other_allele", "p_value" = "pval_nominal", "beta" = "slope")
  pos_to_snp_map <- fread(file.path(proj_dir, "temps/pos_to_rsid.txt")) %>% dplyr::mutate(snp_pos_id = paste0(V1, "_", V2, "_", V3, "_", V4)) %>% dplyr::pull(V5, snp_pos_id)
  gtex_whbl_tab <- fread("/vol/projects/BIIM/resources/GTEx/v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz") %>%
    dplyr::mutate(gene_id = stringr::str_remove_all(gene_id, "\\.[0-9]+$"), variant_id = stringr::str_remove_all(variant_id, "^chr|_b38$")) %>%
    dplyr::mutate(snp_id = pos_to_snp_map[variant_id]) %>%
    dplyr::filter(!is.na(snp_id), snp_id != "") %>%
    tidyr::separate(variant_id, c("snp_chrom", "snp_pos", "assessed_allele", "other_allele"), sep = "_") %>%
    dplyr::select(dplyr::all_of(gtex_whbl_cols)) %>%
    dplyr::filter(!snp_chrom %in% c("X", "Y", "MT")) %>%
    dplyr::mutate(snp_chrom = as.integer(snp_chrom), snp_pos = as.integer(snp_pos))

  comb_tab <- dplyr::left_join(eqtl_tab, eqtl_gen_tab, by = c("gene_id", "snp_id"), suffix = c(".eqtl", "")) %>%
    dplyr::left_join(gtex_whbl_tab, by = c("gene_id", "snp_id"), suffix = c(".eqtlgen", ".gtex"))

  comb_tab %>% data.table::fwrite(save_to)
} else {
  comb_tab <- fread(save_to)
}

# Number of replicated eVariants
plot_tab <- rep_tab %>%
  dplyr::filter(!is.na(p_value.eqtlgen)) %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(
    `eVariants` = n(),
    `Replicated` = sum(eqtlgen.rep),
    `Non-replicated` = `eVariants` - `Replicated`
  ) %>%
  tidyr::pivot_longer(cols = c(`Non-replicated`, `Replicated`), names_to = "Category", values_to = "Counts") %>%
  dplyr::mutate(celltype = factor(celltype, celltype_vec)) %>%
  dplyr::mutate(
    Percent = Counts / `eVariants` * 100,
    label = paste0(Counts, "(", round(Percent, 1), "%)"),
    label = dplyr::if_else(Category == "Non-replicated", "", label)
  ) %>%
  dplyr::mutate(Category = factor(Category, c("Replicated", "Non-replicated")))

p <- plot_tab %>%
  ggplot() +
  geom_bar(aes(x = celltype, y = Percent, fill = Category), stat = "identity", position = "stack") +
  geom_text(aes(x = celltype, y = Percent/2, label = label), position = "stack") +
  scale_fill_npg() +
  labs(x = "Cell type", y = "Percentage") +
  theme_classic() +
  theme(legend.position = "top") +
  coord_flip()

plot_save_to <- file.path(proj_dir, "outputs/pseudo_bulk/replication/normal.replication.by_existing_evariant.bar_plot.pdf")
ggsave(plot_save_to, plot = p, height = 3, width = 5)

# Concordance of eVariants
plot_tab <- rep_tab %>%
  dplyr::filter(!is.na(p_value.eqtlgen)) %>%
  dplyr::mutate(beta.eqtlgen = dplyr::if_else(assessed_allele.eqtl == assessed_allele.eqtlgen, beta.eqtlgen, -beta.eqtlgen)) %>%
  dplyr::mutate(beta.eqtlgen = sign(beta.eqtlgen) * log10(abs(beta.eqtlgen))) %>%
  dplyr::mutate(beta.eqtl = beta.eqtl / beta_se, celltype = factor(celltype, celltype_vec))

p <- plot_tab %>%
  ggplot(aes(x = beta.eqtlgen, y = beta.eqtl)) +
  geom_point(alpha = 0.5, size = 0.5) +
  facet_wrap(~celltype, nrow = 1) +
  geom_vline(xintercept = c(log10(1.96), -log10(1.96)), color = "red", linetype = "dotted") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = c(1.96, -1.96), color = "red", linetype = "dotted") +
  geom_hline(yintercept = 0) +
  labs(x = "-log10(abs(Z-score)) in eQTLGen", y = "Z-score in sc-eQTL") +
  theme_bw()

plot_save_to <- file.path(proj_dir, "outputs/pseudo_bulk/replication/normal.concordance.by_existing_evariant.point_plot.pdf")
ggsave(plot_save_to, plot = p, height = 3, width = 12)




# Plotting
plot_tab <- comb_tab %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(
    n_total = n(),
    n_rep_either = sum(either.replication %in% c("RF", "RS")),
    n_rep_gtex = sum(GTEx.replication %in% c("RF", "RS")), n_rep_eqtlgen = sum(eqtlGen.replication %in% c("RF", "RS")),
    rep_ratio = n_rep_either / n_total * 100, rep_ratio_lab = paste0(n_total, " (" , round(rep_ratio, 1), "%)")
  ) %>%
  dplyr::mutate(celltype = forcats::fct_reorder(celltype, n_total, .desc = TRUE))

plot_tab_long <- plot_tab %>%
  tidyr::pivot_longer(cols = c("n_rep_gtex", "n_rep_eqtlgen"), names_to = "replication", values_to = "n_rep")

p <- ggplot() +
  geom_col(aes(y = celltype, x = n_total), plot_tab, fill = "gray") +
  geom_label(aes(y = celltype, x = 600, label = rep_ratio_lab), plot_tab, position = position_stack(vjust = 0.85), size = 4) +
  geom_col(aes(y = celltype, x = n_rep, fill = replication), plot_tab_long, position = position_stack()) +
  geom_text(aes(y = celltype, x = n_rep, label = n_rep, group = replication), plot_tab_long, position = position_stack(vjust = 0.5), size = 4) +
  labs(x = "Number of eQTL", y = "Cell type", fill = "Replicated in:") +
  # scale_fill_npg(labels = c("GTEx only", "Both", "eQTLGen only")) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 10), legend.position = "top")

plot_save_to <- file.path(proj_dir, "outputs/pseudo_bulk/replication/normal.replication.bar_plot.pdf")
ggsave(plot_save_to, plot = p, height = 3.5, width = 7)
