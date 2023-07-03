#!/usr/bin/env Rscript
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: May 21, 2023
# Updated: Jul 03, 2023

options(stringsAsFactors = FALSE, data.table.verbose = FALSE, error = NULL)
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggsci)
})

proj_dir <- "~/Documents/projects/wp_bcg_eqtl"
time_vec <- c("T0", "T3m")
stim_vec <- c("RPMI", "LPS")
mode_vec <- c("normal", "interaction")
celltype_vec <- c("Monocytes", "CD4T", "CD8T", "NK", "B")
condition_vec <- c("T0_LPS.vs.T0_RPMI", "T3m_LPS.vs.T3m_RPMI", "T3m_RPMI.vs.T0_RPMI") #, "T3m_LPS.vs.T0_RPMI")
col_pattern_vec <- c("normal" = "_Common", "interaction" = ".vs.")
effect_vec <- c("LPS eff. (T0)", "LPS eff. (T3m)", "BCG eff.")
names(effect_vec) <- condition_vec

overwrite <- FALSE
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/replication/normal.eqtlgen.csv")
final_cols <- c("snp_chromosome", "snp_position", "snp_id", "gene_name", "gene_ensembl_id", "celltype", "BCG.assessed_allele", "BCG.beta", "BCG.p_value", "BCG.global_corrected_pValue", "GTEx.assessed_allele", "GTEx.other_allele", "GTEx.p_value", "GTEx.beta", "eqtlGen.assessed_allele", "eqtlGen.other_allele", "eqtlGen.p_value", "eqtlGen.beta")
if (!file.exists(save_to) || overwrite) {
  eqtl_gen_cols <- c(
    "gene_id" = "Gene", "snp_chrom" = "SNPChr", "snp_pos" = "SNPPos", "snp_id" = "SNP", "assessed_allele" = "AssessedAllele",
    "other_allele" = "OtherAllele", "p_value" = "Pvalue", "beta" = "Zscore" # The real meaning of beta is Zscore.
  )

  eqtl_gen_tab <- fread("/vol/projects/BIIM/resources/eqtlGen/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded-pValue0.05.txt") %>%
    dplyr::select(dplyr::all_of(eqtl_gen_cols))

  gtex_whbl_cols <- c("gene_id", "snp_chrom", "snp_pos", "assessed_allele" = "snp_alt_allele", "other_allele" = "snp_ref_allele", "p_value" = "pval_nominal", "beta" = "slope")
  gtex_whbl_tab <- fread("/vol/projects/BIIM/resources/GTEx/v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz") %>%
    dplyr::mutate(gene_id = stringr::str_remove_all(gene_id, "\\.[0-9]+$"), variant_id = stringr::str_remove_all(variant_id, "^chr|_b38$")) %>%
    tidyr::separate(variant_id, c("snp_chrom", "snp_pos", "snp_ref_allele", "snp_alt_allele"), sep = "_") %>%
    dplyr::select(dplyr::all_of(gtex_whbl_cols)) %>%
    dplyr::mutate(snp_id = "")

  sceqtl <- file.path(proj_dir, "outputs/pseudo_bulk/summary_statistic/normal", celltype_vec) %>%
    lapply(function(pp) {
      celltype <- basename(pp)
      cat(celltype, "\n")
      sceqtl_tab <- fread(file.path(pp, "top_qtl_results_all_FDR0.05.txt")) %>%
        dplyr::select(snp_id, gene_name, snp_chromosome, snp_position, p_value, beta, ensembl_gene_id, assessed_allele, global_corrected_pValue) %>%
        dplyr::group_by(snp_id, gene_name) %>%
        dplyr::summarise(replication = {
          per_snp_id <- dplyr::cur_group()$snp_id
          per_snp_chrom <- dplyr::cur_data()$snp_chromosome
          per_snp_pos <- dplyr::cur_data()$snp_position
          per_gene_id <- dplyr::cur_data()$ensembl_gene_id

          r3 <- dplyr::cur_data() %>%
            dplyr::mutate(other_allele = "") %>%
            dplyr::select(beta, assessed_allele, other_allele, p_value, global_corrected_pValue) %>%
            dplyr::mutate(data_source = "BCG")

          r1 <- dplyr::filter(eqtl_gen_tab, snp_chrom == per_snp_chrom, snp_id == per_snp_id, gene_id == per_gene_id)
          if (nrow(r1) > 0) {
            r3 <- dplyr::select(r1, beta, assessed_allele, other_allele, p_value) %>%
              dplyr::mutate(global_corrected_pValue = NA, data_source = "eqtlGen") %>%
              rbind(r3)
          }

          r2 <- dplyr::filter(gtex_whbl_tab, snp_chrom == per_snp_chrom, snp_pos == per_snp_pos, gene_id == per_gene_id)
          if (nrow(r2) > 0) {
            r3 <- dplyr::select(r2, beta, assessed_allele, other_allele, p_value) %>%
              dplyr::mutate(global_corrected_pValue = NA, data_source = "GTEx") %>%
              rbind(r3)
          }

          dplyr::mutate(r3, snp_pos = per_snp_pos, snp_chrom = per_snp_chrom, gene_id = per_gene_id)
      }) %>%
        as.data.table() %>%
        dplyr::mutate(celltype = celltype)
    }) %>%
    Reduce(rbind, .)

  sceqtl %>%
    dplyr::rename(gene_ensembl_id = replication.gene_id, snp_chromosome = replication.snp_chrom, snp_position = replication.snp_pos, data_source = replication.data_source) %>%
    dplyr::relocate(snp_chromosome, snp_position, snp_id, gene_name, gene_ensembl_id, celltype, everything()) %>%
    tidyr::pivot_wider(names_from = data_source, values_from = dplyr::starts_with("replication"), names_sep = ".", names_glue = "{data_source}.{.value}") %>%
    dplyr::rename_with(~stringr::str_remove_all(., "replication\\."), dplyr::contains("replication")) %>%
    dplyr::select(-GTEx.global_corrected_pValue, -eqtlGen.global_corrected_pValue) %>%
    dplyr::arrange(snp_chromosome, snp_position) %>%
    dplyr::select(dplyr::all_of(final_cols)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(GTEx.replicated = {
      !(is.na(GTEx.beta) && is.na(GTEx.p_value)) &&
      ((BCG.assessed_allele == GTEx.assessed_allele && BCG.beta * GTEx.beta > 0) || (BCG.assessed_allele == GTEx.assessed_allele && BCG.beta * GTEx.beta > 0))
    }, eqtlGen.replicated = {
      !(is.na(eqtlGen.beta) && is.na(eqtlGen.p_value)) &&
        ((BCG.assessed_allele == eqtlGen.assessed_allele && BCG.beta * eqtlGen.beta > 0) || (BCG.assessed_allele == eqtlGen.assessed_allele && BCG.beta * eqtlGen.beta > 0))
    }, both.replicated = GTEx.replicated & eqtlGen.replicated, either.replicated = GTEx.replicated | eqtlGen.replicated) %>%
    fwrite(save_to)
} else {
  sceqtl <- fread(save_to)
}

# Plotting
plot_tab <- sceqtl %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(
    n_total = n(),
    n_rep_either = sum(either.replicated), n_rep_both = sum(both.replicated),
    n_rep_gtex = sum(GTEx.replicated), n_rep_eqtlgen = sum(eqtlGen.replicated),
    rep_ratio = n_rep_either / n_total * 100, rep_ratio_lab = paste0(n_total, " (" , round(rep_ratio, 1), "%)")
  ) %>%
  dplyr::mutate(celltype = forcats::fct_reorder(celltype, n_total, .desc = TRUE))

plot_tab_long <- plot_tab %>%
  dplyr::mutate(n_rep_gtex = n_rep_gtex - n_rep_both, n_rep_eqtlgen = n_rep_eqtlgen - n_rep_both) %>%
  tidyr::pivot_longer(cols = c("n_rep_both", "n_rep_gtex", "n_rep_eqtlgen"), names_to = "replication", values_to = "n_rep") %>%
  dplyr::mutate(replication = factor(replication, levels = c("n_rep_gtex", "n_rep_both", "n_rep_eqtlgen")))

p <- ggplot() +
  geom_col(aes(y = celltype, x = n_total), plot_tab, fill = "gray") +
  geom_label(aes(y = celltype, x = 600, label = rep_ratio_lab), plot_tab, position = position_stack(vjust = 0.85), size = 4) +
  geom_col(aes(y = celltype, x = n_rep, fill = replication), plot_tab_long, position = position_stack()) +
  geom_text(aes(y = celltype, x = n_rep, label = n_rep, group = replication), plot_tab_long, position = position_stack(vjust = 0.5), size = 4) +
  labs(x = "Number of eQTL", y = "Cell type", fill = "Replicated in:") +
  scale_fill_npg(labels = c("GTEx only", "Both", "eQTLGen only")) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 10), legend.position = "top")

plot_save_to <- file.path(proj_dir, "outputs/pseudo_bulk/replication/normal.replication.bar_plot.pdf")
ggsave(plot_save_to, plot = p, height = 3.5, width = 7)
