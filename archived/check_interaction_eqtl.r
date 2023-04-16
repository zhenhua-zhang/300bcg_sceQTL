#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(lme4)

#' Loading eQTL results.
load_eqtl_tab <- function(dir, fdr = 0.05, fdr_col = "global_corrected_pValue", other_info = NULL) {
  joinby <- c(
    "snp_id", "ensembl_gene_id", "feature_start", "feature_end",
    "feature_chromosome", "feature_id", "gene_name", "feature_strand",
    "n_samples", "n_e_samples", "snp_chromosome", "snp_position",
    "assessed_allele"
  )

  discard <- c(
    "p_value", "beta", "beta_se", "empirical_feature_p_value", "alpha_param",
    "beta_param", "call_rate", "maf", "hwe_p"
  )

  top_qtl_fpath <- file.path(dir, "top_qtl_results_all_FDR.txt")
  top_qtltab <- data.table::fread(top_qtl_fpath, data.table = FALSE) %>%
    dplyr::select(-dplyr::one_of(discard)) %>%
    (function(dat) dat[dat[, fdr_col] < fdr, ])

  all_qtl_fpath <- file.path(dir, "qtl_results_all.txt")
  all_qtltab <- data.table::fread(all_qtl_fpath, data.table = FALSE)

  cmb_qtltab <- dplyr::right_join(top_qtltab, all_qtltab, by = joinby) %>%
    dplyr::mutate(QTL = stringr::str_c(gene_name, snp_id, sep = "-"))

  if (!is.null(other_info) && !is.null(names(other_info))) {
    for (nm in names(other_info)) { cmb_qtltab[nm] <- other_info[nm] }
  }

  cmb_qtltab %>% dplyr::filter(dplyr::if_any(dplyr::everything(), ~ !is.na(.x)))
}


proj_dir <- "~/Documents/projects/wp_bcg_eqtl"

# Code block to pick up a candidate gene.
if (FALSE) {
  cd4_com_eqtl_file <- "~/Documents/projects/wp_bcg_eqtl/outputs/pseudo_bulk/normal/CD4T"
  cd4_com_tab <- load_eqtl_tab(cd4_com_eqtl_file) %>% select(gene_name, snp_id, p_value, global_corrected_pValue) %>% dplyr::filter(p_value < 0.001)

  cd8_com_eqtl_file <- "~/Documents/projects/wp_bcg_eqtl/outputs/pseudo_bulk/normal/CD8T"
  cd8_com_tab <- load_eqtl_tab(cd8_com_eqtl_file) %>% select(gene_name, snp_id, p_value, global_corrected_pValue) %>% dplyr::filter(p_value < 0.001)

  cd4_iter_eqtl_file <- "~/Documents/projects/wp_bcg_eqtl/outputs/pseudo_bulk/interaction/CD4T/T3m_RPMI.vs.T0_RPMI"
  cd4_iter_tab <- load_eqtl_tab(cd4_iter_eqtl_file) %>% select(gene_name, snp_id, p_value, global_corrected_pValue) %>% dplyr::filter(p_value < 0.001)

  cd8_iter_eqtl_file <- "~/Documents/projects/wp_bcg_eqtl/outputs/pseudo_bulk/interaction/CD8T/T3m_RPMI.vs.T0_RPMI"
  cd8_iter_tab <- load_eqtl_tab(cd8_iter_eqtl_file) %>% select(gene_name, snp_id, p_value, global_corrected_pValue) %>% dplyr::filter(p_value < 0.001)

  cmb_tab <- reduce(list(cd4_com_tab, cd8_com_tab, cd4_iter_tab, cd8_iter_tab), full_join, by = c("gene_name", "snp_id"))

  colnames(cmb_tab) <- c("gene_name", "snp_id",
    paste("cd4_com", c("p_value", "global_corrected_pValue"), sep = "."),
    paste("cd8_com", c("p_value", "global_corrected_pValue"), sep = "."),
    paste("cd4_iter", c("p_value", "global_corrected_pValue"), sep = "."),
    paste("cd8_iter", c("p_value", "global_corrected_pValue"), sep = "."))


  cmb_tab %>% filter(cd4_com.p_value < 5e-2, cd8_com.p_value < 5e-2, cd4_iter.p_value < 5e-2 | cd8_iter.p_value < 5e-2) %>% arrange(cd4_iter.p_value) %>% head()
}

# GSTM3-rs1537236, CD4 common (1.861849e-05), CD8 common (3.370226e-06), CD4 interaction (0.0002773551), CD8 interaction (NA)

geno_tab <- fread("/vol/projects/zzhang/projects/wp_bcg_eqtl/temps/rs1537236.txt") %>%
  pivot_longer(cols = everything(), values_to = "genotypes", names_to = "genotype_id") %>%
  mutate(genotypes = case_when(genotypes == "0|0" ~ 0, genotypes %in% c("1|0", "0|1") ~ 1, genotypes == "1|1" ~ 2))

# geno_tab
#  genotype_id genotypes
#  300BCG044           2
#  300BCG045           2
#  300BCG046           1
#  300BCG047           1
#  300BCG049           0
#  300BCG055           2

# Common eQTL
base_dir <- "/vol/projects/zzhang/projects/wp_bcg_eqtl/inputs/pseudo_bulk"
for (celltype in c("CD8T", "CD4T")) {
  in_dir <- file.path(base_dir, celltype, "T3m_RPMI.vs.T0_RPMI")

  gtop_tab <- fread(file.path(in_dir, "sample_mapping.tsv"), header = FALSE) %>%
    rename(c("genotype_id" = "V1", "sample_id" = "V2"))

  pheno_tab <- fread(file.path(in_dir, "phenotypes.tsv")) %>%
    filter(feature_id == "GSTM3") %>%
    select(-feature_id) %>%
    pivot_longer(everything(), names_to = "sample_id", values_to = "expression")
# pheno_tab example
#     sample_id              expression
#   1 CD8T_T0_RPMI_300BCG044     0.437
#   2 CD8T_T0_RPMI_300BCG045     0.542
#   3 CD8T_T0_RPMI_300BCG046     0.518
#   4 CD8T_T0_RPMI_300BCG047     0.437
#   5 CD8T_T0_RPMI_300BCG049     0.0656

  covar_tab <- fread(file.path(in_dir, "covariates_wpeer.tsv")) # All variables were treated as numeric, which is the same as limix_qtl
# covar_tab example
#                  sample_id time stim age gender
#   1 CD8T_T0_RPMI_300BCG044    0    0  23      1
#   2 CD8T_T0_RPMI_300BCG055    0    0  24      1
#   3 CD8T_T0_RPMI_300BCG061    0    0  23      0
#   4 CD8T_T0_RPMI_300BCG109    0    0  35      0
#   5 CD8T_T0_RPMI_300BCG123    0    0  24      1
#   6 CD8T_T0_RPMI_300BCG214    0    0  19      1

  per_reg_tab <- left_join(
    left_join(pheno_tab, covar_tab, by = "sample_id"),
    left_join(geno_tab, gtop_tab, by = "genotype_id"),
    by = "sample_id") %>%
    as.data.frame()

  # Interaction model
  itm <- lm(expression ~ time + age + gender + genotypes + time:genotypes, per_reg_tab)
  itc <- coef(itm) # coefficients of interaction models.

  common_vars <- c("time", "age", "gender", "genotypes")
  tmp_reg_tab <- per_reg_tab[, common_vars] %>% mutate(time = as.numeric(time), gender = as.numeric(gender))
  it_vals <- fitted(itm) + residuals(itm) - rowSums(tmp_reg_tab * itc[common_vars]) - itc["(Intercept)"]
  per_reg_tab <- mutate(per_reg_tab, exp_corrected = it_vals, interaction = paste(genotypes, "x", if_else(time == 0, "T0", "T3m")))

  gbp <- ggplot(per_reg_tab, aes(x = as.factor(interaction), y = exp_corrected, fill = as.factor(genotypes))) +
    geom_boxplot(width = 0.8) +
    geom_violin(width = 0.6, alpha = 0.5) +
    geom_point(size = 1) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = "Genotype", y = "Expression")

  save_to <- file.path(proj_dir, "temps", paste0(celltype, "-GSTM3-rs1537236-interaction.pdf"))
  ggsave(save_to, plot = gbp, width = 8, height = 5)


  # Common model
  cmm <- lm(expression ~ time + age + gender + genotypes, per_reg_tab)
  cmc <- coef(cmm) # coefficients of interaction models.

  common_vars <- c("time", "age", "gender")
  geno_vals <- fitted(cmm) + residuals(cmm) - rowSums(per_reg_tab[, common_vars] * cmc[common_vars]) - cmc["(Intercept)"]
  per_reg_tab["exp_corrected"] <- geno_vals # Genotypes interacted with time
  gbp <- ggplot(per_reg_tab, aes(x = as.factor(genotypes), y = exp_corrected, fill = as.factor(genotypes))) +
    geom_boxplot(width = 0.8) +
    geom_violin(width = 0.6, alpha = 0.5) +
    geom_point(size = 1) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = "Genotype", y = "Expression")

  save_to <- file.path(proj_dir, "temps", paste0(celltype, "-GSTM3-rs1537236-common.pdf"))
  ggsave(save_to, plot = gbp, width = 5, height = 5)
}
