#!/usr/bin/env Rscript
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang@gmail.com, zhenhua.zhang3@helmholtz-hzi.de
# Created: 2022 Jul 02
# Updated: 2023 Mar 02

# Options, packages used in the analysis
options(stringsAsFactors = FALSE, future.globals.maxSize = 10000 * 1024 ^ 2)
library(magrittr)


#' Load eQTL summary statistics
#'
#'@description A utility function to load summary statistics from the disk
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
    for (nm in names(other_info)) cmb_qtltab[nm] <- other_info[nm]
  }

  cmb_qtltab %>% dplyr::filter(dplyr::if_any(dplyr::everything(), ~ !is.na(.x)))
}


proj_dir <- "~/Documents/projects/wp_bcg_eqtl"

# Two model used in the analysis.
mode_vec <- c("normal", "interaction") # Main cell types.
cell_type_vec <- c("Monocytes", "CD4T", "CD8T", "NK", "B") # , "pDC", "mDC")

# All comparisons used in the analysis.
condition_vec <- c("T0_LPS.vs.T0_RPMI", "T3m_LPS.vs.T3m_RPMI", "T3m_RPMI.vs.T0_RPMI")#, "T3m_LPS.vs.T0_RPMI")

qtl_tab <- list()
for (celltype in cell_type_vec) {
  for (mode in mode_vec) {
    if (mode == "normal") {
      run_id <- paste(celltype, "common", sep = "_")
      indir <- file.path(proj_dir, "outputs/pseudo_bulk/summary_statistic", mode, celltype)
      qtl_tab[[run_id]] <- load_eqtl_tab(indir, fdr = 0.1) %>%
        dplyr::filter(!is.na(global_corrected_pValue), p_value <= 0.05) %>% #, 0.1 <= alpha_param, alpha_param <= 500, 0.1 <= beta_param, beta_param <= 500) %>%
        dplyr::select(snp_id, snp_chromosome, snp_position, maf, feature_id, p_value, beta, beta_se, global_corrected_pValue) %>%
        dplyr::mutate(celltype = celltype, condition = "Common")
    } else {
      for (cond in condition_vec) {
        run_id <- paste(celltype, cond, sep = "_")
        indir <- file.path(proj_dir, "outputs/pseudo_bulk/summary_statistic", mode, celltype, cond)
        qtl_tab[[run_id]] <- load_eqtl_tab(indir, fdr = 0.1) %>%
          dplyr::filter(!is.na(global_corrected_pValue), p_value <= 0.05) %>% #, 0.1 <= alpha_param, alpha_param <= 500, 0.1 <= beta_param, beta_param <= 500) %>%
          dplyr::select(snp_id, snp_chromosome, snp_position, maf, feature_id, p_value, beta, beta_se, global_corrected_pValue) %>%
          dplyr::mutate(celltype = celltype, condition = cond)
      }
    }
  }
}

top_qtl_tab <- Reduce(rbind, qtl_tab) %>%
  dplyr::select(snp_id, feature_id) %>%
  dplyr::distinct() %>%
  dplyr::mutate(QTL = paste0(snp_id, "-", feature_id)) %>%
  dplyr::pull(QTL)

all_qtl_tab <- list()
for (celltype in cell_type_vec) {
  for (mode in mode_vec) {
    if (mode == "normal") {
      run_id <- paste(celltype, "common", sep = "_")
      indir <- file.path(proj_dir, "outputs/pseudo_bulk/summary_statistic", mode, celltype)
      all_qtl_tab[[run_id]] <- load_eqtl_tab(indir, fdr = 0.1) %>%
        # dplyr::filter(0.1 <= alpha_param, alpha_param <= 500, 0.1 <= beta_param, beta_param <= 500) %>%
        dplyr::select(snp_id, snp_chromosome, snp_position, maf, feature_id, p_value, beta, beta_se, global_corrected_pValue) %>%
        dplyr::mutate(celltype = celltype, condition = "Common", QTL = paste0(snp_id, "-", feature_id)) %>%
        dplyr::filter(QTL %in% top_qtl_tab)
    } else {
      for (cond in condition_vec) {
        run_id <- paste(celltype, cond, sep = "_")
        indir <- file.path(proj_dir, "outputs/pseudo_bulk/summary_statistic", mode, celltype, cond)
        all_qtl_tab[[run_id]] <- load_eqtl_tab(indir, fdr = 0.1) %>%
          # dplyr::filter(0.1 <= alpha_param, alpha_param <= 500, 0.1 <= beta_param, beta_param <= 500) %>%
          dplyr::select(snp_id, snp_chromosome, snp_position, maf, feature_id, p_value, beta, beta_se, global_corrected_pValue) %>%
          dplyr::mutate(celltype = celltype, condition = cond, QTL = paste0(snp_id, "-", feature_id)) %>%
          dplyr::filter(QTL %in% top_qtl_tab)
      }
    }
  }
}

Reduce(rbind, all_qtl_tab) %>%
  tidyr::pivot_wider(
    names_from = c("celltype", "condition"),
    values_from = c("p_value", "beta", "beta_se", "global_corrected_pValue"),
    values_fn = function(tab) head(tab, 1)
  ) %>%
  dplyr::arrange(snp_chromosome, snp_position, feature_id) %>%
  as.data.frame() %>%
  data.table::fwrite(file.path(proj_dir, "outputs/pseudo_bulk/overview/filtered.all_top_qtl.FDR0.1.v2.csv"))
