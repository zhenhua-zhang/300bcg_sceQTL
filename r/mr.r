#!/usr/bin/env Rscript

#' Mendelian randomization analysis by R/TwoSampleMR
NULL

# FIXME: While the MR-Egger method is a worthwhile sensitivity analysis for detecting violations of
# the instrumental variable assumptions, there are several reasons why causal estimates from the
# MR-Egger method may be biased and have inflated Type 1 error rates in practice, including violations
# of the InSIDE assumption and the influence of outlying variants.

# Options, packages used in the analysis
library(magrittr)
library(tidyverse)
library(data.table)
library(TwoSampleMR)
options(stringsAsFactors = FALSE, future.globals.maxSize = 10000 * 1024^2)


#' Mendelian randomization analysis
mr_test <- function(hm_dat, min_exp_pval = 5e-8, min_snps = 5, prune = FALSE) {
  # Do prune
  if (prune) hm_dat <- TwoSampleMR::power_prune(hm_dat)

  # Filter out SNPs not suitable for MR test and keep significant exposure SNPs
  if ("clump_keep" %in% colnames(hm_dat)) hm_dat <- dplyr::filter(hm_dat, mr_keep, clump_keep, pval.exposure < min_exp_pval)
  else hm_dat <- dplyr::filter(hm_dat, mr_keep, pval.exposure < min_exp_pval)

  # A list to store results
  mr_res <- list(hm_dat = NULL, mr = NULL, mr_single_snp = NULL, mr_pleio = NULL, mr_hetero = NULL)
  if (nrow(hm_dat)) {
    suppressMessages({
      mr_res[["mr"]] <- TwoSampleMR::mr(hm_dat) # Basic MR analysis
      mr_res[["mr_single_snp"]] <- TwoSampleMR::mr_singlesnp(hm_dat) # Single SNP test
      mr_res[["mr_pleio"]] <- TwoSampleMR::mr_pleiotropy_test(hm_dat) # Horizontal pleiotropy
      mr_res[["mr_hetero"]] <- TwoSampleMR::mr_heterogeneity(hm_dat) # Check heterogeneity
    })
  }

  return(mr_res)
}


#' Plot MR for single SNPs
#'
#' @description Function to plot Figure x
plot_mrres <- function(mrpath, save_to = "./", override = TRUE) {
  hm_tab <- data.table::fread(file.path(mrpath, "harmonized_data.csv"))
  mr_all_tab <- data.table::fread(file.path(mrpath, "mendelian_randomization_test.csv"))

  mr_persnp_tab <- data.table::fread(file.path(mrpath, "mendelian_randomization_test_persnp.csv"))
  mr_persnp_tab %>%
    dplyr::group_by(exposure, outcome) %>%
    dplyr::summarise(mr_figs = {
      per_exposure <- dplyr::cur_group()$exposure
      per_outcome <- dplyr::cur_group()$outcome

      per_mr_persnp <- dplyr::cur_data_all() %>% dplyr::filter(!duplicated(SNP)) %>% as.data.frame()

      all_mr_sig <- per_mr_persnp %>% dplyr::filter(stringr::str_detect(SNP, "^All "), !is.na(p)) %>%
        dplyr::mutate(is_sig = p < 0.05) %>% dplyr::select(is_sig) %>% unlist() %>% sum()
      any_snp_sig <- per_mr_persnp %>% dplyr::filter(!stringr::str_detect(SNP, "^All "), !is.na(p)) %>%
        dplyr::mutate(is_sig = p < 0.05) %>% dplyr::select(is_sig) %>% unlist() %>% sum()

      to_plot <- !is.na(all_mr_sig) && all_mr_sig >= 1 && !is.na(any_snp_sig) && any_snp_sig > 0
      if (to_plot) {
        per_outcome_xx <- stringr::str_remove_all(per_outcome, "[|]+ id:") %>% stringr::str_replace_all("[ ]+", "_")

        # Forest plot
        fig_path <- file.path(save_to, paste("mr_forest", per_exposure, per_outcome_xx, "pdf", sep = "."))
        if (!file.exists(fig_path) || override) {
          p <- TwoSampleMR::mr_forest_plot(per_mr_persnp)
          ggsave(fig_path, plot = p[[1]], width = 7, height = 7)
        }

        # MR plot
        fig_path <- file.path(save_to, paste("mr_scatter", per_exposure, per_outcome_xx, "pdf", sep = "."))
        if (!file.exists(fig_path) || override) {
          per_mr_all <- mr_all_tab %>% dplyr::filter(exposure == per_exposure, outcome == per_outcome) %>% as.data.frame()
          per_hm_tab <- hm_tab %>% dplyr::filter(exposure == per_exposure, outcome == per_outcome, SNP %in% per_mr_persnp$SNP) %>% as.data.frame()

          p <- TwoSampleMR::mr_scatter_plot(per_mr_all, per_hm_tab)
          ggsave(fig_path, plot = p[[1]], width = 7, height = 7)
        }
      }

      NULL
    })

  return(NULL)
}


#
## Main steps.
#
proj_dir <- "/home/zzhang/Documents/projects/wp_bcg_eqtl"
in_dir <- file.path(proj_dir, "outputs/pseudo_bulk")

# Two model used in the analysis.
mode_vec <- c("normal", "interaction")

# Main cell types.
cell_type_vec <- c("Monocytes", "CD4T", "CD8T", "NK", "B") # , "pDC", "mDC")

# All comparisons used in the analysis.
condition_vec <- c("T0_LPS.vs.T0_RPMI", "T3m_LPS.vs.T0_RPMI", "T3m_LPS.vs.T3m_RPMI", "T3m_RPMI.vs.T0_RPMI")


# Per loop per eQTL-mapping-run
for (mode in mode_vec) {
  for (cell_type in cell_type_vec) {
    base_dir <- file.path(proj_dir, "outputs/pseudo_bulk/outcomes", mode)

    if (mode == "normal") {
      wk_dir <- file.path(base_dir, paste(mode, cell_type, sep = "_"))

      hm_dat_path <- file.path(wk_dir, "harmonized_data.csv")
      hm_dat <- data.table::fread(hm_dat_path)

      mr_save_to <- file.path(wk_dir, "mendelian_randomization_test.csv")
      mr_ps_save_to <- file.path(wk_dir, "mendelian_randomization_test_persnp.csv")
      if (!file.exists(mr_save_to) || !file.exists(mr_ps_save_to)) {
        cat("[I]: Estimating causality by TwoSampleMR ...\n")
        mr_results <- mr_test(hm_dat, min_exp_pval = 5e-6)
        data.table::fwrite(mr_results$mr, mr_save_to, verbose = FALSE, showProgress = FALSE)
        data.table::fwrite(mr_results$mr_single_snp, mr_ps_save_to, verbose = FALSE, showProgress = FALSE)
      }

      cat("[I]: Plotting scatter and forest plot for the MR results ...\n")
      mrp_saveto <- file.path(wk_dir, "mr_plots")
      if (!dir.exists(mrp_saveto)) { dir.create(mrp_saveto, recursive = TRUE) }
      .tmp <- plot_mrres(wk_dir, save_to = mrp_saveto)
    } else {
      for (condition in condition_vec) {
        wk_dir <- file.path(base_dir, paste(mode, cell_type, sep = "_"), condition)

        hm_dat_path <- file.path(wk_dir, "harmonized_data.csv")
        hm_dat <- data.table::fread(hm_dat_path)

        mr_save_to <- file.path(wk_dir, "mendelian_randomization_test.csv")
        mr_ps_save_to <- file.path(wk_dir, "mendelian_randomization_test_persnp.csv")
        if (!file.exists(mr_save_to) || !file.exists(mr_ps_save_to)) {
          cat("[I]: Estimating causality by TwoSampleMR ...\n")
          mr_results <- mr_test(hm_dat, min_exp_pval = 5e-6)
          data.table::fwrite(mr_results$mr, mr_save_to, verbose = FALSE, showProgress = FALSE)
          data.table::fwrite(mr_results$mr_single_snp, mr_ps_save_to, verbose = FALSE, showProgress = FALSE)
        }

        mrp_saveto <- file.path(wk_dir, "mr_plots")
        if (!dir.exists(mrp_saveto)) dir.create(mrp_saveto, recursive = TRUE)
        .tmp <- plot_mrres(wk_dir, save_to = mrp_saveto)
      }
    }
  }
}


# Code to plot example MR results in for figures.
if (FALSE) {
  sub_hm_dat <- hm_dat %>% dplyr::filter(exposure == "ADCY3", outcome == "BodyMassIndex", mr_keep, clump_keep)
  mr_methods <- c("mr_weighted_mode", "mr_ivw", "mr_simple_mode")
  mr_res <- TwoSampleMR::mr(sub_hm_dat, method_list = mr_methods) # Basic MR analysis
  mr_single_snp <- TwoSampleMR::mr_singlesnp(sub_hm_dat, all_method = mr_methods) # Single SNP test
  mr_hetero <- TwoSampleMR::mr_heterogeneity(sub_hm_dat) # Check heterogeneity
  mr_pleio <- TwoSampleMR::mr_pleiotropy_test(sub_hm_dat) # Horizontal pleiotropy

  p <- TwoSampleMR::mr_forest_plot(mr_single_snp)
  ggsave("mr_forest_plot-ADCY3-BMI.pdf", plot = p[[1]], width = 6, height = 3)

  p <- TwoSampleMR::mr_scatter_plot(mr_res, sub_hm_dat)
  ggsave("mr_scatter_plot-ADCY3-BMI.pdf", plot = p[[1]], width = 6, height = 6)
}
