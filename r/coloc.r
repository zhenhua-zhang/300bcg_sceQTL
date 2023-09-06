#!/usr/bin/env Rscript
# File: coloc.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Jul 16, 2022
# Updated: Jul 16, 2023

# Colocalization analysis by R/coloc
options(stringsAsFactors = FALSE, data.table.verbose = FALSE, future.globals.maxSize = 10000 * 1024^2)
suppressPackageStartupMessages({
  library(coloc)
  library(magrittr)
  library(tidyverse)
  library(org.Hs.eg.db)
})


#' Fetch data for coloc analysis from the harmonized dataset
fetch_coloc_data <- function(hm_dat, what = "exposure", warn_minp = 5e-5, dtype = "quant") {
  pos_cols <- c("SNP", "pos.exposure")
  src_cols <- c("beta", "se", "eaf", "samplesize")

  cols <- c(paste(src_cols, what, sep = "."), pos_cols)
  names(cols) <- c("beta", "varbeta", "MAF", "N", "snp", "position")

  tab <- hm_dat %>%
    dplyr::select(dplyr::all_of(cols)) %>%
    dplyr::filter(dplyr::if_all(dplyr::everything(), ~ !is.na(.x))) %>%
    dplyr::mutate(type = dtype, varbeta = varbeta ^ 2) %>%
    (function(dat) {
      res <- NULL
      if (nrow(dat)) {
        beta <- dplyr::select(dat, "snp", "beta") %>% tibble::deframe()
        varbeta <- dplyr::select(dat, "snp", "varbeta") %>% tibble::deframe()
        maf <- dplyr::select(dat, "snp", "MAF") %>% tibble::deframe()

        res <- list(beta = beta, varbeta = varbeta, snp = dat$snp, position = dat$position,
                    type = dat$type[1], N = max(dat$N), MAF = maf)
      }
      res
    })

  if (!is.null(tab)) coloc::check_dataset(tab, warn.minp = warn_minp)
  return(tab)
}


#' Colocalization anlaysis
colocalization_test <- function(hm_dat, min_snps = 10, gwas_type = "quant", gwas_pval = 5e-5, ...) {
  if (nrow(hm_dat) < min_snps) return(NA)
  if (all(hm_dat$pval.outcome > gwas_pval)) return(NA)

  qry <- fetch_coloc_data(hm_dat, "exposure", dtype = "quant")
  ref <- fetch_coloc_data(hm_dat, "outcome", dtype = gwas_type)

  cc <- c(nsnps = NA, H0 = NA, H1 = NA, H2 = NA, H3 = NA, H4 = NA)
  if (!is.null(qry) && !is.null(ref)) { sink("/dev/null"); cc <- coloc::coloc.abf(qry, ref, ...)$summary; sink() }

  data.frame(nsnps = cc["nsnps"], H0 = cc["PP.H0.abf"], H1 = cc["PP.H1.abf"], H2 = cc["PP.H2.abf"], H3 = cc["PP.H3.abf"], H4 = cc["PP.H4.abf"])
}


#
## Main steps
#
proj_dir <- "/home/zzhang/Documents/projects/wp_bcg_eqtl"
in_dir <- file.path(proj_dir, "outputs/pseudo_bulk")

# Two model used in the analysis.
mode_vec <- c("normal", "interaction")

# Main cell types.
celltype_vec <- c("Monocytes", "CD4T", "CD8T", "NK", "B") # , "pDC", "mDC")

# All comparisons used in the analysis.
condition_vec <- c("T0_LPS.vs.T0_RPMI", "T3m_LPS.vs.T3m_RPMI", "T3m_RPMI.vs.T0_RPMI")#, "T3m_LPS.vs.T0_RPMI")

# All GWAS summary to be used in the colocalization analysis
gwas_meta <- tibble::tribble(
  ~gwas, ~type,
  "AlzheimerDiseases", "cc",
  "AsthmaDT1", "cc",
  "AsthmaDT2", "cc",
  "AtopicDermatitis", "cc",
  "BladderCancer", "cc",
  "BodyMassIndex", "quant",
  "COVID19Release7", "cc",
  "ColorectalCancer", "cc",
  "CoronaryHeartDisease", "cc",
  "CrohnsDisease", "cc",
  "GoutDisease", "cc",
  "HDLCholesterol", "quant",
  "Height", "quant",
  "InflammatoryBowelDisease", "cc",
  "LDLCholesterol", "quant",
  "LungCancer", "cc",
  "MultipleSclerosis", "cc",
  "ObesityClass1", "cc",
  "ObesityClass2", "cc",
  "ObesityClass3", "cc",
  "OvarianCancer", "cc",
  "ProstateCancer", "cc",
  "Psoriasis", "cc",
  "RheumatoidArthritis", "cc",
  "Schizophrenia", "cc",
  "ThyroidCancer", "cc",
  "Triglycerides", "cc",
  "Type2Diabetes", "cc",
  "UlcerativeColitis", "cc",
)

harm_dir <- file.path(proj_dir, "outputs/pseudo_bulk/harmonization")
coloc_dir <- file.path(proj_dir, "outputs/pseudo_bulk/colocalization")
override <- TRUE
for (per_mode in mode_vec) {
  # per_mode <- "normal"
  in_base_dir <- file.path(harm_dir, per_mode)
  out_base_dir <- file.path(coloc_dir, per_mode)

  for (per_gwas in gwas_meta$gwas) {
    # per_gwas <- "COVID19Release7"
    gwas_type <- gwas_meta$type[gwas_meta$gwas == per_gwas]
    for (per_celltype in celltype_vec) {
      # per_celltype <- "CD4T"
      if (per_mode == "normal") {
        cc_save_to <- file.path(out_base_dir, per_celltype)
        if (!dir.exists(cc_save_to)) dir.create(cc_save_to, recursive = TRUE)

        out_file <- file.path(cc_save_to, paste0(per_gwas, ".csv"))
        if (!file.exists(out_file) || override) {
          cat("[I]: Estimating co-localization by coloc for", per_mode, per_gwas, per_celltype, "...\n")

          infile <- file.path(harm_dir, per_mode, per_celltype, paste0(per_gwas, "_harmonized_data.csv.gz"))
          hm_dat <- data.table::fread(infile, showProgress = FALSE)
          ccres <- dplyr::group_by(hm_dat, exposure, outcome, id.outcome) %>%
            dplyr::summarise(cc_res = colocalization_test(dplyr::cur_data_all(), gwas_type = gwas_type), .groups = "keep") %>%
            data.table::as.data.table() %>%
            dplyr::rename_with(dplyr::starts_with("cc_res."), .fn = ~ stringr::str_remove(.x, "cc_res.")) %>%
            dplyr::filter(dplyr::if_all(dplyr::starts_with("H"), .fns = ~ !is.na(.x)))

          data.table::fwrite(ccres, out_file)
        }
      } else {
        for (per_condition in condition_vec) {
          # per_condition <- "T0_LPS.vs.T0_RPMI"
          cc_save_to <- file.path(out_base_dir, per_celltype, per_condition)
          if (!dir.exists(cc_save_to)) dir.create(cc_save_to, recursive = TRUE)

          out_file <- file.path(cc_save_to, paste0(per_gwas, ".csv"))
          if (!file.exists(out_file) || override) {
            cat("[I]: Estimating co-localization by coloc for", per_mode, per_gwas, per_celltype, per_condition, "...\n")

            infile <- file.path(harm_dir, per_mode, per_celltype, per_condition, paste0(per_gwas, "_harmonized_data.csv.gz"))
            hm_dat <- data.table::fread(infile, showProgress = FALSE)
            ccres <- dplyr::group_by(hm_dat, exposure, outcome, id.outcome) %>%
              dplyr::summarise(cc_res = colocalization_test(dplyr::cur_data_all(), gwas_type = gwas_type), .groups = "keep") %>%
              data.table::as.data.table() %>%
              dplyr::rename_with(dplyr::starts_with("cc_res."), .fn = ~ stringr::str_remove(.x, "cc_res.")) %>%
              dplyr::filter(dplyr::if_all(dplyr::starts_with("H"), .fns = ~ !is.na(.x)))

            data.table::fwrite(ccres,  out_file)
          }
        }
      }
    }
  }
}
