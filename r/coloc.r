#!/usr/bin/env Rscript

# Colocalization analysis by R/coloc
options(stringsAsFactors = FALSE, data.table.verbose = FALSE, data.table.fread.progress = FALSE, future.globals.maxSize = 10000 * 1024^2)
suppressPackageStartupMessages({
  library(coloc)
  library(magrittr)
  library(tidyverse)
  library(org.Hs.eg.db)
})


#' Fetch data for coloc analysis from the harmonized dataset
fetch_coloc_data <- function(hm_dat, what = "exposure", warn_minp = 5e-5) {
  pos_cols <- c("SNP", "pos.exposure")
  src_cols <- c("beta", "se", "eaf", "samplesize")

  cols <- c(paste(src_cols, what, sep = "."), pos_cols)
  names(cols) <- c("beta", "varbeta", "MAF", "N", "snp", "position")

  tab <- hm_dat %>%
    dplyr::select(dplyr::all_of(cols)) %>%
    dplyr::filter(dplyr::if_all(dplyr::everything(), ~ !is.na(.x))) %>%
    dplyr::mutate(type = "quant", varbeta = varbeta ^ 2) %>%
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
colocalization_test <- function(hm_dat, min_snps = 10, ...) {
  ccsum <- c(nsnps = NA, H0 = NA, H1 = NA, H2 = NA, H3 = NA, H4 = NA)

  if (nrow(hm_dat) < min_snps) return(NA)

  qry <- fetch_coloc_data(hm_dat, "exposure")
  ref <- fetch_coloc_data(hm_dat, "outcome")

  if (!is.null(qry) && !is.null(ref)) {
    sink("/dev/null")
    ccsum <- coloc::coloc.abf(qry, ref, ...)$summary
    sink()
  }

  data.frame(
    nsnps = ccsum["nsnps"], H0 = ccsum["PP.H0.abf"], H1 = ccsum["PP.H1.abf"],
    H2 = ccsum["PP.H2.abf"], H3 = ccsum["PP.H3.abf"], H4 = ccsum["PP.H4.abf"]
  )
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
gwas_vec <- c(
  "AlzheimerDiseases", "AsthmaDT1", "AsthmaDT2", "AtopicDermatitis", "BladderCancer",
  "BodyMassIndex", "COVID19Release7", "ColorectalCancer", "CoronaryHeartDisease",
  "CrohnsDisease", "GoutDisease", "HDLCholesterol", "Height", "InflammatoryBowelDisease",
  "LDLCholesterol", "LungCancer", "MultipleSclerosis", "ObesityClass1", "ObesityClass2",
  "ObesityClass3", "OvarianCancer", "ProstateCancer", "Psoriasis", "RheumatoidArthritis",
  "Schizophrenia", "ThyroidCancer", "Triglycerides", "Type2Diabetes", "UlcerativeColitis"
)
gwas_vec <- c("LDLCholesterol", "Triglycerides")

harm_dir <- file.path(proj_dir, "outputs/pseudo_bulk/harmonization")
coloc_dir <- file.path(proj_dir, "outputs/pseudo_bulk/colocalization")
override <- TRUE
for (per_mode in mode_vec) {
  in_base_dir <- file.path(harm_dir, per_mode)
  out_base_dir <- file.path(coloc_dir, per_mode)

  for (per_gwas in gwas_vec) {
    for (per_celltype in celltype_vec) {
      if (per_mode == "normal") {
        cc_save_to <- file.path(out_base_dir, per_celltype)
        if (!dir.exists(cc_save_to)) dir.create(cc_save_to, recursive = TRUE)

        out_file <- file.path(cc_save_to, paste0(per_gwas, ".csv"))
        if (!file.exists(out_file) || override) {
          cat("[I]: Estimating co-localization by coloc for", per_mode, per_gwas, per_celltype, "...\n")

          infile <- file.path(harm_dir, per_mode, per_celltype, paste0(per_gwas, "_harmonized_data.csv.gz"))
          hm_dat <- data.table::fread(infile, showProgress = FALSE)
          ccres <- dplyr::group_by(hm_dat, exposure, outcome, id.outcome) %>%
            dplyr::summarise(cc_res = colocalization_test(dplyr::cur_data_all())) %>%
            data.table::as.data.table() %>%
            dplyr::rename_with(dplyr::starts_with("cc_res."), .fn = ~ stringr::str_remove(.x, "cc_res.")) %>%
            dplyr::filter(dplyr::if_all(dplyr::starts_with("H"), .fns = ~ !is.na(.x)))

          data.table::fwrite(ccres, out_file)
        }
      } else {
        for (per_condition in condition_vec) {
          cc_save_to <- file.path(out_base_dir, per_celltype, per_condition)
          if (!dir.exists(cc_save_to)) dir.create(cc_save_to, recursive = TRUE)

          out_file <- file.path(cc_save_to, paste0(per_gwas, ".csv"))
          if (!file.exists(out_file) || override) {
            cat("[I]: Estimating co-localization by coloc for", per_mode, per_gwas, per_celltype, per_condition, "...\n")

            infile <- file.path(harm_dir, per_mode, per_celltype, per_condition, paste0(per_gwas, "_harmonized_data.csv.gz"))
            hm_dat <- data.table::fread(infile, showProgress = FALSE)
            ccres <- dplyr::group_by(hm_dat, exposure, outcome, id.outcome) %>%
              dplyr::summarise(cc_res = colocalization_test(dplyr::cur_data_all())) %>%
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
