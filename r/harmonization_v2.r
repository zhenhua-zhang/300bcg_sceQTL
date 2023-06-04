#!/usr/bin/env Rscript

# NOTE:
#   1. limix_qtl estimate the additive effects of major alleles and report MAF in the summary statistics.

#' Harmonizing variants of outcome and exposure for MR and COLOC analysis
# quit(save = "no")

# Options, packages used in the analysis
options(stringsAsFactors = FALSE, future.globals.maxSize = 10000 * 1024^2)
suppressPackageStartupMessages({
  library(DBI)
  library(magrittr)
  library(tidyverse)
  library(data.table)
})


#' Load eQTL summary statistics
#'
#' @description A utility function to load summary statistics from the disk
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
  top_qtltab <- data.table::fread(top_qtl_fpath, data.table = FALSE, verbose = FALSE) %>%
    dplyr::select(-dplyr::one_of(discard)) %>%
    (function(dat) dat[dat[, fdr_col] < fdr, ])

  all_qtl_fpath <- file.path(dir, "qtl_results_all.txt")
  all_qtltab <- data.table::fread(all_qtl_fpath, data.table = FALSE, verbose = FALSE)

  cmb_qtltab <- dplyr::right_join(top_qtltab, all_qtltab, by = joinby) %>%
    dplyr::mutate(QTL = stringr::str_c(gene_name, snp_id, sep = "-"))

  if (!is.null(other_info) && !is.null(names(other_info))) {
    for (nm in names(other_info)) {
      cmb_qtltab[nm] <- other_info[nm]
    }
  }

  cmb_qtltab %>% dplyr::filter(dplyr::if_any(dplyr::everything(), ~ !is.na(.x)))
}


#' Create database to store the association summary statistics
#'
#' @description This function create database to store all data.
mkdb_gwas <- function(gwas) {
}

mkdb_eqtl <- function(path) {
}


#' Harmonize variants
#'
#' @description This function try to harmonize variants between two summary statistic and store the
#'              harmonized variants into a sqlite database. If the database has a record, then it
#'              ignore it.
harmonize <- function(eqtl, gwas, dbpath, hmtbl_name = "hmtable", eqtl_p_col = "p_value") {
  con <- DBI::dbConnect(drv = RSQLite::SQLite(), dbname = dbpath)
  DBI::dbListTables(con)
}


proj_dir <- "/home/zzhang/Documents/projects/wp_bcg_eqtl"
in_dir <- file.path(proj_dir, "outputs/pseudo_bulk/summary_statistic")

# Two model used in the analysis.
mode_vec <- c("normal", "interaction")

# Main cell types.
cell_type_vec <- c("Monocytes", "CD4T", "CD8T", "NK", "B")# , "pDC", "mDC")

# All comparisons used in the analysis.
condition_vec <- c("T0_LPS.vs.T0_RPMI", "T3m_LPS.vs.T3m_RPMI", "T3m_RPMI.vs.T0_RPMI")# , "T3m_LPS.vs.T0_RPMI")

#
snp_info_cols <- c("snp_id", "snp_chromosome", "snp_position", "assessed_allele") # , "other_allele")

mode <- "normal"
cell_type <- "B"
condition <- ""
root_dir <- file.path(in_dir, mode, cell_type, condition)
eqtltab <- load_eqtl_tab(root_dir, 0.2, other_info = c("cell_type" = cell_type))

hm_tab <- data.table::fread(file.path(proj_dir, "outputs/pseudo_bulk/harmonization/normal/B.harmonized_data.csv"))

dbcon <- dbConnect(drv = RSQLite::SQLite(), dbname = "snp.sqlite")
snp_info_tab <- eqtltab %>% dplyr::select(dplyr::all_of(snp_info_cols))

dbWriteTable(dbcon, "SNP", snp_info_tab)
dbListTables(dbcon)
dbGetQuery(dbcon, "SELECT * FROM SNP LIMIT 5")
dbDisconnect(dbcon)

dbcon <- dbConnect(drv = RSQLite::SQLite(), dbname = "eqtl.sqlite")
dbWriteTable(dbcon, "normal_B", eqtltab)
dbDisconnect(dbcon)

# 把所有的eQTL作为不同的table写进database
