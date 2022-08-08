#!/usr/bin/env Rscript
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Mar 14, 2022
# Updated: Mar 14, 2022

library(coloc)
library(magrittr)
library(tidyverse)
library(data.table)

wk_dir <- "~/Documents/projects/wp_bcg_eqtl"
tmp_dir <- file.path(wk_dir, "temps")

DEBUG <- TRUE
if (DEBUG) {
  n_rows <- 5000
} else {
  n_rows <- NULL
}


# Load genotype of GRCh38.
tar_cols <- c("ID", "REF", "ALT", "AF")
gntp_file <- file.path(wk_dir, "inputs/genotypes/300BCG_sub40_imp_hg38_ids.vcf.gz")
gntp_dtfm <- fread(gntp_file, tmpdir = tmp_dir, skip = "#CHROM", nrow = n_rows) %>%
  dplyr::mutate(
    AF = str_extract_all(INFO, ";AF=.*?;"),
    AF = as.double(str_remove_all(AF, ";AF=|;$"))
  ) %>%
  dplyr::select(one_of(tar_cols))


# Load eQTL summary statistics
eqtl_col_nm <- c("GENE", "CHROM", "ID", "POS", "PVAL", "ES", "SE")
eqtl_file <- file.path(wk_dir, "inputs/bcg_eqtl/coloc/bcg.lps_pdc_nominal_coloc.txt")
eqtl_dtfm <- fread(eqtl_file, tmpdir = tmp_dir, nrows = n_rows) %>%
  dplyr::select(-V8) %>%
  rlang::set_names(eqtl_col_nm) %>%
  dplyr::left_join(gntp_dtfm, by = "ID") %>%
  dplyr::mutate(LP = -log10(PVAL))


# Load pub GWAS summary statistics
tar_ss <- "ebi-a-GCST007799"
fmt_fields <- c("ES", "SE", "LP", "AF", "ID")
gwas_tar_cols <- c("ID", "REF", "ALT")
gwas_file <- file.path(wk_dir, "inputs/public_gwas/ebi-a-GCST007799.vcf.gz")
gwas_dtfm <- fread(gwas_file, tmpdir = tmp_dir, nrows = n_rows, skip = "#CHROM") %>%
  tidyr::separate(tar_ss, into = fmt_fields, sep = ":") %>%
  dplyr::select(one_of(c(gwas_tar_cols, fmt_fields)))

inner_join(eqtl_dtfm, gwas_dtfm, by = "ID", suffix = c(".eqtl", ".gwas")) %>%
  dplyr::group_by(GENE, ID, CHROM, POS) %>%
  dplyr::summarise(align = {
    cur_data() %>%
      dplyr::mutate(
        pass = case_when(
          (AF.eqtl < 0.5 & AF.gwas < 0.5) | (AF.eqtl > 0.5 & AF.gwas > 0.5) ~ TRUE,
          T ~ FALSE
        )
      )
  })



eqtl_coloc <- list()
gwas_coloc <- list()

eqtl_dtfm %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    summary = {
    }
  )
