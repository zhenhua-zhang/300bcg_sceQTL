#!/usr/bin/env Rscript
# File: harmonization.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Oct 12, 2022
# Updated: Jul 07, 2023


# NOTE:
#   1. limix_qtl estimate the additive effects of major alleles and report MAF in the summary statistics.

# Options beforehand
options(stringsAsFactors = FALSE, data.table.verbose = FALSE, future.globals.maxSize = 10000 * 1024^2)

# Options, packages used in the analysis
suppressPackageStartupMessages({
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
  top_qtltab <- data.table::fread(top_qtl_fpath, data.table = FALSE) %>%
    dplyr::select(-dplyr::one_of(discard)) %>%
    (function(dat) dat[dat[, fdr_col] < fdr, ])

  all_qtl_fpath <- file.path(dir, "qtl_results_all.txt")
  all_qtltab <- data.table::fread(all_qtl_fpath, data.table = FALSE)

  cmb_qtltab <- dplyr::right_join(top_qtltab, all_qtltab, by = joinby) %>%
    dplyr::mutate(QTL = stringr::str_c(gene_name, snp_id, sep = "-"))

  if (!is.null(other_info) && !is.null(names(other_info))) {
    for (nm in names(other_info)) {
      cmb_qtltab[nm] <- other_info[nm]
    }
  }

  cmb_qtltab %>% dplyr::filter(dplyr::if_any(dplyr::everything(), ~ !is.na(.x)))
}


#' Harmonize variants
harmonize_vars <- function(eqtl, gwas, eqtl_pval_col = "p_value", eqtl_pval = 1, exposure_pval = 5e-6,
                           clump = TRUE, clump_kb = 50, clump_r2 = 0.5, clump_bfile = NULL, do_proxy = FALSE,
                           plink_bin = "plink", bcftools_bin = "bcftools", save_to = "./", override = FALSE) {
  # Set up bcftools and plink path
  gwasvcf::set_bcftools(bcftools_bin)
  gwasvcf::set_plink(plink_bin)

  # Sample size
  gwas_path <- lapply(gwas, function(e) e[[1]]) %>% unlist()
  sample_size_dict <- lapply(gwas, function(e) e[[2]]) %>% unlist()

  # Select and rename columns in the eQTL summary statistics.
  tar_col <- c(
    "SNP" = "snp_id", "beta" = "beta", "se" = "beta_se", "effect_allele" = "assessed_allele",
    "eaf" = "maf", "Phenotype" = "gene_name", "chr" = "snp_chromosome", "pos" = "snp_position",
    "samplesize" = "n_samples", "pval" = eqtl_pval_col
  )

  # Load the eQTL summary statistics as the exposure dataset.
  tar_egenes <- dplyr::filter(eqtl, global_corrected_pValue < 0.2) %>% dplyr::pull(feature_id) %>% unique()

  eqtl_tab <- dplyr::filter(eqtl, feature_id %in% tar_egenes) %>%
    dplyr::select(dplyr::all_of(tar_col)) %>%
    dplyr::filter(pval < eqtl_pval, !is.na(pval)) %>%
    dplyr::mutate(eaf = 1 - eaf) %>%
    TwoSampleMR::format_data(type = "exposure")

  # Fetch eQTL SNPs from public available GWAS.
  tar_snps <- eqtl_tab %>% dplyr::filter(pval.exposure < 0.5) %$% SNP %>% unique()
  do_proxy <- ifelse(do_proxy && !is.null(clump_bfile), "yes", "no")
  lapply(gwas_path, function(path, .rsid, .pval) {
    outcome_info <- stringr::str_remove_all(basename(path), ".vcf.gz") %>% stringr::str_split("_", n = 3, simplify = TRUE)
    .outcome <- outcome_info[2]
    .outcome_id <- outcome_info[3]

    hm_tab_save_to <- file.path(save_to, paste0(.outcome, "_harmonized_data.csv"))
    if (file.exists(hm_tab_save_to) || file.exists(paste0(hm_tab_save_to, ".gz")) && (!override)) {
      cat("[W]:", hm_tab_save_to, "exists, skipping it ...\n")

      return(NULL)
    }

    # Fetch SNP information from GWAS summary.
    cat("[I]: Fetching SNPs from GWAS summary", .outcome, .outcome_id, "...\n")
    per_gwas_tab <- gwasvcf::query_gwas(path, rsid = .rsid, proxies = do_proxy, bfile = clump_bfile, threads = 4) %>%
      gwasvcf::query_gwas(pval = .pval) %>%
      gwasglue::gwasvcf_to_TwoSampleMR("outcome") %>%
      dplyr::filter(!is.na(effect_allele.outcome), !is.na(other_allele.outcome), nchar(effect_allele.outcome) == 1, nchar(other_allele.outcome) == 1)

    # Harmonized exposure and outcome dataset.
    cat("[I]: Harmonizing ...\n")
    suppressMessages({ # Here we suppress messages from the package to make the stdout less mess.
      per_hm_tab <- TwoSampleMR::harmonise_data(eqtl_tab, per_gwas_tab) %>%
        dplyr::mutate(samplesize.outcome = dplyr::if_else(is.na(samplesize.outcome), as.integer(sample_size_dict[id.outcome]), as.integer(samplesize.outcome)),
                      outcome = .outcome, id.outcome = .outcome_id) %>%
        dplyr::relocate(exposure, id.exposure, outcome, id.outcome)
    })
    data.table::fwrite(per_hm_tab, hm_tab_save_to)

    NULL
  }, .rsid = tar_snps, .pval = 0.5)
}


#
## Main steps
#
proj_dir <- "/home/zzhang/Documents/projects/wp_bcg_eqtl"
in_dir <- file.path(proj_dir, "outputs/pseudo_bulk/summary_statistic")

# Two model used in the analysis.
mode_vec <- c("normal", "interaction")

# Main cell types.
cell_type_vec <- c("Monocytes", "CD4T", "CD8T", "NK", "B")# , "pDC", "mDC")

# All comparisons used in the analysis.
condition_vec <- c("T0_LPS.vs.T0_RPMI", "T3m_LPS.vs.T3m_RPMI", "T3m_RPMI.vs.T0_RPMI")# , "T3m_LPS.vs.T0_RPMI")

# Reference genotypes panel for local clumping
eur_1kg_geno <- file.path(proj_dir, "inputs/reference/genotypes/GRCh38/EUR")

# Public GWAS list: run 1
pub_gwas_list <- list(
  # # Cancers
  # "ieu-b-4965" = c("2021_ColorectalCancer_ieu-b-4965.vcf.gz", NA),
  # "ieu-b-4809" = c("2021_ProstateCancer_ieu-b-4809.vcf.gz", NA),
  # "ieu-b-4874" = c("2021_BladderCancer_ieu-b-4874.vcf.gz", NA),
  # "ieu-a-1082" = c("2013_ThyroidCancer_ieu-a-1082.vcf.gz", NA),
  # "ieu-b-4963" = c("2017_OvarianCancer_ieu-b-4963.vcf.gz", NA),
  # "ieu-a-966" = c("2014_LungCancer_ieu-a-966.vcf.gz", NA),

  # # Autoimmune disease
  # # "ebi-a-GCST010681" = c("2022_Type1Diabetes_ebi-a-GCST010681.vcf.gz", NA), No GWAS VCF available, TODO: format available SS into GWAS VCF format.
  # "ukb-d-M13_RHEUMA" = c("2018_RheumatoidArthritis_ukb-d-M13_RHEUMA.vcf.gz", NA),
  # "ieu-a-31" = c("2015_InflammatoryBowelDisease_ieu-a-31.vcf.gz", NA),
  # "ukb-b-17670" = c("2019_MultipleSclerosis_ukb-b-17670.vcf.gz", NA),
  # "ieu-a-996" = c("2014_AtopicDermatitis_ieu-a-996.vcf.gz", NA),
  # "ieu-a-32" = c("2015_UlcerativeColitis_ieu-a-32.vcf.gz", NA),
  # "ieu-a-30" = c("2015_CrohnsDisease_ieu-a-30.vcf.gz", NA),
  # "ukb-a-100" = c("2017_Psoriasis_ukb-a-100.vcf.gz", NA),

  # # Infectious disease
  # "ebi-a-GCST010780" = c("2020_COVID19Release4_ebi-a-GCST010780.vcf.gz", NA),
  # "HGI-A2-ALL-eur-leave23andme-20220403" = c("2022_COVID19Release7_HGI-A2-ALL-eur-leave23andme-20220403.vcf.gz", NA),
  "2023_COVID19PlosOne" = c("2023_COVID19PlosOne_Severity.vcf.gz", NA),

  # # Brain disorders
  # "ieu-b-5067" = c("2022_AlzheimerDiseases_ieu-b-5067.vcf.gz", NA),
  # "ieu-b-42" = c("2014_Schizophrenia_ieu-b-42.vcf.gz", NA),

  # # Other genetic-related disease
  # "ebi-a-GCST006867" = c("2018_Type2Diabetes_ebi-a-GCST006867.vcf.gz", NA),
  # "ukb-d-J10_ASTHMA" = c("2018_AsthmaDT1_ukb-d-J10_ASTHMA.vcf.gz", NA),
  # "ieu-a-7" = c("2013_CoronaryHeartDisease_ieu-a-7.vcf.gz", NA),
  # "ukb-a-107" = c("2017_GoutDisease_ukb-a-107.vcf.gz", NA),
  # "ukb-a-255" = c("2017_AsthmaDT2_ukb-a-255.vcf.gz", NA),
  # "ieu-b-90" = c("2013_ObesityClass1_ieu-a-90.vcf.gz", NA),
  # "ieu-b-91" = c("2013_ObesityClass2_ieu-a-91.vcf.gz", NA),
  # "ieu-b-92" = c("2013_ObesityClass3_ieu-a-92.vcf.gz", NA),

  # # Others traits
  # "ieu-b-109" = c("2020_HDLCholesterol_ieu-b-109.vcf.gz", 403943),
  # "ieu-b-110" = c("2020_LDLCholesterol_ieu-b-110.vcf.gz", 440546),
  # "ieu-b-111" = c("2020_Triglycerides_ieu-b-111.vcf.gz", 441016),
  # "ieu-b-40" = c("2018_BodyMassIndex_ieu-b-40.vcf.gz", NA),
  # "ieu-a-89" = c("2014_Height_ieu-a-89.vcf.gz", NA)
) %>%
  lapply(function(e) c(file.path(proj_dir, "inputs/public_gwas", e[1]), e[2]))


# Per loop per eQTL-mapping-run
for (mode in mode_vec) {
  save_to <- file.path(proj_dir, "outputs/pseudo_bulk/harmonization", mode)
  for (cell_type in cell_type_vec) {
    # Load eQTL summary
    qtltab_list <- list()
    if (mode == "normal") {
      root_dir <- file.path(in_dir, mode, cell_type)
      run_id <- cell_type
      qtltab_list[[run_id]] <- load_eqtl_tab(root_dir, 0.2, other_info = c("cell_type" = cell_type))
    } else {
      for (condition in condition_vec) {
        root_dir <- file.path(in_dir, mode, cell_type, condition)
        run_id <- file.path(cell_type, condition)
        qtltab_list[[run_id]] <- load_eqtl_tab(root_dir, 0.2, other_info = c("cell_type" = cell_type, "condition" = condition))
      }
    }

    # Harmonization
    cat("[I]: Start to harmonize summary statistics ...\n")
    .tmp <- lapply(names(qtltab_list), function(pr_run_id) {
      .tmp_tab <- qtltab_list[[pr_run_id]] %>%
        (function(dat) {
          top_egenes <- dplyr::filter(dat, !is.na(global_corrected_pValue)) %>% dplyr::pull(gene_name)
          dplyr::filter(dat, gene_name %in% top_egenes, !stringr::str_starts(gene_name, "HLA"))
        })

      pr_save_to <- file.path(save_to, pr_run_id)
      if (!dir.exists(pr_save_to)) dir.create(pr_save_to, recursive = TRUE)
      harmonize_vars(.tmp_tab, pub_gwas_list, exposure_pval = 5e-5, clump_bfile = eur_1kg_geno, plink_bin = "plink-quiet", save_to = pr_save_to)
      NULL
    })
  }
}
