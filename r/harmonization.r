#!/usr/bin/env Rscript

#' Harmonizing variants of outcome and exposure for MR and COLOC analysis
# quit(save = "no")

# Options, packages used in the analysis
options(stringsAsFactors = FALSE, future.globals.maxSize = 10000 * 1024 ^ 2)
library(magrittr)
library(tidyverse)
library(data.table)


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
    for (nm in names(other_info)) { cmb_qtltab[nm] <- other_info[nm] }
  }

  cmb_qtltab %>% dplyr::filter(dplyr::if_any(dplyr::everything(), ~ !is.na(.x)))
}


#' Harmonize variants
harmonize_vars <- function(eqtl, gwas, alt_af_db, eqtl_pval_col = "p_value", eqtl_pval = 0.1,
                           clump = TRUE, clump_kb = 20, clump_r2 = 0.8, clump_bfile = NULL, do_proxy = TRUE,
                           plink_bin = "plink", bcftools_bin = "bcftools") {
  # Set up bcftools and plink path
  gwasvcf::set_bcftools(bcftools_bin)
  gwasvcf::set_plink(plink_bin)

  # Select and rename columns in the eQTL summary statistics.
  tar_col <- c(
    "SNP" = "snp_id", "beta" = "beta", "se" = "beta_se", "effect_allele" = "assessed_allele", "eaf" = "maf", "Phenotype" = "gene_name",
    "chr" = "snp_chromosome", "pos" = "snp_position", "samplesize" = "n_samples", "pval" = eqtl_pval_col
  )

  # Load the eQTL summary statistics as the exposure dataset.
  eqtl_tab <- dplyr::select(eqtl, dplyr::all_of(tar_col)) %>%
    dplyr::filter(pval < eqtl_pval, !is.na(pval)) %>%
    dplyr::mutate(eaf = alt_af_db[SNP]) %>%
    TwoSampleMR::format_data(type = "exposure")

  # Sample size
  gwas_path <- lapply(gwas, function(e) e[[1]]) %>% unlist()
  sample_size_dict <- lapply(gwas, function(e) e[[2]]) %>% unlist()

  # Fetch public available GWAS.
  tar_snps <- unique(eqtl_tab$SNP)
  do_proxy <- ifelse(do_proxy && !is.null(clump_bfile), "no", "yes")
  gwas_tab <- lapply(gwas_path, function(path, .rsid, .pval) {
    outcome_info <- basename(path) %>% stringr::str_remove_all(".vcf.gz") %>% stringr::str_split("_", n = 3, simplify = TRUE)
    gwasvcf::query_gwas(path, rsid = .rsid, proxies = do_proxy, bfile = proxy_ld_file, threads = 4) %>% 
      gwasvcf::query_gwas(pval = .pval) %>%
      gwasglue::gwasvcf_to_TwoSampleMR("outcome") %>%
      dplyr::filter(!is.na(effect_allele.outcome), !is.na(other_allele.outcome), nchar(effect_allele.outcome) == 1, nchar(other_allele.outcome) == 1) %>%
      dplyr::mutate(outcome = outcome_info[2], id.outcome = outcome_info[3])
  }, .rsid = tar_snps, .pval = 0.1) %>%
  Reduce(rbind, .)

  # Harmonized exposure and outcome dataset.
  suppressMessages({ # Here we suppress messages from the package to make the stdout less mess.
    hm_tab <- TwoSampleMR::harmonise_data(eqtl_tab, gwas_tab) %>%
      dplyr::mutate(samplesize.outcome = dplyr::if_else(is.na(samplesize.outcome), as.integer(sample_size_dict[id.outcome]), as.integer(samplesize.outcome)))
  })

  # Do clumping locally.
  if (clump && !is.null(clump_bfile)) {
    clump_snps <- dplyr::select(hm_tab, SNP, pval.exposure, exposure, outcome, mr_keep) %>%
      dplyr::filter(mr_keep, pval.exposure < 5e-6, pval.exposure < 0.05) %>%
      dplyr::group_by(exposure, outcome)
    clump_cols <- c("rsid" = "SNP", "pval" = "pval.exposure", "id" = "exposure")
    hm_tab <- dplyr::summarise(clump_snps, clumped = {
      if (n() >= 3)
        tryCatch({
          dplyr::cur_data_all() %>%
            dplyr::select(dplyr::all_of(clump_cols)) %>%
            ieugwasr::ld_clump(plink_bin = plink_bin, bfile = clump_bfile, clump_kb = clump_kb, clump_r2 = clump_r2) %>%
            dplyr::select(dplyr::all_of(c("SNP" = "rsid"))) %>%
            dplyr::mutate(clump_keep = TRUE)
        }, error = function(e) data.frame(SNP = NULL, clump_keep = NULL))
      else
        dplyr::select(dplyr::cur_data(), SNP) %>% dplyr::mutate(clump_keep = TRUE)
    }) %>%
      data.table::as.data.table() %>%
      dplyr::rename_with(.cols = dplyr::starts_with("clumped."), .fn = ~ stringr::str_remove_all(.x, "clumped.")) %>%
      dplyr::right_join(hm_tab, by = c("outcome", "exposure", "SNP")) %>%
      dplyr::mutate(clump_keep = dplyr::if_else(is.na(clump_keep), FALSE, clump_keep))
  }

  hm_tab
}


#
## Main steps
#
proj_dir <- "/home/zzhang/Documents/projects/wp_bcg_eqtl"
in_dir <- file.path(proj_dir, "outputs/pseudo_bulk")

# Two model used in the analysis.
mode_vec <- c("normal", "interaction")

# Main cell types.
cell_type_vec <- c("Monocytes", "CD4T", "CD8T", "NK", "B") # , "pDC", "mDC")

# All comparisons used in the analysis.
condition_vec <- c("T0_LPS.vs.T0_RPMI", "T3m_LPS.vs.T0_RPMI", "T3m_LPS.vs.T3m_RPMI", "T3m_RPMI.vs.T0_RPMI")

# Reference genotypes panel for local clumping
eur_1kg_geno <- file.path(proj_dir, "inputs/reference/genotypes/EUR")

# Allele frequecies of effect allele used in sceQTL mapping.
afdb_path <- file.path(proj_dir, "inputs/genotypes/300BCG_sub40_imp_hg38_ref_allele_frequency.csv")
afdb <- data.table::fread(afdb_path) %>% tibble::deframe()

# Public GWAS list
pub_gwas_list <- list(
  # Cancers
  "ieu-b-4965" = c("2021_ColorectalCancer_ieu-b-4965.vcf.gz", NA),
  "ieu-b-4809" = c("2021_ProstateCancer_ieu-b-4809.vcf.gz", NA),
  "ieu-b-4874" = c("2021_BladderCancer_ieu-b-4874.vcf.gz", NA),
  "ieu-a-1082" = c("2013_ThyroidCancer_ieu-a-1082.vcf.gz", NA),
  "ieu-b-4963" = c("2017_OvarianCancer_ieu-b-4963.vcf.gz", NA),
  "ieu-a-966" = c("2014_LungCancer_ieu-a-966.vcf.gz", NA),

  # Autoimmune disease
  # "ebi-a-GCST010681" = c("2022_Type1Diabetes_ebi-a-GCST010681.vcf.gz", NA), No GWAS VCF available, TODO: format available SS into GWAS VCF format.
  "ukb-d-M13_RHEUMA" = c("2018_RheumatoidArthritis_ukb-d-M13_RHEUMA.vcf.gz", NA),
  "ieu-a-31" = c("2015_InflammatoryBowelDisease_ieu-a-31.vcf.gz", NA),
  "ukb-b-17670" = c("2019_MultipleSclerosis_ukb-b-17670.vcf.gz", NA),
  "ieu-a-996" = c("2014_AtopicDermatitis_ieu-a-996.vcf.gz", NA),
  "ieu-a-32" = c("2015_UlcerativeColitis_ieu-a-32.vcf.gz", NA),
  "ieu-a-30" = c("2015_CrohnsDisease_ieu-a-30.vcf.gz", NA),
  "ukb-a-100" = c("2017_Psoriasis_ukb-a-100.vcf.gz", NA),

  # Infectious disease
  "ebi-a-GCST010780" = c("2020_COVID19Release4_ebi-a-GCST010780.vcf.gz", NA),

  # Brain disorders
  "ieu-b-5067" = c("2022_AlzheimerDiseases_ieu-b-5067.vcf.gz", NA),
  "ieu-b-42" = c("2014_Schizophrenia_ieu-b-42.vcf.gz", NA),

  # Other genetic-related disease
  "ebi-a-GCST006867" = c("2018_Type2Diabetes_ebi-a-GCST006867.vcf.gz", NA),
  "ukb-d-J10_ASTHMA" = c("2018_Asthma_ukb-d-J10_ASTHMA.vcf.gz", NA),
  "ieu-a-7" = c("2013_CoronaryHeartDisease_ieu-a-7.vcf.gz", NA),
  "ukb-a-107" = c("2017_GoutDisease_ukb-a-107.vcf.gz", NA),

  # Others traits
  "ieu-a-89" = c("2014_Height_ieu-a-89.vcf.gz", NA),
  "ieu-b-40" = c("2018_BodyMassIndex_ieu-b-40.vcf.gz", NA),
  "ieu-b-111" = c("2020_Triglycerides_ieu-b-111.vcf.gz", 441016),
  "ieu-b-109" = c("2020_HDLCholesterol_ieu-b-109.vcf.gz", 403943),
  "ieu-b-110" = c("2020_LDLCholesterol_ieu-b-110.vcf.gz", 440546)
) %>%
lapply(function(e) c(file.path(proj_dir, "inputs/public_gwas", e[1]), e[2]))


# Per loop per eQTL-mapping-run
override <- FALSE
for (mode in mode_vec) {
  qtltab_list <- list()
  for (cell_type in cell_type_vec) {
    if (mode == "normal") {
      root_dir <- file.path(in_dir, mode, cell_type)
      run_id <- paste(mode, cell_type, sep = "_")

      if (file.exists(file.path(root_dir, "top_qtl_results_all_FDR.txt"))) {
        cat("[I]: Loading results from", root_dir, "...\n")
        qtltab_list[[run_id]] <- load_eqtl_tab(root_dir, 0.1, other_info = c("cell_type" = cell_type))
      } else {
        cat("[W]: Not top QTL results available for", run_id, "Skipping ...\n")
      }
    } else {
      for (condition in condition_vec) {
        root_dir <- file.path(in_dir, mode, cell_type, condition)
        run_id <- paste(mode, cell_type, condition, sep = "_")

        if (file.exists(file.path(root_dir, "top_qtl_results_all_FDR.txt"))) {
          cat("[I]: Loading results from", root_dir, "...\n")
          qtltab_list[[run_id]] <- load_eqtl_tab(root_dir, 0.1, other_info = c("cell_type" = cell_type, "condition" = condition))
        } else {
          cat("[W]: Not top QTL results available for", run_id, "Skipping ...\n")
        }
      }
    }
  }

  save_to <- file.path(proj_dir, "outputs/pseudo_bulk/outcomes", mode)

  # Harmonization
  .tmp <- lapply(names(qtltab_list), function(pr_run_id) {
    tmp_tab <- qtltab_list[[pr_run_id]] %>%
      (function(dat) {
        top_egenes <- dplyr::filter(dat, !is.na(global_corrected_pValue)) %>% dplyr::select(gene_name) %>% unlist() %>% as.vector()
        dplyr::filter(dat, gene_name %in% top_egenes, p_value < 0.1, !stringr::str_starts(gene_name, "HLA"))
      })

    pr_save_to <- file.path(save_to, pr_run_id)
    if (!dir.exists(pr_save_to)) dir.create(pr_save_to, recursive = TRUE)

    hm_save_to <- file.path(pr_save_to, "harmonized_data.csv")
    if (!file.exists(hm_save_to) || override) {
      cat("[I]: Harmonizing variants ...\n")
      hm_dat <- harmonize_vars(tmp_tab, pub_gwas_list, afdb, clump_bfile = eur_1kg_geno, plink_bin = "plink-quiet")
      data.table::fwrite(hm_dat, hm_save_to, verbose = FALSE, showProgress = FALSE)
    } else {
      cat("[W]: Found", hm_save_to, "Skipping ...\n")
    }
    NULL
  })
}
