#!/usr/bin/env Rscript
library(magrittr)
library(tidyverse)


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


#
## Main steps
#
proj_dir <- "/home/zzhang/Documents/projects/wp_bcg_eqtl"
in_dir <- file.path(proj_dir, "outputs/pseudo_bulk")

# Two model used in the analysis.
mode_vec <- c("normal", "interaction")

# Main cell types.
cell_type_vec <- c("Monocytes", "CD4T", "CD8T", "NK", "B")# , "pDC", "mDC")

# All comparisons used in the analysis.
condition_vec <- c("T0_LPS.vs.T0_RPMI", "T3m_LPS.vs.T3m_RPMI", "T3m_RPMI.vs.T0_RPMI")# , "T3m_LPS.vs.T0_RPMI")

# Reference genotypes panel for local clumping
eur_1kg_geno <- file.path(proj_dir, "inputs/reference/genotypes/GRCh38/EUR")

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
  "ukb-d-J10_ASTHMA" = c("2018_AsthmaDT1_ukb-d-J10_ASTHMA.vcf.gz", NA),
  "ukb-a-255" = c("2017_AsthmaDT2_ukb-a-255.vcf.gz", NA),
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


plink_bin <- "plink-quiet"
gwasvcf::set_plink(plink_bin)
gwasvcf::set_bcftools("bcftools")

mode <- "normal"
cell_type <- "Monocytes"

clump_kb <- 50
clump_r2 <- 0.2
clump <- TRUE
do_proxy <- "no"
eqtl_pval <- 1e-1
eqtl_pval_col <- "p_value"
clump_bfile <- file.path(proj_dir, "inputs/reference/genotypes/GRCh38/EUR")

# GWAS summary
gwas_path <- lapply(pub_gwas_list, function(e) e[[1]]) %>% unlist()
sample_size_dict <- lapply(pub_gwas_list, function(e) e[[2]]) %>% unlist()

# eQTL summary
eqtl <- load_eqtl_tab(file.path(in_dir, mode, cell_type), 0.2, other_info = c("cell_type" = cell_type))

# Select and rename columns in the eQTL summary statistics.
tar_col <- c(
  "SNP" = "snp_id", "beta" = "beta", "se" = "beta_se", "effect_allele" = "assessed_allele",
  "eaf" = "maf", "Phenotype" = "gene_name", "chr" = "snp_chromosome", "pos" = "snp_position",
  "samplesize" = "n_samples", "pval" = eqtl_pval_col
)

tar_egenes <- dplyr::filter(eqtl, global_corrected_pValue < 0.2) %$% feature_id
eqtl_tab <- dplyr::filter(eqtl, feature_id %in% tar_egenes) %>%
  dplyr::select(dplyr::all_of(tar_col)) %>%
  dplyr::filter(pval < eqtl_pval, !is.na(pval)) %>%
  dplyr::mutate(eaf = 1 - eaf) %>%
  TwoSampleMR::format_data(type = "exposure")

tar_snps <- eqtl_tab %>% dplyr::filter(pval.exposure < 0.1) %$% SNP %>% unique()
hm_tab <- lapply(gwas_path, function(path, .rsid, .pval) {
  outcome_info <- stringr::str_remove_all(basename(path), ".vcf.gz") %>% stringr::str_split("_", n = 3, simplify = TRUE)
  outcome <- outcome_info[3]
  outcome_id <- outcome_info[2]

  # Fetch SNP information from GWAS summary.
  cat("[I]: Fetching SNPs from GWAS summary", outcome, outcome_id, "...\n")
  gwas_tab <- gwasvcf::query_gwas(path, rsid = .rsid, proxies = do_proxy, bfile = clump_bfile, threads = 4) %>%
    gwasvcf::query_gwas(pval = .pval) %>%
    gwasglue::gwasvcf_to_TwoSampleMR("outcome") %>%
    dplyr::filter(!is.na(effect_allele.outcome), !is.na(other_allele.outcome), nchar(effect_allele.outcome) == 1, nchar(other_allele.outcome) == 1) %>%
    dplyr::mutate(outcome = outcome, id.outcome = outcome_id)

  # Harmonized exposure and outcome dataset.
  cat("[I]: Harmonizing ...\n")
  suppressMessages({ # Here we suppress messages from the package to make the stdout less mess.
    TwoSampleMR::harmonise_data(eqtl_tab, gwas_tab) %>%
      dplyr::mutate(samplesize.outcome = dplyr::if_else(is.na(samplesize.outcome), as.integer(sample_size_dict[id.outcome]), as.integer(samplesize.outcome)))
  })
}, .rsid = tar_snps, .pval = 0.1) %>%
  Reduce(rbind, .)


# Do clumping locally.
clump_snps <- dplyr::select(hm_tab, SNP, pval.exposure, exposure, id.exposure, outcome, id.outcome, mr_keep) %>%
  dplyr::filter(mr_keep, pval.exposure < 5e-6) %>%
  dplyr::group_by(exposure, id.exposure, outcome, id.outcome)
clump_cols <- c("rsid" = "SNP", "pval" = "pval.exposure", "id" = "exposure")

dplyr::summarise(clump_snps, clumped = {
  if (dplyr::n() >= 3) {
    tryCatch(
      {
        dplyr::cur_data_all() %>%
          dplyr::select(dplyr::all_of(clump_cols)) %>%
          ieugwasr::ld_clump(clump_kb, clump_r2, bfile = clump_bfile, plink_bin = plink_bin) %>%
          dplyr::select(dplyr::all_of(c("SNP" = "rsid"))) %>%
          dplyr::mutate(clump_keep = TRUE)
      },
      error = function(e) {
        print(e)
        data.frame(SNP = NULL, clump_keep = NULL)
      }
    )
  } else {
    dplyr::select(dplyr::cur_data(), SNP) %>% dplyr::mutate(clump_keep = TRUE)
  }
}) %>%
  data.table::as.data.table() %>%
  dplyr::rename_with(.cols = dplyr::starts_with("clumped."), .fn = ~ stringr::str_remove_all(.x, "clumped.")) %>%
  dplyr::right_join(hm_tab, by = c("id.exposure", "exposure", "outcome", "id.outcome", "SNP")) %>%
  dplyr::mutate(clump_keep = dplyr::if_else(is.na(clump_keep), FALSE, clump_keep)) %>%
  data.table::fwrite(paste0(mode, ".", cell_type, ".harmonized.csv"))
