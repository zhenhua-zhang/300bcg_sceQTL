#!/usr/bin/env Rscript
# Author : Zhenhua Zhang
# E-mail : zhenhua.zhang217@gmail.com
# Created: Feb 15, 2022
# Updated: Apr 18, 2022

# SCeQTL using a zero-inflated negative binomial regression.

options(datatable.verbose = FALSE, stringsAsFactors = FALSE)

library(data.table)
library(tidyverse)
library(magrittr)
library(optparse)
library(glmmTMB)

# Readme
readme <- function() {
  cat(
    "1. Spearman's rank correlation testing.",
    "2. Sex, age, six genotype PCs, and two PEER factors.",
    "3. Correlate expression residuals to genotypes (0/1/2) per individual.",
    "4. Expression residuals are calculate by regression co-factors.",
    sep = "\n"
  )
}

proj_dir <- "~/Documents/projects/wp_bcg_eqtl"

#' Get CLI options
get_opts <- function() {
  opt <- optparse::OptionParser()
  opt <- optparse::add_option(opt, c("-g", "--gntp-dir"),
    help = "A directory including genotype files named by gene symbols."
  )
  opt <- optparse::add_option(opt, c("-p", "--phtp-file"),
    help = "Path to phenotype file. Columns are genes and rows are cells."
  )
  opt <- optparse::add_option(opt, c("-c", "--cvrt-file"),
    help = "Path to covariate file. Columns are covariates and rows are sample."
  )

  opt <- optparse::add_option(opt, c("-r", "--rand-seed"),
    help = "Random seed."
  )
  opt <- optparse::add_option(opt, c("-n", "--n-cpus"),
    help = "Nr. of parallel jobs."
  )
  opt <- optparse::add_option(opt, c("-o", "--out-dir"),
    help = "Directory for outputs."
  )

  optparse::parse_args2(opt)
}


#
## Main steps
#
opts <- get_opts()

# Load genotypes
gntp_dir <- opts$gntp_dir
if (is.null(gntp_dir)) {
  stop("Required genotype directory by -g/--gntp is missing!")
}
gntp_files <- list.files(gntp_dir, pattern = "*.csv", full.names = TRUE)

# Load expression matrix
phtp_file <- opts$phtp_file
if (is.null(phtp_file)) {
  stop("Required phenotype files by -p/--phtp is missing!")
}
phtp_tab <- fread(phtp_file) %>% dplyr::rename("PatientID" = "ids")


# Reproducibility
rand_seed <- ifelse(is.null(opts$rand_seed), 1024, opts$rand_seed)
set.seed(rand_seed)

# IO
out_dir <- ifelse(is.null(opts$out_dir), "./sceqtls", opts$out_dir)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load covariates
cvrt_file <- opts$cvrt
if (!is.null(cvrt_file))
  cvrt_tab <- fread(cvrt_file) %>% dplyr::rename("PatientID" = "ids")


n_cpus <- ifelse(is.null(opts$n_cpus), 1, as.integer(opts$n_cpus))
# Integrate genotypes, phenotypes, and covariates. Perform regression analysis.
lapply(gntp_files, function(per_fp, phtp, cvrt) {
  gene_id <- per_fp %>%
    basename() %>%
    str_remove_all(".csv")

  fread(per_fp) %>%
    dplyr::group_by(GeneID, SNPs) %>%
    dplyr::summarise(est = {
      geneid <- cur_group()$GeneID %>% as.character()
      snpid <- cur_group()$SNPs %>% as.character()

      reg_tab <- NULL
      if (geneid %in% colnames(phtp)) {
        sub_phtp <- phtp %>% select(one_of(geneid, "PatientID"))
        sub_gntp <- cur_data() %>%
          pivot_longer(everything(), names_to = "PatientID", values_to = snpid)

        fml <- formula(
          paste(geneid, "~", snpid, "+ (1 | PatientID) + age + gender")
        )
        zifml <- formula(paste("~", snpid, "+ age + gender"))

        reg_tab <- sub_gntp %>%
          inner_join(sub_phtp, by = "PatientID") %>%
          inner_join(cvrt, by = "PatientID") %>%
          glmmTMB(fml,
            data = ., ziformula = ~1, family = poisson,
            control = glmmTMBControl(parallel = n_cpus)
          ) %>%
          summary() %>%
          coef() %$% cond %>%
          as.data.frame() %>%
          filter(str_detect(rownames(.), "rs")) %>%
          set_names(c("Est", "StdErr", "ZValue", "PValue"))
      }

      reg_tab
    }) %>%
    as.data.table() %>%
    rename_with(.fn = function(e) str_remove_all(e, "est.")) %>%
    fwrite(paste0(out_dir, gene_id, out_suf, ".csv"))

  NULL
}, phtp = phtp_tab, cvrt = cvrt_tab)

# vim: set ai tw=500:
