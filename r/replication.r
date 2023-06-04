#!/usr/bin/env Rscript
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: May 21, 2023
# Updated: May 21, 2023

options(stringsAsFactors = FALSE, data.table.verbose = FALSE, error = NULL)
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
})

proj_dir <- "~/Documents/projects/wp_bcg_eqtl"
time_vec <- c("T0", "T3m")
stim_vec <- c("RPMI", "LPS")
mode_vec <- c("normal", "interaction")
celltype_vec <- c("Monocytes", "CD4T", "CD8T", "NK", "B")
condition_vec <- c("T0_LPS.vs.T0_RPMI", "T3m_LPS.vs.T3m_RPMI", "T3m_RPMI.vs.T0_RPMI") #, "T3m_LPS.vs.T0_RPMI")
col_pattern_vec <- c("normal" = "_Common", "interaction" = ".vs.")
effect_vec <- c("LPS eff. (T0)", "LPS eff. (T3m)", "BCG eff.")
names(effect_vec) <- condition_vec

eqtl_gen_tab <- fread("/vol/projects/BIIM/resources/eqtlGen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz") %>%
  dplyr::select(SNP, GeneSymbol, FDR, AssessedAllele, Zscore)

sceqtl <- file.path(proj_dir, "outputs/pseudo_bulk/summary_statistic/normal", celltype_vec) %>%
  lapply(function(pp) {
    celltype <- basename(pp)
    sceqtl_tab <- fread(file.path(pp, "top_qtl_results_all_FDR.txt")) %>%
      dplyr::select(snp_id, p_value, beta, beta_se, gene_name, assessed_allele, global_corrected_pValue) %>%
      dplyr::group_by(snp_id, gene_name) %>%
      dplyr::summarise(replication = {
        per_gene <- cur_group()$gene_name
        per_snp <- cur_group()$snp_id

        dplyr::filter(eqtl_gen_tab, SNP == per_snp, GeneSymbol == per_gene) %>%
          cbind(dplyr::cur_data())
    }) %>%
      as.data.table() %>%
      dplyr::mutate(celltype = celltype)
  }) %>%
  Reduce(rbind, .) %>%
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/replication/normal.eqtlgen.csv")
fwrite(sceqtl, save_to)

sceqtl %>%
  dplyr::filter(!is.na(replication.FDR), replication.global_corrected_pValue < 0.1) %>%
  dplyr::mutate(is_relicated = dplyr::case_when(
    replication.AssessedAllele == replication.assessed_allele & (replication.beta * replication.Zscore > 0) ~ TRUE,
    replication.AssessedAllele != replication.assessed_allele & (replication.beta * replication.Zscore < 0) ~ TRUE,
    TRUE ~ FALSE
  )) %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(n_replicated = sum(is_relicated))

# # A tibble: 5 x 2, FDR < 0.05
#   celltype  n_replicated
#   <chr>            <int>
# 1 B                  113
# 2 CD4T               218
# 3 CD8T               151
# 4 Monocytes          149
# 5 NK                 132


# # A tibble: 5 x 2, FDR < 0.1
#   celltype  n_replicated
#   <chr>            <int>
# 1 B                  170
# 2 CD4T               293
# 3 CD8T               228
# 4 Monocytes          203
# 5 NK                 173
