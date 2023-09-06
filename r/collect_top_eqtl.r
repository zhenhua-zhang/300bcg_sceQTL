#!/usr/bin/env Rscript
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang@gmail.com, zhenhua.zhang3@helmholtz-hzi.de
# Created: 2022 Jul 02
# Updated: 2023 Mar 02

# Options, packages used in the analysis
options(stringsAsFactors = FALSE, future.globals.maxSize = 10000 * 1024 ^ 2)
library(magrittr)
library(qvalue)

#' Load eQTL summary statistics
#'
#'@description A utility function to load summary statistics from the disk
load_eqtl_tab <- function(indir, fdr = 0.05, max_p_val = NULL, gene_bl = NULL) {
  # Gene black list by default, if the parameter is NULL
  if (is.null(gene_bl)) {
    gene_bl <- c(
      "RPLP0", "RPLP1", "RPLP2", "RPL3", "RPL3L", "RPL5", "RPL6", "RPL7", "RPL7A", "RPL7L1", "RPL8", "RPL9", "RPL10",
      "RPL10A", "RPL10L", "RPL11", "RPL12", "RPL13", "RPL13A", "RPL14", "RPL15", "RPL17", "RPL18", "RPL18A", "RPL19",
      "RPL21", "RPL22", "RPL23", "RPL23A", "RPL24", "RPL26", "RPL26L1", "RPL27", "RPL27A", "RPL28", "RPL29", "RPL30",
      "RPL31", "RPL32", "RPL34", "RPL35", "RPL35A", "RPL36", "RPL36A", "RPL36AL", "RPL37", "RPL37A", "RPL38", "RPL39",
      "RPL39L", "UBA52", "RPL41",
      "RPSA", "RPS2", "RPS3", "RPS3A", "RPS4X", "RPS4Y1", "RPS4Y2", "RPS5", "RPS6", "RPS7", "RPS8", "RPS9", "RPS10",
      "RPS11", "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A", "RPS16", "RPS17", "RPS18", "RPS19", "RPS20", "RPS21",
      "RPS23", "RPS24", "RPS25", "RPS26", "RPS27", "RPS27A", "RPS27L", "RPS28", "RPS29", "FAU",
      "MRPL1", "MRPL2", "MRPL3", "MRPL4", "MRPL9", "MRPL10", "MRPL11", "MRPL12", "MRPL13", "MRPL14", "MRPL15", "MRPL16",
      "MRPL17", "MRPL18", "MRPL19", "MRPL20", "MRPL21", "MRPL22", "MRPL23", "MRPL24", "MRPL27", "MRPL28", "MRPL30",
      "MRPL32", "MRPL33", "MRPL34", "MRPL35", "MRPL36", "MRPL37", "MRPL38", "MRPL39", "MRPL40", "MRPL41", "MRPL42",
      "MRPL43", "MRPL44", "MRPL45", "MRPL46", "MRPL47", "MRPL48", "MRPL49", "MRPL50", "MRPL51", "MRPL52", "MRPL53",
      "MRPL54", "MRPL55", "MRPL58", "MRPL57", "GADD45GIP1", "MRPS30", "MRPS18A", "MRPS2", "MRPS5", "MRPS6", "MRPS7",
      "MRPS9", "MRPS10", "MRPS11", "MRPS12", "MRPS14", "MRPS15", "MRPS16", "MRPS17", "MRPS18C", "MRPS21", "MRPS22",
      "MRPS23", "MRPS24", "MRPS25", "MRPS26", "MRPS27", "MRPS28", "DAP3", "MRPS31", "MRPS33", "MRPS34", "MRPS35",
      "CHCHD1", "AURKAIP1", "PTCD3", "MRPS18B",
      "HLA-A", "HLA-B", "HLA-C", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1", "HLA-DPA2", "HLA-DPA3",
      "HLA-DPB1", "HLA-DPB2", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DQB3", "HLA-DRA", "HLA-DRB1",
      "HLA-DRB2", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-DRB6", "HLA-DRB7", "HLA-DRB8", "HLA-DRB9", "HLA-E",
      "HLA-F", "HLA-G", "HLA-H", "HLA-J", "HLA-K", "HLA-L", "HLA-N", "HLA-P", "HLA-S", "HLA-T", "HLA-U", "HLA-V",
      "HLA-W", "HLA-X", "HLA-Y", "HLA-Z", "HLA-F-AS1"
    )
  }

  # All associations
  all_qtl_tab <- file.path(indir, "qtl_results_all.txt") %>%
    data.table::fread(showProgress = FALSE) %>%
    dplyr::filter(!feature_id %in% gene_bl)

  # eGene, defined using the lowest p-value per feature if the q-value is <= fdr
  egene_tab <- dplyr::group_by(all_qtl_tab, feature_id) %>%
    dplyr::arrange(empirical_feature_p_value) %>%
    dplyr::slice(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(q_value = qvalue::qvalue_truncp(empirical_feature_p_value)$qvalues) %>%
    dplyr::filter(q_value <= fdr) %>%
    dplyr::mutate(QTL = paste(feature_id, snp_id, sep = "-"))

  # eGene q-values
  egene_qval_vec <- dplyr::pull(egene_tab, q_value, QTL)

  # Max p-value to define a significant eVariant
  max_p_by_fdr <- dplyr::select(egene_tab, q_value, empirical_feature_p_value) %>%
    dplyr::slice_max(order_by = q_value) %>%
    dplyr::pull(empirical_feature_p_value) %>%
    head(1)

  # eVariant. A much looser definition can be applied by setting max_p_val
  if (is.null(max_p_val)) {
    evariant_tab <- dplyr::filter(all_qtl_tab, empirical_feature_p_value <= max_p_by_fdr)
  } else {
    cat("[W]: Variants with p-value <=", max_p_val, "are also retained as eVariants.\n")
    evariant_tab <- dplyr::filter(all_qtl_tab, empirical_feature_p_value <= max_p_by_fdr | p_value <= max_p_val)
  }

  dplyr::mutate(evariant_tab, QTL = paste(feature_id, snp_id, sep = "-"), q_value = egene_qval_vec[QTL])
}


proj_dir <- "~/Documents/projects/wp_bcg_eqtl"

mode_vec <- c("normal", "interaction")# , "per_condition") # Two model used in the analysis.
celltype_vec <- c("Monocytes", "CD4T", "CD8T", "NK", "B") # , "pDC", "mDC") # Main cell types.
condition_vec <- c("T0_LPS", "T0_RPMI", "T3m_LPS", "T3m_RPMI") # Four conditions.
comparison_vec <- c("T0_LPS.vs.T0_RPMI", "T3m_LPS.vs.T3m_RPMI", "T3m_RPMI.vs.T0_RPMI")#, "T3m_LPS.vs.T0_RPMI") # Three comparisons.

qtl_tab <- list()
tar_cols<- c(
  "snp_id", "snp_chromosome", "snp_position", "assessed_allele", "maf", "feature_id", "ensembl_gene_id", "p_value",
  "beta", "beta_se", "empirical_feature_p_value", "q_value", "QTL"
)
for (pct in celltype_vec) {
  for (pmd in mode_vec) {
    tar_condition <- if (pmd == "normal") {"common"} else if (pmd == "interaction") {comparison_vec} else if (pmd == "per_condition") {condition_vec}

    for (pcd in tar_condition) {
      runid <- paste(pct, pmd, sep = "_")
      indir <- file.path(proj_dir, "outputs/pseudo_bulk/summary_statistic", pmd, pct)
      if (pmd %in% c("interaction", "per_condition")) {
        indir <- file.path(indir, pcd)
        runid <- paste(runid, pcd, sep = "_")
      }
      qtl_tab[[runid]] <- load_eqtl_tab(indir, 0.05, 5e-5) %>% dplyr::select(dplyr::one_of(tar_cols)) %>% dplyr::mutate(celltype = pct, condition = pcd)
    }
  }
}

Reduce(rbind, qtl_tab) %>%
  dplyr::arrange(snp_chromosome, snp_position, feature_id) %>%
  as.data.frame() %>%
  data.table::fwrite(file.path(proj_dir, "outputs/pseudo_bulk/overview/filtered.eGene_eVariants.FDR0.05.csv"))
