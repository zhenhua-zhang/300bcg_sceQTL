#!/usr/bin/env Rscript
# File: atac_seq.r
# Author: zhenhua zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Aug 30, 2023
# Updated: Aug 30, 2023

options(stringsAsFactors = FALSE, datatable.showProgress = FALSE, datatable.verbose = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(DESeq2)
  library(Seurat)
})

# Analysis steps:
# 1. Normalization (library size)
# 2. Peak-to-gene links, ATAC-seq peak (bulk or sorted cells) to gene expression linkage (pseudo-bulk, including CD55 and TFs).
# 3. acQTL in PBMC and sorted cells (monocytes, CD8+ T, NK cells).

proj_dir <- "~/Documents/projects/wp_bcg_eqtl"
time_vec <- c("T0", "T3m")
stim_vec <- c("RPMI", "LPS")
mode_vec <- c("normal", "interaction")
celltype_vec <- c("Monocytes", "CD4T", "CD8T", "NK", "B")
condition_vec <- c("T0_LPS.vs.T0_RPMI", "T3m_LPS.vs.T3m_RPMI", "T3m_RPMI.vs.T0_RPMI") #, "T3m_LPS.vs.T0_RPMI")
col_pattern_vec <- c("normal" = "Common", "interaction" = ".vs.")
effect_vec <- c("LPS eff. (T0)", "LPS eff. (T3m)", "BCG eff.")
names(effect_vec) <- condition_vec

# PBMC RNA-seq
so_path <- file.path(proj_dir, "inputs/sc_rnaseq/bcg4-0712.rds")
pbmc <- readRDS(so_path)
DefaultAssay(pbmc) <- "RNA"
pbmc <- pbmc[, (!pbmc@meta.data$clusters1 %in% c("HSP(T)"))] # Removing HSP(T) cells
new_cell_lvl <- c("Monocytes", "CD4+ T", "CD8+ T", "NK", "B", "Platelet", "mDC", "pDC", "Undefined")
pbmc@meta.data$clusters1 <- factor(pbmc@meta.data$clusters1, levels = new_cell_lvl)
pbmc@meta.data$ts <- factor(pbmc@meta.data$ts, levels = c("T0_RPMI", "T0_LPS", "T3m_RPMI", "T3m_LPS"))
Idents(pbmc) <- "clusters1"

group_idx <- c("clusters1", "time", "stim", "ids")
tar_cells <- pbmc@meta.data$clusters1 %in% c("Monocytes", "CD4+ T", "CD8+ T", "NK", "B")
tar_features <- c("CD55", "SPI1", "STAT1", "ESR1", "CEBPA", "KLF5", "TBX3", "HLF", "FOXP1", "RARA", "CTCF")
# Gene expression
exp_mat <- pbmc[tar_features, tar_cells] %>%
  AverageExpression(assays = "RNA", group.by = group_idx) %>%
  as.data.frame() %>%
  dplyr::mutate(feature_id = rownames(.)) %>%
  dplyr::rename_with(starts_with("RNA."), .fn = ~ str_remove_all(.x, "RNA|\\.")) %>%
  dplyr::relocate(feature_id) %>%
  tidyr::pivot_longer(-feature_id, names_to = "sample", values_to = "expression") %>%
  tidyr::separate(sample, into = c("celltype", "time", "stim", "ids"), sep = "_") %>%
  dplyr::filter(time == "T0", stim == "RPMI") %>%
  dplyr::select(-time, -stim)

# Genotype
gt_per_ind <- fread(file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/tar_snps.genotypes.csv")) %>%
  dplyr::rename_with(.fn = ~ stringr::str_remove_all(.x, "\\[[0-9]{1,2}\\]|:GT$|^# ")) %>%
  tidyr::pivot_longer(dplyr::starts_with("300"), values_to = "GT_01", names_to = "ids") %>%
  dplyr::mutate(GT_RA = dplyr::case_when(GT_01 %in% c("0|0") ~ paste0(REF, REF), GT_01 %in% c("1|0", "0|1") ~ paste0(REF, ALT), GT_01 %in% c("1|1") ~ paste0(ALT, ALT))) %>%
  dplyr::select(ID, CHROM, POS, REF, ALT, GT_RA, ids) %>%
  tidyr::pivot_wider(names_from = c(ID, CHROM, POS, REF, ALT), values_from = GT_RA)

# Meta data
meta_data <- pbmc@meta.data %>% as.data.frame() %>% dplyr::select(ids, age, gender) %>% dplyr::distinct()

# PBMC ATAC-seq readcounts
pbmc_rcmat <- file.path(proj_dir, "inputs/atac_seq", paste0("quantification_filtered_PBMC.csv.gz")) %>% fread()
pbmc_coldata <- colnames(pbmc_rcmat) %>%
  purrr::discard(~.x == "ID") %>%
  tibble::as_tibble() %>%
  tidyr::separate(value, into = c("donor", "timepoint", "source"), sep = "_", remove = FALSE) %>%
  dplyr::filter(timepoint == "V1") %>%
  (function(tab) {tab <- as.data.frame(tab); rownames(tab) <- tab$value; tab$value <- NULL; as.matrix(tab)})

pbmc_rcmat <- pbmc_rcmat %>% 
  dplyr::select(dplyr::all_of(c("ID", rownames(pbmc_coldata)))) %>%
  (function(tab) {tab <- as.data.frame(tab); rownames(tab) <- tab$ID; tab$ID <- NULL; as.matrix(tab)})
pbmc_dds <- DESeqDataSetFromMatrix(countData = pbmc_rcmat, colData = pbmc_coldata, design = ~ 1)
pbmc_rcmat_norm <- pbmc_dds[rowSums(counts(pbmc_dds)) >= 10, ] %>% normTransform() %>% assay() %>% as.data.frame() %>% dplyr::mutate(ID = rownames(.)) %>% dplyr::relocate(ID)

# Peak annotation
pbmc_pkmat <- file.path(proj_dir, "inputs/atac_seq/peaks_filtered.csv.gz") %>% fread() %>% dplyr::filter(gene_name %in% c("CD55")) 

# Joining data for regression
pbmc_regtab <- dplyr::inner_join(pbmc_pkmat, pbmc_rcmat_norm, by = c("peak_id" = "ID")) %>%
  tidyr::pivot_longer(dplyr::starts_with("300BCG"), names_to = "sample", values_to = "readdepth") %>%
  tidyr::separate(sample, into = c("ids"), extra = "drop", sep = "_") %>%
  dplyr::inner_join(exp_mat, by = "ids") %>%
  dplyr::inner_join(gt_per_ind, by = "ids") %>%
  dplyr::inner_join(meta_data, by = "ids")

# Correlation by linear regression
pbmc_regres <- pbmc_regtab %>%
  dplyr::filter(feature_id == "CD55") %>%
  dplyr::group_by(peak_id, feature_id, celltype) %>%
  dplyr::summarize(reg = {
    data <- dplyr::cur_data()
    m <- lm(expression ~ readdepth + age + gender, data = data)

    summary(m)$coefficients %>%
      as.data.frame() %>%
      dplyr::mutate(var = rownames(.)) %>%
      dplyr::select(var, estimate = Estimate, stderr = `Std. Error`, t_stat = `t value`, p_value = `Pr(>|t|)`)
  }) %>%
  tidyr::unnest(reg)

pbmc_regres %>% dplyr::filter(var == "readdepth", p_value < 0.05) %>% as.data.frame



# # Monocytes ATAC-seq readcounts, No overlapped samples
# mono_rctab <- file.path(proj_dir, "inputs/atac_seq", paste0("quantification_filtered_cd8t_monocyte_nkcell.csv.gz")) %>%
#   fread() %>%
#   dplyr::select(ID, dplyr::contains("V1_monocyte"))
# 
# mono_coldata <- colnames(mono_rctab) %>%
#   purrr::discard(~.x == "ID") %>%
#   tibble::as_tibble() %>%
#   tidyr::separate(value, into = c("donor", "timepoint", "source"), sep = "_", remove = FALSE) %>%
#   dplyr::filter(timepoint == "V1") %>%
#   (function(tab) {tab <- as.data.frame(tab); rownames(tab) <- tab$value; tab$value <- NULL; as.matrix(tab)})
# 
# mono_rcmat <- mono_rctab %>% 
#   dplyr::select(dplyr::all_of(c("ID", rownames(mono_coldata)))) %>%
#   (function(tab) {tab <- as.data.frame(tab); rownames(tab) <- tab$ID; tab$ID <- NULL; as.matrix(tab)})
# 
# mono_dds <- DESeqDataSetFromMatrix(countData = mono_rcmat, colData = mono_coldata, design = ~ 1)
# mono_rcmat_norm <- mono_dds[rowSums(counts(mono_dds)) >= 10, ] %>% normTransform() %>% assay() %>% as.data.frame() %>% dplyr::mutate(ID = rownames(.)) %>% dplyr::relocate(ID)
# 
# mono_pkmat <- file.path(proj_dir, "inputs/atac_seq/peaks_filtered_monocyte.csv.gz") %>% fread() %>% dplyr::filter(gene_name %in% c("CD55")) 
# 
# mono_regtab <- dplyr::inner_join(mono_pkmat, mono_rcmat_norm, by = c("peak_id" = "ID")) %>%
#   tidyr::pivot_longer(dplyr::starts_with("300BCG"), names_to = "sample", values_to = "readdepth") %>%
#   tidyr::separate(sample, into = c("ids"), extra = "drop", sep = "_") %>%
#   dplyr::inner_join(exp_mat, by = "ids") %>%
#   dplyr::inner_join(gt_per_ind, by = "ids") %>%
#   dplyr::inner_join(meta_data, by = "ids")
# 
# mono_regres <- mono_regtab %>%
#   dplyr::filter(feature_id == "CD55") %>%
#   dplyr::group_by(peak_id, feature_id) %>%
#   dplyr::summarize(reg = {
#     data <- dplyr::cur_data()
#     m <- lm(expression ~ readdepth + age + gender, data = data)
# 
#     summary(m)$coefficients %>%
#       as.data.frame() %>%
#       dplyr::mutate(var = rownames(.)) %>%
#       dplyr::select(var, estimate = Estimate, stderr = `Std. Error`, t_stat = `t value`, p_value = `Pr(>|t|)`)
#   }) %>%
#   tidyr::unnest(reg)
# 
# mono_regres %>% dplyr::filter(var == "readdepth", p_value < 0.05) %>% as.data.frame
