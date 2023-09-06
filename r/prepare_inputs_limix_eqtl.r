#!/usr/bin/env Rscript
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: 2022 May 15
# Updated: 2022 May 20
# File: prepare_inputs_limix_eqtl.r
options(error = NULL)

library(peer)
library(Seurat)
library(biomaRt)
library(magrittr)
library(tidyverse)
library(data.table)



#' Add PEERs to the covariates matrix
#'
#' @param cov_path Path to the covariates table.
#' @param exp_path Path to the expression table. Row is gene, colum is sample.
#' @param n_peers Nr. of PEERs speculated for the model. Default is 10.
#' @param n_cov_peers Nr. of PEERs to be added to cov matrix. Default is 5.
#' @return The covariate matrix with PEERs.
add_peer <- function(cov_path, exp_path, n_peers = 10, n_cov_peers = 5, n_iters = 5000) {
  stopifnot(n_peers > n_cov_peers)

  # Expression matrix
  exp_tab <- data.table::fread(exp_path, header = TRUE)
  exp_mat <- exp_tab %>% dplyr::select(-feature_id) %>% as.matrix() %>% t()

  # Covariates matrix
  cov_tab <- data.table::fread(cov_path, header = TRUE)
  cov_mat <- cov_tab %>% dplyr::select(-sample_id) %>% as.matrix()

  # The PEER model
  model <- peer::PEER()

  # Customized parameters
  ## Maximum number of interations to estimate the parameters.
  peer::PEER_setNmax_iterations(model, n_iters)

  ## Limiting values (tolerances)
  # peer::PEER_setTolerance(model, 1)
  # peer::PEER_setVarTolerance(model, 0.1)

  ## Prior parameters
  # peer::PEER_setPriorAlpha(model, 0.001, 0.1)
  # peer::PEER_setPriorEps(model, 0.1, 10.)

  # set the observed data.
  peer::PEER_setPhenoMean(model, exp_mat)

  # Speculate a 10 hidden confounders.
  peer::PEER_setNk(model, n_peers)
  peer::PEER_getNk(model)

  # Add covariates
  peer::PEER_setCovariates(model, cov_mat)

  # Perform the inference
  peer::PEER_update(model)

  # Add the mean expression as additional factor (covariate)
  peer::PEER_setAdd_mean(model, TRUE)

  # Get statistic for the posterior matrix
  fct_mat <- peer::PEER_getX(model)
  .wgt_mat <- peer::PEER_getW(model)
  .res_mat <- peer::PEER_getResiduals(model)
  .prc_mat <- peer::PEER_getAlpha(model)

  # Save the factors to the disk
  cov_tab %>%
    base::cbind(fct_mat[, 1:n_cov_peers]) %>%
    dplyr::rename_with(
      dplyr::starts_with("V"),
      .fn = ~ stringr::str_replace(.x, "V", "PEER_")
    )
}



# The mode to prepare phenotypes, covariants, and sample map.
# mode <- "pseudo_time"
mode <- "pseudo_bulk"

# Chunk size used to split the genome
chunk_size <- 1000000

# Maximum zeros ratios
max_zero_ratio <- 0.3

# Project main folder
proj_dir <- "~/Documents/projects/wp_bcg_eqtl"



#
## Prepare genotypes
#
cat(
  "Check more information at",
  file.path(proj_dir, "scripts/bash/prepare_kinship_matrix.sh"), "\n"
)


# Load scRNA-seq results
pbmc_file <- file.path(proj_dir, "inputs/sc_rnaseq/bcg4-0712.rds")
pbmc <- readRDS(pbmc_file)

DefaultAssay(pbmc) <- "RNA"

meta_data <- pbmc@meta.data

# Samples per condition per stimulation per cell type
meta_data %>%
  dplyr::select(time, stim, clusters1, ids) %>%
  dplyr::filter(!clusters1 %in% c("HSP(T)", "Platelet", "Undefined", "mDC", "pDC")) %>%
  dplyr::distinct() %>%
  dplyr::group_by(clusters1, time, stim) %>%
  dplyr::summarise(n = n()) %>%
  print(n = 40)


#
## Prepare annotation file
#
tar_features <- pbmc %>% rownames() %>% unique()

save_to <- file.path(proj_dir, "inputs/annotations/annotations_hg38.tsv")
if (!file.exists(save_to)) {
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  feature_tab <- getBM(
    attributes = c(
      "external_gene_name", "chromosome_name", "transcript_start",
      "transcript_end", "ensembl_gene_id", "external_gene_name", "strand"
    ),
    filters = c(
      "chromosome_name", "transcript_is_canonical", "external_gene_name"
    ),
    values = list(as.character(1:22), TRUE, tar_features),
    mart = ensembl
  ) %>%
    dplyr::rename(
      "feature_id" = "external_gene_name", "chromosome" = "chromosome_name",
      "start" = "transcript_start", "end" = "transcript_end",
      "gene_name" = "external_gene_name.1", "feature_strand" = "strand"
    ) %>%
    dplyr::mutate(feature_strand = if_else(feature_strand == -1, "-", "+")) %>%
    dplyr::filter(feature_id != "", !duplicated(feature_id))

  feature_tab %>% fwrite(save_to, sep = "\t")
} else {
  feature_tab <- fread(save_to)
}



#
## Prepare phenotypes, covariates, and sample map by cell type, pseudo-bulk
#
common_features <- feature_tab$feature_id %>% unique()

# These cell types will be removed.
# TODO: HSP(T) need more information
celltype_blacklist <- c("HSP(T)", "Platelet", "mDC", "pDC", "Undefined")

if (mode == "single_cell") {
  stop("Single-cell mode is not ready!!")
  ## Single-cell mode
  meta_data %>%
    dplyr::mutate(cell_barcode = rownames(.)) %>%
    dplyr::group_by(clusters1, time, stim) %>%
    dplyr::summarise(emat = {
      group_names <- dplyr::cur_group()
      cell_type <- str_remove_all(group_names$clusters1, " |\\+")
      time_point <- group_names$time
      stimulation <- group_names$stim

      save_to <- paste0(cell_type, "_", time_point, "_", stimulation, ".tsv")

      # 1. Phenotypes, expression matrix by cells.
      tar_cb <- cur_data()$cell_barcode
      exp_mat <- pbmc[common_features, tar_cb]@assays$RNA@data %>%
        as.data.frame() %>%
        dplyr::mutate(feature_id = rownames(.)) %>%
        dplyr::relocate(feature_id) %>%
        fwrite(
          file.path(proj_dir, "inputs", mode, "phenotypes", save_to),
          sep = "\t"
        )

      # 2. Covariates, age and gender are included one per cell.
      cov_mat <- cur_data() %>%
        dplyr::select(cell_barcode, age, gender) %>%
        dplyr::mutate(gender = if_else(gender == "f", 1, 0)) %>%
        dplyr::rename("sample_id" = "cell_barcode") %>%
        fwrite(
          file.path(proj_dir, "inputs", mode, "covariates", save_to),
          sep = "\t"
        )

      # 3. Sample mapping, genotype_id -> cell barcode.
      smp_map <- cur_data() %>%
        dplyr::select(ids, cell_barcode) %>%
        fwrite(
          file.path(proj_dir, "inputs", mode, "sample_mapping", save_to),
          col.names = FALSE, sep = "\t"
        )
    })
} else {
  if (mode != "pseudo_bulk") {
    cat("Unknown mode ", mode, ". Using 'pseudo-bulk' by default.\n", sep = "")
  }

  group_idx <- c("clusters1", "time", "stim", "ids")
  tar_cells <- !meta_data$clusters1 %in% celltype_blacklist

  # 1. Phenotypes, average expression matrix per individual
  exp_mat <- pbmc[common_features, tar_cells] %>%
    AverageExpression(assays = "RNA", group.by = group_idx) %>%
    as.data.frame() %>%
    dplyr::mutate(feature_id = rownames(.)) %>%
    dplyr::rename_with(starts_with("RNA."), .fn = ~ str_remove_all(.x, "RNA|\\.")) %>%
    dplyr::mutate(zero_ratio = rowSums(. == 0) / dim(.)[2]) %>%
    dplyr::filter(zero_ratio < max_zero_ratio) %>%
    dplyr::relocate(feature_id)

  # 2. Covariates, age and gender.
  cov_mat <- meta_data %>%
    as.data.table() %>%
    dplyr::filter(!clusters1 %in% celltype_blacklist) %>%
    dplyr::select(one_of(group_idx), age, gender) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      clusters1 = str_remove_all(clusters1, "[ +]"),
      gender = if_else(gender == "m", 0, 1),
      sample_id = str_c(clusters1, time, stim, ids, sep = "_"),
      time = if_else(time == "T0", 0, 1),
      stim = if_else(stim == "RPMI", 0, 1),
    ) %>%
    dplyr::select(-c(clusters1, ids)) %>%
    dplyr::relocate(sample_id)
  cov_mat %>% head()

  # 3. Sample mapping. genotype_id -> sample_id
  smp_map <- meta_data %>%
    as.data.table() %>%
    dplyr::filter(!clusters1 %in% celltype_blacklist) %>%
    dplyr::select(one_of(group_idx)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      clusters1 = str_remove_all(clusters1, "[ +]"),
      sample_id = str_c(clusters1, time, stim, ids, sep = "_")
    ) %>%
    dplyr::select(ids, sample_id)
  smp_map %>% head()

  # Non-celltype phenotypes, covariants, and sample mapping.
  save_to <- file.path(proj_dir, "inputs", mode, "all")
  if (!dir.exists(save_to)) dir.create(save_to, recursive = TRUE)

  fwrite(exp_mat, file.path(save_to, "phenotypes.tsv"), sep = "\t")
  fwrite(cov_mat, file.path(save_to, "covariates.tsv"), sep = "\t")
  fwrite(smp_map, file.path(save_to, "sample_mapping.tsv"), col.names = FALSE, sep = "\t")


  # Per cell type phenotypes, covariates, and sample mapping.
  tar_celltype <- exp_mat %>%
    colnames() %>%
    purrr::discard(~ .x %in% c("feature_id", "zero_ratio")) %>%
    str_split("_", 2, simplify = TRUE) %>%
    as.data.frame() %$%
    V1 %>%
    unique()

  . <- tar_celltype %>%
    purrr::map(.f = function(x) {
      save_to <- file.path(proj_dir, "inputs", mode, x)
      if (!dir.exists(save_to)) dir.create(save_to, recursive = TRUE)

      dplyr::select(exp_mat, feature_id, starts_with(x)) %>%
        fwrite(file.path(save_to, "phenotypes.tsv"), sep = "\t")

      dplyr::filter(cov_mat, str_starts(sample_id, x)) %>%
        fwrite(file.path(save_to, "covariates.tsv"), sep = "\t")

      dplyr::filter(smp_map, str_starts(sample_id, x)) %>%
        fwrite(file.path(save_to, "sample_mapping.tsv"), col.names = FALSE, sep = "\t")
    })
}



#
## Add PEERs to the covariate matrix.
#
for (ct in c("Monocytes", "CD4T", "CD8T", "NK", "B", "all")) {
  base_path <- file.path(proj_dir, "inputs/pseudo_bulk", ct)
  exp_path <- file.path(base_path, "phenotypes.tsv")
  cov_path <- file.path(base_path, "covariates.tsv")

  cov_wp_tab <- add_peer(cov_path, exp_path)
  cov_wp_tab %>%
    data.table::fwrite(
      file.path(base_path, "covariates_wpeer.tsv"), sep = "\t"
    )
}



#
## Prepare chunk file, based on the annotation file
#
# version 1
feature_tab %>%
  dplyr::select(chromosome, start, end) %>%
  dplyr::group_by(chromosome) %>%
  dplyr::summarise(chunk = {
    min_base <- cur_data()$start %>% min()
    max_base <- cur_data()$end %>% max()

    n_chunks <- as.integer((max_base - min_base) / chunk_size)
    cut(start, n_chunks)
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(chunk = str_remove_all(chunk, "\\(|\\]")) %>%
  tidyr::separate(chunk, into = c("start", "end"), sep = ",") %>%
  dplyr::mutate(
    start = as.integer(start),
    start = ifelse(start >= 0, start, 0),
    end = as.integer(end),
    chunk = paste0(
      chromosome, ":", format(start, trim = TRUE, scientific = FALSE), "-",
      format(end, trim = TRUE, scientific = FALSE)
    )
  ) %>%
  dplyr::arrange(chromosome, start, end) %>%
  dplyr::select(chunk) %>%
  unique() %>%
  fwrite(file.path(proj_dir, "inputs/chunks_file.txt"), col.names = FALSE)

# version 2
feature_tab %>%
  dplyr::select(chromosome, start, end) %>%
  dplyr::arrange(chromosome, start, end) %>%
  dplyr::mutate(index = dplyr::row_number()) %>%
  dplyr::group_by(chromosome) %>%
  dplyr::summarise(chunk = {
    min_base <- cur_data()$start %>% min()
    max_base <- cur_data()$end %>% max()

    n_chunks <- as.integer((max_base - min_base) / chunk_size)
    cut(start, n_chunks)
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(chunk = str_remove_all(chunk, "\\(|\\]")) %>%
  tidyr::separate(chunk, into = c("start", "end"), sep = ",") %>%
  dplyr::mutate(
    start = as.integer(start),
    start = ifelse(start >= 0, start, 0),
    end = as.integer(end),
    chunk = paste0(
      chromosome, ":", format(start, trim = TRUE, scientific = FALSE), "-",
      format(end, trim = TRUE, scientific = FALSE)
    )
  ) %>%
  dplyr::arrange(chromosome, start, end) %>%
  dplyr::select(chunk) %>%
  unique() %>%
  fwrite(file.path(proj_dir, "inputs/chunks_file_v2.txt"), col.names = FALSE)

# vim: set ai tw=200:
