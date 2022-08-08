#!/usr/bin/env Rscript
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Apr 15, 2022
# Updated: Apr 15, 2022

options(datatable.verbose = FALSE, stringsAsFactors = FALSE)

library(data.table)
library(tidyverse)
library(magrittr)
library(biomaRt)
library(Seurat)


proj_dir <- file.path("~/Documents/projects/wp_bcg_eqtl")
obj <- readRDS(file.path(proj_dir, "inputs/sc_rnaseq/bcg4-0712.rds"))

# Filter cells>
DefaultAssay(obj) <- "RNA"
obj <- subset(obj,
  subset = percent.mt < 10 & nFeature_RNA <= 2.5e3 & nCount_RNA >= 2.5e2
)

# Get meta data.
tar_meta <- c("ids", "age", "gender", "time", "stim", "batch", "pool")
meta_tab <- obj@meta.data %>%
  dplyr::select(one_of(tar_meta)) %>%
  unique() %>%
  as.data.table()

meta_tab_file <- file.path(proj_dir, "outputs/sc_eqtl_lmm/meta_info.csv")
if (!file.exists(meta_tab_file)) {
  meta_tab %>% fwrite(meta_tab_file, quote = FALSE)
}

# Rename cells idents
cm <- c(
  "Monocytes" = "Monocytes", "CD4+ T" = "CD4pT", "CD8+ T" = "CD8pT",
  "B" = "B", "NK" = "NK", "HSP(T)" = "HSP_T", "Platelet" = "Platelet",
  "Undefined" = "Undefined", "pDC" = "pDC", "mDC" = "mDC"
)
clusters2 <- cm[obj@meta.data$clusters1]
names(clusters2) <- colnames(obj)
obj <- AddMetaData(obj, clusters2, col.name = "clusters2")

tar_gene <- obj@meta.data %>%
  dplyr::filter(!(clusters2 %in% c("Undefined", "HSP_T", "Platelet"))) %>%
  dplyr::mutate(barcodes = rownames(.)) %>%
  dplyr::group_by(time, stim, clusters2) %>%
  dplyr::summarise(tar_gene = {
    p_time <- cur_group()$time
    p_stim <- cur_group()$stim
    p_cell <- cur_group()$clusters2
    p_ids <- cur_data()$ids

    save_to <- file.path(
      proj_dir, "outputs/sc_eqtl_lmm/phenotypes",
      paste0(p_cell, "_", p_stim, "_", p_time, ".csv")
    )

    if (!file.exists(save_to)) {
      tar_cells <- cur_data()$barcodes
      tar_genes <- obj[, tar_cells] %>%
        GetAssayData(slot = "counts") %>%
        is_weakly_greater_than(3) %>%
        rowSums() %>%
        divide_by(n()) %>%
        is_greater_than(0.2)

      obj[tar_genes, tar_cells] %>%
        GetAssayData(slot = "counts") %>%
        as.data.frame() %>%
        t() %>%
        as.data.frame() %>%
        dplyr::mutate(time = p_time, stim = p_stim, ids = p_ids) %>%
        fwrite(save_to, quote = FALSE)
    }

    names(tar_genes)[tar_genes]
  }) %$%
  tar_gene %>%
  unique()


# Gene coordination.
# Apr 15, 2022, current using GRCh38. Only autosomes are included for now.
ensembl <- useEnsembl("ensembl", "hsapiens_gene_ensembl")
gene_tab <- getBM(
  c("chromosome_name", "start_position", "end_position", "hgnc_symbol"),
  c("chromosome_name", "hgnc_symbol"),
  list(1:22, tar_gene),
  ensembl
)

# Genotypes
gntp_file <- file.path(
  proj_dir, "inputs/genotypes/300BCG_sub40_imp_hg38_ids.vcf.gz"
)
gntp <- fread(gntp_file, skip = "#CHROM", tmpdir = dirname(gntp_file)) %>%
  dplyr::select(!one_of("QUAL", "FILTER", "INFO", "FORMAT")) %>%
  dplyr::rename("CHROM" = "#CHROM") %>%
  dplyr::mutate(
    CHROM = as.integer(str_remove(CHROM, "chr")),
    POS = as.integer(POS),
    # keep dosage
    dplyr::across(
      dplyr::starts_with("300BCG"),
      .fns = ~ str_split(.x, ":", simplify = TRUE)[, 2] %>% as.numeric()
    )
  ) %>%
  dplyr::rowwise() %>%
  dplyr::filter(min(table(round(c_across(starts_with("300BCG"))))) >= 3)

gene_tab %>%
  dplyr::group_by(chromosome_name) %>%
  summarise(n_genes = n())

# Genotypes (SNPs) per gene per chromosome.
gene_tab %>%
  dplyr::group_by(chromosome_name, hgnc_symbol) %>%
  dplyr::mutate(
    start_position = as.integer(start_position),
    end_position = as.integer(end_position)
  ) %>%
  summarise(test = {
    cgroup <- cur_group()
    chrom <- cgroup$chromosome_name
    geneid <- cgroup$hgnc_symbol

    chrom_dir <- file.path(proj_dir, "outputs/sc_eqtl_lmm/genotypes", chrom)
    if (!dir.exists(chrom_dir)) dir.create(chrom_dir, recursive = TRUE)

    sub_gntp_file <- file.path(
      proj_dir, "outputs/sc_eqtl_lmm/genotypes", chrom, paste0(geneid, ".csv")
    )
    start <- cur_data()$start_position[1] - 2.5e5
    end <- cur_data()$end_position[1] + 2.5e5

    if (!file.exists(sub_gntp_file)) {
      gntp %>%
        dplyr::filter(CHROM == chrom & start <= POS & POS <= end) %>%
        dplyr::select(-c(CHROM, POS, REF, ALT)) %>%
        dplyr::rename("SNPs" = "ID") %>%
        dplyr::mutate(GeneID = geneid) %>%
        (function(d) if (nrow(d)) fwrite(d))
    }

    sub_gntp_file
  })

# vim: set ai tw=200:
