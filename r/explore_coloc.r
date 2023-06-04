#!/usr/bin/env Rscript
library(pheatmap)
library(tidyverse)
library(RColorBrewer)


proj_dir <- "~/Documents/projects/wp_bcg_eqtl"
work_dir <- file.path(proj_dir, "outputs/pseudo_bulk/outcomes/normal")

coloc_tab <- list.files(work_dir, pattern = "normal_*", full.names = TRUE) %>%
  file.path("colocalization_test.csv") %>%
  lapply(FUN = function(fpath) {
    ct <- stringr::str_split(fpath, pattern = "/", simplify = TRUE)[1, 11] %>%
    stringr::str_split(pattern = "_", simplify = TRUE)
    tab <- data.table::fread(fpath) %>% dplyr::mutate(celltype = ct[1, 2])
  }) %>%
  Reduce(rbind, .)

coloc_tab <- dplyr::mutate(coloc_tab, outcome = paste0(outcome, " (", id.outcome, ")")) %>%
  dplyr::filter(H4 > 0.5) %>%
  dplyr::select(-c(id.outcome, nsnps, H0, H1, H2, H3)) %>%
  tidyr::pivot_wider(names_from = "outcome", values_from = "H4") %>%
  dplyr::mutate(exposure = paste0(celltype, "_", exposure)) %>%
  replace(is.na(.), 0)

coloc_mat <- dplyr::select(coloc_tab, -celltype) %>%
  (function(tab) {
    row_idx <- tab$exposure
    tab <- dplyr::select(tab, -exposure) %>% as.matrix()
    rownames(tab) <- row_idx
    tab
  })


celltype <- dplyr::select(coloc_tab, celltype) %>% dplyr::mutate(celltype = as.factor(celltype)) %>% as.data.frame()
rownames(celltype) <- rownames(coloc_mat)

save_to <- file.path(proj_dir, "/outputs/pseudo_bulk/outcomes/normal/coloc_summary.pdf")

breaks_list <- seq(0, 20, by = 1)

pdf(save_to)
pheatmap(coloc_mat, annotation_row = celltype, show_rownames = FALSE, clustering_method = "complete",
         color = colorRampPalette(brewer.pal(n = 7, name = "Purples"))(length(breaks_list)))
dev.off()
