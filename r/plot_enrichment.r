#!/usr/bin/env Rscript
library(tidyverse)
library(pheatmap)
library(gplots)
library(RColorBrewer)

proj_dir <- "~/Documents/projects/wp_bcg_eqtl"
work_dir <- file.path(proj_dir, "outputs/pseudo_bulk/outcomes/annotation_enrichment")

enrtab <- list.files(work_dir, pattern = "*.csv", recursive = TRUE, full.names = TRUE) %>%
  lapply(function(e) {
    qtl_from <- stringr::str_split(basename(e), pattern = "\\.", simplify = TRUE)[1]
    data.table::fread(e) %>% dplyr::mutate(qtl_from = qtl_from)
  }) %>%
  Reduce(rbind, .)

# Roadmap annotations
annodat <- enrtab %>% dplyr::filter(data_source == "roadmap")
annoplot <- ggplot(annodat, aes(x = annotype, y = celltype)) +
  geom_tile(aes(fill = log2od), stat = "identity") +
  facet_wrap(~qtl_from, ncol = 1) +
  labs(x = "Annotation type", y = "Reference cell type") +
  scale_fill_gradient2(name = "Log2(odds-ratio)", low = "blue", mid = "white", high = "red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

saveto <- file.path(work_dir, "normal/figure/roadmap_annot.pdf")
ggsave(saveto, plot = annoplot, height = 10)

# 300BCG ATAC-seq annotations
bcgdat <- dplyr::filter(enrtab, data_source == "300BCG_ATACseq")
annoplot <- ggplot(bcgdat, aes(x = annotype, y = celltype)) +
  geom_tile(aes(fill = log2od), stat = "identity") +
  facet_wrap(~qtl_from, nrow = 1) +
  labs(x = "Annotation type", y = "Reference cell type") +
  scale_fill_gradient2(name = "Log2(odds-ratio)", low = "blue", mid = "white", high = "red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
saveto <- file.path(work_dir, "normal/figure/300bcg_atacseq_annot.pdf")
ggsave(saveto, plot = annoplot, width = 7, height = 3)

# GWAS SNPs
gwastab <- enrtab %>% dplyr::filter(data_source == "ieuopengwas") %>%
  dplyr::select(c("Trait" = "celltype", "QTL" = "qtl_from", "Log2OD" = "log2od")) %>%
  tidyr::pivot_wider(id_cols = QTL, values_from = Log2OD, names_from = Trait) %>%
  replace(is.na(.), 0) %>%
  as.data.frame() %>%
  (function(dat) {rownames(dat) <- dat$QTL; return(dat[-1])})

blist <- seq(-3.5, 3.5, by = 0.25)
cmap <- colorRampPalette(c("navy", "white", "red"))(length(blist))
saveto <- file.path(work_dir, "normal/figure/ieuopengwas_gwas.pdf")
pdf(saveto, width = 7, height = 3.5)
pheatmap(gwastab, color = cmap, breaks = blist, cluster_rows = FALSE)
dev.off()
