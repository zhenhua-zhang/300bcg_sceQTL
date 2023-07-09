# File: misc.r
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Jun 10, 2023
# Updated: Jun 10, 2023

options(stringsAsFactors = FALSE, data.table.verbose = FALSE, error = NULL)
suppressPackageStartupMessages({
  library(DropletUtils)
  library(data.table)
  library(tidyverse)
  library(Seurat)
})

# Subset of PBMC
proj_dir <- "~/Documents/projects/wp_bcg_eqtl"
pbmc <- readRDS(file.path(proj_dir, "/inputs/sc_rnaseq/bcg4-0712.rds"))
tar_ids <- pbmc@meta.data %>%
  dplyr::select(ids, age, gender, time, stim) %>%
  dplyr::filter(stim == "RPMI", time == "T0") %>%
  dplyr::group_by(ids, age, gender) %>%
  dplyr::summarize(n = n()) %>%
  dplyr::filter(20 <= age, age <= 25, n > 1600) %>%
  dplyr::pull(ids)

tar_cells <- pbmc@meta.data %>%
  dplyr::filter(ids %in% tar_ids) %>%
  rownames()

sub_pbmc <- pbmc[, tar_cells]
write10xCounts(sub_pbmc@assays$RNA@counts, path = file.path(proj_dir, "misc/5samples_300BCG_baseline"), version = "3")

sub_pbmc@meta.data %>%
  dplyr::mutate(cell_barcodes = rownames(.)) %>%
  dplyr::select(orig.ident, nCount_RNA, nFeature_RNA, cell_barcodes, status, ids, age, gender, batch, pool) %>%
  fwrite(file.path(proj_dir, "temps/5samples_300BCG_baseline/meta.csv"))


# Check the summary statistics of COVID-19 GWAS from https://doi.org/10.1371/journal.pone.0279356
sumstat_file <- "/vol/projects/zzhang/projects/wp_bcg_eqtl/inputs/COVID19_summary_statistics/Severe_SNPs.txt"
sumstat_tab <- fread(sumstat_file)

# P-values unadjusted
plot_y_lab <- "GWAS -log10(p-value), raw"

# P-values adjusted by SPA
p <- sumstat_tab %>%
  dplyr::filter(pvalue_severe < 5e-2, Chrom != "X") %>%
  dplyr::mutate(Chrom = factor(Chrom, levels = 1:22)) %>%
  ggplot() +
  geom_point(aes(x = Position, y = -log10(pvalue_severe), group = Chrom, color = Chrom), position = "stack", size = 0.2) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", linewidth = 0.5, color = "red") +
  geom_hline(yintercept = -log10(5e-7), linetype = "dashed", linewidth = 0.5, color = "blue") +
  geom_hline(yintercept = -log10(5e-6), linetype = "dashed", linewidth = 0.5, color = "black") +
  facet_grid(~ Chrom, scales = "free_x", space="free", switch = "x") +
  labs(x = "Chromosome", y = plot_y_lab) +
  lims(y = c(-log10(5e-2), NA)) +
  theme_classic() +
  theme(
    panel.spacing = unit(0, "mm"), legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    strip.text.x.bottom = element_text(angle = 90, hjust = 0), strip.background = element_blank(), strip.placement = "outside"
  )
plot_save_to <- file.path(proj_dir, "temps/COVID19_GWAS_raw.pdf")
ggsave(plot_save_to, plot = p, width = 12, height = 4)
