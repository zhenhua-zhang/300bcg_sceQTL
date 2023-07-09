#!/usr/bin/env Rscript

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)

  library(ggsci)
  library(ggrepel)
  library(webr) # To plot PieDonut chart

  library(Seurat)
})

cat("Prepare supplementary data ...\n")

# Project directory
proj_dir <- "~/Documents/projects/wp_bcg_eqtl"
save_dir <- file.path(proj_dir, "outputs/pseudo_bulk/overview")

# Individual characteristics
pbmc_obj_file <- file.path(proj_dir, "inputs/sc_rnaseq/bcg4-0712.rds")
pbmc <- readRDS(pbmc_obj_file)
meta_info <- pbmc@meta.data %>%
  dplyr::select(age, gender, ids) %>%
  dplyr::distinct() %>%
  as.data.table()

# Gender pie plot
p_gender <- meta_info %>%
  dplyr::select(ids, gender) %>%
  dplyr::group_by(gender) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(n) %>%
  dplyr::mutate(pos = cumsum(n) - 0.5 * n) %>%
  ggplot() +
  geom_bar(aes(x = 2, y = n, fill = gender), stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(x = 2, y = pos, label = n), color = "white") +
  scale_fill_npg() +
  labs(fill = "Gender") +
  xlim(0.5, 2.5) +
  theme_void() +
  theme(legend.position = "none")

# Age distribution
p_age <- meta_info %>%
  dplyr::select(ids, age, gender) %>%
  dplyr::mutate(gender = dplyr::if_else(gender == "f", "Female", "Male")) %>%
  ggplot() +
  geom_violin(aes(x = gender, y = age, fill = gender), draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_point(aes(x = gender, y = age), position = position_jitter(width = 0.2)) +
  scale_fill_npg() +
  labs(x = NULL, y = "Age", fill = "Gender") +
  theme_classic()

p_demo <- p_gender + p_age
save_to <- file.path(save_dir, "cohort_demographical_plot.pdf")
ggsave(save_to, p_demo, width = 6, height = 3)


# Number of cells
save_to <- file.path(save_dir, "cell_count_plot.pdf")
pdf(save_to, width = 6, height = 6)
pbmc@meta.data %>%
  dplyr::filter(clusters1 != "HSP(T)") %>%
  dplyr::group_by(ts, clusters1) %>%
  dplyr::summarise(n = n()) %>%
  PieDonut(aes(pies = ts, donuts = clusters1, count = n), r0 = 0.2, r1 = 0.8, labelpositionThreshold = 0.04, showPieName = FALSE)
dev.off()

# eQTL replication, including a figure and a table.
cat("Check ./replication.r for more details... \n")
