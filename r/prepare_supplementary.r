#!/usr/bin/env Rscript

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)

  library(ggsci)
  library(ggrepel)

  library(Seurat)
})

cat("Prepare supplementary data ...\n")

# Project directory
proj_dir <- "~/Documents/projects/wp_bcg_eqtl"
save_dir <- file.path(proj_dir, "outputs/pseudo_bulk/overview")

# Individual characteristics
pbmc_obj_file <- file.path(proj_dir, "inputs/sc_rnaseq/bcg4-0712.rds")
pbmc <- readRDS(pbmc_obj_file)

# eQTL replication, including a figure and a table.
cat("Check ./replication.r for more details... \n")
