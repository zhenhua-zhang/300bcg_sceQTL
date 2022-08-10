#!/usr/bin/env Rscript
library(Seurat)
library(magrittr)
library(tidyverse)
library(data.table)
library(RColorBrewer)


#' Draw a plot to show the cell propoptions
#'
draw_cellprop <- function(obj, discard_cells = c("Undefined", "Platelet"), width = 6, height = 5, save_to = "./") {
  cpp_tab <- obj@meta.data %>%
    dplyr::filter(status == "singlet") %>%
    dplyr::group_by(time, stim, ids) %>%
    dplyr::summarise(cpp = {
      n_pcond <- dplyr::n()
      n_pcell_pcond <- dplyr::cur_data()$clusters1 %>% table()

      data.frame(n_pcell_pcond / n_pcond)
    }) %>%
    as.data.table() %>%
    dplyr::rename("cell_type" = "cpp..", "cell_prop" = "cpp.Freq") %>%
    dplyr::mutate(stim = factor(stim, levels = c("RPMI", "LPS"))) %>%
    dplyr::filter(!cell_type %in% discard_cells)

  g_cpp <- cpp_tab %>%
    ggplot(aes(x = time, y = cell_prop)) +
    geom_boxplot(aes(fill = stim), color = "black", outlier.color = "white") +
    geom_point(aes(color = stim), position = position_dodge2(0.75), size = 0.5) +
    facet_wrap(~ reorder(cell_type, -cell_prop), nrow = 1) +
    scale_fill_manual(values = c("white", "white")) +
    scale_color_manual(values = c("darkblue", "darkred")) +
    theme_classic()

  ggsave(file.path(save_to, "cell_prop.pdf"), g_cpp, width = width, height = height)
}



# Define the number of colors you want
nb_cols <- 10
mycolors <- colorRampPalette(brewer.pal(8, "Blues"))(nb_cols)

proj_dir <- "~/Documents/projects/wp_bcg_eqtl"
so_path <- file.path(proj_dir, "inputs/sc_rnaseq/bcg4-0712.rds")
pbmc <- readRDS(so_path)

Idents(pbmc) <- "clusters1"


# Draw the UMAP
g_umap <- DimPlot(pbmc, reduction = "umap", raster = FALSE)
ggsave(file.path(proj_dir, "sc_rnaseq_umap.pdf"), g_umap, width = 6, height = 6)


# Draw the cell-proportion boxplot
draw_cellprop(pbmc, c("Undefined", "Platelet"), width = 7)
