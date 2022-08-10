#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(magrittr)

#' Plot heatmap for functional enrichment by GoShifter
#'
#' @description Function to plot Figure x
plot_goshifter <- function(info_list, save_to = "GoShifter_annotation_enrichment.png", width = 8, height = 4, nrow = 1) {
  cmb_tab <- info_list %>%
    lapply(function(e) {
      tmp_tab <- fread(e["path"])

      obs_enrichment <- tmp_tab$enrichment[1]
      extreme_obs <- sum(tmp_tab$enrichment >= obs_enrichment) - 1
      p_val <- extreme_obs / nrow(tmp_tab)

      data.frame(
        P_val_perm = p_val, Enrichment = obs_enrichment, Cell_type = e["cell_type"],
        Annotation = e["anno_type"], Condition = e["condition"], Background = e["background"]
      )
    }) %>%
    Reduce(rbind, .) %>%
    as.data.table() %>%
    dplyr::mutate(P_val_perm_adj = p.adjust(P_val_perm), Significant = P_val_perm_adj < 0.05)

  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(ggplot2::aes(x = Cell_type, y = Annotation, fill = Enrichment), cmb_tab) +
    ggplot2::geom_text(ggplot2::aes(x = Cell_type, y = Annotation, label = round(P_val_perm, digits = 3)), cmb_tab) +
    ggplot2::geom_point(ggplot2::aes(x = Cell_type, y = Annotation, shape = Significant), (cmb_tab %>% dplyr::filter(Significant))) +
    ggplot2::scale_fill_gradient(low = "white", high = "darkred") +
    ggplot2::facet_wrap(~Background + Condition, nrow = nrow) +
    ggplot2::theme_classic()

  ggplot2::ggsave(save_to, plot = p, width = width, height = height)
}

proj_dir <- "~/Documents/projects/wp_bcg_eqtl"
wkdir <- file.path(proj_dir, "outputs/pseudo_bulk/outcomes/annotation_enrichment")

# Normal model results
enrich_out_list <- list.files(file.path(wkdir, "normal"), pattern = "*.nperm2000.enrich", recursive = TRUE, full.name = TRUE) %>%
  as.list() %>%
  lapply(function(e) {
     cell_type <- e %>% dirname() %>% basename()
     info <- e %>% basename() %>% stringr::str_split("\\.", simplify = TRUE)

     c(path = e, cell_type = cell_type, condition = info[3], anno_type = info[1], background = info[2])
  })
plot_goshifter(enrich_out_list, "GoShifter_annotation_enrichment-normal.png", 12, nrow = 1)


enrich_out_list <- list.files(file.path(wkdir, "interaction"), pattern = "*.enrich", recursive = TRUE, full.name = TRUE) %>%
  as.list() %>%
  lapply(function(e) {
     condition <- e %>% dirname() %>% basename()
     cell_type <- e %>% dirname() %>% dirname() %>% basename()
     info <- e %>% basename() %>% stringr::str_split("\\.", simplify = TRUE)

     c(path = e, cell_type = cell_type, condition = condition, anno_type = info[1], background = info[2])
  })
plot_goshifter(enrich_out_list, "GoShifter_annotation_enrichment-interaction.png", 14, 8)
