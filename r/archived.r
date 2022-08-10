# FIXME: Using MASH to estimate the shared signals
# cat("[I]: Estimating correlations for", mode, "eQTL ...\n")
# qtlcor_tab <- make_cor_tab(qtltab_list, save_to = save_to) %>%
#   dplyr::mutate(x = stringr::str_remove(x, paste0(mode, "_")), y = stringr::str_remove(y, paste0(mode, "_")))
# cat("[I]: Drawing heatmap for the correlation of", mode, "eQTL ...\n")
# draw_heatmap(qtlcor_tab, save_to, 3.25, 4.5)

#' Combine the QTL tables for a heatmap plot
#'
#'@description Function to plot Figure x
make_cor_tab <- function(qtltabs, overwrite = TRUE, save_to = "./") {
  if (is.null(names(qtltabs))) stop("The qtltabs should be a named list")

  tab_path <- file.path(save_to, "shared_QTL_correlation_table.csv")
  if (file.exists(tab_path) && !overwrite) {
    cat("[I]: Found", tab_path, "so loading it from the disk.\n")
    cor_tab <- data.table::fread(tab_path)
  } else {
    cor_tab <- names(qtltabs) %>%
      combn(2) %>%
      t() %>%
      as.data.frame() %>%
      apply(1, function(pair, .qtltabs) {
        tab_x <- pair[1]
        tab_y <- pair[2]

        cor_res <- .qtltabs[[tab_x]] %>%
          dplyr::inner_join(.qtltabs[[tab_y]], by = c("snp_id", "ensembl_gene_id")) %>%
          dplyr::select(dplyr::matches("beta.[xy]")) %>%
          (function(df) cor.test(df$beta.x, df$beta.y))

        data.frame(x = tab_x, y = tab_y, p_value = cor_res$p.value, correlation = cor_res$estimate)
      }, .qtltabs = qtltabs) %>%
      Reduce(rbind, .) %>%
      data.table::as.data.table() %>%
      (function(df) {
        comp <- c(df$x, df$y) %>% unique()
        ext_df <- data.frame(x = comp, y = comp, p_value = 0, correlation = 1)
        rbind(df, ext_df)
      })
    cor_tab %>% data.table::fwrite(tab_path, quote = FALSE)
  }

  cor_tab
}


#' Plot a heatmap of share eQTL
#'
#' @description Function to plot Figure x
draw_heatmap <- function(cortab, save_to = "./", height = 7, width = 7, override = TRUE) {
  fig_path <- file.path(save_to, "shared_QTL_correlation_heatmap.png")
  if (!file.exists(fig_path) || override) {
    ghmap <- cortab %>%
      dplyr::mutate(x = factor(x, unique(.$x)), y = factor(y, unique(.$x))) %>%
      ggplot() +
      geom_tile(aes(x = x, y = y, fill = correlation)) +
      scale_fill_gradient2(low = "blue3", mid = "white", high = "brown4") +
      labs(x = NULL, y = NULL) +
      theme_classic() +
      theme(axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1))

    ggplot2::ggsave(fig_path, plot = ghmap, width = width, height = height)
  } else {
    cat("[W]: Found", fig_path, "Skipping it ...\n")
  }
}
