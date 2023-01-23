#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, data.table.verbose = FALSE)
suppressPackageStartupMessages({
  # Functional annotation, enrichment, background databases, etc.
  library(msigdbr)
  library(ReactomePA) # reactome.db_1.81.0
  library(clusterProfiler)
  library(DOSE)
  library(enrichplot)
  library(ComplexHeatmap)

  # Misc, color, data manipulation, etc.
  library(circlize)
  library(magrittr)
  library(tidyverse)
  library(data.table)
})


#
## Variables
#
proj_dir <- "~/Documents/projects/wp_bcg_eqtl"

mode_vec <- c("normal", "interaction")
celltype_vec <- c("Monocytes", "CD4T", "CD8T", "NK", "B")
condition_vec <- c("T0_LPS.vs.T0_RPMI", "T3m_LPS.vs.T3m_RPMI", "T3m_RPMI.vs.T0_RPMI") #, "T3m_LPS.vs.T0_RPMI")
col_pattern_vec <- c("normal" = "_Common", "interaction" = ".vs.")

token <- "filtered" # or "unfiltered"
all_eqtl_file <- file.path(proj_dir, "outputs/pseudo_bulk/overview", paste0(token, ".all_top_qtl.FDR0.1.csv"))
eqtl_tab <- data.table::fread(all_eqtl_file)


#
## Parameters
#
fdr_max <- 0.1
min_isize <- 10
override <- TRUE


#
## Enrichment analysis
#
for (per_celltype in celltype_vec) {
  pval_col <- paste0("global_corrected_pValue_", per_celltype, "_Common")
  beta_col <- paste0("beta_", per_celltype, "_Common")
  gene_info <- eqtl_tab %>%
    dplyr::filter(!!rlang::sym(pval_col) < 0.05) %>%
    dplyr::select(dplyr::one_of(c("feature_id", beta_col))) %>%
    tibble::deframe() %>%
    (function(vec) {
      clusterProfiler::bitr(names(vec), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") %>% dplyr::mutate(Beta = vec[SYMBOL])
    }) %>%
    dplyr::filter(!stringr::str_starts(SYMBOL, "HLA-"))

  # GSEA
  cat("Gene set enrichment analysis for", per_celltype, "...\n")
  gsea_res_all <- NULL
  for (per_gscat in c("H", paste0("C", 1:8))) {
    save_to <- paste(per_celltype, per_gscat, "gseaplot.pdf", sep = "-")
    if (!file.exists(save_to) || override) {
      gsea_res <- gene_info %>%
        dplyr::select(ENTREZID, Beta) %>%
        dplyr::mutate(ENTREZID = as.integer(ENTREZID), Beta = abs(Beta)) %>%
        dplyr::arrange(dplyr::desc(Beta)) %>%
        tibble::deframe() %>%
        (function(vec) {
          msigdbr(species = "Homo sapiens") %>%
            dplyr::filter(gs_cat == per_gscat) %>%
            dplyr::select(gs_name, entrez_gene) %>%
            clusterProfiler::GSEA(vec, TERM2GENE = ., pvalueCutoff = 1.0, minGSSize = 10, maxGSSize = 8000)
        })

      gsea_res_all <- gsea_res@result %>%
        dplyr::mutate(
          Celltype = per_celltype,
          Category = per_gscat,
          Genesymbol = lapply(core_enrichment,
            function(a) {
              stringr::str_split(a, pattern = "/", simplify = TRUE) %>%
              as.integer() %>%
              clusterProfiler::bitr(fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Hs.eg.db") %>%
              dplyr::pull(SYMBOL) %>%
              paste(collapse = "/")
        })) %>%
        rbind(gsea_res_all, .)

      if (nrow(gsea_res@result) < 1) next

      p <- gseaplot2(gsea_res, geneSetID = 1)
      ggsave(save_to, width = 7, height = 4.5)
    }
  }
  gsea_res_all %>% data.table::fwrite(paste0(per_celltype, "_gsea.csv"))

  # Over-representation analysis
  cat("Over representation analysis for", per_celltype, "...\n")
  list(
    enrichGO(gene_info$ENTREZID, pvalueCutoff = 1, readable = TRUE, OrgDb = "org.Hs.eg.db") %>%
      (function(res) res@result %>% dplyr::mutate(Celltype = per_celltype, Category = "GO")),
    enrichPathway(gene_info$ENTREZID, pvalueCutoff = 1, readable = TRUE) %>%
      (function(res) res@result %>% dplyr::mutate(Celltype = per_celltype, Category = "ReactomePathway")),
    enrichKEGG(gene_info$ENTREZID, pvalueCutoff = 1, organism = "hsa") %>%
      (function(res) res@result %>% dplyr::mutate(Celltype = per_celltype, Category = "KEGG")),
    enrichWP(gene_info$ENTREZID, pvalueCutoff = 1, organism = "Homo sapiens") %>%
      (function(res) res@result %>% dplyr::mutate(Celltype = per_celltype, Category = "WikiPathway")),
    enrichDO(gene_info$ENTREZID, pvalueCutoff = 1, readable = TRUE) %>%
      (function(res) res@result %>% dplyr::mutate(Celltype = per_celltype, Category = "Disease"))
  ) %>%
    Reduce(rbind, .) %>%
    data.table::fwrite(paste0(per_celltype, "_enrichment.csv"))
}


# Plot the enrichment analysis
enr_tab <- paste0(celltype_vec, "_enrichment.csv") %>%
  lapply(function(p) data.table::fread(p)) %>%
  Reduce(rbind, .) %>%
  dplyr::mutate(
    Log2OddRatio = purrr::map2(GeneRatio, BgRatio, ~ log2(eval(parse(text = .x)) / eval(parse(text = .y)))),
    Log2OddRatio = unlist(Log2OddRatio)
  )

col_fun_prop <- colorRamp2(c(min(abs(enr_tab$Log2OddRatio)), max(abs(enr_tab$Log2OddRatio))), c("gray95", "#2166AC"))
for (pvalmax in c(0.01, 0.05, 0.1)) {
  hm_tab <- enr_tab %>%
    dplyr::filter(p.adjust < pvalmax, Log2OddRatio >= 1 | Log2OddRatio <= -1) %>%
    tidyr::pivot_wider(id_cols = c(ID, Category, Description), names_from = Celltype, values_from = Log2OddRatio) %>%
    dplyr::arrange(Category, ID) %>%
    dplyr::mutate(dplyr::across(-c(ID, Category, Description), ~ dplyr::if_else(is.na(.x), 0, .x)))

  hm_plot_list <- list()
  # for (pcat in c("KEGG", "WikiPathway", "ReactomePathway", "Disease", "GO")) {
  for (pcat in c("KEGG", "WikiPathway", "ReactomePathway", "Disease", "GO")) {
    # for (pcat in c("WikiPathway", "Disease", "GO")) {
    hm_plot_list[[pcat]] <- hm_tab %>%
      dplyr::filter(Category == pcat) %>%
      (function(tab) {
        mat <- dplyr::select(tab, -c(ID, Description, Category)) %>%
          as.matrix() %>%
          t() # %>% (function(e) -log10(e))
        # colnames(mat) <- tab$Description
        mat
      }) %>%
      Heatmap(name = pcat, column_title = pcat, col = col_fun_prop, show_heatmap_legend = FALSE)
  }

  pdf(paste0("enrichment_heatmap_padj", pvalmax, ".pdf"), width = 10, height = 2.5)
  hm_plot <- Reduce(`+`, hm_plot_list)
  draw(hm_plot)

  lgd <- Legend(col_fun = col_fun_prop, title = "Log2(OddRatio)", direction = "horizontal")
  draw(lgd, x = unit(0.89, "npc"), y = unit(0.875, "npc"), just = c("left"))
  dev.off()
}
