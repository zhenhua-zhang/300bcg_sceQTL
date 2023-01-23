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
  library(ggrepel)
  library(ggpubr)
  library(ggsci)
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

#
## Plot number of eQTLs
#
col_dict <- c("B" = "#E64B35", "CD4T" = "#4DBBD5", "CD8T" = "#00A087", "NK" = "#3C5488", "Monocytes" = "#F39B7F")
plot_size <- list("normal" = c("width" = 7, "height" = 3), "interaction" = c("width" = 15, "height" = 4))
for (mode in mode_vec) {
  run_pattern <- col_pattern_vec[mode] %>% as.vector()
  genesets <- eqtl_tab %>%
    dplyr::select(c(feature_id, snp_id, dplyr::contains(run_pattern) & dplyr::contains("global_corrected_pValue"))) %>%
    tidyr::pivot_longer(cols = dplyr::contains("global_"), names_to = "Celltype", values_to = "global_corrected_pValue") %>%
    dplyr::filter(global_corrected_pValue < fdr_max) %>%
    dplyr::mutate(Celltype = stringr::str_remove_all(Celltype, "global_corrected_pValue_|_Common")) %>%
    dplyr::group_by(Celltype) %>%
    dplyr::summarize(GeneSets = {
      geneset <- list(unique(feature_id))
      names(geneset) <- cur_group()[[1]]
      geneset
    })

  mat <- make_comb_mat(genesets$GeneSets)
  mat <- mat[comb_size(mat) >= min_isize]

  row_annotations <- stringr::str_split(ComplexHeatmap::set_name(mat), pattern = "_|\\.", n = 2, simplify = TRUE)

  celltypes <- row_annotations[, 1]
  celltypes_col <- col_dict[celltypes]

  row_labels <- row_annotations[, ifelse(mode == "normal", 1, 2)]
  comb_col_vec <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4")[comb_degree(mat)]

  pw <- plot_size[[mode]]["width"]
  ph <- plot_size[[mode]]["height"]

  pdf(paste0(token, ".", mode, ".eGene_number_upset_FDR", fdr_max, ".pdf"), width = pw, height = ph)
  us <- UpSet(mat,
    comb_col = comb_col_vec,
    comb_order = order(comb_size(mat)),
    top_annotation = upset_top_annotation(mat, add_numbers = TRUE, annotation_name_rot = 90, numbers_rot = 90),
    right_annotation = upset_right_annotation(mat, gp = gpar(fill = "gray", col = "gray25"), add_numbers = TRUE),
    left_annotation = rowAnnotation(Celltype = celltypes, col = list(Celltype = celltypes_col)),
    row_labels = row_labels,
    row_names_gp = gpar(fontsize = 9),
    row_names_max_width = max_text_width(rownames(mat), gp = gpar(fontsize = 9))
  )
  draw(us)
  dev.off()
}



#
## Plot shared signals across different summary statistics
#
cat("check mashr.r for more details\n")



#
## Enrichment analysis
#
for (per_celltype in celltype_vec) {
  pval_col <- paste0("global_corrected_pValue_", per_celltype, "_Common")
  beta_col <- paste0("beta_", per_celltype, "_Common")
  gene_info <- eqtl_tab %>%
    dplyr::filter(!!rlang::sym(pval_col) < 0.1) %>%
    dplyr::select(dplyr::one_of(c("feature_id", beta_col))) %>%
    tibble::deframe() %>%
    abs() %>%
    sort(decreasing = TRUE) %>%
    (function(vec) {
      clusterProfiler::bitr(names(vec), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") %>% dplyr::mutate(Beta = vec[SYMBOL])
    }) %>%
    dplyr::filter(!stringr::str_starts(SYMBOL, "HLA-"))

  # GSEA
  cat("Gene set enrichment analysis for", per_celltype, "...\n")
  gsea_res_all <- NULL
  for (per_gscat in c("H", paste0("C", 1:8))) {
    save_to <- paste(per_celltype, per_gscat, "gseaplot.pdf", sep = "-")
    gsea_res <- gene_info %>%
      dplyr::select(ENTREZID, Beta) %>%
      dplyr::mutate(ENTREZID = as.integer(ENTREZID)) %>%
      tibble::deframe() %>%
      (function(vec) {
        msigdbr(species = "Homo sapiens") %>%
          dplyr::filter(gs_cat == per_gscat) %>%
          dplyr::select(gs_name, entrez_gene) %>%
          clusterProfiler::GSEA(vec, TERM2GENE = ., pvalueCutoff = 1.0, minGSSize = 10, maxGSSize = 8000)
      })
    gsea_res_all <- gsea_res@result %>% dplyr::mutate(Celltype = per_celltype, Category = per_gscat) %>% rbind(gsea_res_all, .)

    if (nrow(gsea_res@result) < 1) next
    p <- gseaplot2(gsea_res, geneSetID = 1)
    ggsave(save_to, width = 7, height = 4.5)
  }
  gsea_res_all %>% data.table::fwrite(paste0(per_celltype, "_gsea.csv"))

  # Over-representation analysis
  cat("Over representation analysis for", per_celltype, "...\n")
  list(
    enrichGO(gene_info$ENTREZID, ont = "BP", pvalueCutoff = 1, readable = TRUE, OrgDb = "org.Hs.eg.db") %>%
      (function(res) {
        res@result %>% dplyr::mutate(Celltype = per_celltype, Category = "GO")
      }),
    enrichPathway(gene_info$ENTREZID, pvalueCutoff = 1, readable = TRUE) %>%
      (function(res) {
        res@result %>% dplyr::mutate(Celltype = per_celltype, Category = "ReactomePathway")
      }),
    enrichKEGG(gene_info$ENTREZID, pvalueCutoff = 1, organism = "hsa") %>%
      (function(res) {
        res@result %>% dplyr::mutate(Celltype = per_celltype, Category = "KEGG")
      }),
    enrichWP(gene_info$ENTREZID, pvalueCutoff = 1, organism = "Homo sapiens") %>%
      (function(res) {
        res@result %>% dplyr::mutate(Celltype = per_celltype, Category = "WikiPathway")
      }),
    enrichDO(gene_info$ENTREZID, pvalueCutoff = 1, readable = TRUE) %>%
      (function(res) {
        res@result %>% dplyr::mutate(Celltype = per_celltype, Category = "Disease")
      })
  ) %>%
    Reduce(rbind, .) %>%
    data.table::fwrite(paste0(per_celltype, "_enrichment.csv"))
}


#
## Functional enrichment analysis
#
enr_tab <- paste0(celltype_vec, "_enrichment.csv") %>%
  lapply(function(p) data.table::fread(p)) %>%
  Reduce(rbind, .) %>%
  dplyr::mutate(
    Log2OddRatio = purrr::map2(GeneRatio, BgRatio, ~ log2(eval(parse(text = .x)) / eval(parse(text = .y)))),
    Log2OddRatio = unlist(Log2OddRatio)
  )

col_fun_prop <- colorRamp2(c(min(abs(enr_tab$Log2OddRatio)), max(abs(enr_tab$Log2OddRatio))), c("gray95", "#00AFBB"))
for (pvalmax in c(1e-3, 5e-3, 1e-2, 5e-2, 1e-1)) {
  hm_tab <- enr_tab %>%
    dplyr::filter(pvalue < pvalmax) %>%
    dplyr::mutate(p.adjust = -log10(p.adjust)) %>%
    tidyr::pivot_wider(id_cols = c(ID, Category, Description), names_from = Celltype, values_from = p.adjust) %>%
    dplyr::arrange(Category, ID) %>%
    dplyr::mutate(dplyr::across(-c(ID, Category, Description), ~ dplyr::if_else(is.na(.x), 0, .x)))

  hm_plot_list <- list()
  for (pcat in c("KEGG", "WikiPathway", "ReactomePathway", "GO")) {
    # for (pcat in c("WikiPathway", "Disease", "GO")) {
    hm_plot_list[[pcat]] <- hm_tab %>%
      dplyr::filter(Category == pcat) %>%
      (function(tab) {
        mat <- dplyr::select(tab, -c(ID, Description, Category)) %>% as.matrix() %>% t()
        # rownames(mat) <- tab$Description

        mat
      }) %>%
      Heatmap(
        name = pcat, column_title = pcat, col = col_fun_prop, show_heatmap_legend = FALSE,
        row_names_max_width = max_text_width(rownames(.), gp = gpar(fontsize = 7))
      )
  }

  pdf(paste0("enrichment_heatmap_pval", pvalmax, ".pdf"), width = 10, height = 2.5)
  hm_plot <- Reduce(`+`, hm_plot_list)
  draw(hm_plot)

  lgd <- Legend(col_fun = col_fun_prop, title = "-Log10(p.adjust)", direction = "horizontal")
  draw(lgd, x = unit(0.875, "npc"), y = unit(0.875, "npc"), just = c("left"))
  dev.off()
}


#
## eQTL rank plot
#
celltype <- "B"
highlight <- ""
tar_pval_col <- paste0("p_value_", celltype, "_Common")
tar_beta_col <- paste0("beta_", celltype, "_Common")

qtlrank_tab <- eqtl_tab %>%
  dplyr::filter(!!rlang::sym(tar_pval_col) < 0.05, !stringr::str_detect(QTL, "HLA-")) %>%
  dplyr::group_by(feature_id) %>%
  dplyr::slice_min(!!rlang::sym(tar_pval_col), n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::select(dplyr::all_of(c("p_value" = tar_pval_col, "beta" = tar_beta_col, "feature_id"))) %>%
  dplyr::arrange(p_value) %>%
  dplyr::mutate(rank = dplyr::row_number(), p_value = -log10(p_value))

myplot <- ggplot() +
  geom_point(aes(x = rank, y = p_value), qtlrank_tab, shape = 21, color = "gray", fill = "gray", size = 0.75) +
  geom_point(aes(x = rank, y = p_value), head(qtlrank_tab, 10), shape = 21, color = "black", fill = "purple", size = 3, alpha = 0.5) +
  geom_text_repel(
    aes(rank, p_value, label = feature_id), dplyr::filter(qtlrank_tab, feature_id != highlight) %>% head(9),
    min.segment.length = 0, max.overlaps = 20, box.padding = 0.5
  ) +
  geom_label_repel(
    aes(rank, p_value, label = feature_id), dplyr::filter(qtlrank_tab, feature_id == highlight), color = "red", box.padding = 7
  ) +
  labs(x = "Rank by p-value", y = "-Log10(p-value)") +
  theme_classic()
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/overview", paste0("eqtl_rank.", celltype, ".pdf"))
ggsave(save_to, plot = myplot, width = 6, height = 4)


#
## Choose an example eQTL, excepting ADCY3-rs11687089
#
tar_qtl <- eqtl_tab %>%
  dplyr::select(c(feature_id, snp_id, dplyr::contains("_Common") & dplyr::contains("global_corrected_pValue"))) %>%
  tidyr::pivot_longer(cols = dplyr::contains("global_"), names_to = "Celltype", values_to = "global_corrected_pValue") %>%
  dplyr::filter(global_corrected_pValue < fdr_max & !is.na(global_corrected_pValue)) %>%
  dplyr::mutate(Celltype = stringr::str_remove_all(Celltype, "global_corrected_pValue_|_Common")) %>%
  dplyr::group_by(feature_id, snp_id) %>%
  dplyr::summarize(
    QTL = paste(snp_id, feature_id, sep = "-"),
    Celltype_specific = n() == 1,
    TCell_specific = sum(cur_data()$Celltype %in% c("CD4T", "CD8T")) == 2
  )

# ADCY3-rs11687089
# DEAF1-rs10794331
# GFM1-rs1367313
# PPFIBP2-rs4758197
# HAGH-rs344352
t_cell_specific <- tar_qtl %>% dplyr::filter(TCell_specific) %$% QTL
eqtl_tab %>%
  dplyr::select(feature_id, snp_id, QTL, dplyr::contains("_Common")) %>%
  dplyr::filter(QTL %in% t_cell_specific) %>%
  data.table::fwrite("TCell_specific.eQTL.csv")

eqtl_tab %>%
  dplyr::filter(
    is.na(global_corrected_pValue_Monocytes_Common),
    !is.na(global_corrected_pValue_Monocytes_T0_LPS.vs.T0_RPMI),
    dplyr::if_any(dplyr::contains("_Common") & dplyr::contains("_pValue"), ~ !is.na(.x))
  ) %>%
dplyr::select(-dplyr::contains("beta"), -dplyr::contains("global")) %>%
as.data.frame()


#
## Plot number of coloc loci
#
min_h4 <- 0.5

# Number of colocalizations.
mode <- "normal"
coloc_files <- list.files(file.path(proj_dir, "outputs/pseudo_bulk/coloc", mode),
  pattern = "colocalization_test.csv", recursive = TRUE, full.names = TRUE) %>%
  purrr::discard(~ stringr::str_detect(.x, "archive"))

mat <- coloc_files %>%
  lapply(function(fp) {
    if (mode == "normal") {
      celltype <- fp %>% dirname() %>% basename()
      condition <- ""
    } else {
      celltype <- fp %>% dirname() %>% dirname() %>% basename()
      condition <- fp %>% dirname() %>% basename()
    }

    data.table::fread(fp) %>%
      dplyr::mutate(Celltype = celltype, Condition = condition) %>%
      dplyr::filter(H4 >= min_h4) %>%
      dplyr::group_by(outcome, Celltype, Condition) %>%
      dplyr::summarise(n_genes = n(), gene_list = paste(exposure, collapse = "|"))
  }) %>%
  Reduce(rbind, .) %>%
  tidyr::pivot_wider(id_cols = outcome, names_from = c(Celltype, Condition), values_from = n_genes) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(dplyr::across(-outcome, ~ dplyr::if_else(is.na(.x), as.integer(0), .x))) %>%
  (function(m) {
    tmp <- m %>% dplyr::select(-outcome) %>% as.matrix()
    rownames(tmp) <- m$outcome
    tmp
  })


col_fun_prop <- colorRamp2(c(0, max(mat)), c("gray95", "#2166AC"))
pdf(paste0(mode, ".colocalization.per_celltype.pdf"), width = 8, height = 6)
hm <- Heatmap(t(mat),
  col = col_fun_prop, name = "Nr. of \ncolocs",
  top_annotation = columnAnnotation(Counts = anno_barplot(rowSums(mat)), show_annotation_name = FALSE),
  right_annotation = rowAnnotation(Counts = anno_barplot(colSums(mat)), show_annotation_name = FALSE),
  column_names_rot = 45, column_names_side = "top", column_dend_side = "bottom"
)
draw(hm)
dev.off()


#
## Plot interaction eQTL
#
tribble(
  ~tar_snp, ~tar_feature, ~tar_baseline, ~tar_treatment,
  "rs11731570", "EXOSC9", "T0_RPMI", "T0_LPS",
  "rs11731570", "EXOSC9", "T0_RPMI", "T3m_RPMI",
  "rs11731570", "EXOSC9", "T3m_RPMI", "T3m_LPS",
  "rs11687089", "ADCY3", "T0_RPMI", "T0_LPS",
  "rs11687089", "ADCY3", "T0_RPMI", "T3m_RPMI",
  "rs11687089", "ADCY3", "T3m_RPMI", "T3m_LPS",
  "rs67260737", "NUMA1", "T0_RPMI", "T0_LPS",
  "rs67260737", "NUMA1", "T0_RPMI", "T3m_RPMI",
  "rs67260737", "NUMA1", "T3m_RPMI", "T0_RPMI",
) %>%
apply(1, function(vec) {
  per_snp <- vec[1]
  per_feature <- vec[2]
  per_baseline <- vec[3]
  per_treatment <- vec[4]

  qtl_info_tab <- data.table::fread(paste0(per_feature, "-", per_snp, ".QTL_info.csv")) %>%
    dplyr::mutate(time = dplyr::if_else(time == 0, "T0", "T3m"),
                  stim = dplyr::if_else(stim == 0, "RPMI", "LPS"),
                  patient = stringr::str_extract(V1, "300BCG[0-9]+"),
                  celltype = factor(celltype, celltype_vec)
    ) %>%
    dplyr::select(-V1) %>%
    tidyr::pivot_wider(names_from = c(time, stim), values_from = !!rlang::sym(per_feature))

  g <- ggplot(data = qtl_info_tab, mapping = aes_string(x = per_baseline, y = per_treatment, color = per_snp)) +
    geom_point() +
    geom_abline(linetype = "dashed") +
    stat_smooth(method = "lm") +
    facet_wrap(~celltype, scales = "free", nrow = 1) +
    theme_classic()
  save_to <- paste0("interaction.", per_feature, "-", per_snp, ".", per_treatment, ".", per_baseline, ".pdf")
  ggsave(save_to, width = 12, height = 2.5)
})


#
## Plot boxplot, common effect
#
per_feature <- "ADCY3"
per_snp <- "rs11687089"

per_feature <- "DEAF1"
per_snp <- "rs10794331"

per_feature <- "EXOSC9"
per_snp <- "rs11731570"
per_celltype <- "CD4T"

per_feature <- "NUMA1"
per_snp <- "rs67260737"
per_celltype <- "Monocytes"

per_feature <- "SUCO"
per_snp <- "rs1011731"
per_celltype <- "CD8T"

per_feature <- "COMMD1"
per_snp <- "rs7562347"
per_celltype <- "NK"

qtl_info_tab <- data.table::fread(paste0(per_feature, "-", per_snp, ".QTL_info.csv")) %>%
  dplyr::mutate(
    time = dplyr::if_else(time == 0, "T0", "T3m"), stim = dplyr::if_else(stim == 0, "RPMI", "LPS"),
    patient = stringr::str_extract(V1, "300BCG[0-9]+"), celltype = factor(celltype, celltype_vec)
  ) %>%
  dplyr::arrange(dplyr::desc(!!rlang::sym(per_snp)))


y_lab <- paste0("Expression: ", per_feature)
# Common effects
plot <- qtl_info_tab %>%
  ggboxplot(x = per_snp, y = per_feature, fill = per_snp, palette = "npg", xlab = FALSE, ylab = y_lab) +
  geom_pwc(label = "p.adj.signif", hide.ns = TRUE, vjust = 0.5, tip.length = 0) +
  theme_classic2()
plot <- facet(plot, nrow = 1, facet.by = "celltype")
ggsave(paste0("normal.", per_feature, "-", per_snp, ".boxplot.pdf"), plot = plot, height = 3, width = 8)

# LPS effect per time point
time_point <- c("T3m", "T0")
for (per_tp in time_point) {
  plot <- qtl_info_tab %>%
    dplyr::select(time, stim, celltype, patient, one_of(per_snp, per_feature)) %>%
    dplyr::filter(time == per_tp, celltype == per_celltype) %>%
    # tidyr::pivot_wider(names_from = stim, values_from = !!rlang::sym(per_feature)) %>%
    ggboxplot(x = per_snp, y = per_feature, fill = per_snp, palette = "npg", xlab = FALSE, ylab = y_lab) +
    # ggpaired(cond1 = "RPMI", cond2 = "LPS", fill = per_snp, palette = "npg", xlab = FALSE, ylab = y_lab) +
    geom_pwc(label = "p.adj.signif", hide.ns = TRUE, vjust = 0.5, tip.length = 0) +
    theme_classic2()
  plot <- facet(plot, nrow = 1, facet.by = "stim")
  save_to <- paste0("interaction.", per_feature, "-", per_snp, ".", per_tp, ".", per_celltype, ".per_snp.pdf")
  ggsave(save_to, plot = plot, width = 4, height = 3)
}

# BCG effects only RPMI
stimulation <- "RPMI"
plot <- qtl_info_tab %>%
  dplyr::select(time, stim, celltype, patient, one_of(per_snp, per_feature)) %>%
  dplyr::filter(stim == stimulation, celltype == per_celltype) %>%
  # tidyr::pivot_wider(names_from = time, values_from = !!rlang::sym(per_feature)) %>%
  ggboxplot(x = per_snp, y = per_feature, fill = per_snp, palette = "npg", xlab = FALSE, ylab = y_lab) +
  # ggpaired(cond1 = "T0", cond2 = "T3m", fill = per_snp, palette = "npg", xlab = FALSE, ylab = y_lab) +
  geom_pwc(label = "p.adj.signif", hide.ns = TRUE, vjust = 0.5, tip.length = 0) +
  theme_classic2()
plot <- facet(plot, nrow = 1, facet.by = "time")
save_to <- paste0("interaction.", per_feature, "-", per_snp, ".", stimulation, ".", per_celltype, ".per_snp.pdf")
ggsave(save_to, plot = plot, width = 4, height = 3)
