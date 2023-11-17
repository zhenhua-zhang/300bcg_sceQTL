#!/usr/bin/env Rscript
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com, zhenhua.zhang3@helmholtz-hzi.de
# Created: 2022 Oct 12
# Updated: Sep 05, 2023

options(stringsAsFactors = FALSE, data.table.verbose = FALSE, datatable.showProgress = FALSE, digits = 2)
suppressPackageStartupMessages({
  # Single-cell ATAC-seq
  # library(ArchR)

  # Single-cell RNA-seq
  library(SeuratObject)
  library(Seurat)
  library(Matrix)

  # Trackplot
  library(Gviz)
  library(GenomicRanges)

  # Correlation analysis
  library(lme4)
  library(lmerTest)

  # Co-expression analysis
  library(CSCORE)
  library(WGCNA)

  # Functional enrichment
  library(DOSE)
  library(msigdbr)
  library(ReactomePA) # reactome.db_1.81.0
  library(clusterProfiler)

  # Data manipulation
  library(pzfx)
  library(readxl)
  library(magrittr)
  library(tidyverse)
  library(data.table)

  # Plotting
  # library(circlize)
  # library(ggalluvial)
  # library(ggfortify)
  # library(ggsignif)
  # library(ggdist)
  library(enrichplot)
  library(colorRamp2)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(patchwork)
  library(ggbeeswarm)
  library(ggbreak)
  library(ggforce)
  library(ggrepel)
  library(ggpubr)
  library(ggh4x)
  library(ggsci)
  library(webr) # To plot PieDonut chart
})


#
## Variables
#
proj_dir <- "~/Documents/projects/wp_bcg_eqtl"
time_vec <- c("T0", "T3m")
stim_vec <- c("RPMI", "LPS")
time_stim_vec <- c("T0_RPMI", "T0_LPS", "T3m_RPMI", "T3m_LPS")
mode_vec <- c("normal", "interaction")
celltype_vec <- c("Monocytes", "CD4T", "CD8T", "NK", "B")
condition_vec <- c("T0_LPS.vs.T0_RPMI", "T3m_LPS.vs.T3m_RPMI", "T3m_RPMI.vs.T0_RPMI") #, "T3m_LPS.vs.T0_RPMI")
col_pattern_vec <- c("normal" = "Common", "interaction" = ".vs.")
effect_vec <- c("LPS eff. (T0)", "LPS eff. (T3m)", "BCG eff.")
names(effect_vec) <- condition_vec


#
## Parameters
#
fdr_max <- 0.05
min_h4 <- 0.3

# Define the number of colors
nb_cols <- 10
mycolors <- colorRampPalette(brewer.pal(8, "Blues"))(nb_cols)

#
## Load single-cell RNA-seq results
#
pbmc_300bcg <- readRDS(file.path(proj_dir, "inputs/sc_rnaseq/bcg4-0712.rds"))
pbmc_300bcg <- pbmc_300bcg[, (!pbmc_300bcg@meta.data$clusters1 %in% c("HSP(T)"))] # Removing HSP(T) cells
new_cell_lvl <- c("Monocytes", "CD4+ T", "CD8+ T", "NK", "B", "Platelet", "mDC", "pDC", "Undefined")
pbmc_300bcg@meta.data$clusters1 <- factor(pbmc_300bcg@meta.data$clusters1, levels = new_cell_lvl)
pbmc_300bcg@meta.data$ts <- factor(pbmc_300bcg@meta.data$ts, levels = time_stim_vec)
Idents(pbmc_300bcg) <- "clusters1"
DefaultAssay(pbmc_300bcg) <- "RNA"

#
## Plot single-cell UMAP
#
g_umap <- DimPlot(pbmc_300bcg, reduction = "umap") + NoAxes(keep.text = TRUE)
ggsave(file.path(proj_dir, "outputs/pseudo_bulk/overview/sc_rnaseq_umap.pdf"), g_umap, width = 7, height = 6)
g_umap_pcond <- DimPlot(pbmc_300bcg, reduction = "umap", split.by = "ts", ncol = 2) + NoAxes(keep.text = TRUE)
ggsave(file.path(proj_dir, "outputs/pseudo_bulk/overview/sc_rnaseq_umap.per_condition.pdf"), g_umap_pcond, width = 8, height = 8)


#
## Draw feature map of cell type marker genes
#
DefaultAssay(pbmc_300bcg) <- "RNA"
ct_markers <- c(
  "CTSS", "FCN1", "NEAT1", "LYZ", "PSAP", "S100A9", "AIF1", # "MNDA", "TYROBP", # Monocytes
  "IL7R", "MAL", "LTB", "LDHB", "TPT1", "TRAC", "TMSB10", # "CD3D", "CD4", "CD3G", # CD4 T
  "CD8B", "CD8A", "CD3D", "TMSB10", "HCST", "CD3G", "LINC02446", "CTSW", # "CD3E", "TRAC", # CD8 T
  "NKG7", "KLRD1", "TYROBP", "GNLY", "FCER1G", "PRF1", "CD247", "KLRF1", # "CST7", "GZMB", # NK
  "CD79A", "RALGPS2", "CD79B", "MS4A1", "BANK1", "CD74", "TNFRSF13C", "HLA-DQA1", # "IGHM", "MEF2C", # B
  "CD74", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "CCDC88A", "HLA-DRA", "HLA-DMA", "CST3" #, "HLA-DQB1", "HLA-DRB1" # DC
)
g_ct_marker <- FeaturePlot(pbmc_300bcg, features = ct_markers, raster = TRUE) & NoLegend() & FontSize(main = 12) & NoAxes()
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/overview/sc.celltype_markers.featureplot.pdf")
ggsave(save_to, g_ct_marker)

DefaultAssay(pbmc_300bcg) <- "integrated"
tar_celltype <- c("Monocytes", "CD4+ T", "CD8+ T", "NK", "B", "mDC", "pDC")
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/overview/sc.celltype_markers.dotplot.horizontal.pdf")
g_ct_marker <- DotPlot(pbmc_300bcg, assay = "RNA", features = unique(ct_markers), cols = "RdBu", idents = tar_celltype) + RotatedAxis() + coord_flip()
ggsave(save_to, g_ct_marker, width = 12, height = 4)

save_to <- file.path(proj_dir, "outputs/pseudo_bulk/overview/sc.celltype_markers.dotplot.vertical.pdf")
g_ct_marker <- DotPlot(pbmc_300bcg, assay = "RNA", features = unique(ct_markers), cols = "RdBu", idents = tar_celltype) + RotatedAxis() + coord_flip()
ggsave(save_to, g_ct_marker, width = 6, height = 10)


#
## Draw the cell-proportion boxplot
#
cpp_tab <- pbmc_300bcg@meta.data %>%
  dplyr::filter(status == "singlet") %>%
  dplyr::group_by(time, stim, ids) %>%
  dplyr::summarise(cpp = { data.frame(table(dplyr::cur_data()$clusters1) / dplyr::n()) }) %>%
  as.data.table() %>%
  dplyr::rename("cell_type" = "cpp..", "cell_proportion" = "cpp.Freq") %>%
  dplyr::filter(!cell_type %in% c("Undefined", "Platelet")) %>%
  dplyr::mutate(stim = factor(stim, levels = c("RPMI", "LPS")), cell_type = forcats::fct_reorder(cell_type, cpp))

g_cpp <- ggboxplot(
  cpp_tab, x = "time", y = "cell_proportion", fill = "stim", color = "black", size = 0.5,
  facet.by = "cell_type", add = "ggscatter", palette = "npg", width = 0.6, xlab = "Timepoint",
  ylab = "Cell proportion"
) %>%
facet(facet.by = "cell_type", nrow = 1) +
scale_y_break(c(0.63, 0.73), scale = 0.05) +
scale_y_break(c(0.49, 0.62), scale = 0.05, expand = expansion(add = 0.05)) +
scale_y_continuous(breaks = seq(0, 10) / 10.0) +
theme_classic() +
theme(axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank(), legend.position = "top")

save_to <- file.path(proj_dir, "outputs/pseudo_bulk/overview/cell_prop.pdf")
pdf(save_to, width = 7.5, height = 3.5, onefile = FALSE)
print(g_cpp)
dev.off()


#
## Overiview of the cohort
#
meta_info <- pbmc_300bcg@meta.data %>% dplyr::select(age, gender, ids) %>% dplyr::distinct() %>% as.data.table()

# Gender pie plot
plot_tab <- dplyr::select(meta_info, ids, gender) %>% dplyr::group_by(gender) %>% dplyr::summarise(n = n()) %>% dplyr::ungroup() %>% dplyr::arrange(n) %>% dplyr::mutate(pos = cumsum(n) - 0.5 * n)
p_gender <- ggplot(data = plot_tab) +
  geom_bar(aes(x = 2, y = n, fill = gender), stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(x = 2, y = pos, label = n), color = "white") +
  scale_fill_npg() +
  labs(fill = "Gender") +
  xlim(0.5, 2.5) +
  theme_void() +
  theme(legend.position = "none")

# Age distribution
plot_tab <- dplyr::select(meta_info, ids, age, gender) %>% dplyr::mutate(gender = dplyr::if_else(gender == "f", "Female", "Male"))
p_age <- ggplot(data = plot_tab) +
  geom_violin(aes(x = gender, y = age, fill = gender), draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_point(aes(x = gender, y = age), position = position_jitter(width = 0.2)) +
  scale_fill_npg() +
  labs(x = NULL, y = "Age", fill = "Gender") +
  theme_classic()

p_demo <- p_gender + p_age
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/overview/cohort_demographical_plot.pdf")
ggsave(save_to, p_demo, width = 6, height = 3)


#
## Number of cells per condition
#
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/overview/cell_count_plot.pdf")
tar_clusters <- c("Monocytes", "CD4+ T", "CD8+ T", "NK", "B")
pdf(save_to, width = 6, height = 6)
pbmc_300bcg@meta.data %>%
  dplyr::filter(clusters1 != "HSP(T)") %>%
  dplyr::mutate(clusters1 = as.character(clusters1), clusters1 = dplyr::if_else(clusters1 %in% tar_clusters, clusters1, "Others")) %>%
  dplyr::group_by(ts, clusters1) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(clusters1 = forcats::fct_reorder(clusters1, n)) %>%
  PieDonut(aes(pies = ts, donuts = clusters1, count = n), r0 = 0.2, r1 = 0.7, ratioByGroup = TRUE, showPieName = FALSE, showRatioThreshold = 0.02)
dev.off()


#
## Plot number of expressed genes
#
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/overview/average_expression.per_condition.all_genes.csv")
if (!file.exists(save_to)) {
  exp_cts_tab <- AverageExpression(pbmc_300bcg, assays = "RNA", group.by = c("clusters1", "ids", "ts")) %>%
    as.data.frame() %>%
    dplyr::mutate(gene_symbol = rownames(.)) %>%
    tidyr::pivot_longer(-gene_symbol) %>%
    tidyr::separate(name, into = c("celltype", "ids", "time", "stim"), sep = "_") %>%
    tidyr::pivot_wider(names_from = ids, values_from = value) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(n_non_zero = sum(c_across(!c(gene_symbol, celltype, time, stim)) > 0, na.rm = TRUE), p_non_zero = n_non_zero / as.numeric(38.0))
  fwrite(exp_cts_tab, save_to)
} else {
  exp_cts_tab <- fread(save_to)
}

plot_tab <- exp_cts_tab %>%
  dplyr::group_by(celltype, time, stim) %>%
  dplyr::summarise(
    `05%` = (cur_data()$p_non_zero > 0.05) %>% sum(),
    `10%` = (cur_data()$p_non_zero > 0.1) %>% sum(),
    `20%` = (cur_data()$p_non_zero > 0.2) %>% sum(),
    `30%` = (cur_data()$p_non_zero > 0.30) %>% sum(),
    `50%` = (cur_data()$p_non_zero > 0.50) %>% sum(),
    `70%` = (cur_data()$p_non_zero > 0.70) %>% sum(),
  ) %>%
  tidyr::pivot_longer(-c(celltype, time, stim)) %>%
  dplyr::mutate(celltype = stringr::str_remove_all(celltype, "RNA\\.|\\.\\."), ts = paste0(time, "_", stim)) %>%
  dplyr::filter(!celltype %in% c("Undefined")) %>%
  dplyr::mutate(celltype = factor(celltype, levels = c("CD4T", "CD8T", "Monocytes", "NK", "B", "mDC", "pDC", "Platelet"))) %>%
  dplyr::mutate(ts = factor(ts, levels = c("T0_RPMI", "T0_LPS", "T3m_RPMI", "T3m_LPS")))

p <- plot_tab %>%
  ggplot() +
  geom_col(aes(x = name, y = value, fill = celltype), position = "dodge2") +
  facet_wrap(~ts, ncol = 2) +
  labs(x = "Fraction of expressed genes", y = "Number of genes", fill = "Cell type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/overview/n_expressed_genes.per_condition.pdf")
ggsave(save_to, p, width = 10, height = 4)



eqtl_tab <- fread(file.path(proj_dir, "outputs/pseudo_bulk/overview/filtered.eGene_eVariants.FDR0.05.v3.csv"))
#
## Plot number of eQTLs
#
col_dict <- c("B" = "#E64B35", "CD4T" = "#4DBBD5", "CD8T" = "#00A087", "NK" = "#3C5488", "Monocytes" = "#F39B7F")
plot_size <- list("normal" = c("width" = 7, "height" = 3), "interaction" = c("width" = 15, "height" = 4))
for (mode in mode_vec) {
  min_isize <- ifelse(mode == "normal", 5, 20)
  run_pattern <- col_pattern_vec[mode] %>% as.vector()
  genesets <- eqtl_tab %>%
    dplyr::filter(stringr::str_detect(condition, run_pattern)) %>%
    dplyr::filter(q_value < fdr_max & !is.na(q_value)) %>%
    dplyr::group_by(celltype, condition) %>%
    dplyr::summarize(GeneSets = {
      geneset <- list(unique(feature_id))
      if (mode == "normal")
        names(geneset) <- cur_group()$celltype
      else
        names(geneset) <- cur_group() %>% unlist() %>% paste0(collapse = "-")
      geneset
    })

  mat <- make_comb_mat(genesets$GeneSets)
  mat <- mat[comb_size(mat) >= min_isize]

  row_annotations <- stringr::str_split(ComplexHeatmap::set_name(mat), pattern = "-|\\.", n = 2, simplify = TRUE)

  celltypes <- row_annotations[, 1]
  celltypes_col <- col_dict[celltypes]

  row_labels <- row_annotations[, ifelse(mode == "normal", 1, 2)]
  comb_col_vec <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4")[comb_degree(mat)]

  pw <- plot_size[[mode]]["width"]
  ph <- plot_size[[mode]]["height"]

  save_to <- file.path(proj_dir, "outputs/pseudo_bulk/overview", paste0(mode, ".eGene_number.FDR", fdr_max, ".upset_plot.v2.pdf"))
  pdf(save_to, width = pw, height = ph)
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
## Enrichment analysis. TODO: Create a supplementary table
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
    (function(vec) {clusterProfiler::bitr(names(vec), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") %>% dplyr::mutate(Beta = vec[SYMBOL])}) %>%
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
      (function(res) { res@result %>% dplyr::mutate(Celltype = per_celltype, Category = "GO") }),
    enrichPathway(gene_info$ENTREZID, pvalueCutoff = 1, readable = TRUE) %>%
      (function(res) { res@result %>% dplyr::mutate(Celltype = per_celltype, Category = "ReactomePathway") }),
    enrichKEGG(gene_info$ENTREZID, pvalueCutoff = 1, organism = "hsa") %>%
      (function(res) { res@result %>% dplyr::mutate(Celltype = per_celltype, Category = "KEGG") }),
    enrichWP(gene_info$ENTREZID, pvalueCutoff = 1, organism = "Homo sapiens") %>%
      (function(res) { res@result %>% dplyr::mutate(Celltype = per_celltype, Category = "WikiPathway") }),
    enrichDO(gene_info$ENTREZID, pvalueCutoff = 1, readable = TRUE) %>%
      (function(res) { res@result %>% dplyr::mutate(Celltype = per_celltype, Category = "Disease") })
  ) %>%
    Reduce(rbind, .) %>%
    data.table::fwrite(paste0(per_celltype, "_enrichment.csv"))
}

enr_tab <- paste0(celltype_vec, "_enrichment.csv") %>%
  lapply(function(p) data.table::fread(p)) %>%
  Reduce(rbind, .) %>%
  dplyr::mutate(Log2OddRatio = purrr::map2(GeneRatio, BgRatio, ~ log2(eval(parse(text = .x)) / eval(parse(text = .y)))), Log2OddRatio = unlist(Log2OddRatio))

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
      (function(tab) {dplyr::select(tab, -c(ID, Description, Category)) %>% as.matrix() %>% t()}) %>%
      Heatmap(name = pcat, column_title = pcat, col = col_fun_prop, show_heatmap_legend = FALSE, row_names_max_width = max_text_width(rownames(.), gp = gpar(fontsize = 7)))
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
highlight <- c("CD55", "ADCY3")
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
  geom_text_repel(aes(rank, p_value, label = feature_id), dplyr::filter(qtlrank_tab, !feature_id %in% highlight) %>% head(9), min.segment.length = 0, max.overlaps = 20, box.padding = 0.5) +
  geom_label_repel(aes(rank, p_value, label = feature_id), dplyr::filter(qtlrank_tab, feature_id %in% highlight), color = "red", box.padding = 7) +
  labs(x = "Rank by p-value", y = "-Log10(p-value)") +
  theme_classic()
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/overview", paste0("eqtl_rank.", celltype, ".pdf"))
ggsave(save_to, plot = myplot, width = 6, height = 4)


#
## Plot number of coloc loci
#
gwas_abbr <- c(
  "AsthmaDT1" = "ASDT1", "CoronaryHeartDisease" = "CHD", "ObesityClass1" = "OB1", "ObesityClass2" = "OB2",
  "ObesityClass3" = "OB3", "ThyroidCancer" = "TC", "AtopicDermatitis" = "AD", "Height" = "Height", "LungCancer" = "LC",
  "Schizophrenia" = "SZ", "Urate" = "UT", "CrohnsDisease" = "CrD", "InflammatoryBowelDisease" = "IBD",
  "SystemicLupusErythematosus" = "SLE", "UlcerativeColitis" = "UC", "AsthmaDT2" = "AS", "GoutDisease" = "GD",
  "OvarianCancer" = "OC", "Psoriasis" = "PR", "BodyMassIndex" = "BMI", "DiastolicBloodPressure" = "DBP",
  "RheumatoidArthritis" = "RA", "SystolicBloodPressure" = "SBP", "Type2Diabetes" = "T2D", "MultipleSclerosis" = "MS",
  "COVID19Release4" = "COVID19_R4", "HDLCholesterol" = "HDL", "LDLCholesterol" = "LDL", "Triglycerides" = "TRI",
  "BladderCancer" = "BC", "ColorectalCancer" = "CC", "ProstateCancer" = "PC", "AlzheimerDiseases" = "AzD",
  "COVID19Release7" = "COVID19"
)

# Number of independent loci per GWAS.
independent_loci <- fread(file.path(proj_dir, "outputs/public_gwas/all_of_indenpendent_loci_per_gwas.txt")) %>%
  dplyr::mutate(GWAS = stringr::str_replace_all(GWAS, "^(.*?)_(.*?)_(.*)", "\\2"))

plot_tab <- independent_loci %>%
  dplyr::group_by(GWAS) %>%
  dplyr::summarise(nr_independent_loci = n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(nr_independent_loci) %>%
  dplyr::mutate(
    GWAS_abbr = gwas_abbr[GWAS],
    Category = dplyr::case_when(
      GWAS_abbr %in% c("HDL", "Height", "BMI", "TRI", "LDL") ~ "General traits",
      GWAS_abbr %in% c("AD", "SLE", "GD", "PR", "RA", "MS") ~ "Autoimmune diseases",
      GWAS_abbr %in% c("PC", "TC", "LC", "BC", "CC", "OC") ~ "Cancer",
      GWAS_abbr %in% c("OB1", "OB2", "OB3") ~ "Obesity",
      GWAS_abbr %in% c("IBD", "UC", "CrD") ~ "IBD",
      GWAS_abbr %in% c("COVID19") ~ "Infectious diseases",
      T ~ "Others"),
    bar_label = paste0(GWAS_abbr, " (", nr_independent_loci, ")"),
    nr_independent_loci = dplyr::if_else(GWAS_abbr == "MS", as.integer(100), nr_independent_loci), # MS has 101 independent loci, just for ordering purpose.
  ) %>%
  dplyr::filter(!GWAS_abbr %in% c("COVID19_R4", "SBP", "DBP", "UT", "ASDT1")) %>%
  dplyr::mutate(
    GWAS_abbr = forcats::fct_reorder(GWAS_abbr, nr_independent_loci),
    angle = 90 - 360 * (1:n() - 0.5) / n(),
    hjust = ifelse(angle <= -90, 1, 0),
    angle = ifelse(angle <= -90, angle + 180, angle)
  )

p <- plot_tab %>%
  ggplot(aes(x = GWAS_abbr, y = nr_independent_loci, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = bar_label, hjust = hjust), angle = plot_tab$angle, size = 3) +
  scale_y_log10(limits = c(1, 20000)) +
  coord_polar(start = 0) +
  theme_void()
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/colocalization/nr_independent_loci_per_gwas.pdf")
ggsave(save_to, plot = p, width = 6, height = 6)


# Number of colocalizations.
ps_list <- data.frame("normal" = c(width = 8, height = 4), "interaction" = c(width = 9, height = 5))
for (per_mode in mode_vec) {
  coloc_files <- list.files(file.path(proj_dir, "outputs/pseudo_bulk/colocalization", per_mode), pattern = "*.csv", recursive = TRUE, full.names = TRUE)
  coloc_mat <- lapply(coloc_files, function(fp) {
    if (per_mode == "normal") {
      celltype <- fp %>% dirname() %>% basename()
      condition <- "common"
    } else if (per_mode == "interaction") {
      celltype <- fp %>% dirname() %>% dirname() %>% basename()
      condition <- fp %>% dirname() %>% basename()
    } else { stop(paste("Unknown mode", per_mode)) }

    data.table::fread(fp) %>%
      dplyr::mutate(Celltype = celltype, Condition = condition, outcome = as.character(outcome)) %>%
    dplyr::group_by(outcome, Celltype, Condition) %>%
    dplyr::filter(H4 >= min_h4)
  }) %>%
  Reduce(rbind, .)

  overlapped_loci <- dplyr::group_by(coloc_mat, outcome) %>%
    dplyr::summarise(nr_coloc = {
      per_outcome <- dplyr::cur_group()$outcome
      .indep_loci_pg <- dplyr::filter(independent_loci, GWAS == per_outcome) 
      .indep_loci_nr_pg <- nrow(.indep_loci_pg)
      .clumpped_snps_pg <- stringr::str_extract_all(.indep_loci_pg$SP2, "rs[0-9]+") %>% Reduce(c, .)
      cat("[I]:", per_outcome, "has", n(), "loci colocalized with eQTL.", "The GWAS has", length(.clumpped_snps_pg), "SNPs belongs to", .indep_loci_nr_pg, "independent loci. \n")

      dplyr::cur_data() %>%
        dplyr::rowwise() %>%
        dplyr::mutate(nr_overlapped_loci = sum(dplyr::filter(eqtl_tab, feature_id == exposure, celltype == celltype, condition == condition)$snp_id %in% .clumpped_snps_pg)) %>%
        tidyr::unnest(cols = c(nr_overlapped_loci)) %>%
        dplyr::mutate(nr_indep_loci = .indep_loci_nr_pg)
    }) %>%
    tidyr::unnest(cols = c(nr_coloc)) %>%
    dplyr::mutate(abbrev = gwas_abbr[outcome]) %>%
    dplyr::filter(!abbrev %in% c("COVID19_R4", "SBP", "DBP", "UT", "ASDT1"))

  prop_mat <- dplyr::group_by(overlapped_loci, abbrev, Celltype, Condition, nr_indep_loci) %>%
    dplyr::summarise(nr_overlaps = sum(nr_overlapped_loci>0)) %>%
    dplyr::mutate(coloc_prop = nr_overlaps / nr_indep_loci * 100) %>%
    tidyr::pivot_wider(id_cols = abbrev, names_from = c(Celltype, Condition), values_from = coloc_prop) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dplyr::across(-abbrev, ~ dplyr::if_else(is.na(.x), as.double(0), .x))) %>%
    dplyr::rename_with(~stringr::str_remove(.x, pattern = "_common$")) %>%
    (function(m) {tmp <- dplyr::select(m, -abbrev) %>% as.matrix(); rownames(tmp) <- m$abbrev; tmp})

  count_mat <- tidyr::pivot_wider(overlapped_loci, id_cols = abbrev, names_from = c(Celltype, Condition), values_from = nr_indep_loci, values_fn = length) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dplyr::across(-abbrev, ~ dplyr::if_else(is.na(.x), as.integer(0), .x))) %>%
    dplyr::rename_with(~stringr::str_remove(.x, pattern = "_common$")) %>%
    (function(m) {tmp <- m %>% dplyr::select(-abbrev) %>% as.matrix(); rownames(tmp) <- m$abbrev; tmp})

  save_to <- file.path(proj_dir, "outputs/pseudo_bulk/colocalization", paste0(per_mode, ".colocalization.per_celltype.v3.pdf"))
  col_fun_prop <- colorRamp2(c(0.0, max(prop_mat) + 0.15), c("gray95", "#2166AC"))

  pdf(save_to, width = ps_list["width", per_mode], height = ps_list["height", per_mode])
  row_ann <- rowAnnotation(Counts = anno_barplot(colSums(count_mat)), show_annotation_name = FALSE)
  col_ann <- columnAnnotation(Counts = anno_barplot(rowSums(count_mat)), show_annotation_name = FALSE)
  hm <- Heatmap(t(prop_mat),
    col = col_fun_prop, heatmap_legend_param = list(title = "Perc. of colocs", direction = "horizontal", legend_width = unit(6, "cm")),
    top_annotation = col_ann, right_annotation = row_ann,
    column_names_rot = 45, column_names_side = "top", column_dend_side = "bottom",
    row_names_max_width = max_text_width(colnames(prop_mat))
  )
  draw(hm, heatmap_legend_side = "top")
  dev.off()
}


#
## Example eQTL, G and GxE
#
work_dir <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL")
example_eqtl <- tibble::tribble(
  ~snp, ~feature,
  "rs11687089", "ADCY3", # Good common eQTL across different conditions.
  "rs2564978", "CD55", # Good example contrast pattern across cell types (i.e., rs2564978-C is associated higher expression in Monocytes but lower expression in CD4T/CD8T).
  "rs11080327", "SLFN5", # Also good ieQTL, but not cell-type specific. GETex (WB, 1.5e-49)
  # "rs883416", "SLFN5", # Also good ieQTL, but not cell-type specific. GETex (WB, 1.5e-49), replacint rs11080327.
)

apply(example_eqtl, 1, function(e) {
  per_snp <- e["snp"]
  per_feature <- e["feature"]
  # per_snp <- "rs2564978"; per_feature <- "CD55"
  gt_col <- paste0(per_snp, "_GT")

  per_infile <- file.path(work_dir, paste0(per_feature, "-", per_snp, ".QTL_info.csv"))
  eqtl_tab <- data.table::fread(per_infile) %>%
    dplyr::mutate(
      time = factor(dplyr::if_else(time == 0, "T0", "T3m"), levels = time_vec),
      stim = factor(dplyr::if_else(stim == 0, "RPMI", "LPS"), levels = stim_vec),
      donor_id = stringr::str_extract(V1, "300BCG[0-9]+"),
      celltype = factor(celltype, celltype_vec)) %>%
    dplyr::select(-V1) %>%
    as.data.frame() %>%
    (function(tab) {
      new_level <- NULL
      for (pg in unique(tab[, gt_col])) {
        pg_1 <- substr(pg, 1, 1)
        pg_2 <- substr(pg, 2, 2)
        if (pg_1 != pg_2) new_level <- paste0(c(pg_1, pg_1, pg_2), c(pg_1, pg_2, pg_2))
      }
      tab[gt_col] <- factor(tab[, gt_col], levels = new_level)

      tab
    }) %>%
    as.data.frame()

  y_lab <- paste0("Avg. exp.: ", per_feature)

  # Common effects
  plot <- ggboxplot(eqtl_tab, x = gt_col, y = per_feature, fill = gt_col, palette = "lancet", xlab = FALSE, ylab = y_lab, facet.by = "celltype")
  if (per_feature == "KANSL1") {
    plot <- plot + geom_pwc(label = "p.adj.signif", hide.ns = TRUE, vjust = 0.5, tip.length = 0, y.position = 7, step.increase = 0.03) +
      scale_y_break(c(10, 27), ticklabels = c(27.5)) +
      lims(y = c(0, 28))
  } else {
    plot <- plot + geom_pwc(label = "p.adj.signif", hide.ns = TRUE, vjust = 0.5, tip.length = 0)
  }
  plot <- facet(plot, facet.by = "celltype", nrow = 1)
  save_to <- file.path(work_dir, paste0("normal.", per_feature, "-", per_snp, ".all_cell_type.boxplot.pdf"))
  ggsave(save_to, plot = plot, width = 6, height = 4)

  # Per condition
  plot_tab <- dplyr::mutate(eqtl_tab, condition = paste0(time, "_", stim), condition = factor(condition, c("T0_RPMI", "T0_LPS", "T3m_RPMI", "T3m_LPS")))
  plot <- ggboxplot(plot_tab, x = gt_col, y = per_feature, fill = gt_col, palette = "lancet", xlab = FALSE, ylab = y_lab, facet.by = c("condition", "celltype")) +
    geom_pwc(label = "p", tip.length = 0, bracket.nudge.y = -0.3)
  plot <- facet(plot, facet.by = c("condition", "celltype"), nrow = 5)
  save_to <- file.path(work_dir, paste0("per_condition.", per_feature, "-", per_snp, ".all_cell_type.boxplot.pdf"))
  ggsave(save_to, plot = plot, width = 8, height = 7)

  # Interaction effects
  plot_tab <- lapply(condition_vec, function(pp) {
      per_cmp_vec <- stringr::str_split(pp, pattern = ".vs.|_", simplify = TRUE)
      per_time_vec <- union(per_cmp_vec[1], per_cmp_vec[3])
      per_stim_vec <- union(per_cmp_vec[2], per_cmp_vec[4])
      compby_time <- ifelse(length(per_time_vec) == 2, TRUE, FALSE)
      dplyr::filter(eqtl_tab, time %in% per_time_vec, stim %in% per_stim_vec) %>%
        dplyr::mutate(
          comparison = pp, Effect = effect_vec[pp],
          Treatment = dplyr::case_when(stim == "LPS" | (stim == "RPMI" & time == "T3m" & Effect == "BCG eff.") ~ "Stim", TRUE ~ "BL")
      )
    }) %>%
    Reduce(rbind, .) %>%
    dplyr::mutate(Treatment = factor(Treatment, levels = c("BL", "Stim")))

  lapply(celltype_vec, function(per_celltype, .plot_tab) {
    plot <- dplyr::filter(.plot_tab, celltype == per_celltype) %>%
      dplyr::mutate(cross_comp = paste0(Treatment, !!rlang::sym(gt_col))) %>%
      (function(tab) {
        y_lab <- paste("Avg. exp.:", per_feature, "in", per_celltype)

        plot_dot <- ggline(tab, x = "Treatment", y = per_feature, group = gt_col, color = gt_col, palette = "jama", ylab = y_lab, add = c("mean_ci", "jitter"), facet.by = "Effect", xlab = FALSE, legend.title = per_snp) %>%
          facet(facet.by = "Effect", nrow = 1, strip.position = "top") %>%
          ggpar(legend = "right")

        plot_box <- ggboxplot(tab, x = gt_col, y = per_feature, fill = "Treatment", palette = "lancet", ylab = y_lab, xlab = FALSE, facet.by = "Effect", legend.title = "Condition") %>%
          facet(facet.by = "Effect", nrow = 1) %>%
          ggpar(legend = "right") +
          theme(strip.background = element_blank(), strip.text = element_blank())

        plot_dot / plot_box
      })

    save_to <- file.path(work_dir, paste0("interaction.", per_feature, "-", per_snp, ".", per_celltype, ".box_and_dot.pdf"))
    ggsave(save_to, plot = plot, width = 7, height = 5)
  }, .plot_tab = plot_tab)

  NA
})


#
## Per genotype expression of ADCY3/CD55/SLFN5
#
gt_per_ind <- fread(file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/tar_snps.genotypes.csv")) %>%
  dplyr::rename_with(.fn = ~ stringr::str_remove_all(.x, "\\[[0-9]{1,2}\\]|:GT$|^# ")) %>%
  tidyr::pivot_longer(dplyr::starts_with("300"), values_to = "GT_01", names_to = "ids") %>%
  dplyr::mutate(GT_RA = dplyr::case_when(GT_01 %in% c("0|0") ~ paste0(REF, REF), GT_01 %in% c("1|0", "0|1") ~ paste0(REF, ALT), GT_01 %in% c("1|1") ~ paste0(ALT, ALT))) %>%
  dplyr::select(ID, GT_RA, ids) %>%
  tidyr::pivot_wider(names_from = ID, values_from = GT_RA)


#
## Allele-specific TF binding affinities
#
# abs_fc: fold-change of allele-specific binding affinity calculated by -log2(alt/ref). Positive values mean alt allele has higher binding affinity.
# motif_fc: fold-change of log2-ratio between motif P-values for the ref and alt alleles. Positive values indicate Alt-ASBs (preferred binding to the alt allele).
tar_celltype <- "CD4+ T"; tar_feature <- "SLFN5"; tar_snp <- "rs11080327"
tar_celltype <- "Monocytes"; tar_feature <- "CD55"; tar_snp <- "rs2564978"
all_tf_features <- list.files(file.path(proj_dir, "inputs//ADASTRA/TF")) %>% stringr::str_remove("_HUMAN.tsv") %>% unique()

## Allele-specific TF binding affinities, CL
tar_tissue <- c(
  "CD14+ monocytes", "GM10847 female B-cells Lymphoblastoid Cell Lines",
  "GM12878 female B-cells lymphoblastoid cell line", "GM12892 female B-cells lymphoblastoid cell line",
  "GM18951 B-cells Lymphoblastoid Cell Lines", "HUES64 embryonic stem cells",
  "ID00014 lymphoblastoid cell DiGeorge syndrome", "ID00015 lymphoblastoid cell",
  "KMS-11 Plasma cell myeloma", "MOLM-13 acute monocytic leukemia AML Homo sapiens",
  "MOLM14 Adult acute myeloid leukemia", "MV-4-11 Childhood acute monocytic leukemia",
  "OCI-LY3 diffuse large B-cell lymphoma", "OCI-LY7 diffuse large B-cell lymphoma",
  "SET2 Adult acute megakaryoblastic leukemia Homo Sapiens", "activated primary CD4+ T cells from blood",
  "adult erythroblasts", "iPS-derived dermal fibroblasts",
  "monocyte-derived macrophages from peripheral blood", "monocytes from periferal blood",
  "naive IgD+ B cells", "neural progenitors", "neutrophils from periferal blood", "primary B-CLL cells"
)

astfb_cl_tab <- data.table::fread(file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/ASTFB", paste0(tar_snp, ".astfb_cl.txt")))
tar_astfb_cl <- astfb_cl_tab %>%
  dplyr::filter(!is.na(es_mean_ref), !is.na(es_mean_alt)) %>%
  dplyr::mutate(asb_fc = -log2(fdrp_bh_alt / fdrp_bh_ref)) %>%
  dplyr::arrange(asb_fc) %>%
  dplyr::mutate(
    CL = stringr::str_replace_all(CL, "_+", " ") %>% stringr::str_remove_all(" ?.tsv"),
    asb_sig = dplyr::if_else(fdrp_bh_alt < 0.05 | fdrp_bh_ref < 0.05, "Yes", "No") %>% factor(levels = c("Yes", "No")),
  ) %>%
  dplyr::filter(CL %in% tar_tissue) %>%
  dplyr::mutate(CL = stringr::str_trunc(CL, width = 40), CL = forcats::fct_reorder(CL, asb_fc))

asb_cl_bar_plot <- ggplot(tar_astfb_cl, aes(x = asb_fc, y = CL, fill = asb_sig)) +
  geom_bar(stat = "identity") +
  scale_fill_npg() +
  labs(y = NULL, x = "Allelic binding FC (ALT/REF)", fill = "FDR < 0.05") +
  theme_classic() +
  theme(axis.ticks.y.left = element_blank(), axis.line.y.left = element_blank(), legend.position = "top")

save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/ASTFB", paste0(tar_snp, ".astfb_cl.v1.pdf"))
ggsave(save_to, plot = asb_cl_bar_plot, width = 6.5, height = 5)


# Allele-specific TF binding affinities, TF
DefaultAssay(pbmc_300bcg) <- "RNA"
tar_bcg_features <- rownames(pbmc_300bcg)

astfb_tf_tab <- data.table::fread(file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/ASTFB", paste0(tar_snp, ".astfb_tf.txt")))
tar_astfb <- astfb_tf_tab %>%
  dplyr::filter(TF %in% tar_bcg_features) %>%
  dplyr::mutate(asb_fc = -log2(fdrp_bh_alt / fdrp_bh_ref)) %>%
  dplyr::arrange(asb_fc) %>%
  dplyr::mutate(
    asb_sig = dplyr::if_else(fdrp_bh_alt < 0.05 | fdrp_bh_ref < 0.05, "Yes", "No"),
    asb_sig = factor(asb_sig, levels = c("Yes", "No")),
    motif_sig = dplyr::if_else(motif_log_palt > -log10(0.05) | motif_log_pref > -log10(0.05), "Yes", "No"),
    motif_sig = factor(motif_sig, levels = c("Yes", "No"))
  ) %>%
  dplyr::filter(asb_sig == "Yes" | motif_sig == "Yes") %>%
  dplyr::mutate(TF = forcats::fct_reorder(TF, asb_fc))

asb_bar_plot <- ggplot(tar_astfb, aes(x = asb_fc, y = TF, fill = asb_sig)) +
  geom_bar(stat = "identity") +
  scale_fill_npg() +
  labs(y = NULL, x = "Allelic binding FC (ALT/REF)") +
  theme_classic() +
  theme(axis.ticks.y.left = element_blank(), axis.line.y.left = element_blank(), legend.position = "none")

motif_bar_plot <- ggplot(tar_astfb, aes(x = motif_fc, y = TF, fill = motif_sig)) +
  geom_bar(stat = "identity") +
  scale_fill_npg() +
  labs(y = NULL, x = "Motif FC", fill = "FDR < 0.05") +
  theme_classic() +
  theme(axis.ticks.y.left = element_blank(), axis.text.y.left = element_blank(), axis.line.y.left = element_blank(), legend.position = "top")

asb_tf_bar_plot <- asb_bar_plot + motif_bar_plot
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/ASTFB", paste0(tar_snp, ".astfb_tf.pdf"))
ggsave(save_to, plot = asb_tf_bar_plot, width = 6.5, height = 3.5)

## Expression profiles on target features and prioritized TFs
DefaultAssay(pbmc_300bcg) <- "RNA"
coexp_features <- c(tar_feature, rev(levels(tar_astfb$TF)))

p <- DotPlot(pbmc_300bcg, features = rev(coexp_features)) + RotatedAxis() + coord_flip()
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.dot_plot.300BCG.pdf"))
ggsave(save_to, plot = p, width = 6, height = 4.5)

for (per_ts in time_stim_vec) {
  cells <- pbmc_300bcg@meta.data %>% as.data.frame() %>% dplyr::filter(ts == per_ts) %>% rownames()
  p <- DotPlot(pbmc_300bcg[, cells], features = rev(coexp_features)) + RotatedAxis() + coord_flip()
  save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.dot_plot.300BCG.", per_ts, ".pdf"))
  ggsave(save_to, plot = p, width = 6, height = 5)
}


## Co-expression of TF and target gene, linear mixed model, considering time and stim as random effect.
coexp_tab <- lapply(celltype_vec, function(cc) {
  avgexp_tab <- fread(file.path(proj_dir, "inputs/pseudo_bulk", cc, "phenotypes.tsv")) %>%
    dplyr::filter(feature_id %in% coexp_features) %>%
    tidyr::pivot_longer(-feature_id, names_to = "sample_id", values_to = "avgexp") %>%
    tidyr::pivot_wider(names_from = feature_id, values_from = avgexp)
  covar_tab <- fread(file.path(proj_dir, "inputs/pseudo_bulk", cc, "covariates_wpeer.tsv"))
  dplyr::inner_join(avgexp_tab, covar_tab, by = "sample_id") %>%
    tidyr::separate(sample_id, into = c("celltype", "x", "xx", "ids"), sep = "_", extra = "drop") %>%
    dplyr::select(-x, -xx) %>%
    dplyr::mutate(stim = factor(stim), time = factor(time), gender = factor(gender))
}) %>%
  Reduce(rbind, .) %>%
  dplyr::left_join(gt_per_ind, by = "ids") %>%
  tidyr::pivot_longer(-c(celltype, ids, !!tar_feature, time, stim, age, gender, dplyr::starts_with("rs"), dplyr::starts_with("PEER_")))

for (per_ct in celltype_vec) {
  for (per_tf in unique(coexp_tab$name)) {
    cat("---------", per_tf, "---------\n")
    m <- dplyr::filter(coexp_tab, celltype == per_ct, name == per_tf) %>%
      lmer(CD55 ~ value + (value | stim) + (value | time) + (1 | time) + (1 | stim) + age + gender, data = .)

    summary(m) %>% print()
  }
}


## Co-expression of CD55/SLFN4 and TFs by CSCORE
estimate_coexp <- function(obj, cells, features, fig_save_to = NULL, do_cluster = "none", annotate = TRUE, label_list = NULL, fdr = 5e-2, min_features = 20, ...) {
  shared_features <- intersect(features, rownames(obj))
  if (length(shared_features) != length(features)) { warning("Some features are not in the object"); features <- shared_features }

  # Conexpression
  coexp <- obj[, cells] %>% CSCORE(genes = features)

  # Obtain estimtations
  est <- as.data.frame(coexp$est) %>% dplyr::mutate(tf_name_x = rownames(.)) %>% tidyr::pivot_longer(-tf_name_x, names_to = "tf_name_y", values_to = "correlation")

  # Adjust p-values
  p_val_adj <- coexp$p_value %>%
    (function(mat) {
      mat[upper.tri(mat)] <- p.adjust(mat[upper.tri(mat)], method = "BH")
      mat[lower.tri(mat)] <- p.adjust(mat[lower.tri(mat)], method = "BH")
      as.data.frame(mat)
    }) %>%
    dplyr::mutate(tf_name_x = rownames(.)) %>%
    tidyr::pivot_longer(-tf_name_x, names_to = "tf_name_y", values_to = "p_value_adj")
  p_val <- as.data.frame(coexp$p_value) %>%
    dplyr::mutate(tf_name_x = rownames(.)) %>%
    tidyr::pivot_longer(-tf_name_x, names_to = "tf_name_y", values_to = "p_value") %>%
    dplyr::inner_join(p_val_adj, by = c("tf_name_x", "tf_name_y"))

  # Obtain statistic
  test_stat <- as.data.frame(coexp$test_stat) %>%
    dplyr::mutate(tf_name_x = rownames(.)) %>%
    tidyr::pivot_longer(-tf_name_x, names_to = "tf_name_y", values_to = "statistic")

  coexp_tab <- dplyr::left_join(est, p_val, by = c("tf_name_x", "tf_name_y")) %>%
    dplyr::left_join(test_stat, by = c("tf_name_x", "tf_name_y")) %>%
    dplyr::mutate(label = ifelse(p_value_adj < 0.05 & abs(correlation) >= 5e-2 & correlation < 0.999, format(round(correlation, 2)), ""))

  if (!is.null(fig_save_to)) {
    features_order_x <- features_order_y <- features
    if (do_cluster %in% c("x", "y", "both")) {
      coexp_hc <- hclust(as.dist(coexp$est))
      if (do_cluster %in% c("x", "both")) features_order_x <- features[coexp_hc$order]
      if (do_cluster %in% c("y", "both")) features_order_y <- features[coexp_hc$order]
    }

    plot_tab <- dplyr::mutate(coexp_tab, tf_name_y = factor(tf_name_y, levels = rev(features_order_y)), tf_name_x = factor(tf_name_x, levels = features_order_x))
    p <- ggplot(plot_tab) + geom_tile(aes(x = tf_name_x, y = tf_name_y, fill = correlation))
    if (length(features) <= min_features && annotate) {p <- p + geom_text(aes(x = tf_name_x, y = tf_name_y, label = label))}
    p <- p + scale_x_discrete(position = "top") +
      scale_fill_gradient2(low = "darkblue", mid = "gray98", high = "darkred", midpoint = 0, limits = c(-1, 1)) +
      labs(x = "", y = "", fill = "Co-expression") +
      theme_classic() +
      theme(axis.line = element_blank(), legend.position = "top")
    if (length(features) <= min_features) {
      p <- p + theme(axis.text.x.top = element_text(angle = 45, hjust = 0))
    } else if (length(label_list) > 0) {
      p <- p + scale_y_discrete(breaks = label_list, position = "right") + scale_x_discrete(breaks = NULL)
    } else {
      p <- p + theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "top")
    }
    ggsave(fig_save_to, p, ...)
  }
  invisible(coexp_tab)
}

## 300BCG
potential_tf_from_lieature <- c("CD55", "IL6", "IL1B", "IFNG", "TNFRSF1A", "TNFRSF1B", "NOD2", "LAMTOR5")
# potential_tf_from_lieature %in% rownames(pbmc_300bcg)

# pbmc_300BCG <- readRDS(file.path(proj_dir, "inputs/sc_rnaseq/300BCG.rds"))
# coexp_features <- c("SLFN5", "STAT1", "IRF4", "HIC1", "GATA1", "CTCF", "PRDM1", "PBX1", "PATZ1", "GATA2", "SP3", "SPI1")
cells <- pbmc_300bcg@meta.data %>% as.data.frame() %>% dplyr::filter(clusters1 %in% c("Monocytes")) %>% rownames()
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.coexpression.single_cell.300bcg.monocyte.pdf"))
estimate_coexp(pbmc_300bcg, cells, coexp_features, save_to, width = 6, height = 5)
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.from_literature.coexpression.single_cell.300bcg.monocyte.pdf"))
estimate_coexp(pbmc_300bcg, cells, potential_tf_from_lieature, save_to, width = 6, height = 5)

cells <- pbmc_300bcg@meta.data %>% as.data.frame() %>% dplyr::filter(clusters1 %in% c("CD4+ T")) %>% rownames()
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.coexpression.single_cell.300bcg.CD4T.pdf"))
estimate_coexp(pbmc_300bcg, cells, coexp_features, save_to, width = 6, height = 5)

cells <- pbmc_300bcg@meta.data %>% as.data.frame() %>% dplyr::filter(clusters1 %in% c("CD8+ T")) %>% rownames()
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.coexpression.single_cell.300bcg.CD8T.pdf"))
estimate_coexp(pbmc_300bcg, cells, coexp_features, save_to, width = 6, height = 5)

cells <- pbmc_300bcg@meta.data %>% as.data.frame() %>% dplyr::filter(clusters1 %in% c("NK")) %>% rownames()
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.coexpression.single_cell.300bcg.NK.pdf"))
estimate_coexp(pbmc_300bcg, cells, coexp_features, save_to, width = 6, height = 5)

cells <- pbmc_300bcg@meta.data %>% as.data.frame() %>% dplyr::filter(clusters1 %in% c("B")) %>% rownames()
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.coexpression.single_cell.300bcg.B.pdf"))
estimate_coexp(pbmc_300bcg, cells, coexp_features, save_to, width = 6, height = 5)

cells <- pbmc_300bcg@meta.data %>% as.data.frame() %>% dplyr::filter(ts == "T0_RPMI", clusters1 %in% c("Monocytes")) %>% rownames()
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.coexpression.single_cell.300bcg.monocyte.T0_RPMI.pdf"))
estimate_coexp(pbmc_300bcg, cells, coexp_features, save_to, width = 6, height = 5)
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.from_literature.coexpression.single_cell.300bcg.monocyte.T0_RPMI.pdf"))
estimate_coexp(pbmc_300bcg, cells, potential_tf_from_lieature, save_to, width = 6, height = 5)

cells <- pbmc_300bcg@meta.data %>% as.data.frame() %>% dplyr::filter(ts == "T0_LPS", clusters1 %in% c("Monocytes")) %>% rownames()
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.coexpression.single_cell.300bcg.monocyte.T0_LPS.pdf"))
estimate_coexp(pbmc_300bcg, cells, coexp_features, save_to, width = 6, height = 5)
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.from_literature.coexpression.single_cell.300bcg.monocyte.T0_LPS.pdf"))
estimate_coexp(pbmc_300bcg, cells, potential_tf_from_lieature, save_to, width = 6, height = 5)

cells <- pbmc_300bcg@meta.data %>% as.data.frame() %>% dplyr::filter(ts == "T3m_RPMI", clusters1 %in% c("Monocytes")) %>% rownames()
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.coexpression.single_cell.300bcg.monocyte.T3m_RPMI.pdf"))
estimate_coexp(pbmc_300bcg, cells, coexp_features, save_to, width = 6, height = 5)
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.from_literature.coexpression.single_cell.300bcg.monocyte.T3m_RPMI.pdf"))
estimate_coexp(pbmc_300bcg, cells, potential_tf_from_lieature, save_to, width = 6, height = 5)

cells <- pbmc_300bcg@meta.data %>% as.data.frame() %>% dplyr::filter(ts == "T3m_LPS", clusters1 %in% c("Monocytes")) %>% rownames()
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.coexpression.single_cell.300bcg.monocyte.T3m_LPS.pdf"))
estimate_coexp(pbmc_300bcg, cells, coexp_features, save_to, width = 6, height = 5)
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.from_literature.coexpression.single_cell.300bcg.monocyte.T3m_LPS.pdf"))
estimate_coexp(pbmc_300bcg, cells, potential_tf_from_lieature, save_to, width = 6, height = 5)

## Berlin cohort
pbmc_berlin <- readRDS(file.path(proj_dir, "inputs/sc_rnaseq/berlin.rds"))
cells <- pbmc_berlin@meta.data %>% as.data.frame() %>% dplyr::filter(Condition == "Control", celltypeL0 %in% c("CD14+ Monocytes", "CD16+ Monocytes")) %>% rownames()
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.coexpression.single_cell.berlin.monocyte.pdf"))
estimate_coexp(pbmc_berlin, cells, coexp_features, save_to, width = 6, height = 5)
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.from_literature.coexpression.single_cell.berlin.monocyte.pdf"))
estimate_coexp(pbmc_berlin, cells, potential_tf_from_lieature, save_to, width = 6, height = 5)

cells <- pbmc_berlin@meta.data %>% as.data.frame() %>% dplyr::filter(Condition == "Control", celltypeL0 %in% c("B")) %>% rownames()
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.coexpression.single_cell.berlin.b_cell.pdf"))
estimate_coexp(pbmc_berlin, cells, coexp_features, save_to, width = 6, height = 5)
 
## Bonn cohort
pbmc_bonn <- readRDS(file.path(proj_dir, "inputs/sc_rnaseq/bonn.rds"))
cells <- pbmc_bonn@meta.data %>% as.data.frame() %>% dplyr::filter(Condition == "Control", celltypeL0 %in% c("CD14+ Monocytes", "CD16+ Monocytes")) %>% rownames()
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.coexpression.single_cell.bonn.monocyte.pdf"))
estimate_coexp(pbmc_bonn, cells, coexp_features, save_to, width = 6, height = 5)
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.from_literature.coexpression.single_cell.bonn.monocyte.pdf"))
estimate_coexp(pbmc_bonn, cells, potential_tf_from_lieature, save_to, width = 6, height = 5)

cells <- pbmc_bonn@meta.data %>% as.data.frame() %>% dplyr::filter(Condition == "Control", celltypeL0 %in% c("B")) %>% rownames()
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.coexpression.single_cell.bonn.b_cell.pdf"))
estimate_coexp(pbmc_bonn, cells, coexp_features, save_to, width = 6, height = 5)

## MHH50 cohort
pbmc_mhh50 <- readRDS(file.path(proj_dir, "outputs/COVID_MHH50/PBMC_scRNAseq.rds"))
cells <- pbmc_mhh50@meta.data %>% as.data.frame() %>% dplyr::filter(Disease == "Convalescent", celltypeL1 %in% c("merged.cMonos")) %>% rownames()
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.coexpression.single_cell.mhh50.monocyte.pdf"))
estimate_coexp(pbmc_mhh50, cells, coexp_features, save_to, width = 6, height = 5)

cells <- pbmc_mhh50@meta.data %>% as.data.frame() %>% dplyr::filter(Disease == "Convalescent", celltypeL1 %in% c("B")) %>% rownames()
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(tar_feature, "_TF.coexpression.single_cell.mhh50.b_cell.pdf"))
estimate_coexp(pbmc_mhh50, cells, coexp_features, save_to, "both", width = 6, height = 5)


#
## Estimate the gene expression of target genes in COVID-19 context
#
pbmc_file <- file.path(proj_dir, "outputs/COVID_MHH50/PBMC_scRNAseq.rds")
if (file.exists(pbmc_file)) {
  pbmc_mhh50 <- readRDS(pbmc_file)
} else {
  pbmc_mhh50 <- readRDS("/vol/projects/BIIM/Covid_50MHH/NubesSubmission/scRNA/PBMC_scRNAseq.rds")
  pbmc_mhh50$Disease <- c("Active" = "Hospitalized", "Convalescent" = "Convalescent")[pbmc_mhh50$Disease] %>% factor(levels = c("Convalescent", "Hospitalized"))

  # bcftools query -H -f '%CHROM,%POS,%ID,%REF,%ALT[,%DS]\n' -i 'ID=="rs11080327"' /vol/projects/BIIM/Covid_50MHH/ASoC/outputs/genotypes/variants_fl_addFI.vcf.gz > ~/Documents/projects/wp_bcg_eqtl/outputs/pseudo_bulk/example_eQTL/rs11080327_mhh50_genotype_DS.csv
  map_tab <- fread("/vol/projects/zzhang/projects/wp_covid19_mhh50/inputs/idmapping/id_mapping.txt")

  # Adding rs11080327
  gt_order <- c("GG", "GA", "AA")
  gt_per_ind <- fread(file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/rs11080327_mhh50_genotype_GT.csv")) %>%
    dplyr::rename_with(.fn = ~ stringr::str_remove_all(.x, "\\[[0-9]{1,2}\\]|:GT$|^# ")) %>%
    tidyr::pivot_longer(-c(CHROM, POS, ID, REF, ALT), values_to = "GT_01", names_to = "ids") %>%
    dplyr::mutate(GT_RA = dplyr::case_when(GT_01 %in% c("0|0") ~ paste0(REF, REF), GT_01 %in% c("1|0", "0|1") ~ paste0(REF, ALT), GT_01 %in% c("1|1") ~ paste0(ALT, ALT))) %>%
    dplyr::select(ID, CHROM, POS, REF, ALT, GT_RA, ids) %>%
    tidyr::pivot_wider(names_from = c(ID, CHROM, POS, REF, ALT), values_from = GT_RA) %>%
    dplyr::left_join(map_tab, by = c("ids" = "genoID")) %>%
    dplyr::select(patientID, dplyr::starts_with("rs"))

  ds_per_ind <- fread(file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/rs11080327_mhh50_genotype_DS.csv")) %>%
    dplyr::rename_with(.fn = ~ stringr::str_remove_all(.x, "\\[[0-9]{1,2}\\]|:DS$|^# ")) %>%
    tidyr::pivot_longer(-c(CHROM, POS, ID, REF, ALT), values_to = "GT_DS", names_to = "ids") %>%
    dplyr::select(ID, CHROM, POS, REF, ALT, GT_DS, ids) %>%
    tidyr::pivot_wider(names_from = c(ID, CHROM, POS, REF, ALT), values_from = GT_DS) %>%
    dplyr::left_join(map_tab, by = c("ids" = "genoID")) %>%
    dplyr::select(patientID, dplyr::starts_with("rs"))

  new_meta_data <- pbmc_mhh50@meta.data %>%
    dplyr::left_join(gt_per_ind, by = c("patient" = "patientID")) %>%
    dplyr::left_join(ds_per_ind, by = c("patient" = "patientID"), suffix = c(".gt", ".ds"))

  # Genotype. Hospitalized and convalescent
  pbmc_mhh50$rs11080327_gt_pc_2c <- new_meta_data %>%
    dplyr::mutate(gt_pc = paste(Disease, rs11080327_17_35244527_G_A.gt, sep = "_")) %>%
    dplyr::pull(gt_pc, cellbarcodes) %>%
    factor(levels = c(paste0("Hospitalized_", gt_order), paste0("Convalescent_", gt_order)))

  # Genotype. Mild, severe, and convalescent
  pbmc_mhh50$rs11080327_gt_pc_3c <- new_meta_data %>%
    dplyr::mutate(gt_pc = paste(Severity, rs11080327_17_35244527_G_A.gt, sep = "_")) %>%
    dplyr::pull(gt_pc, cellbarcodes) %>%
    factor(levels = c(paste0("severe_", gt_order), paste0("mild_", gt_order), paste0("post_", gt_order)))

  # Dosage
  pbmc_mhh50$rs11080327_ds_pc <- new_meta_data %>% dplyr::pull(rs11080327_17_35244527_G_A.ds, cellbarcodes)


  # Adding rs2564978, TODO: redundant code, merging them in to the above section.
  gt_order <- c("TT", "TC", "CC")
  gt_per_ind <- fread(file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/rs2564978_mhh50_genotype_GT.csv")) %>%
    dplyr::rename_with(.fn = ~ stringr::str_remove_all(.x, "\\[[0-9]{1,2}\\]|:GT$|^# ")) %>%
    tidyr::pivot_longer(-c(CHROM, POS, ID, REF, ALT), values_to = "GT_01", names_to = "ids") %>%
    dplyr::mutate(GT_RA = dplyr::case_when(GT_01 %in% c("0|0") ~ paste0(REF, REF), GT_01 %in% c("1|0", "0|1") ~ paste0(REF, ALT), GT_01 %in% c("1|1") ~ paste0(ALT, ALT))) %>%
    dplyr::select(ID, CHROM, POS, REF, ALT, GT_RA, ids) %>%
    tidyr::pivot_wider(names_from = c(ID, CHROM, POS, REF, ALT), values_from = GT_RA) %>%
    dplyr::left_join(map_tab, by = c("ids" = "genoID")) %>%
    dplyr::select(patientID, dplyr::starts_with("rs"))

  ds_per_ind <- fread(file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/rs2564978_mhh50_genotype_DS.csv")) %>%
    dplyr::rename_with(.fn = ~ stringr::str_remove_all(.x, "\\[[0-9]{1,2}\\]|:DS$|^# ")) %>%
    tidyr::pivot_longer(-c(CHROM, POS, ID, REF, ALT), values_to = "GT_DS", names_to = "ids") %>%
    dplyr::select(ID, CHROM, POS, REF, ALT, GT_DS, ids) %>%
    tidyr::pivot_wider(names_from = c(ID, CHROM, POS, REF, ALT), values_from = GT_DS) %>%
    dplyr::left_join(map_tab, by = c("ids" = "genoID")) %>%
    dplyr::select(patientID, dplyr::starts_with("rs"))

  new_meta_data <- pbmc_mhh50@meta.data %>%
    dplyr::left_join(gt_per_ind, by = c("patient" = "patientID")) %>%
    dplyr::left_join(ds_per_ind, by = c("patient" = "patientID"), suffix = c(".gt", ".ds"))

  # Genotype. Hospitalized and convalescent
  pbmc_mhh50$rs2564978_gt_pc_2c <- new_meta_data %>%
    dplyr::mutate(gt_pc = paste(Disease, rs2564978_1_207321071_T_C.gt, sep = "_")) %>%
    dplyr::pull(gt_pc, cellbarcodes) %>%
    factor(levels = c(paste0("Hospitalized_", gt_order), paste0("Convalescent_", gt_order)))

  # Genotype. Mild, severe, and convalescent
  pbmc_mhh50$rs2564978_gt_pc_3c <- new_meta_data %>%
    dplyr::mutate(gt_pc = paste(Severity, rs2564978_1_207321071_T_C.gt, sep = "_")) %>%
    dplyr::pull(gt_pc, cellbarcodes) %>%
    factor(levels = c(paste0("severe_", gt_order), paste0("mild_", gt_order), paste0("post_", gt_order)))

  # Dosage
  pbmc_mhh50$rs2564978_ds_pc <- new_meta_data %>% dplyr::pull(rs2564978_1_207321071_T_C.ds, cellbarcodes)

  saveRDS(pbmc_mhh50, pbmc_file)
}


# Differential expression analysis between conditions, i.e., hospitalized vs convalescent
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/MHH50.two_classes.SLFN5.de_gene.by_condition.csv")
if (!file.exists(save_to)) {
  mhh50_deg <- NULL
  for (pct in unique(pbmc_mhh50$celltypeL1)) {
    if (pct %in% c("Plasmablast")) next

    mhh50_subset <- pbmc_mhh50["SLFN5", pbmc_mhh50$celltypeL1 == pct]
    Idents(mhh50_subset) <- "Disease"

    mhh50_deg <- FindMarkers(mhh50_subset, ident.1 = "Hospitalized", ident.2 = "Convalescent", logfc.threshold = 0.01) %>%
      dplyr::mutate(Gene = rownames(.), Comparison = "hospitalized.vs.convalescent", Cohort = "MHH50", Celltype = pct) %>%
      rbind(mhh50_deg)
  }
  mhh50_deg %>% dplyr::mutate(p_val_adj = p.adjust(p_val)) %>% data.table::fwrite(save_to)

  Idents(pbmc_mhh50) <- "celltypeL1"

  # Per condition
  p <- DotPlot(pbmc_mhh50, features = "SLFN5", group.by = "Disease", cols = "RdBu", idents = c("CD4.T")) +
    scale_size(range = c(6, 12), breaks = c(32, 32.5, 33)) +
    labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p_save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/MHH50.two_classes.SLFN5.de_gene.CD4T.dot_plot.pdf")
  ggsave(p_save_to, plot = p, width = 4.25, height = 3.25)

  p <- VlnPlot(pbmc_mhh50, features = "SLFN5", split.by = "Disease", idents = c("CD4.T"), pt.size = 0) +
    labs(x = NULL, y = NULL, title = NULL) +
    scale_fill_manual(values = c("Convalescent" = "#1f78b4", "Hospitalized" = "#e31a1c")) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    coord_flip()
  p_save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/MHH50.two_classes.SLFN5.de_gene.CD4T.violin_plot.pdf")
  ggsave(p_save_to, plot = p, width = 5, height = 3.5)

  p <- FeaturePlot(pbmc_mhh50, features = "SLFN5", split.by = "Disease", order = TRUE) + patchwork::plot_layout(ncol = 1, nrow = 2)
  p_save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/MHH50.two_classes.SLFN5.de_gene.CD4T.feature_plot.pdf")
  ggsave(p_save_to, plot = p, width = 4, height = 8)

  # Per genotype per condition
  p <- DotPlot(pbmc_mhh50, features = "SLFN5", group.by = "gt_pc_2c", cols = "RdBu", idents = c("CD4.T")) +
    scale_size(range = c(6, 12), breaks = c(30, 34, 38)) +
    labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p_save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/MHH50.two_classes.SLFN5.de_gene.per_gt.CD4T.dot_plot.pdf")
  ggsave(p_save_to, plot = p, width = 4.5, height = 3)

  p <- VlnPlot(pbmc_mhh50, features = "SLFN5", split.by = "gt_pc_2c", idents = c("CD4.T"), pt.size = 0) +
    labs(x = NULL, y = NULL, title = NULL) +
    scale_fill_manual(values = c("Convalescent_GG" = "white", "Convalescent_GA" = "gray", "Convalescent_AA" = "black",
                                 "Hospitalized_GG" = "white", "Hospitalized_GA" = "gray", "Hospitalized_AA" = "black")) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    coord_flip()
  p_save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/MHH50.two_classes.SLFN5.de_gene.per_gt.CD4T.violin_plot.pdf")
  ggsave(p_save_to, plot = p, width = 5.5, height = 4)

  p <- FeaturePlot(pbmc_mhh50, features = "SLFN5", split.by = "gt_pc_2c", order = TRUE) + patchwork::plot_layout(ncol = 3, nrow = 2)
  p_save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/MHH50.two_classes.SLFN5.de_gene.per_gt.CD4T.feature_plot.pdf")
  ggsave(p_save_to, plot = p, width = 4, height = 8)
}



# Differential expression genes between groups of patients with different genotypes.
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/SLFN5.deg_by_sevrity.MHH50.three_classes.csv")
if (!file.exists(save_to)) {
  mhh50_deg <- NULL
  for (pct in unique(pbmc_mhh50$celltypeL1)) {
    if (pct %in% c("Plasmablast", "ncMono", "pDC", "mDC")) next

    mhh50_subset <- pbmc_mhh50["SLFN5", pbmc_mhh50$celltypeL1 == pct]
    Idents(mhh50_subset) <- "Severity"

    cat("Estimating DE", pct, "\n")
    pairs <- combn(c("severe", "mild", "post"), 2)
    for (pp in 1:ncol(pairs)) {
      ident_1 <- pairs[1, pp]
      ident_2 <- pairs[2, pp]
      mhh50_deg <- FindMarkers(mhh50_subset, ident.1 = ident_1, ident.2 = ident_2, logfc.threshold = 0.01) %>%
        dplyr::mutate(Gene = rownames(.), Comparison = paste0(ident_1, ".vs.", ident_2), Cohort = "MHH50", Celltype = pct) %>%
        rbind(mhh50_deg)
    }
  }
  mhh50_slfn_deg <- mhh50_slfn_deg %>% as.data.table() %>% dplyr::mutate(p_val_adj = p.adjust(p_val), Celltype = dplyr::if_else(Celltype == "merged.cMonos", "Mono", Celltype))
  fwrite(mhh50_slfn_deg, save_to)
}

# Used in the results.
mhh50_slfn_deg <- NULL
for (pct in unique(pbmc_mhh50$celltypeL1)) {
  if (pct %in% c("Plasmablast", "ncMono", "pDC", "mDC")) next

  for (psv in c("Hospitalized", "Convalescent")) {
    mhh50_subset <- pbmc_mhh50[, pbmc_mhh50$celltypeL1 == pct & pbmc_mhh50$Disease == psv]

    Idents(mhh50_subset) <- "gt_pc_2c"
    cat("Estimating DE", pct, psv, "\n")

    pairs <- combn(sort(unique(mhh50_subset$gt_pc_2c)), 2)
    for (pp in 1:ncol(pairs)){
      ident_1 <- pairs[1, pp]
      ident_2 <- pairs[2, pp]
      mhh50_slfn_deg <- FindMarkers(mhh50_subset, ident.1 = ident_1, ident.2 = ident_2, features = c("SLFN5"), logfc.threshold = 0, test.use = "LR", latent.vars = c("gender", "Age")) %>%
        dplyr::mutate(Gene = "SLFN5", Comparison = paste0(ident_1, ".vs.", ident_2), Cohort = "MHH50", Celltype = pct, Test = "LR_adj_age_sex") %>%
        rbind(mhh50_slfn_deg)
    }
  }
}
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/SLFN5.deg_by_genotype.MHH50.two_classes.csv")
mhh50_slfn_deg <- mhh50_slfn_deg %>% as.data.table() %>% dplyr::mutate(p_val_adj = p.adjust(p_val), Celltype = dplyr::if_else(Celltype == "merged.cMonos", "Mono", Celltype))
fwrite(mhh50_slfn_deg, save_to)


# Differential expression analysis between conditions.
mhh50_slfn_deg <- NULL
for (pct in unique(pbmc_mhh50$celltypeL1)) {
  if (pct %in% c("Plasmablast", "ncMono", "pDC", "mDC")) next

  for (psv in c("severe", "mild", "post")) {
    mhh50_subset <- pbmc_mhh50[, pbmc_mhh50$celltypeL1 == pct & pbmc_mhh50$Severity == psv]

    Idents(mhh50_subset) <- "gt_pc_3c"
    cat("Estimating DE", pct, psv, "\n")

    pairs <- combn(sort(unique(mhh50_subset$gt_pc_3c)), 2)
    for (pp in 1:ncol(pairs)){
      ident_1 <- pairs[1, pp]
      ident_2 <- pairs[2, pp]
      mhh50_slfn_deg <- FindMarkers(mhh50_subset, ident.1 = ident_1, ident.2 = ident_2, features = c("SLFN5"), logfc.threshold = 0, test.use = "LR", latent.vars = c("gender", "Age")) %>%
        dplyr::mutate(Gene = "SLFN5", Comparison = paste0(ident_1, ".vs.", ident_2), Cohort = "MHH50", Celltype = pct, Test = "LR_adj_age_sex") %>%
        rbind(mhh50_slfn_deg)
    }
  }
}
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/SLFN5.deg_by_genotype.MHH50.three_classes.csv")
mhh50_slfn_deg <- mhh50_slfn_deg %>% as.data.table() %>% dplyr::mutate(p_val_adj = p.adjust(p_val), Celltype = dplyr::if_else(Celltype == "merged.cMonos", "Mono", Celltype))
fwrite(mhh50_slfn_deg, save_to)

rownames(mhh50_deg) <- NULL
mhh50_deg <- mhh50_deg %>% dplyr::mutate(p_val_adj = p.adjust(p_val), Celltype = dplyr::if_else(Celltype == "merged.cMonos", "Mono", Celltype))
mhh50_deg %>% data.table::fwrite(file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/SLFN5.differential_expression.MHH50.three_classes.csv"))

p <- mhh50_deg %>%
  dplyr::filter(Comparison != "severe.vs.mild") %>%
  dplyr::mutate(Celltype = forcats::fct_reorder(Celltype, avg_log2FC)) %>%
  ggplot() +
  geom_line(aes(x = avg_log2FC, y = Celltype, group = Celltype), linewidth = 7, alpha = 0.3, lineend = "round") +
  geom_point(aes(x = avg_log2FC, y = Celltype, color = Comparison, size = -log10(p_val_adj))) +
  scale_size_continuous(range = c(1, 5), breaks = c(5, 15, 25)) +
  scale_color_jama() +
  labs(x = "Average log2(fold-change)", y = NULL) +
  theme_classic()
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/SLFN5.differential_expression.MHH50.three_classes.pdf")
ggsave(save_to, plot = p, height = 2.5, width = 5)


group_by_col <- "gt_pc_3c"
dp_1 <- DotPlot(pbmc_mhh50, feature = "SLFN5", idents = c("CD4.T"), group.by = group_by_col) +
  labs(x = "CD4.T") +
  labs(y = NULL, title = NULL) +
  scale_color_gradient(low = "gray95", high = "blue", breaks = c(-1.5, 0, 1.5)) +
  scale_size(range = c(1, 8)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/SLFN5_expression_per_genotype.CD4T.three_classes.dot_plot.pdf"), plot = dp_1, height = 3.5, width = 4.5)

dp_2 <- DotPlot(pbmc_mhh50, feature = "SLFN5", idents = c("CD8.T"), group.by = group_by_col) +
  labs(x = "CD8.T") +
  labs(y = NULL, title = NULL) +
  scale_color_gradient(low = "gray95", high = "blue", breaks = c(-1.5, 0, 1.5)) +
  scale_size(range = c(1, 8)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/SLFN5_expression_per_genotype.CD8T.three_classes.dot_plot.pdf"), plot = dp_2, height = 3.5, width = 4.5)

# dp <- (dp_1 + dp_2 + plot_layout(width = c(11, 12))) &
#   labs(y = NULL, title = NULL) &
#   scale_color_gradient2(low = "#00A087", mid = "gray75", high = "#DC0000") & #, limits = c(-1, 0.25), breaks = c(-1, 0, 0.25)) &
#   scale_size(range = c(1, 8))
# ggsave(file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/SLFN5_expression_per_genotype.three_classes.dot_plot.pdf"), plot = dp, height = 3.5, width = 6)

fp <- FeaturePlot(pbmc_mhh50, feature = "SLFN5", split.by = group_by_col, ncol = 3, order = TRUE)
ggsave(file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/SLFN5_expression_per_genotype.three_classes.feature_plot.pdf"), plot = fp, width = 9)

vp <- VlnPlot(pbmc_mhh50, feature = "SLFN5", idents = c("CD4.T", "CD8.T"), pt.size = 0, split.by = group_by_col) +
  labs(x = NULL, title = NULL) +
  theme(legend.position = "top")
ggsave(file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/SLFN5_expression_per_genotype.three_classes.violin_plot.pdf"), plot = vp, height = 3)


# eQTL effect pseudobulk
avg_exp_tab <- pbmc_mhh50["SLFN5", ] %>%
  AverageExpression(assays = "RNA", group.by = c("patient", "gt_pc_2c", "celltypeL1")) %>%
  as.data.frame() %>%
  tidyr::pivot_longer(cols = tidyr::everything()) %>%
  tidyr::separate(name, into = c("patient", "condition", "genotype", "celltype"), sep = "_") %>%
  dplyr::filter(celltype %in% c("CD4.T", "CD8.T", "B", "NK", "merged.cMonos")) %>%
  dplyr::mutate(genotype = factor(genotype, levels = c("GG", "GA", "AA"))) %>%
  dplyr::mutate(celltype = dplyr::if_else(celltype == "merged.cMonos", "cMono", celltype)) %>%
  dplyr::mutate(celltype = factor(celltype, levels = c("cMono", "CD4.T", "CD8.T", "NK", "B"))) %>%
  dplyr::mutate(condition = dplyr::if_else(condition == "Hospitalized", "Hospitalized", condition)) %>%
  dplyr::mutate(patient = stringr::str_remove(patient, "RNA."))

qtl_effect <- pbmc_mhh50@meta.data %>%
  dplyr::select(patient, Age, gender, dosage_per_condition) %>%
  dplyr::mutate(patient = as.character(patient)) %>%
  dplyr::distinct() %>%
  dplyr::right_join(avg_exp_tab) %>%
  dplyr::group_by(celltype, condition) %>%
  dplyr::summarise(regression = {
    tab <- dplyr::cur_data()
    lm(value ~ dosage_per_condition + Age + gender, data = tab) %>%
      summary() %$% coefficients %>%
      as.data.frame() %>%
      tibble::rownames_to_column("param") %>%
      dplyr::rename(estimate = Estimate, std_error = `Std. Error`, t_value = `t value`, p_value = `Pr(>|t|)`) %>%
      dplyr::filter(param != "(Intercept)")
  }) %>%
  tidyr::unnest(regression)

bp <- (ggboxplot(
        avg_exp_tab, x = "condition", y = "value", fill = "genotype", size = 0.5, facet.by = "celltype",
        add = "ggscatter", palette = "lacet", width = 0.6, xlab = "Genotypes", ylab = "Avg. Expression: SLFN5") +
      geom_pwc(aes(group = genotype), p.adjust.method = "none", hide.ns = "p")
  ) %>%
  facet(facet.by = "celltype", nrow = 1)
file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/SLFN5_expression_per_genotype.box_plot.pdf") %>%
  ggsave(plot = bp, height = 5, width = 8)

eqtl_effect <- pbmc_mhh50@meta.data %>%
  dplyr::mutate(cell_barcode = rownames(.)) %>%
  dplyr::filter(celltypeL1 %in% c("merged.cMonos", "CD4.T", "CD8.T", "NK", "B")) %>%
  dplyr::group_by(celltypeL1, Severity) %>%
  dplyr::summarise(n = n()) %>%
  as.data.frame() %>%
  apply(1, function(vec) {
    per_cluster <- vec[1]
    per_condition <- vec[2]

    tar_cells <- pbmc_mhh50@meta.data %>% as.data.frame() %>% dplyr::filter(celltypeL1 == per_cluster, Severity == per_condition) %>% rownames()
    sub_pbmc <- pbmc_mhh50["SLFN5", tar_cells]

    meta_tab <- sub_pbmc@meta.data %>% as.data.frame()
    expr_tab <- sub_pbmc@assays$RNA@data %>% as.data.frame() %>% t()
    regr_tab <- cbind(meta_tab, expr_tab)

    fml <- formula(stringr::str_glue("SLFN5 ~ dosage_per_condition + (dosage_per_condition | patient) + gender + Age"))
    tryCatch(
      expr = summary(lmer(fml, regr_tab))$coefficients %>% as.data.frame() %>% dplyr::mutate(param = rownames(.)),
      error = function(e) {print(e); NULL}
    ) %>%
      dplyr::mutate(celltype = per_cluster, condition = per_condition)
  }) %>%
  Reduce(rbind, .) %>%
  dplyr::select(dplyr::all_of(c("beta" = "Estimate", "beta_se" = "Std. Error", "t_statistic" = "t value", "p_value" = "Pr(>|t|)", "param", "celltype", "condition"))) %>%
  dplyr::mutate(p_value_adj_fdr = p.adjust(p_value, method = "fdr"))

eqtl_effect %>% dplyr::filter(param %in% c("dosage_per_condition")) %>% as.data.table
eqtl_effect %>% fwrite(file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/SLFN5_eQTL_effect.csv"), quote = FALSE, na = "NA")


#
## Estimate the chromatin accessibility, have to activate the ArchR conda environment.
#
setwd(file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/"))
archr_proj_path <- file.path(proj_dir, "outputs/COVID_MHH50/PBMC_scATACseq")
proj <- ArchR::loadArchRProject(archr_proj_path)

if (FALSE) {
  proj <- addPeak2GeneLinks(proj, maxDist = 1000000)
  saveArchRProject(proj, archr_proj_path)
}

feature_list <- tribble(
  ~feature_id, ~chrom, ~start, ~end,
  # "CD55", "chr1", 207321507, 207321508,
  # "ADCY3", "chr2", 24919838, 24919839,
  "SLFN5", "chr17", 35243071, 35266360,
)

subproj <- getCellColData(proj) %>%
  as.data.frame() %>%
  dplyr::filter(Clusters2 %in% c("cMono", "CD4 T", "CD8 T", "NK", "B")) %>%
  rownames() %>%
  subsetCells(proj, cellNames = .)
new_group <- getCellColData(subproj)$cell.cond %>%
  stringr::str_replace("mild|severe", "Hospitalized") %>%
  stringr::str_replace("post", "Convalescent")
subproj <- addCellColData(subproj, data = new_group, name = "cell.condition_2", cells = rownames(getCellColData(subproj)))

# Subset the p2g links that linked to SLFN5 TSS
tar_feature <- "SLFN5"
p2g_links <- getPeak2GeneLinks(proj, resolution = 1, FDRCutOff = 0.05, corCutOff = 0.35)

set.seed(31415926)
# CD4T
p <- plotBrowserTrack(
  ArchRProj = subproj, groupBy = "cell.condition_2", useGroups = c("CD4 T.Hospitalized", "CD4 T.Convalescent"),
  geneSymbol = tar_feature, upstream = 3500, downstream = 3500, sizes = c(5, 0.5, 1.5, 1), loops = p2g_links
)
save_to <- paste0("Plot-Tracks-", tar_feature, "-with-Peak2GeneLinks.per_condition.CD4T.v3.pdf")
plotPDF(plotList = p, name = save_to, ArchRProj = proj, addDOC = FALSE, width = 6, height = 4)

# All cell types
p <- plotBrowserTrack(
  ArchRProj = subproj, groupBy = "cell.condition_2",
  geneSymbol = tar_feature, upstream = 3500, downstream = 3500, sizes = c(5, 0.5, 1.5, 1), loops = p2g_links
)
save_to <- paste0("Plot-Tracks-", tar_feature, "-with-Peak2GeneLinks.per_condition.v3.pdf")
plotPDF(plotList = p, name = save_to, ArchRProj = proj, addDOC = FALSE, width = 6, height = 4)


# Peaks
getPeakSet(proj) %>%
  (function(tab) { tab$cell_type <- names(tab); tab }) %>%
  as.data.table() %>%
  data.table::fwrite(file.path(proj_dir, "outputs/COVID_MHH50/PBMC_scATACseq/Plots/peaks.csv"))

# Peak annotations
getPeakAnnotation(proj) %>%
  (function(tab) { tab$cell_type <- names(tab); tab }) %>%
  as.data.table() %>%
  data.table::fwrite(file.path(proj_dir, "outputs/COVID_MHH50/PBMC_scATACseq/Plots/peaks.annotated.csv"))

# Peak to gene links
getPeak2GeneLinks(proj, resolution = 1, FDRCutOff = 0.05, corCutOff = 0.1) %>%
  (function(p2glist) p2glist$Peak2GeneLinks %>% as.data.frame()) %>%
  dplyr::rename(peakStart = "start", geneTSS = "end") %>%
  data.table::fwrite(file.path(proj_dir, "outputs/COVID_MHH50/PBMC_scATACseq/Plots/peak2gene_links.csv"))


# Peaks in SLFN5 locus
tar_range <- GRanges(seqnames = "chr17", ranges = IRanges(start = 35228071, end = 35266360))
save_to <- file.path(proj_dir, "outputs/COVID_MHH50/PBMC_scATACseq/Plots/peaks.annotated.SLFN5.csv")
if (!file.exists(save_to)) {
  all_peak_matches <- getMatches(proj)
  tar_peak_matches <- all_peak_matches[overlapsAny(all_peak_matches, tar_range) & rownames(all_peak_matches) %in% c("CD4 T", "CD8 T", "cMono", "NK", "B"), ]
  annotated_peaks <- tar_peak_matches@assays$data %>%
    as.data.frame() %>%
    dplyr::mutate(index = dplyr::row_number()) %>%
    tidyr::pivot_longer(-c(group, group_name, index)) %>%
    dplyr::filter(value) %>%
    dplyr::arrange(name, index) %>%
    dplyr::group_by(name) %>%
    dplyr::reframe(tp = as.data.frame(rowData(tar_peak_matches)[dplyr::cur_data()$index, ])) %>%
    tidyr::unnest(tp) %>%
    as.data.frame() %>%
    dplyr::group_by(idx) %>%
    dplyr::mutate(tf_motifs = paste(name, collapse = "|")) %>%
    dplyr::select(-name) %>%
    dplyr::distinct()

  data.table::fwrite(annotated_peaks, save_to)
}

# Coexpression between TFs and SLFN5
# The target TFs are those that prioritized using ASB database above.
tar_feature <- "SLFN5"

save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/SLFN5_TF.coexpression.all_conditions.single_cell.MHH50.csv")
if (!file.exists(save_to)) {
cor_test_tab <- pbmc_mhh50@meta.data %>%
  dplyr::filter(celltypeL1 %in% c("merged.cMonos", "CD4.T", "CD8.T", "NK", "B")) %>%
  dplyr::mutate(cell_barcode = rownames(.)) %>%
  dplyr::group_by(celltypeL1, Disease) %>%
  dplyr::summarise(n = n()) %>%
  as.data.frame() %>%
  apply(1, function(vec) {
    per_cluster <- vec[1]
    per_cond <- vec[2]
    cat(per_cluster, per_cond, "\n")

    tar_cells <- pbmc_mhh50@meta.data %>% as.data.frame() %>% dplyr::filter(celltypeL1 == per_cluster, Disease == per_cond) %>% rownames()
    sub_pbmc <- pbmc_mhh50[coexp_features, tar_cells]

    keta_tab <- sub_pbmc@meta.data %>% as.data.frame()
    expr_tab <- sub_pbmc@assays$RNA@data %>% as.data.frame() %>% t()
    regr_tab <- cbind(meta_tab, expr_tab)

    regr_res <- NULL
    for (per_tf in coexp_features) {
      if (per_tf == tar_feature) next

      fml <- formula(stringr::str_glue("{tar_feature} ~ {per_tf} + ({per_tf} | patient) + Age + gender"))
      regr_res <- tryCatch(
        expr = summary(lmer(fml, regr_tab))$coefficients %>% as.data.frame() %>% dplyr::mutate(param = rownames(.)) %>% rbind(regr_res),
        error = function(e) NULL
      )
    }
    rownames(regr_res) <- NULL
    dplyr::mutate(regr_res, celltype = per_cluster, disease = per_cond)
  }) %>%
  Reduce(rbind, .) %>%
  dplyr::select(dplyr::all_of(c("beta" = "Estimate", "beta_se" = "Std. Error", "t_statistic" = "t value", "p_value" = "Pr(>|t|)", "param", "celltype", "disease"))) %>%
  dplyr::filter(param %in% coexp_features) %>%
  dplyr::mutate(p_value_adj_fdr = p.adjust(p_value, method = "fdr"))

  fwrite(cor_test_tab, save_to)
} else {
  cor_test_tab <- fread(save_to)
}


if (!file.exists(save_to)) {
  cor_test_tab <- pbmc_mhh50@meta.data %>%
    dplyr::mutate(cell_barcode = rownames(.)) %>%
    dplyr::filter(celltypeL1 %in% c("merged.cMonos", "CD4.T", "CD8.T", "NK", "B")) %>%
    dplyr::group_by(celltypeL1, Disease) %>%
    dplyr::summarise(n = n()) %>%
    as.data.frame() %>%
    apply(1, function(vec) {
      per_cluster <- vec[1]
      per_condition <- vec[2]

      tar_cells <- pbmc_mhh50@meta.data %>% as.data.frame() %>% dplyr::filter(celltypeL1 == per_cluster, Disease == per_condition) %>% rownames()
      sub_pbmc <- pbmc_mhh50[coexp_features, tar_cells]

      meta_tab <- sub_pbmc@meta.data %>% as.data.frame()
      expr_tab <- sub_pbmc@assays$RNA@data %>% as.data.frame() %>% t()
      regr_tab <- cbind(meta_tab, expr_tab)

      regr_res <- NULL
      for (per_tf in coexp_features) {
        if (per_tf == "SLFN5") next
        if (!(per_tf %in% colnames(regr_tab))) next

        fml <- formula(stringr::str_glue("SLFN5 ~ {per_tf} + ({per_tf} | patient) + gender"))
        regr_res <- tryCatch(
          expr = summary(lmer(fml, regr_tab))$coefficients %>% as.data.frame() %>% dplyr::mutate(param = rownames(.)) %>% rbind(regr_res),
          error = function(e) {print(e); NULL}
        )
      }
      if (is.null(regr_res)) {
        return(NA)
      } else {
        rownames(regr_res) <- NULL
        dplyr::mutate(regr_res, celltype = per_cluster, condition = per_condition)
      }
    }) %>%
    Reduce(rbind, .) %>%
    dplyr::select(dplyr::all_of(c("beta" = "Estimate", "beta_se" = "Std. Error", "t_statistic" = "t value", "p_value" = "Pr(>|t|)", "param", "celltype", "condition"))) %>%
    dplyr::filter(param %in% coexp_features) %>%
    dplyr::mutate(p_value_adj_fdr = p.adjust(p_value, method = "fdr"))

  fwrite(cor_test_tab, save_to)
}

plot_tab <- cor_test_tab %>%
  dplyr::filter(celltype == "CD4.T") %>%
  dplyr::mutate(p_value_adj_fdr_log10 = -log10(p_value_adj_fdr), p_value_log10 = -log10(p_value), Correlation = dplyr::if_else(t_statistic > 0, "Pos", "Neg"), Correlation = factor(Correlation, levels = c("Neg", "Pos"))) %>%
  dplyr::arrange(-p_value_adj_fdr_log10)

p <- ggplot() +
  geom_point(data = dplyr::filter(plot_tab, p_value_adj_fdr <= 0.1), mapping = aes(x = t_statistic, y = p_value_log10, fill = Correlation), alpha = 0.75, size = 2.5, shape = 21) +
  geom_point(data = dplyr::filter(plot_tab, p_value_adj_fdr > 0.1), mapping = aes(x = t_statistic, y = p_value_log10, fill = Correlation), size = 1) +
  geom_label_repel(data = dplyr::filter(plot_tab, p_value_adj_fdr <= 0.1), mapping = aes(x = t_statistic, y = p_value_log10, label = param), alpha = .75, min.segment.length = 0, max.overlaps = 30) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  scale_fill_lancet() +
  theme_classic() +
  theme(legend.position = "top", legend.title = element_text(size = 12), axis.text = element_text(size = 12)) +
  labs(y = "-Log10(p-value)", x = "Z-score") +
  xlim(c(-6, 6)) +
  facet_wrap(~condition)

save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/SLFN5_TF.coexpression.CD4T.all_conditions.single_cell.MHH50.pdf")
ggsave(save_to, plot = p, width = 7, height = 5)


# CD55 labwork
inhibitor <- c("BL", "BL", "BL", "DB1976", "DB1976", "DB1976", "DB1976", "DB1976", "DB1976", "DB2313", "DB2313", "DB2313", "DB2313", "DB2313", "DB2313", "Ruxo", "Ruxo", "Ruxo")
dosage <- c(0, 0, 0, 0.5, 0.5, 0.5, 2, 2, 2, 0.5, 0.5, 0.5, 2, 2, 2, 0, 0, 0)
mock <- c(1.000, 1.000, 1.000, 1.330, 0.756, 0.677, 0.786, 0.920, 0.694, 1.345, 1.149, 0.797, 1.532, 1.053, 0.897, 2.913, 1.220, 1.052)
lps <- c(2.096, 7.743, 1.946, 4.884, 3.397, 1.597, 6.689, 4.103, 2.820, 3.053, 1.006, 1.877, 3.382, 2.553, 1.294, 4.709, 2.852, 2.331)

cd55_labwork_tab <- pzfx::read_pzfx(file.path(proj_dir, "inputs/Labwork/PU1 inhibitors Hua.pzfx"))
cd55_labwork_tab %>% tidyr::pivot_longer(-c(ROWTITLE)) %>% dplyr::pull(name)

plot_tab <- data.frame(inhibitor = inhibitor, dosage = dosage, Mock = mock, LPS = lps) %>%
  dplyr::mutate(dosage = factor(dosage, levels = c("0", "0.5", "2"))) %>%
  tidyr::pivot_longer(cols = c("Mock", "LPS"), names_to = "condition", values_to = "level") %>%
  dplyr::mutate(condition = factor(condition, c("Mock", "LPS"))) %>%
  dplyr::filter(inhibitor %in% c("BL", "DB2313"))

p <- plot_tab %>%
  ggplot(mapping = aes(x = dosage, y = level)) +
  geom_boxplot() +
  geom_quasirandom(aes(color = condition), size = 3, alpha = 0.5) +
  scale_color_jama() +
  facet_nested(~ condition + inhibitor, scale = "free_x", space="free", nest_line = element_line()) +
  labs(color = "Condition", x = "Inhibitor dosage", y = "Fold-change") +
  theme_classic() +
  theme(strip.background = element_blank())

ggsave(file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/CD55_labwork_v2.pdf"), plot = p, width = 4.5, height = 3.5)

res_tab <- file.path(proj_dir, "inputs/Labwork/PU1 Inhibitorn5_7 .xlsx")
cd55_rtpcr_tab <- excel_sheets(res_tab) %>%
  lapply(function(x, .res_tab) read_excel(.res_tab, sheet = x) %>% dplyr::mutate(batch = x, well_index = rep(1:24 + 1:24%%2, 3), repeat_index = rep(c(1, 2), as.integer(nrow(.)/2))), .res_tab = res_tab) %>%
  dplyr::bind_rows() %>%
  dplyr::select(sample_id = `Sample Name`, batch, well_index, repeat_index, `gene` = `Target Name`, `ddCT`) %>%
  dplyr::mutate(sample_id = dplyr::if_else(sample_id == "1976 low", "1976 low mock", sample_id)) %>%
  dplyr::mutate(sample_id = stringr::str_replace(sample_id, "2312", "2313")) %>%
  dplyr::mutate(treatment = dplyr::if_else(stringr::str_detect(sample_id, "LPS"), "LPS", "Mock")) %>%
  dplyr::mutate(treatment = factor(treatment, levels = c("Mock", "LPS"))) %>%
  dplyr::mutate(inhibitor = stringr::str_remove_all(sample_id, " (mock|high|low|LPS)")) %>%
  dplyr::mutate(inhibitor = factor(inhibitor, levels = c("DMSO", "1976", "2313", "Ruxo"))) %>%
  dplyr::mutate(dosage = dplyr::case_when(stringr::str_detect(sample_id, "low") ~ "0.5", stringr::str_detect(sample_id, "high") ~ "2", stringr::str_detect(sample_id, "Ruxo") ~ "1", TRUE ~ "0")) %>%
  dplyr::filter(gene %in% c("CD55")) %>%
  dplyr::group_by(sample_id, batch, gene, well_index, inhibitor, treatment) %>%
  dplyr::summarise(`log2(2-ddCt)` = mean(ddCT), dosage = head(dosage, 1))

p <- ggplot(cd55_rtpcr_tab, aes(x = dosage, y = `log2(2-ddCt)`)) +
  geom_boxplot() +
  geom_dots(side = "both", binwidth = unit(2, "mm")) +
  facet_nested(~treatment + inhibitor, scale = "free", space = "free", nest_line = element_line(linetype = 2)) +
  theme_classic() +
  theme(strip.background = element_blank(), ggh4x.facet.nestline = element_line(colour = "black"))
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/CD55_rtpcr.all_inhibitors.pdf")
ggsave(save_to, plot = p, width = 5.5, height = 3.5)

p <- cd55_rtpcr_tab %>% dplyr::filter(treatment == "Mock") %>%
  ggplot(aes(x = dosage, y = `log2(2-ddCt)`)) +
  geom_boxplot() +
  geom_dots(side = "both", binwidth = unit(2, "mm")) +
  facet_nested(~treatment + inhibitor, scale = "free_x", space = "free_x", nest_line = element_line(linetype = 2)) +
  theme_classic() +
  theme(strip.background = element_blank(), ggh4x.facet.nestline = element_line(colour = "black"))
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/CD55_rtpcr.all_inhibitors.mock.pdf")
ggsave(save_to, plot = p, width = 3, height = 3.5)

p <- cd55_rtpcr_tab %>% dplyr::filter(treatment == "LPS") %>%
  ggplot(aes(x = dosage, y = `log2(2-ddCt)`)) +
  geom_boxplot() +
  geom_dots(side = "both", binwidth = unit(2, "mm")) +
  facet_nested(~treatment + inhibitor, scale = "free_x", space = "free_x", nest_line = element_line(linetype = 2)) +
  theme_classic() +
  theme(strip.background = element_blank(), ggh4x.facet.nestline = element_line(colour = "black"))
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/CD55_rtpcr.all_inhibitors.LPS.pdf")
ggsave(save_to, plot = p, width = 3, height = 3.5)


p <- cd55_rtpcr_tab %>% dplyr::filter(treatment == "LPS", inhibitor %in% c("DMSO", "2313", "Ruxo")) %>%
  ggplot(aes(x = dosage, y = `log2(2-ddCt)`)) +
  geom_boxplot() +
  geom_dots(side = "both", binwidth = unit(2, "mm")) +
  facet_nested(~treatment + inhibitor, scale = "free_x", space = "free_x", nest_line = element_line(linetype = 2)) +
  theme_classic() +
  theme(strip.background = element_blank(), ggh4x.facet.nestline = element_line(colour = "black"))
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/CD55_rtpcr.2313_Ruxo.LPS.pdf")
ggsave(save_to, plot = p, width = 2.5, height = 3.5)

p <- cd55_rtpcr_tab %>% dplyr::filter(treatment == "Mock", inhibitor %in% c("DMSO", "2313", "Ruxo")) %>%
  ggplot(aes(x = dosage, y = `log2(2-ddCt)`)) +
  geom_boxplot() +
  geom_dots(side = "both", binwidth = unit(2, "mm")) +
  facet_nested(~treatment + inhibitor, scale = "free_x", space = "free_x", nest_line = element_line(linetype = 2)) +
  theme_classic() +
  theme(strip.background = element_blank(), ggh4x.facet.nestline = element_line(colour = "black"))
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/CD55_rtpcr.2313_Ruxo.mock.pdf")
ggsave(save_to, plot = p, width = 2.5, height = 3.5)

# Replication in long-covid cohort
eqtl_data <- fread(file.path(proj_dir, "outputs/pseudo_bulk/replication/longcovid_from_Qiuyao/slfn5_eqtl.txt")) %>%
  tidyr::pivot_longer(cols = dplyr::starts_with("AveExp_"), names_to = "Celltype", values_to = "AveExp") %>%
  dplyr::mutate(Celltype = stringr::str_remove_all(Celltype, "AveExp_")) %>%
  dplyr::select(IID = combo, Age, Gender, Timepoint, Stimulation, Disease = condition, Celltype, Genotype, Dosage = DS, SLFN5_AveExp = AveExp) %>%
  dplyr::mutate(Celltype = c("B" = "B", "cd4" = "CD4T", "cd8" = "CD8T", "mo" = "Monocytes", "nk" = "NK")[Celltype])

caqtl_data <- fread(file.path(proj_dir, "outputs/pseudo_bulk/replication/longcovid_from_Qiuyao/slfn5_caqtl.txt")) %>%
  dplyr::select(IID = combo, Timepoint = timepoint, Stimulation = stim, Celltype = celltype, Peak_AveCA = chr17.35241744.35244833)

plot_tab <- dplyr::left_join(eqtl_data, caqtl_data, by = c("IID", "Timepoint", "Stimulation", "Celltype")) %>%
  dplyr::filter(Disease %in% c("Severe/ICU", "Severe", "severe"))

set.seed(31415)
p_eqtl <- plot_tab %>%
  ggplot(mapping = aes(x = Genotype, y = SLFN5_AveExp)) +
  geom_boxplot(aes(fill = Genotype), outlier.color = "white", outlier.size = 0) +
  geom_point(color = "black", position = position_jitter(width = 0.15)) +
  facet_nested(~ Celltype, scale = "free_x", space="free", nest_line = element_line()) +
  labs(fill = "Genotype", x = NULL, y = "SLFN5 expression") +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text = element_text(size = 13), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 12))

p_caqtl <- plot_tab %>%
  ggplot(mapping = aes(x = Genotype, y = Peak_AveCA)) +
  geom_boxplot(aes(fill = Genotype), outlier.color = "white", outlier.size = 0) +
  geom_point(color = "black", position = position_jitter(width = 0.15)) +
  facet_nested(~ Celltype, scale = "free_x", space="free", nest_line = element_line()) +
  labs(fill = "Genotype", x = NULL, y = "Chrom. acc. at rs11080327") +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text = element_text(size = 13), strip.text = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 12), axis.ticks.x = element_blank())

p_ca2gex <- plot_tab %>%
  ggplot(mapping = aes(y = SLFN5_AveExp, x = Peak_AveCA, color = Genotype)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_nested(~ Celltype, scale = "free_x", nest_line = element_line()) +
  labs(fill = "Genotype", y = "SLFN5 Average expression", x = "Chrom. acc. at rs11080327") +
  theme_classic() +
  theme(strip.background = element_blank(), strip.text = element_blank(), legend.position = "none")

save_to <- file.path(proj_dir, "outputs/pseudo_bulk/replication/slfn5_eqtl_caqtl.severe_icu.pdf")
p <- p_eqtl / p_caqtl
ggsave(save_to, plot = p, width = 9, height = 5)



# ------- Code to be reviewed -------
if (FALSE) {
  new_meta <- pbmc_300bcg@meta.data %>%
    as.data.frame() %>%
    dplyr::mutate(cell_barcode = rownames(.)) %>%
    dplyr::left_join(gt_per_ind, by = "ids") %>%
    (function(tab) { rownames(tab) <- tab$cell_barcode; tab }) %>%
    dplyr::select(dplyr::starts_with("rs"))

  pbmc_300bcg <- AddMetaData(pbmc_300bcg, new_meta, colnames(new_meta))

  cond_tab <- tibble::tribble(
    ~celltype, ~time, ~stim,
    "Monocytes", "T0", "RPMI",
    "Monocytes", "T0", "LPS",
    "Monocytes", "T3m", "RPMI",
    "Monocytes", "T3m", "LPS",
  )

  fplist <- cond_tab %>%
    apply(1, function(vec) {
      per_celltype <- vec["celltype"]
      per_time <- vec["time"]
      per_stim <- vec["stim"]

      condition <- paste(per_time, per_stim)
      tar_cells <- pbmc_300bcg@meta.data %>% dplyr::filter(time == per_time, stim == per_stim) %>% rownames()
      count <- 0
      FeaturePlot(pbmc_300bcg, feature = c("CD55"), cells = tar_cells, split.by = "rs2564978_chr1_207321071_T_C", combine = FALSE) %>%
        purrr::map2(seq_along(.), function(pp, ii) {
          pp <- pp + theme_classic() + labs(title = NULL, x = NULL, y = NULL, subtitle = NULL)
          if (ii == 1) {
            pp + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y.right = element_blank(), axis.ticks.y.right = element_blank(), axis.line.y.right = element_blank())
          } else if (ii == 2) {
            pp + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
          } else { # ii is 3
            pp + theme(axis.title.y = element_blank(), axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank(), axis.line.y.left = element_blank())
          }
        }) %>%
        Reduce(`|`, .)
    }) %>%
    Reduce(`/`, .)

  vllist <- cond_tab %>%
    apply(1, function(vec) {
      per_celltype <- vec["celltype"]
      per_time <- vec["time"]
      per_stim <- vec["stim"]

      condition <- paste(per_time, per_stim)
      tar_cells <- pbmc_300bcg@meta.data %>% dplyr::filter(clusters1 == per_celltype, time == per_time, stim == per_stim) %>% rownames()
      VlnPlot(pbmc_300bcg[, tar_cells], cols = "white", feature = c("CD55"), pt.size = 0, split.by = "rs2564978_chr1_207321071_T_C") &
        NoAxes() & NoLegend() & theme(plot.margin = margin(0, 0, 0, 0, "mm")) & labs(title = NULL) & ylab(NULL)
    }) %>%
    Reduce(`/`, .)

  plot <- (vllist | fplist) + plot_layout(width = c(1, 4))
  save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/CD55_expression.per_genotype.pdf")
  ggsave(save_to, plot = plot, width = 10, height = 10)
}
