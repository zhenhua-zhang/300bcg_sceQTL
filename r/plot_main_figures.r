#!/usr/bin/env Rscript
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com, zhenhua.zhang3@helmholtz-hzi.de
# Created: 2022 Oct 12
# Updated: 2023 May 15

options(stringsAsFactors = FALSE, data.table.verbose = FALSE)
suppressPackageStartupMessages({
  # Functional annotation, enrichment, background databases, etc.
  library(msigdbr)
  library(ReactomePA) # reactome.db_1.81.0
  library(clusterProfiler)
  library(DOSE)
  library(enrichplot)
  library(ComplexHeatmap)
  library(circlize)
  library(SeuratObject)
  library(Seurat)
  library(Matrix)
  library(ArchR)
  library(lme4)
  library(lmerTest)
  library(GenomicRanges)

  # Misc, color, data manipulation, etc.
  library(RColorBrewer)
  library(magrittr)
  library(tidyverse)
  library(data.table)
  library(patchwork)
  library(ggbreak)
  library(ggforce)
  library(ggrepel)
  library(ggpubr)
  library(ggh4x)
  library(ggsci)
  library(gtools)
})


#
## Variables
#
proj_dir <- "~/Documents/projects/wp_bcg_eqtl"
time_vec <- c("T0", "T3m")
stim_vec <- c("RPMI", "LPS")
mode_vec <- c("normal", "interaction")
celltype_vec <- c("Monocytes", "CD4T", "CD8T", "NK", "B")
condition_vec <- c("T0_LPS.vs.T0_RPMI", "T3m_LPS.vs.T3m_RPMI", "T3m_RPMI.vs.T0_RPMI") #, "T3m_LPS.vs.T0_RPMI")
col_pattern_vec <- c("normal" = "_Common", "interaction" = ".vs.")
effect_vec <- c("LPS eff. (T0)", "LPS eff. (T3m)", "BCG eff.")
names(effect_vec) <- condition_vec


#
## Parameters
#
fdr_max <- 0.1
min_isize <- 10


#
## Plot single-cell UMAP, cell type markers, and cell proportion
#
# Define the number of colors you want
nb_cols <- 10
mycolors <- colorRampPalette(brewer.pal(8, "Blues"))(nb_cols)

proj_dir <- "~/Documents/projects/wp_bcg_eqtl"
so_path <- file.path(proj_dir, "inputs/sc_rnaseq/bcg4-0712.rds")
pbmc <- readRDS(so_path)

pbmc <- pbmc[, (!pbmc@meta.data$clusters1 %in% c("HSP(T)"))]

new_cell_lvl <- c("Monocytes", "CD4+ T", "CD8+ T", "NK", "B", "Platelet", "mDC", "pDC", "Undefined")
pbmc@meta.data$clusters1 <- factor(pbmc@meta.data$clusters1, levels = new_cell_lvl)
pbmc@meta.data$ts <- factor(pbmc@meta.data$ts, levels = c("T0_RPMI", "T0_LPS", "T3m_RPMI", "T3m_LPS"))
Idents(pbmc) <- "clusters1"
DefaultAssay(pbmc) <- "RNA"


#
## Draw the UMAP
#
g_umap <- DimPlot(pbmc, reduction = "umap") + NoAxes(keep.text = TRUE)
ggsave(file.path(proj_dir, "outputs/pseudo_bulk/overview/sc_rnaseq_umap.pdf"), g_umap, width = 7, height = 6)

# Figure Sx
g_umap_pcond <- DimPlot(pbmc, reduction = "umap", split.by = "ts", ncol = 2) + NoAxes(keep.text = TRUE)
ggsave(file.path(proj_dir, "outputs/pseudo_bulk/overview/sc_rnaseq_umap.per_condition.pdf"), g_umap_pcond, width = 8, height = 8)


#
## Draw feature map of cell type marker genes
#
DefaultAssay(pbmc) <- "RNA"
ct_markers <- c(
  "CTSS", "FCN1", "NEAT1", "LYZ", "PSAP", "S100A9", "AIF1", # "MNDA", "TYROBP", # Monocytes
  "IL7R", "MAL", "LTB", "LDHB", "TPT1", "TRAC", "TMSB10", # "CD3D", "CD4", "CD3G", # CD4 T
  "CD8B", "CD8A", "CD3D", "TMSB10", "HCST", "CD3G", "LINC02446", "CTSW", # "CD3E", "TRAC", # CD8 T
  "NKG7", "KLRD1", "TYROBP", "GNLY", "FCER1G", "PRF1", "CD247", "KLRF1", # "CST7", "GZMB", # NK
  "CD79A", "RALGPS2", "CD79B", "MS4A1", "BANK1", "CD74", "TNFRSF13C", "HLA-DQA1", # "IGHM", "MEF2C", # B
  "CD74", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "CCDC88A", "HLA-DRA", "HLA-DMA", "CST3" #, "HLA-DQB1", "HLA-DRB1" # DC
)
g_ct_marker <- FeaturePlot(pbmc, features = ct_markers, raster = TRUE) & NoLegend() & FontSize(main = 12) & NoAxes()
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/overview/sc.celltype_markers.featureplot.pdf")
ggsave(save_to, g_ct_marker)

DefaultAssay(pbmc) <- "integrated"
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/overview/sc.celltype_markers.dotplot.horizontal.pdf")
g_ct_marker <- DotPlot(pbmc, assay = "RNA", features = unique(ct_markers), cols = "RdBu",
                       idents = c("Monocytes", "CD4+ T", "CD8+ T", "NK", "B", "mDC", "pDC")) +
                RotatedAxis() + coord_flip()
ggsave(save_to, g_ct_marker, width = 12, height = 4)

save_to <- file.path(proj_dir, "outputs/pseudo_bulk/overview/sc.celltype_markers.dotplot.vertical.pdf")
g_ct_marker <- DotPlot(pbmc, assay = "RNA", features = unique(ct_markers), cols = "RdBu",
                       idents = c("Monocytes", "CD4+ T", "CD8+ T", "NK", "B", "mDC", "pDC")) +
                RotatedAxis() + coord_flip()
ggsave(save_to, g_ct_marker, width = 6, height = 10)


#
## Draw the cell-proportion boxplot, Figure 1x
#
cpp_tab <- pbmc@meta.data %>%
  dplyr::filter(status == "singlet") %>%
  dplyr::group_by(time, stim, ids) %>%
  dplyr::summarise(cpp = {
    n_pcond <- dplyr::n()
    n_pcell_pcond <- dplyr::cur_data()$clusters1 %>% table()
    data.frame(n_pcell_pcond / n_pcond)
  }) %>%
  as.data.table() %>%
  dplyr::rename("cell_type" = "cpp..", "cell_proportion" = "cpp.Freq") %>%
  dplyr::filter(!cell_type %in% c("Undefined", "Platelet")) %>%
  dplyr::mutate(
    stim = factor(stim, levels = c("RPMI", "LPS")),
    cell_type = factor(cell_type, levels = c("CD4+ T", "CD8+ T", "Monocytes", "NK", "B", "pDC", "mDC"))
  )

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

# ADCY3-rs11687089, DEAF1-rs10794331, GFM1-rs1367313, PPFIBP2-rs4758197, HAGH-rs344352
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
min_h4 <- 0.3
nr_indep_loci_per_gwas <- tibble::tribble(
  ~GWAS, ~nr_independent_loci, ~abbrev,
  # "AsthmaDT1", 13, "ASDT1",
  "CoronaryHeartDisease", 55, "CHD",
  "ObesityClass1", 24, "OB1",
  "ObesityClass2", 13, "OB2",
  "ObesityClass3", 3, "OB3",
  "ThyroidCancer", 741, "TC",
  "AtopicDermatitis", 19, "AD",
  "Height", 1044, "HT",
  "LungCancer", 13, "LC",
  "Schizophrenia", 152, "SZ",
  "Urate", 6, "UT",
  "CrohnsDisease", 91, "CrD",
  "InflammatoryBowelDisease", 109, "IBD",
  "SystemicLupusErythematosus", 145, "SLE",
  "UlcerativeColitis", 62, "UC",
  "AsthmaDT2", 42, "AS",
  "GoutDisease", 38, "GD",
  "OvarianCancer", 2, "OC",
  "Psoriasis", 116, "PR",
  "BodyMassIndex", 1457, "BMI",
  "DiastolicBloodPressure", 1350, "DBP",
  "RheumatoidArthritis", 50, "RA",
  "SystolicBloodPressure", 1279, "SBP",
  "Type2Diabetes", 169, "T2D",
  "MultipleSclerosis", 47, "MS",
  "COVID19Release4", 4, "COVID19_R4",
  "HDLCholesterol", 1198, "HDL",
  "LDLCholesterol", 587, "LDL",
  "Triglycerides", 1028, "TRI",
  "BladderCancer", 1, "BC",
  "ColorectalCancer", 9, "CC",
  "ProstateCancer", 66, "PC",
  "AlzheimerDiseases", 15, "AzD",
  "COVID19Release7", 52, "COVID19"#_R7"
)

# Number of colocalizations.
mode <- "normal"
coloc_files <- list.files(file.path(proj_dir, "outputs/pseudo_bulk/colocalization", mode), pattern = "*.csv", recursive = TRUE, full.names = TRUE)

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
      dplyr::mutate(Celltype = celltype, Condition = condition, outcome = as.character(outcome)) %>%
      dplyr::group_by(outcome, Celltype, Condition) %>%
      dplyr::filter(H4 >= min_h4) %>%
      dplyr::summarise(n_genes = n(), gene_list = paste(exposure, collapse = "|"), .groups = "keep")
  }) %>%
  Reduce(rbind, .) %>%
  dplyr::inner_join(nr_indep_loci_per_gwas, by = c("outcome" = "GWAS")) %>%
  dplyr::mutate(coloc_prop = n_genes / nr_independent_loci)

prop_mat <- mat %>%
  tidyr::pivot_wider(id_cols = abbrev, names_from = c(Celltype, Condition), values_from = coloc_prop) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(dplyr::across(-abbrev, ~ dplyr::if_else(is.na(.x), as.double(0), .x))) %>%
  dplyr::rename_with(~stringr::str_remove(.x, pattern = "_$")) %>%
  (function(m) {
    tmp <- m %>% dplyr::select(-abbrev) %>% as.matrix()
    rownames(tmp) <- m$abbrev
    tmp
  })

count_mat <- mat %>%
  tidyr::pivot_wider(id_cols = abbrev, names_from = c(Celltype, Condition), values_from = n_genes) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(dplyr::across(-abbrev, ~ dplyr::if_else(is.na(.x), as.integer(0), .x))) %>%
  dplyr::rename_with(~stringr::str_remove(.x, pattern = "_$")) %>%
  (function(m) {
    tmp <- m %>% dplyr::select(-abbrev) %>% as.matrix()
    rownames(tmp) <- m$abbrev
    tmp
  })


col_fun_prop <- colorRamp2(c(0, max(prop_mat)), c("gray95", "#2166AC"))
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/colocalization", paste0(mode, ".colocalization.per_celltype.pdf"))
pdf(save_to, width = 8, height = 3)
hm <- Heatmap(t(prop_mat),
  col = col_fun_prop, name = "Prop. of \ncolocs",
  top_annotation = columnAnnotation(Counts = anno_barplot(rowSums(count_mat)), show_annotation_name = FALSE),
  right_annotation = rowAnnotation(Counts = anno_barplot(colSums(count_mat)), show_annotation_name = FALSE),
  column_names_rot = 45, column_names_side = "top", column_dend_side = "bottom"
)
draw(hm)
dev.off()


#
## Example eQTL, G and GxE
#
work_dir <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL")
example_eqtl <- tibble::tribble(
  ~snp, ~feature,
  "rs11687089", "ADCY3", # Good common eQTL across different conditions.
  "rs2564978", "CD55", # Good example contrast pattern across cell types (i.e., rs2564978-C is associated higher expression in Monocytes but lower expression in CD4T/CD8T).
  "rs11080327", "SLFN5", # Also good ieQTL, but not cell-type specific. GETex (WB, 1.5e-49)
)

apply(example_eqtl, 1, function(e) {
  per_snp <- e["snp"]
  per_feature <- e["feature"]
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
        if (pg_1 != pg_2)
          new_level <- paste0(c(pg_1, pg_1, pg_2), c(pg_1, pg_2, pg_2))
      }
      tab[gt_col] <- factor(tab[, gt_col], levels = new_level)

      tab
    }) %>%
    as.data.frame()

  y_lab <- paste0("Avg. expression: ", per_feature)

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
  ggsave(save_to, plot = plot, width = 8, height = 4)

  # Per condition
  plot_tab <- dplyr::mutate(eqtl_tab, condition = paste0(time, "_", stim), condition = factor(condition, c("T0_RPMI", "T0_LPS", "T3m_RPMI", "T3m_LPS")))
  plot <- ggboxplot(plot_tab, x = gt_col, y = per_feature, fill = gt_col, palette = "lacet", xlab = FALSE, ylab = y_lab, facet.by = c("condition", "celltype")) +
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
        dplyr::mutate(comparison = pp, Effect = effect_vec[pp], Treatment = dplyr::case_when(stim == "LPS" | (stim == "RPMI" & time == "T3m" & Effect == "BCG eff.") ~ "Exposed", TRUE ~ "Baseline"))
    }) %>%
    Reduce(rbind, .) %>%
    dplyr::mutate(Treatment = factor(Treatment, levels = c("Baseline", "Exposed")))

  lapply(celltype_vec, function(per_celltype, .plot_tab) {
    plot <- dplyr::filter(.plot_tab, celltype == per_celltype) %>%
      dplyr::mutate(cross_comp = paste0(Treatment, !!rlang::sym(gt_col))) %>%
      (function(tab) {
        y_lab <- paste("Avg. expression:", per_feature, "in", per_celltype)
        plot_box <- ggboxplot(tab, x = gt_col, y = per_feature, fill = "Treatment", palette = "lancet", ylab = y_lab, facet.by = "Effect") %>%
          facet(facet.by = "Effect", nrow = 3) +
          theme(strip.background = element_blank(), strip.text = element_blank())

        plot_dot <- ggline(tab, x = "Treatment", y = per_feature, group = gt_col, color = gt_col, palette = "npg", add = c("mean_ci", "jitter"), facet.by = "Effect", ylab = FALSE) %>%
          facet(facet.by = "Effect", nrow = 3, strip.position = "right") +
          theme(axis.text.y.left = element_blank())

        plot_box & plot_dot
      })

    save_to <- file.path(work_dir, paste0("interaction.", per_feature, "-", per_snp, ".", per_celltype, ".box_and_dot.pdf"))
    ggsave(save_to, plot = plot, width = 6, height = 7)
  }, .plot_tab = plot_tab)

  NULL
})


#
## Per genotype expression of ADCY3/CD55/SLFN5
#
genotype_tab <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/tar_snps.genotypes.csv") %>%
  fread() %>%
  dplyr::rename_with(.fn = ~ stringr::str_remove_all(.x, "\\[[0-9]{1,2}\\]|:GT$|^# "))

gt_per_ind <- genotype_tab %>%
  tidyr::pivot_longer(dplyr::starts_with("300"), values_to = "GT_01", names_to = "ids") %>%
  dplyr::mutate(GT_RA = dplyr::case_when(
    GT_01 %in% c("0|0") ~ paste0(REF, REF),
    GT_01 %in% c("1|0", "0|1") ~ paste0(REF, ALT),
    GT_01 %in% c("1|1") ~ paste0(ALT, ALT)
  )) %>%
  dplyr::select(ID, CHROM, POS, REF, ALT, GT_RA, ids) %>%
  tidyr::pivot_wider(names_from = c(ID, CHROM, POS, REF, ALT), values_from = GT_RA)

new_meta <- pbmc@meta.data %>%
  as.data.frame() %>%
  dplyr::mutate(cell_barcode = rownames(.)) %>%
  dplyr::left_join(gt_per_ind, by = "ids") %>%
  (function(tab) { rownames(tab) <- tab$cell_barcode; tab }) %>%
  dplyr::select(dplyr::starts_with("rs"))

pbmc <- AddMetaData(pbmc, new_meta, colnames(new_meta))

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
    tar_cells <- pbmc@meta.data %>% dplyr::filter(time == per_time, stim == per_stim) %>% rownames()
    count <- 0
    FeaturePlot(pbmc, feature = c("CD55"), cells = tar_cells, split.by = "rs2564978_chr1_207321071_T_C", combine = FALSE) %>%
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
    tar_cells <- pbmc@meta.data %>% dplyr::filter(clusters1 == per_celltype, time == per_time, stim == per_stim) %>% rownames()
    VlnPlot(pbmc[, tar_cells], cols = "white", feature = c("CD55"), pt.size = 0, split.by = "rs2564978_chr1_207321071_T_C") &
      NoAxes() & NoLegend() & theme(plot.margin = margin(0, 0, 0, 0, "mm")) & labs(title = NULL) & ylab(NULL)
  }) %>%
  Reduce(`/`, .)

plot <- (vllist | fplist) + plot_layout(width = c(1, 4))
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/CD55_expression.per_genotype.pdf")
ggsave(save_to, plot = plot, width = 10, height = 10)

pbmc@meta.data$rs2564978_chr1_207321071_T_C <- factor(pbmc@meta.data$rs2564978_chr1_207321071_T_C, levels = c("TT", "TC", "CC"))
g_cd55_exp <- FeaturePlot(pbmc, feature = c("CD55"), split.by = "rs2564978_chr1_207321071_T_C") & NoLegend() & NoAxes()
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/CD55_expression.per_genotype.feature_plot.pdf")
ggsave(save_to, plot = g_cd55_exp, width = 7, height = 3)

g_cd55_exp <- VlnPlot(pbmc, feature = c("CD55"), pt.size = 0, idents = c("Monocytes", "CD4+ T", "CD8+ T", "NK", "B"), split.by = "rs2564978_chr1_207321071_T_C") & labs(x = NULL)
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/CD55_expression.per_genotype.violin_plot.pdf")
ggsave(save_to, plot = g_cd55_exp, width = 7, height = 3)

g_spi1_exp <- FeaturePlot(pbmc, feature = c("SPI1"), split.by = "rs2564978_chr1_207321071_T_C") & NoLegend() & NoAxes()
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/SPI1_expression.per_genotype.feature_plot.pdf")
ggsave(save_to, plot = g_spi1_exp, width = 7, height = 3)

g_spi1_exp <- VlnPlot(pbmc, feature = c("SPI1"), pt.size = 0, idents = c("Monocytes", "CD4+ T", "CD8+ T", "NK", "B"), split.by = "rs2564978_chr1_207321071_T_C") & labs(x = NULL)
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/SPI1_expression.per_genotype.violin_plot.pdf")
ggsave(save_to, plot = g_spi1_exp, width = 7, height = 3)

pbmc@meta.data$rs11687089_chr2_24860057_T_C <- factor(pbmc@meta.data$rs11687089_chr2_24860057_T_C, levels = c("TT", "TC", "CC"))
g_adcy3_exp <- FeaturePlot(pbmc, feature = c("ADCY3"), split.by = "rs11687089_chr2_24860057_T_C") & NoLegend() & NoAxes()
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/ADCY3_expression.per_genotype.feature_plot.pdf")
ggsave(save_to, plot = g_adcy3_exp, width = 7.5, height = 3)

g_adcy3_exp <- VlnPlot(pbmc, feature = c("ADCY3"), pt.size = 0, idents = c("Monocytes", "CD4+ T", "CD8+ T", "NK", "B"), split.by = "rs11687089_chr2_24860057_T_C")
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/ADCY3_expression.per_genotype.violin_plot.pdf")
ggsave(save_to, plot = g_adcy3_exp, width = 7.5, height = 3)

pbmc@meta.data$rs11080327_chr17_35244527_G_A <- factor(pbmc@meta.data$rs11080327_chr17_35244527_G_A, levels = c("GG", "GA", "AA"))
g_slfn5_exp <- FeaturePlot(pbmc, feature = c("SLFN5"), split.by = "rs11080327_chr17_35244527_G_A") & NoLegend() & NoAxes()
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/SLFN5_expression.per_genotype.feature_plot.pdf")
ggsave(save_to, plot = g_slfn5_exp, width = 7.5, height = 3)

g_slfn5_exp <- VlnPlot(pbmc, feature = c("SLFN5"), pt.size = 0, idents = c("Monocytes", "CD4+ T", "CD8+ T", "NK", "B"), split.by = "rs11080327_chr17_35244527_G_A")
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/SLFN5_expression.per_genotype.violin_plot.pdf")
ggsave(save_to, plot = g_slfn5_exp, width = 7.5, height = 3)


#
##  Correlation between BMI and ADCY3
#
work_dir <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL")
meta_info <- file.path(proj_dir, "inputs/phenotypes/meta_information.with_bmi.csv") %>% data.table::fread()

example_qtl <- tibble::tribble(
  ~feature, ~snp, ~celltype, ~dw, ~up,
  "ADCY3", "rs11687089", "Monocytes", NA, NA,
  "ADCY3", "rs11687089", "CD4T", NA, NA,
  "ADCY3", "rs11687089", "CD8T", NA, NA,
  "ADCY3", "rs11687089", "NK", NA, NA,
  "ADCY3", "rs11687089", "B", NA, NA
)

for (per_snp in c("rs11687089")) {
  for (per_feature in c("ADCY3")) {
    per_infile <- file.path(work_dir, paste0(per_feature, "-", per_snp, ".QTL_info.csv"))
    qtl_info_tab <- data.table::fread(per_infile) %>%
      dplyr::mutate(
        time = dplyr::if_else(time == 0, "T0", "T3m"),
        time = factor(time, levels = time_vec),
        stim = dplyr::if_else(stim == 0, "RPMI", "LPS"),
        stim = factor(stim, levels = stim_vec),
        donor_id = stringr::str_extract(V1, "300BCG[0-9]+"), celltype = factor(celltype, celltype_vec)
      ) %>%
      dplyr::arrange(dplyr::desc(!!rlang::sym(paste0(per_snp, "_GT")))) %>%
      dplyr::left_join(meta_info, by = c("donor_id" = "DonorID")) %>%
      dplyr::mutate(gender = factor(gender, levels = c(0, 1)))

    m <- lmer(BMI ~ rs11687089_GT + ADCY3 + (ADCY3 | celltype) + (ADCY3 | donor_id) + age + gender, qtl_info_tab)
    m <- lmer(BMI ~ rs11687089_DS + ADCY3 + (ADCY3 | celltype) + (ADCY3 | donor_id) + age + gender, qtl_info_tab)

    spearman_rho_tab <- expand.grid(c("T0", "T3m"), c("RPMI", "LPS"), celltype_vec) %>%
      as.data.frame() %>%
      apply(1, function(line, .qtab, .tar_feature) {
        per_time <- line["Var1"]
        per_stim <- line["Var2"]
        per_celltype <- line["Var3"]

        per_tab <- dplyr::filter(.qtab, celltype == per_celltype, time == per_time, stim == per_stim) %>%
          as.data.frame()

        ctres <- cor.test(per_tab[, "BMI"], per_tab[, .tar_feature], method = "spearman")
        m <- lm(BMI ~ ADCY3 + age + gender, per_tab)
        cat(paste(per_stim, per_time, per_celltype, "\n"))
        print(summary(m))
        data.frame(celltype = per_celltype, time = per_time, stim = per_stim, p_value = ctres$p.value, spearman_rho = ctres$estimate)
      }, .qtab = qtl_info_tab, .tar_feature = per_feature) %>%
      Reduce(rbind, .) %>%
      dplyr::mutate(feature = per_feature, snp = per_snp)
  }
}


tmp <- apply(example_qtl, 1, function(line, .mtab) {
  per_snp <- line["snp"]
  per_feature <- line["feature"]
  per_celltype <- line["celltype"]
  per_dw <- as.numeric(line["dw"])
  per_up <- as.numeric(line["up"])

  per_infile <- file.path(work_dir, paste0(per_feature, "-", per_snp, ".QTL_info.csv"))
  qtl_info_tab <- data.table::fread(per_infile) %>%
    dplyr::mutate(
      time = dplyr::if_else(time == 0, "T0", "T3m"),
      time = factor(time, levels = time_vec),
      stim = dplyr::if_else(stim == 0, "RPMI", "LPS"),
      stim = factor(stim, levels = stim_vec),
      donor_id = stringr::str_extract(V1, "300BCG[0-9]+"), celltype = factor(celltype, celltype_vec)
    ) %>%
    dplyr::arrange(dplyr::desc(!!rlang::sym(per_snp))) %>%
    dplyr::left_join(.mtab, by = "donor_id") %>%
    as.data.frame()

  for (per_time in c("T0", "T3m")) {
    for (per_stim in c("RPMI", "LPS")) {
      work_tab <- qtl_info_tab %>%
        dplyr::select(time, stim, celltype, donor_id, Age, Sex, one_of(per_snp, per_feature, "BMI")) %>%
        dplyr::filter(celltype == per_celltype, time == per_time, stim == per_stim)

      ctt <- cor.test(work_tab[, "BMI"], work_tab[, per_feature], method = "spearman")
      data.frame(p_value = ctt$p.value, spearman_rho = ctt$estimate)
    }
  }

  ctt
}, .mtab = meta_info)


#
## Allele-specific TF binding affinities
#
astfb_tab <- data.table::fread(file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/ASTFB/rs2564978.astfb.txt"))
tar_bcg_features <- rownames(pbmc)

tar_astfb <- astfb_tab %>%
  dplyr::filter(!is.na(motif_log_pref), !is.na(motif_log_palt), TF %in% tar_bcg_features) %>%
  dplyr::mutate(asb_fc = -log2(fdrp_bh_alt / fdrp_bh_ref)) %>%
  dplyr::arrange(asb_fc) %>%
  dplyr::mutate(
    asb_sig = dplyr::if_else(fdrp_bh_alt < 0.05 | fdrp_bh_ref < 0.05, "Yes", "No"),
    asb_sig = factor(asb_sig, levels = c("Yes", "No")),
    motif_sig = dplyr::if_else(motif_log_palt > -log10(0.05) | motif_log_pref > -log10(0.05), "Yes", "No"),
    motif_sig = factor(motif_sig, levels = c("Yes", "No")),
    TF = factor(TF, .$TF)
  )
tf_order <- tar_astfb$TF

asb_bar_plot <- ggplot(tar_astfb, aes(x = asb_fc, y = TF, fill = asb_sig)) +
  geom_bar(stat = "identity") +
  scale_fill_npg() +
  labs(y = NULL, x = "Allelic binding FC") +
  theme_classic() +
  theme(
    axis.ticks.y.left = element_blank(), axis.text.y.left = element_blank(),
    axis.line.y.left = element_blank(), legend.position = "none"
  )

motif_bar_plot <- ggplot(tar_astfb, aes(x = motif_fc, y = TF, fill = motif_sig)) +
  geom_bar(stat = "identity") +
  scale_fill_npg(name = "Signif.") +
  labs(y = NULL, x = "Motif FC") +
  theme_classic() +
  theme(
    axis.ticks.y.left = element_blank(), axis.text.y.left = element_blank(),
    axis.line.y.left = element_blank(), legend.position = "top"
  )

motif_pos <- ggplot(tar_astfb, aes(x = motif_pos - 25, y = TF, xend = motif_pos, yend = TF, color = motif_orient)) +
  geom_segment(arrow = arrow(length = unit(0.00, "cm")), linewidth = 3) +
  geom_vline(xintercept = 0, color = "black", linewidth = 1) +
  scale_color_jama(name = "Motif strand") +
  labs(x = "Genome coord.", y = NULL) +
  theme_classic() +
  theme(
    axis.ticks.x.bottom = element_blank(), axis.text.x.bottom = element_blank(),
    legend.position = "top"
  )

bar_plots <- motif_pos + asb_bar_plot + motif_bar_plot
save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/ASTFB/rs2564978.pdf")
ggsave(save_to, width = 6, height = 3.25)


#
## Co-expression of TF and CD55, pseudo-bulk
#
group_idx <- c("clusters1", "time", "stim", "ids")
coexp_features <- c(as.character(tar_astfb$TF), "CD55")

if (use_pseudo_bulk) {
  meta_tab <- pbmc@meta.data %>%
    dplyr::select(ids, age, gender, time, stim) %>%
    dplyr::distinct() %>%
    as.data.table()
  coexp_tab <- pbmc[coexp_features, ] %>%
    AverageExpression(assays = "RNA", group.by = group_idx) %>%
    as.data.frame() %>%
    dplyr::mutate(feature_id = rownames(.)) %>%
    dplyr::rename_with(starts_with("RNA."), .fn = ~ str_remove_all(.x, "RNA|\\.")) %>%
    tidyr::pivot_longer(-feature_id, names_to = "samples", values_to = "average_expression") %>%
    tidyr::separate(samples, into = c("celltype", "time", "stim", "ids"), sep = "_")
  merged_tab <- dplyr::inner_join(meta_tab, coexp_tab, by = c("ids", "time", "stim"))

  cor_test_tab <- merged_tab %>%
    dplyr::select(celltype, time, stim) %>%
    dplyr::distinct() %>%
    dplyr::filter(celltype %in% celltype_vec) %>%
    apply(1, function(vec) {
      per_celltype <- vec[1]
      per_time <- vec[2]
      per_stim <- vec[3]

      per_tmp_tab <- merged_tab %>%
        dplyr::filter(celltype == per_celltype, time == per_time, stim == per_stim) %>%
        tidyr::pivot_wider(names_from = feature_id, values_from = average_expression)

      res_tab <- NULL
      for (per_tf in coexp_features) {
        if (per_tf == "CD55") next
        fml <- stringr::str_glue("CD55 ~ {per_tf} + age + gender")
        res_tab <- summary(lm(fml, per_tmp_tab))$coefficients %>% as.data.frame() %>% dplyr::mutate(param = rownames(.)) %>% rbind(res_tab)
      }
      rownames(res_tab) <- NULL
      dplyr::mutate(res_tab, celltype = per_celltype, time = per_time, stim = per_stim)
    }) %>%
    Reduce(rbind, .) %>%
    dplyr::select(dplyr::all_of(c("beta" = "Estimate", "beta_se" = "Std. Error", "t_statistic" = "t value", "p_value" = "Pr(>|t|)", "param", "celltype", "time", "stim"))) %>%
    dplyr::filter(param %in% coexp_features) %>%
    dplyr::mutate(p_value_adj_fdr = p.adjust(p_value, method = "fdr"))

  p <- cor_test_tab %>%
    dplyr::filter(celltype == "Monocytes") %>%
    dplyr::mutate(
      p_value_label = dplyr::if_else(p_value_adj_fdr < 0.05 & p_value_adj_fdr > 0, base::format(p_value_adj_fdr, digit = 1, scientific = FALSE), ""),
      p_value_adj_fdr_log10 = -log10(p_value_adj_fdr),
      condition = paste(time, stim),
      condition = factor(condition, levels = c("T0 RPMI", "T0 LPS", "T3m RPMI", "T3m LPS"))
    ) %>%
    ggplot() +
    geom_point(aes(x = celltype, y = param, color = beta, size = p_value_adj_fdr_log10)) +
    scale_size(name = "-log10(p_adj_fdr)", range = c(1, 12)) +
    scale_color_gradient2(name = "Spearman's Rho", low = "#3C5488FF", mid = "gray95", high = "#DC0000FF") +
    theme_bw() +
    theme(
      legend.position = "top", legend.title = element_text(size = 14),
      axis.line.x.top = element_blank(), axis.line.y.left = element_blank(),
      axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45, size = 14), axis.text.y = element_text(size = 14),
      strip.background = element_blank(), strip.text = element_text(size = 15)
    ) +
    labs(x = NULL, y = NULL) +
    facet_wrap(~ condition, nrow = 1)

  save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/CD55_TF.coexpression.Monocytes.all_conditions.single_cell.pdf")
  ggsave(save_to, plot = p, width = 7, height = 5.5)
} else {
  cor_test_tab <- pbmc@meta.data %>%
    dplyr::mutate(cell_barcode = rownames(.)) %>%
    dplyr::group_by(clusters1, time, stim) %>%
    dplyr::summarise(n = n()) %>%
    as.data.frame() %>%
    apply(1, function(vec) {
      per_cluster <- vec[1]
      per_time <- vec[2]
      per_stim <- vec[3]

      tar_cells <- pbmc@meta.data %>% as.data.frame() %>% dplyr::filter(clusters1 == per_cluster, time == per_time, stim == per_stim) %>% rownames()
      sub_pbmc <- pbmc[coexp_features, tar_cells]

      meta_tab <- sub_pbmc@meta.data %>% as.data.frame()
      expr_tab <- sub_pbmc@assays$RNA@data %>% as.data.frame() %>% t()
      regr_tab <- cbind(meta_tab, expr_tab)

      regr_res <- NULL
      for (per_tf in coexp_features) {
        if (per_tf == "CD55") next

        fml <- formula(stringr::str_glue("CD55 ~ {per_tf} + ({per_tf} | ids) + age + gender"))
        regr_res <- tryCatch(
          expr = summary(lmer(fml, regr_tab))$coefficients %>% as.data.frame() %>% dplyr::mutate(param = rownames(.)) %>% rbind(regr_res),
          error = function(e) NULL
        )
      }
      rownames(regr_res) <- NULL
      dplyr::mutate(regr_res, celltype = per_cluster, time = per_time, stim = per_stim)
    }) %>%
    Reduce(rbind, .) %>%
    dplyr::select(dplyr::all_of(c("beta" = "Estimate", "beta_se" = "Std. Error", "t_statistic" = "t value", "p_value" = "Pr(>|t|)", "param", "celltype", "time", "stim"))) %>%
    dplyr::filter(param %in% coexp_features) %>%
    dplyr::mutate(p_value_adj_fdr = p.adjust(p_value, method = "fdr"))

  save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/CD55_TF.coexpression.all_conditions.single_cell.csv")
  fwrite(cor_test_tab, save_to)
  p <- cor_test_tab %>%
    dplyr::filter(celltype == "Monocytes") %>%
    dplyr::mutate(
      p_value_label = dplyr::if_else(p_value_adj_fdr < 0.05 & p_value_adj_fdr > 0, "#", ""),
      p_value_adj_fdr_log10 = -log10(p_value_adj_fdr),
      p_value_log10 = -log10(p_value),
      condition = paste(time, stim),
      condition = factor(condition, levels = c("T0 RPMI", "T0 LPS", "T3m RPMI", "T3m LPS")),
      param = factor(param, tf_order)
    ) %>%
    ggplot() +
    geom_point(aes(x = condition, y = param, color = t_statistic, size = p_value_log10)) +
    geom_text(aes(x = condition, y = param, label = p_value_label), size = 5) +
    scale_size(name = "-log10(p_adj)", range = c(1, 12)) +
    scale_color_gradient2(name = "Z-score", low = "#3C5488FF", mid = "gray95", high = "#DC0000FF") +
    theme_bw() +
    theme(
      legend.position = "right", legend.title = element_text(size = 12),
      axis.line.x.top = element_blank(), axis.line.y.left = element_blank(),
      axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45, size = 12), axis.text.y = element_text(size = 12),
    ) +
    labs(x = NULL, y = NULL)

  save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/CD55_TF.coexpression.Monocytes.all_conditions.single_cell.pdf")
  ggsave(save_to, plot = p, width = 5, height = 4.5)
}


#
## Estimate the gene expression of target genes in COVID-19 context
#
pbmc_file <- "/vol/projects/BIIM/Covid_50MHH/NubesSubmission/scRNA/PBMC_scRNAseq.rds"
mhh50 <- readRDS(pbmc_file)

map_tab <- fread("/vol/projects/zzhang/projects/wp_covid19_mhh50/inputs/idmapping/id_mapping.txt")
genotype_tab <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/rs11080327_mhh50_genotype.csv") %>%
  fread() %>%
  dplyr::rename_with(.fn = ~ stringr::str_remove_all(.x, "\\[[0-9]{1,2}\\]|:GT$|^# "))
gt_per_ind <- genotype_tab %>%
  tidyr::pivot_longer(-c(CHROM, POS, ID, REF, ALT), values_to = "GT_01", names_to = "ids") %>%
  dplyr::mutate(GT_RA = dplyr::case_when(
    GT_01 %in% c("0|0") ~ paste0(REF, REF),
    GT_01 %in% c("1|0", "0|1") ~ paste0(REF, ALT),
    GT_01 %in% c("1|1") ~ paste0(ALT, ALT)
  )) %>%
  dplyr::select(ID, CHROM, POS, REF, ALT, GT_RA, ids) %>%
  tidyr::pivot_wider(names_from = c(ID, CHROM, POS, REF, ALT), values_from = GT_RA) %>%
  dplyr::left_join(map_tab, by = c("ids" = "genoID")) %>%
  dplyr::select(patientID, dplyr::starts_with("rs"))

new_meta_data <- mhh50@meta.data %>%
  dplyr::left_join(gt_per_ind, by = c("patient" = "patientID")) %>%
  dplyr::mutate(genotype_per_condition = paste(Disease, rs11080327_17_35244527_G_A, sep = "_")) %>%
  dplyr::pull(genotype_per_condition, cellbarcodes)

mhh50 <- AddMetaData(mhh50, new_meta_data, "genotype_per_condition")

mhh50_deg <- NULL
for (pct in unique(mhh50$celltypeL1)) {
  if (pct %in% c("Plasmablast")) next

  mhh50_subset <- mhh50["SLFN5", mhh50$celltypeL1 == pct]
  Idents(mhh50_subset) <- "Disease"

  mhh50_deg <- FindMarkers(mhh50_subset, ident.1 = "Active", ident.2 = "Convalescent", logfc.threshold = 0.1) %>%
    dplyr::mutate(Gene = rownames(.), Comparison = "hospitalized.vs.convalescent", Cohort = "MHH50", Celltype = pct) %>%
    rbind(mhh50_deg)
}
mhh50_deg %>%
  dplyr::mutate(p_val_adj = p.adjust(p_val)) %>%
  data.table::fwrite(file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/SLFN5.differential_expression.MHH50.csv"))

fp <- FeaturePlot(mhh50, feature = "SLFN5", split = "genotype_per_condition")
ggsave(file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/SLFN5_expression_per_genotype.feature_plot.pdf"), plot = fp)
vp <- VlnPlot(mhh50, feature = "SLFN5", idents = "CD4.T", split.by = "genotype_per_condition")
ggsave(file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/SLFN5_expression_per_genotype.violin_plot.pdf"), plot = vp)

avg_exp_tab <- mhh50["SLFN5", ] %>%
  AverageExpression(assays = "RNA", group.by = c("patient", "genotype_per_condition", "celltypeL1")) %>%
  as.data.frame() %>%
  tidyr::pivot_longer(cols = tidyr::everything()) %>%
  tidyr::separate(name, into = c("patient", "condition", "genotype", "celltype"), sep = "_") %>%
  dplyr::filter(celltype %in% c("CD4.T", "CD8.T", "B", "NK", "merged.cMonos")) %>%
  dplyr::mutate(genotype = factor(genotype, levels = c("GG", "GA", "AA"))) %>%
  dplyr::mutate(celltype = dplyr::if_else(celltype == "merged.cMonos", "cMono", celltype)) %>%
  dplyr::mutate(celltype = factor(celltype, levels = c("cMono", "CD4.T", "CD8.T", "NK", "B"))) %>%
  dplyr::mutate(condition = dplyr::if_else(condition == "Active", "Hospitalized", condition))

bp <- ggboxplot(avg_exp_tab, x = "genotype", y = "value", fill = "condition", size = 0.5,
                 facet.by = "celltype", add = "ggscatter", palette = "lacet", width = 0.6, xlab = "Genotypes",
                 ylab = "Avg. Expression: SLFN5") %>%
  facet(facet.by = "celltype", nrow = 1)
ggsave(file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/SLFN5_expression_per_genotype.box_plot.pdf"), plot = bp, height = 5)


#
## Estimate the chromatin accessibility
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
  # "KANSL1", "chr17", 46029916, 46225389,
  "SLFN5", "chr17", 35243071, 35266360,
  # "ADCY3", "chr2", 24919838, 24919839,
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
p2g_links <- getPeak2GeneLinks(proj, resolution = 1, FDRCutOff = 0.05, corCutOff = 0.1)

p <- plotBrowserTrack(ArchRProj = subproj, groupBy = "cell.condition_2", geneSymbol = tar_feature, upstream = 4000, sizes = c(7, 0.5, 2, 1), ylim = c(0, 0.999), downstream = 4000, loops = p2g_links, facetbaseSize = 11, baseSize = 11)
save_to <- paste0("Plot-Tracks-", tar_feature, "-with-Peak2GeneLinks.per_condition.v3.pdf")
plotPDF(plotList = p, name = save_to, ArchRProj = proj, addDOC = FALSE, width = 8, height = 7)

p <- plotBrowserTrack(ArchRProj = subproj, groupBy = "Clusters2", geneSymbol = tar_feature, upstream = -150000, sizes = c(3.5, 0.35, 1.45, 0.55), ylim = c(0, 0.9875), downstream = 150000, loops = p2g_links)
save_to <- paste0("Plot-Tracks-", tar_feature, "-with-Peak2GeneLinks.per_celltype.v3.pdf")
plotPDF(plotList = p, name = save_to, ArchRProj = proj, addDOC = FALSE, width = 8, height = 4)


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

data.table::fwrite(annotated_peaks, file.path(proj_dir, "outputs/COVID_MHH50/PBMC_scATACseq/Plots/peaks.annotated.SLFN5.csv"))

# Coexpression between TFs and SLFN5
coexp_features <- dplyr::pull(annotated_peaks, tf_motifs) %>%
  stringr::str_split(pattern = "\\|", simplify = TRUE) %>%
  stringr::str_extract(pattern = "([0-9A-Z]+)") %>%
  purrr::discard(~is.na(.x)) %>%
  unique() %>%
  c("SLFN5")

cor_test_tab <- mhh50@meta.data %>%
  dplyr::mutate(cell_barcode = rownames(.)) %>%
  dplyr::filter(celltypeL1 %in% c("merged.cMonos", "CD4.T", "CD8.T", "NK", "B")) %>%
  dplyr::group_by(celltypeL1, Disease) %>%
  dplyr::summarise(n = n()) %>%
  as.data.frame() %>%
  apply(1, function(vec) {
    per_cluster <- vec[1]
    per_condition <- vec[2]

    tar_cells <- mhh50@meta.data %>% as.data.frame() %>% dplyr::filter(celltypeL1 == per_cluster, Disease == per_condition) %>% rownames()
    sub_pbmc <- mhh50[coexp_features, tar_cells]

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

save_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/SLFN5_TF.coexpression.all_conditions.single_cell.MHH50.csv")
fwrite(cor_test_tab, save_to)

plot_tab <- cor_test_tab %>%
  dplyr::filter(celltype == "CD4.T") %>%
  dplyr::mutate(
    p_value_adj_fdr_log10 = -log10(p_value_adj_fdr),
    p_value_log10 = -log10(p_value),
    Correlation = dplyr::if_else(t_statistic > 0, "Pos", "Neg"),
    Correlation = factor(Correlation, levels = c("Neg", "Pos")),
    condition = dplyr::if_else(condition == "Active", "Hospitalized", "Convalescent")
  ) %>%
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
