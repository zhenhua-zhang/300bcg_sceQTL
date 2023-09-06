#!/usr/bin/env Rscript
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang3@helmohltz-hzi.de, zhenhua.zhang217@gmail.com
# Created: Oct 22, 2022
# Updated: Jul 07, 2023

#' Code to generate a track plot by
options(stringsAsFactors = FALSE, datatable.showProgress = FALSE, data.table.verbose = FALSE, Gviz.ucscUrl =  "http://genome-euro.ucsc.edu/cgi-bin/")
suppressPackageStartupMessages({
  library(Gviz)
  library(biomaRt)
  library(magrittr)
  library(tidyverse)
  library(data.table)
  library(AnnotationHub)
  library(GenomicRanges)
  library(GenomicInteractions)
})


proj_dir <- "~/Documents/projects/wp_bcg_eqtl"

# Abbreviation of conditions
condabbr <- c("T0_LPS.vs.T0_RPMI" = "L0-R0", "T3m_LPS.vs.T3m_RPMI" = "L3-R3", "T3m_RPMI.vs.T0_RPMI" = "R3-R0")

# The genome coordination version
genome <- "hg38"

# eQTL summary statistics
base_dir <- file.path(proj_dir, "outputs/pseudo_bulk")
marker_snp <- "rs62063281"; marker_feature <- "ADCY3"
marker_snp <- "rs2564978"; marker_feature <- "CD55"
marker_snp <- "rs11080327"; marker_feature <- "SLFN5"
marker_snp <- "rs883416"; marker_feature <- "SLFN5"; marker_snp_pos <- 35243422

sumstat_query <- tibble::tribble(
  ~feature_id, ~model, ~cell_type, ~condition, ~gwas_id, ~start, ~end, ~skip_atac, ~atac_order,
  "CD55", "normal", "Monocytes", "", "", 0, 0, FALSE, 2,
  # "CD55", "normal", "CD4T", "", "", 0, 0, FALSE, 1,
  # "CD55", "normal", "CD8T", "", "", 0, 0, FALSE, 3,
  # "CD55", "normal", "NK", "", "", 0, 0, FALSE, 4,
  "CD55", "normal", "B", "", "", 0, 0, FALSE, 5,
  # "SLFN5", "interaction", "CD8T", "T0_LPS.vs.T0_RPMI", "", 0, 0, FALSE,
  # "SLFN5", "interaction", "CD4T", "T0_LPS.vs.T0_RPMI", "", 0, 0, FALSE,
  # "SLFN5", "interaction", "CD4T", "T0_LPS.vs.T0_RPMI", "COVID19PlosOne", 0, 0, FALSE,
  # "SLFN5", "interaction", "CD4T", "T0_LPS.vs.T0_RPMI", "COVID19Release7", 0, 0, FALSE,
) %>%
dplyr::filter(feature_id == marker_feature)

sumstat <- apply(sumstat_query, 1, function(vec, .base_dir, .marker_snp) {
  per_feature_id <- vec["feature_id"]
  per_cell_type <- vec["cell_type"]
  per_condition <- vec["condition"]
  per_gwas_id <- vec["gwas_id"]
  per_model <- vec["model"]

  if (per_gwas_id == "") {
    tar_cols <- c(feature_id = "feature_id", snp_id = "snp_id", chrom = "snp_chromosome", pos = "snp_position", log10_pval = "p_value")
    tab <- file.path(.base_dir, "summary_statistic", per_model, per_cell_type, per_condition, "qtl_results_all.txt") %>%
      data.table::fread() %>%
      dplyr::select(feature_id, snp_id, snp_chromosome, snp_position, p_value) %>%
      dplyr::filter(feature_id == per_feature_id) %>%
      dplyr::mutate(p_value = -log10(p_value)) %>%
      dplyr::select(dplyr::all_of(tar_cols))

    per_condition_abbr <- ifelse(per_condition %in% condabbr, condabbr[per_condition], "")
    per_cell_type_abbr <- ifelse(per_cell_type == "Monocytes", "Mono.", per_cell_type)
    track_title <- paste(per_feature_id, per_cell_type_abbr, per_condition_abbr)
  } else {
    per_gwas_file <- paste0(per_gwas_id, "_harmonized_data.csv.gz")
    tar_cols <- c(feature_id = "feature_id", snp_id = "SNP", chrom = "chr.outcome", pos = "pos.outcome", log10_pval = "p_value")

    tab <- file.path(.base_dir, "harmonization", per_model, per_cell_type, per_condition, per_gwas_file) %>%
      data.table::fread() %>%
      dplyr::select(exposure, SNP, chr.outcome, pos.outcome, pval.outcome) %>%
      dplyr::filter(exposure == per_feature_id) %>%
      dplyr::mutate(feature_id = per_gwas_id, p_value = -log10(pval.outcome)) %>%
      dplyr::select(dplyr::all_of(tar_cols))

    track_title <- per_gwas_id
  }
  tab <- dplyr::mutate(tab, chrom = stringr::str_remove_all(chrom, "chr"))

  ylim_vec <- c(0, 9)
  if (.marker_snp %in% tab$snp_id)
    marker_snp_track <- dplyr::filter(tab, snp_id == .marker_snp) %>%
      dplyr::select(-c(feature_id, snp_id)) %>%
      makeGRangesFromDataFrame(seqnames.field = "chrom", start.field = "pos", end.field = "pos", keep.extra.columns = TRUE) %>%
      DataTrack(name = track_title, ylim = ylim_vec, cex.title = 1, cex.axis = 0.7, col = "red", cex = 1.5, pch = 18, col.axis = "black")
  else
    marker_snp_track <- head(tab, 1) %>%
      dplyr::select(-c(feature_id, snp_id)) %>%
      makeGRangesFromDataFrame(seqnames.field = "chrom", start.field = "pos", end.field = "pos", keep.extra.columns = TRUE) %>%
      DataTrack(name = track_title, ylim = ylim_vec, cex.title = 1, cex.axis = 0.7, col.axis = "black")

  all_snps_track <- dplyr::filter(tab, snp_id != .marker_snp) %>%
    dplyr::select(-c(feature_id, snp_id)) %>%
    makeGRangesFromDataFrame(seqnames.field = "chrom", start.field = "pos", end.field = "pos", keep.extra.columns = TRUE) %>%
    DataTrack(name = track_title, ylim = ylim_vec, cex.title = 1, cex.axis = 0.7, col.title = "black", background.title = "white", col.axis = "black")
  track <- OverlayTrack(trackList = list(all_snps_track, marker_snp_track))

  tmp <- list(list(tab, track))
  names(tmp) <- paste(per_feature_id, per_model, per_cell_type, per_condition, per_gwas_id, sep = "|")

  tmp
}, .base_dir = base_dir, .marker_snp = marker_snp) %>%
  unlist(recursive = FALSE)

sumstat_tk <- lapply(sumstat, function(ii) ii[[2]])

# Determine the plotting region based on summary stattistics.
sumstat_tabs <- lapply(sumstat, function(ii) ii[[1]])
plot_range <- lapply(sumstat, function(ii) data.frame(seqnames = head(ii[[1]]$chrom, 1), start = max(ii[[1]]$pos), end = min(ii[[1]]$pos))) %>%
  Reduce(rbind, .) %>%
  (function(tab) c(head(tab$seqnames, 1), min(tab$start, tab$end), max(tab$start, tab$end), sep = "-"))
plot_chrom <- plot_range[1]
plot_start <- as.integer(plot_range[2])
plot_end <- as.integer(plot_range[3])

# Chromosome ideogram track
ideogram_tk <- IdeogramTrack(genome = "hg38", chromosome = plot_chrom, background.title = "white")

# Genome axi track
genome_axis_tk <- GenomeAxisTrack(background.title = "white")

# GeneHancer track
genehancer_file <- file.path(proj_dir, "inputs/annotations/GeneHancer/GeneHancer.CD55.chr1_205376872_209305771.grch38.csv")
genehancer_file <- file.path(proj_dir, "inputs/annotations/GeneHancer/GeneHancer.SLFN5.chr17_33729164_36787563.grch38.csv")
genehancer_tab <- data.table::fread(genehancer_file) %>%
  dplyr::filter(geneName %in% c("CD55", "KANSL1", "SLFN5")) %>%
  dplyr::mutate(score = dplyr::if_else(score > 25, as.integer(22.5), score))
anchor_enhancer <- genehancer_tab %>%
  dplyr::select(geneHancerChrom:geneHancerEnd) %>%
  makeGRangesFromDataFrame(seqnames.field = "geneHancerChrom", start.field = "geneHancerStart", end.field = "geneHancerEnd")
anchor_gene <- genehancer_tab %>%
  dplyr::select(geneChrom:geneEnd, geneStrand) %>%
  makeGRangesFromDataFrame(seqnames.field = "geneChrom", start.field = "geneStart", end.field = "geneEnd", strand.field = "geneStrand")

genehancer_tk <- GenomicInteractions(anchor_enhancer, anchor_gene, counts = genehancer_tab$score) %>%
  InteractionTrack(name = "GeneHancer", chromosome = paste0("chr", plot_chrom))
displayPars(genehancer_tk) <- list(
  plot.outside = TRUE, col.outside = "gray", col.interactions = "darkblue", col.anchors.line = "darkblue",
  interaction.dimension = "height", interaction.measure = "counts", anchor.height = 0.1,
  background.title = "gray95", cex.title = 1.0, col.title = "black", rotation.title = 90, cex.axis = 0.7, col.axis = "black"
)

# ATAC-seq annotation track
celltype_vec <- sumstat_query$cell_type %>% unique()
atacseq_tk <- sumstat_query %>%
  dplyr::filter(gwas_id == "", !skip_atac) %>%
  dplyr::arrange(atac_order) %>%
  apply(1, function(vec, .sumstat_tabs, .marker_snp, .plot_chrom, .plot_start, .plot_end) {
    per_feature_id <- vec["feature_id"]
    per_cell_type <- vec["cell_type"]
    per_condition <- vec["condition"]
    per_gwas_id <- vec["gwas_id"]

    per_model <- vec["model"]
    if (per_cell_type == "Monocytes") {
      anno_title <- "Mono."
      anno_path <- file.path(proj_dir, "inputs/atac_seq/peaks_filtered_monocyte.csv.gz")
    } else if (per_cell_type == "NK") {
      anno_title <- "NK"
      anno_path <- file.path(proj_dir, "inputs/atac_seq/peaks_filtered_nkcell.csv.gz")
    } else if (per_cell_type == "CD8T") {
      anno_title <- "CD8T"
      anno_path <- file.path(proj_dir, "inputs/atac_seq/peaks_filtered_cd8t.csv.gz")
    } else {
      anno_title <- "PBMC"
      anno_path <- file.path(proj_dir, "inputs/atac_seq/peaks_filtered.csv.gz")
    }

    anno_range <- data.table::fread(anno_path) %>%
      dplyr::mutate(chr = as.integer(stringr::str_remove_all(chr, "chr"))) %>%
      dplyr::filter(chr == .plot_chrom, (.plot_start <= start & start <= .plot_end) | (.plot_start <= end & end <= .plot_end)) %>%
      makeGRangesFromDataFrame(TRUE, seqnames.field = "chr", start.field = "start", end.field = "end")

  # Annotations containing a SNP used in the association analysis.
    sumstat_idx <- paste(per_feature_id, per_model, per_cell_type, per_condition, per_gwas_id, sep = "|")
    eqtl_range <- .sumstat_tabs[[sumstat_idx]] %>%
      dplyr::select(-c(feature_id, snp_id)) %>%
      dplyr::filter(log10_pval >= -log10(0.05)) %>%
      dplyr::mutate(chrom = as.integer(stringr::str_remove_all(chrom, "chr"))) %>%
      makeGRangesFromDataFrame(seqnames.field = "chrom", start.field = "pos", end.field = "pos")
    anno_eqtl_range <- anno_range[overlapsAny(anno_range, eqtl_range, maxgap = 100), ]
    anno_eqtl_tk <- AnnotationTrack(
      start = start(anno_eqtl_range), width = width(anno_eqtl_range), chromosome = .plot_chrom,
      strand = strand(anno_eqtl_range), id = anno_eqtl_range$gene_name, genome = genome, name = anno_title,
      col.title = "black", background.title = "gray95", background.panel = "white", cex.title = 1,
      cex.axis = 0.7, stacking = "dense"
    )

  # Annotation containing the marker SNP by .marker_snp
    marker_range <- .sumstat_tabs[[sumstat_idx]] %>%
      dplyr::filter(snp_id == .marker_snp) %>%
      dplyr::select(-c(feature_id, snp_id)) %>%
      dplyr::mutate(chrom = as.integer(stringr::str_remove_all(chrom, "chr"))) %>%
      makeGRangesFromDataFrame(seqnames.field = "chrom", start.field = "pos", end.field = "pos")
    marker_range <- anno_range[overlapsAny(anno_range, marker_range, maxgap = 100), ]
    if (length(marker_range) > 0) {
      cat("Mark SNP overlaps ATAC-seq peaks!\n")
      anno_marker_tk <- AnnotationTrack(
        start = start(marker_range), width = width(marker_range), chromosome = .plot_chrom,
        strand = strand(marker_range), genome = genome, name = anno_title, col = "red", fill = "red",
        col.title = "black", background.title = "gray95", background.panel = "white", cex.title = 1,
        cex.axis = 0.7, stacking = "dense"
      )

      OverlayTrack(trackList = list(anno_eqtl_tk, anno_marker_tk))
    } else {
      anno_eqtl_tk
    }
}, .sumstat_tabs = sumstat_tabs, .marker_snp = marker_snp, .plot_chrom = plot_chrom, .plot_start = plot_start, .plot_end = plot_end)


# Feature structure track.
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_pos_tk <- BiomartGeneRegionTrack(
  genome = genome, chromosome = plot_chrom, start = plot_start, end = plot_end,
  filters = list(biotype = "protein_coding", transcript_is_canonical = TRUE), name = "Gene",
  transcriptAnnotation = "symbol", col.title = "black", background.title = "gray95", background.panel = "white", fontsize.group = 20, biomart = ensembl,
  cex.title = 1, cex.axis = 0.7, rotation.title = 90
)


token <- "normal."
# token <- "interaction."
track_list <- c(ideogram_tk, genome_axis_tk, sumstat_tk, atacseq_tk, genehancer_tk, gene_pos_tk)
track_height_list <- c(.7, 1, rep(3.25, length(sumstat_tk)), rep(1.1, length(atacseq_tk)), 3, 3.25)
save_plot_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL", paste0(token, marker_feature, "-", marker_snp, ".track_plot.eQTL_effect.pdf"))
pdf(save_plot_to, width = 6, height = 9)
plotTracks(track_list, sizes = track_height_list, from = plot_start + 300000, to = plot_end - 300000)
dev.off()


#
## The public data resutls was not included in the main figures
#

# Zoom in to CD55
# Histone marker track
hismod_tk <- tibble::tribble(
  ~histone_marker, ~cell_line, ~file,
  "H3k27ac", "GM12878", "wgEncodeBroadHistoneGm12878H3k27acStdSig.bigWig",
  "H3k27ac", "H1-hESC", "wgEncodeBroadHistoneH1hescH3k27acStdSig.bigWig",
  "H3k27ac", "HSMM", "wgEncodeBroadHistoneHsmmH3k27acStdSig.bigWig",
  "H3k27ac", "HUVEC", "wgEncodeBroadHistoneHuvecH3k27acStdSig.bigWig",
  "H3k27ac", "K562", "wgEncodeBroadHistoneK562H3k27acStdSig.bigWig",
  "H3k27ac", "NHEK", "wgEncodeBroadHistoneNhekH3k27acStdSig.bigWig",
  "H3k27ac", "NHLF", "wgEncodeBroadHistoneNhlfH3k27acStdSig.bigWig",
) %>%
apply(1, function(vec) {
  per_bwf <- file.path(proj_dir, "inputs/annotations/UCSC_GB/bbi/wgEncodeReg/wgEncodeRegMarkH3k27ac", vec["file"])
  DataTrack(
    per_bwf, type = "horizon", genome = "hg38", chromosome = paste0("chr", plot_chrom), from = plot_start, to = plot_end, name = vec["cell_line"],
    fill = "darkblue",
    cex.title = 0.7, cex.axis = 0.7, col.title = "black", col.axis = "black", background.title = "gray95", rotation.title = 0, lwd = 1.5
  )
})

# Methylation track
meth_meta_tab <- file.path(proj_dir, "inputs/methylation/GSE184269/meta_infomation.csv") %>% data.table::fread()
meth_info_tab <- file.path(proj_dir, "inputs/methylation/GSE184269/GSE184269_Matrix_processed_FINAL.grch38.csv.gz") %>% data.table::fread()
meth_tab <- meth_info_tab %>%
  dplyr::filter(seqnames == paste0("chr", plot_chrom), plot_start <= start & start <= plot_end) %>%
  tidyr::pivot_longer(-c(seqnames, start, end, width, strand, cpg_id), values_to = "methylation_level", names_to = "sample_id") %>%
  dplyr::mutate(sample_id = stringr::str_remove(sample_id, pattern = "^X")) %>%
  dplyr::left_join(meth_meta_tab, by = "sample_id") %>%
  dplyr::filter(tissue %in% c("Monocytes", "Naive B", "Naive CD4", "Naive CD8", "Natural Killer cells")) %>%
  dplyr::mutate(tissue = dplyr::case_when(
    tissue == "Monocytes" ~ "Mono.", tissue == "Naive B" ~ "B", tissue == "Naive CD8" ~ "CD8T",
    tissue == "Naive CD4" ~ "CD4T", tissue == "Natural Killer cells" ~ "NK"
  )) %>%
  dplyr::select(-c(sample_id, age, sex, strand, width)) %>%
  tidyr::pivot_wider(names_from = c(tissue, donor_id), values_from = methylation_level)

meth_groups <- colnames(meth_tab) %>%
  purrr::discard(~.x %in% c("seqnames", "start", "end", "strand", "cpg_id")) %>%
  stringr::str_remove_all(pattern = "_GT[0-9]{3,3}")
meth_ranges <- makeGRangesFromDataFrame(meth_tab, TRUE)
methylation_tk <- meth_ranges[overlapsAny(meth_ranges, anchor_enhancer), ] %>%
  DataTrack(
    name = "Avg. Met. Lev.", groups = meth_groups, type = c("a"), cex.title = 1, cex.axis = 0.7, col.title = "black",
    col.axis = "black", background.title = "gray95", lwd = 1.5
  )

# SNP pos track
marker_snp_range <- sumstat_tabs %>%
  Reduce(rbind, .) %>%
  dplyr::filter(snp_id == "rs2564978") %>%
  dplyr::slice_max(log10_pval) %>%
  dplyr::mutate(chrom = paste0("chr", chrom), start = pos - 20, end = pos + 20) %>%
  makeGRangesFromDataFrame()
marker_snp_track <- marker_snp_range %>%
  AnnotationTrack(name = "SNP", fill = "red", col = "red", cex.title = 1, cex.axis = 0.7, col.title = "black", background.title = "gray95", rotation.title = 0)

meth_around_marker_snp <- meth_ranges[overlapsAny(meth_ranges, marker_snp_range, maxgap = 15000), ] %>%
  as.data.frame() %>%
  tidyr::pivot_longer(-c(seqnames, start, end, width, strand, cpg_id), values_to = "methylation_level", names_to = "sample") %>%
  dplyr::mutate(celltype = stringr::str_remove_all(sample, pattern = "_GT[0-9]{3,3}")) %>%
  dplyr::group_by(cpg_id) %>%
  dplyr::summarise(test = {
    aov_md <- dplyr::cur_data() %>% aov(methylation_level ~ celltype, data = .)
    summary(aov_md)[[1]][, "Pr(>F)"]
  }, .groups = "keep")

# Plot all tracks
track_list <- c(ideogram_tk, genome_axis_tk, marker_snp_track, gene_pos_tk, hismod_tk, genehancer_tk, methylation_tk)
track_height_list <- c(.6, 0.8, .5, .3, rep(.45, length(hismod_tk)), 3, 6)
save_plot_to <- file.path(proj_dir, "outputs/pseudo_bulk/example_eQTL/normal.CD55-rs2564978.public_data.track_plot.pdf")
pdf(save_plot_to, width = 7, height = 9)
plotTracks(track_list, sizes = track_height_list, from = plot_start + 955000, to = plot_end - 1010000)
dev.off()
