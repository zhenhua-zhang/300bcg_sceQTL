#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, data.table.verbose = FALSE)
suppressPackageStartupMessages({
  library(Gviz)
  library(biomaRt)
  library(magrittr)
  library(tidyverse)
  library(data.table)
  library(GenomicRanges)
})


#' Load eQTL summary statistics
#'
#' @description A utility function to load summary statistics from the disk
load_eqtl_tab <- function(dir, fdr = 0.05, fdr_col = "global_corrected_pValue", other_info = NULL) {
  joinby <- c(
    "snp_id", "ensembl_gene_id", "feature_start", "feature_end",
    "feature_chromosome", "feature_id", "gene_name", "feature_strand",
    "n_samples", "n_e_samples", "snp_chromosome", "snp_position",
    "assessed_allele"
  )

  discard <- c(
    "p_value", "beta", "beta_se", "empirical_feature_p_value", "alpha_param",
    "beta_param", "call_rate", "maf", "hwe_p"
  )

  top_qtl_fpath <- file.path(dir, "top_qtl_results_all_FDR.txt")
  top_qtltab <- data.table::fread(top_qtl_fpath, data.table = FALSE) %>%
    dplyr::select(-dplyr::one_of(discard)) %>%
    (function(dat) dat[dat[, fdr_col] < fdr, ])

  all_qtl_fpath <- file.path(dir, "qtl_results_all.txt")
  all_qtltab <- data.table::fread(all_qtl_fpath, data.table = FALSE)

  cmb_qtltab <- dplyr::right_join(top_qtltab, all_qtltab, by = joinby) %>%
    dplyr::mutate(QTL = stringr::str_c(gene_name, snp_id, sep = "-"))

  if (!is.null(other_info) && !is.null(names(other_info))) {
    for (nm in names(other_info)) {
      cmb_qtltab[nm] <- other_info[nm]
    }
  }

  cmb_qtltab %>% dplyr::filter(dplyr::if_any(dplyr::everything(), ~ !is.na(.x)))
}

proj_dir <- "~/Documents/projects/wp_bcg_eqtl"
model <- "normal"
model <- "interaction"

tar_gene <- "ADCY3"; tar_celltype <- "Monocytes"; tar_gwas <- c("BodyMassIndex", "ObesityClass1", "ObesityClass2")
tar_gene <- "NUMA1"; tar_celltype <- "Monocytes"; tar_gwas <- c("Type2Diabetes"); condition <- "T3m_RPMI.vs.T0_RPMI"
tar_gene <- "TDRKH"; tar_celltype <- "Monocytes"; tar_gwas <- c("Psoriasis"); condition <- "T3m_LPS.vs.T3m_RPMI"
tar_gene <- "PIAS3"; tar_celltype <- "Monocytes"; tar_gwas <- c("GoutDisease"); condition <- "T0_LPS.vs.T0_RPMI"
tar_gene <- "AGFG1"; tar_celltype <- "Monocytes"; tar_gwas <- c("ThyroidCancer"); condition <- "T3m_LPS.vs.T3m_RPMI"
tar_gene <- "ZNF100"; tar_celltype <- "Monocytes"; tar_gwas <- c("ThyroidCancer"); condition <- "T3m_LPS.vs.T3m_RPMI"
tar_gene <- "SLC9A1"; tar_celltype <- "Monocytes"; tar_gwas <- c("Type2Diabetes"); condition <- "T3m_RPMI.vs.T0_RPMI"
tar_gene <- "UHRF1BP1"; tar_celltype <- "Monocytes"; tar_gwas <- c("Schizophrenia"); condition <- "T3m_LPS.vs.T3m_RPMI"
tar_gene <- "SLC25A45"; tar_celltype <- "Monocytes"; tar_gwas <- c("ThyroidCancer"); condition <- "T3m_LPS.vs.T3m_RPMI"

tar_gene <- "POLR1D"; tar_celltype <- "CD4T"; tar_gwas <- c("CoronaryHeartDisease"); condition <- "T3m_RPMI.vs.T0_RPMI"
tar_gene <- "SUCO"; tar_celltype <- "CD8T"; tar_gwas <- "Type2Diabetes"; condition <- "T3m_RPMI.vs.T0_RPMI"
tar_gene <- "ITGA10"; tar_celltype <- "CD8T"; tar_gwas <- "GoutDisease"; condition <- "T3m_RPMI.vs.T0_RPMI"
tar_gene <- "RARS2"; tar_celltype <- "CD8T"; tar_gwas <- "ThyroidCancer"; condition <- "T3m_LPS.vs.T3m_RPMI"
tar_gene <- "COMMD1"; tar_celltype <- "NK"; tar_gwas <- c("UlcerativeColitis", "InflammatoryBowelDisease"); condition <- "T0_LPS.vs.T0_RPMI"
tar_gene <- "SRSF12"; tar_celltype <- "B"; tar_gwas <- "CrohnsDisease"; condition <- "T3m_RPMI.vs.T0_RPMI"


# The genome coordination version
genome <- "hg38"


# Colocalization track
coloc_track_list <- list()

## Load the harmonized variants between eQTL and GWAS summaries
hm_tab_path <- file.path(proj_dir, "outputs/pseudo_bulk/harmonization", model, tar_celltype, condition, "harmonized_data.csv")
hm_tab <- data.table::fread(hm_tab_path) %>%
  dplyr::filter(exposure == tar_gene, outcome %in% tar_gwas) %>%
  dplyr::select(dplyr::one_of(c("SNP", "outcome", "chr.exposure", "pos.exposure", "pval.outcome", "pval.exposure"))) %>%
  dplyr::mutate(pval.outcome = -log10(pval.outcome), pval.exposure = -log10(pval.exposure))

## Decide the region of the whole plot.
region_chr <- hm_tab$chr.exposure %>% head(1)
region_start <- hm_tab$pos.exposure %>% min()
region_end <- hm_tab$pos.exposure %>% max()

## The top eQTL SNP will be hightlighted in the track plot
top_snp <- hm_tab %>%
  dplyr::select(pval.exposure, SNP) %>%
  dplyr::slice_max(pval.exposure, n = 1) %>%
  dplyr::pull(SNP) %>%
  head(1)

## Dot plots of eQTL
eqtl_tab <- hm_tab %>% dplyr::select(-c(outcome, pval.outcome)) %>% dplyr::distinct()
coloc_ylim_eqtl <- max(eqtl_tab$pval.exposure) + 1
coloc_top_eqtl <- eqtl_tab %>%
  dplyr::filter(SNP %in% top_snp) %>%
  makeGRangesFromDataFrame(seqnames.field = "chr.exposure", start.field = "pos.exposure", end.field = "pos.exposure", keep.extra.columns = TRUE) %>%
  DataTrack(name = tar_gene, ylim = c(0, coloc_ylim_eqtl), cex.title = 1, cex.axis = 0.7, col = "red", cex = 1.5, pch = 18)

eqtl_range <- eqtl_tab %>%
  makeGRangesFromDataFrame(seqnames.field = "chr.exposure", start.field = "pos.exposure", end.field = "pos.exposure", keep.extra.columns = TRUE)
coloc_all_eqtl <- DataTrack(eqtl_range, name = tar_gene, ylim = c(0, coloc_ylim_eqtl), cex.title = 1, cex.axis = 0.7)
coloc_track_list[[tar_gene]] <- OverlayTrack(trackList = list(coloc_all_eqtl, coloc_top_eqtl), background.title = "gray50")

for (per_gwas in tar_gwas) {
  ## Dot plots of GWAS
  gwas_tab <- hm_tab %>% dplyr::filter(outcome == per_gwas) %>% dplyr::select(-c(outcome, pval.exposure))
  coloc_gwas_ylim <- max(gwas_tab$pval.outcome) + 1
  coloc_top_gwas <- gwas_tab %>%
    dplyr::filter(SNP %in% top_snp) %>%
    makeGRangesFromDataFrame(seqnames.field = "chr.exposure", start.field = "pos.exposure", end.field = "pos.exposure", keep.extra.columns = TRUE) %>%
    DataTrack(name = per_gwas, ylim = c(0, coloc_gwas_ylim), cex.title = 1, cex.axis = 0.7, col = "red", cex = 1.5, pch = 18)
  coloc_all_gwas <- gwas_tab %>%
    makeGRangesFromDataFrame(seqnames.field = "chr.exposure", start.field = "pos.exposure", end.field = "pos.exposure", keep.extra.columns = TRUE) %>%
    DataTrack(name = per_gwas, ylim = c(0, coloc_gwas_ylim), cex.title = 1, cex.axis = 0.7)
  coloc_track_list[[per_gwas]] <- OverlayTrack(trackList = list(coloc_all_gwas, coloc_top_gwas), background.title = "gray50")
}


# Annotation track
anno_title <- "ATACseq (Mono.)"
anno_path <- file.path(proj_dir, "inputs/atac_seq/peaks_filtered_monocyte.csv.gz")
anno_tab <- data.table::fread(anno_path)
anno_range <- anno_tab %>%
  dplyr::mutate(chr = as.integer(stringr::str_remove_all(chr, "chr"))) %>%
  dplyr::filter(chr == region_chr, (region_start <= start & start <= region_end) | (region_start <= end & end <= region_end),
                !reg_feature %in% c("reg_NONE")) %>%
  makeGRangesFromDataFrame(seqnames.field = "chr", start.field = "start", end.field = "end", keep.extra.columns = TRUE)
anno_range <- anno_range[overlapsAny(anno_range, eqtl_range), ]

annot_tk <- AnnotationTrack(
    start = start(anno_range), width = width(anno_range), chromosome = region_chr,
    strand = strand(anno_range), id = anno_range$gene_name, genome = genome, name = anno_title,
    background.title = "gray50", background.panel = "white", cex.title = 1, cex.axis = 0.7
)


# Chromosome ideogram track
ideogram_tk <- IdeogramTrack(genome = genome, chromosome = region_chr)


# Genome axi track
genome_axis_tk <- GenomeAxisTrack()


# Feature structure track.
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_pos_track <- BiomartGeneRegionTrack(
  genome = genome, chromosome = region_chr, start = region_start, end = region_end,
  filters = list(biotype = "protein_coding", transcript_is_canonical = TRUE), name = "Ensembl gene",
  transcriptAnnotation = "symbol", background.title = "gray50", background.panel = "white", fontsize.group = 20, biomart = ensembl,
  cex.title = 1, cex.axis = 0.7
)


# Plot all tracks
track_list <- c(
  list(idoogram = ideogram_tk, genome_axis = genome_axis_tk),
  coloc_track_list,
  list(atac_seq = annot_tk, gene_pos = gene_pos_track)
)
track_hight <- c(list(0.5, 1), rep(4, length(coloc_track_list)), 2, 4)

# Save the plot
file_name <- ifelse(model == "normal", paste(tar_gene, "trackplot.pdf", sep = "."), paste(tar_gene, condition, "trackplot.pdf", sep = "."))
save_to <- file.path(proj_dir, "outputs/pseudo_bulk", "example_eQTL", model, file_name)
pdf(save_to, width = 5, height = 12)
plotTracks(track_list, sizes = track_hight, from = region_start - 10000, to = region_end + 10000)
dev.off()
