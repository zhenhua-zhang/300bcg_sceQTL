#!/usr/bin/env Rscript
library(magrittr)
# library(GenomicRanges)

# NOTE:
# 1. Make sure the chromsome names are consistent across all input files, i.e., with prefix "chr".
# 2. The top SNPs of each QTL were not clumpped, meaning their pired-wise LD R^2 could be high.

#' Load eQTL summary statistics
#'
#'@description A utility function to load summary statistics from the disk
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
    for (nm in names(other_info)) { cmb_qtltab[nm] <- other_info[nm] }
  }

  cmb_qtltab %>% dplyr::filter(dplyr::if_any(dplyr::everything(), ~ !is.na(.x)))
}


#' Load annotation files
#'
#' @example The example name of a annotation file is TYPE.CELLTYPE.SOURCE[.txt]
#' The file name has three required parts, if any of them is unknown, one can use any placeholder for it. E.g., Enhancer.NULL.ATACseq.txt
#' The example input content could be:
#' CHROM  START  STOP  METAINFO
#' 1      1000   1500  DPP9
load_ann_tab <- function(dir) {
  ann_db_path <- list.files(dir, recursive = TRUE, full.names = TRUE) %>%
    lapply(function(e) {
      if (!file.exists(e)) cat("[W]: the file does not exist, skip it", e, "...\n")
      extra_info <- basename(e) %>% stringr::str_split(pattern = "\\.", simplify = TRUE)

      ranges <- data.table::fread(e) %>%
        (function(tab) {
          GenomicRanges::GRanges(
            seqnames = tab$Chrom,
            ranges = IRanges::IRanges(tab$Start, tab$Stop),
            strand = S4Vectors::Rle(BiocGenerics::strand(c("*")), nrow(e)),
            metainfo = tab$Metainfo
          )
        })

      list(ranges = ranges, annotype = extra_info[1], celltype = extra_info[2], data_source = extra_info[3])
    })
}


#' Estimate enrichment
#'
#'@description Estimate the enrichment of identified QTL SNPs in genomic features annotations.
#'
estimate_enrichment <- function(opts) {
  qtl_tab <- load_eqtl_tab(opts$qtl_db) %>%
    dplyr::select(snp_chromosome, snp_position, snp_id, gene_name, global_corrected_pValue)

  # Population ranges
  all_snp_ranges <- qtl_tab %>%
    dplyr::group_by(snp_chromosome, snp_position, snp_id) %>%
    dplyr::summarize(
      gene_name = paste0(gene_name, collapse = "|"),
      is_top_snp = any(!is.na(global_corrected_pValue))
    ) %>%
    (function(e) {
      GenomicRanges::GRanges(seqnames = e$snp_chromosome,
        ranges = IRanges::IRanges(e$snp_position, e$snp_position),
        strand = S4Vectors::Rle(BiocGenerics::strand(c("*")), nrow(e)),
        snp_id = e$snp_id,
        gene_name = e$gene_name,
        is_top_snp = e$is_top_snp
      )
    })

  # Observed ranges
  top_snp_ranges <- all_snp_ranges[all_snp_ranges$is_top_snp, ]

  # any(stringr::str_split(SNPS, pattern == "\\|", simplify = TRUE) %in% total_snps)
  # Haplotype blocks
  total_snps <- qtl_tab$snp_id %>% unique()
  hb_tab <- data.table::fread(opts$blocks)
  hb_ranges <- hb_tab %>%
    dplyr::filter(
      as.integer(opts$min_block_size) <= KB, KB <= as.integer(opts$max_block_size)
    ) %>%
    tidyr::separate_rows(SNPS) %>%
    dplyr::filter(SNPS %in% total_snps) %>%
    dplyr::group_by(CHR, BP1, BP2, KB) %>%
    dplyr::summarize(SNPS = paste0(SNPS, collapse = "|")) %>%
    dplyr::ungroup() %>%
    (function(e) {
      GenomicRanges::GRanges(seqnames = e$CHR,
        ranges = IRanges::IRanges(e$BP1, e$BP2),
        strand = S4Vectors::Rle(BiocGenerics::strand(c("*")), nrow(e)),
        span_kbp = e$KB,
        SNP = e$SNPS
      )
    })

  # Population and trials blocks.
  all_blocks <- hb_ranges[IRanges::overlapsAny(hb_ranges, all_snp_ranges), ]
  qtl_blocks <- hb_ranges[IRanges::overlapsAny(hb_ranges, top_snp_ranges), ]
  population_size <- length(all_blocks)
  trial_size <- length(qtl_blocks)

  # Estimate enrichment per genomic feature annotation
  anno_ranges <- load_ann_tab(opts$ann_db)
  lapply(anno_ranges, function(per_range) {
    min_ov_size <- ifelse(per_range$annotype %in% c("GWAS", "gwas"), 1, opts$min_ov_size)
    nr_suc_exp <- IRanges::overlapsAny(all_blocks, per_range$ranges, minoverlap = min_ov_size) %>% sum()
    nr_suc_obs <- IRanges::overlapsAny(qtl_blocks, per_range$ranges, minoverlap = min_ov_size) %>% sum()

    estres <- data.frame(
      annotype = per_range$annotype, celltype = per_range$celltype, data_source = per_range$data_source,
      log2od = NA, p_val = NA, obs_prob = NA, exp_prob = NA,
      expectation = nr_suc_exp, population = population_size, observation = nr_suc_obs, trials = trial_size
    )

    if (nr_suc_exp > 0 && nr_suc_obs > 0) {
      pp_test <- prop.test(nr_suc_obs, trial_size, nr_suc_exp / population_size)
      log2_od <- log2(pp_test$estimate / pp_test$null.value)

      estres$log2od <- log2_od
      estres$p_val <- pp_test$p.value
      estres$obs_prob <- pp_test$estimate
      estres$exp_prob <- pp_test$null.value
    }

    estres
  }) %>%
  Reduce(rbind, .)
}


# CLI options
p <- optparse::OptionParser(description = "Genomic annotation enrichment")
p <- optparse::add_option(p, c("-s", "--qtl-db"), help = "The directory to load eQTL results. Required")
p <- optparse::add_option(p, c("-b", "--blocks"), help = "The haplotype blocks. Required")
p <- optparse::add_option(p, c("-a", "--ann-db"), help = "Annotation database in bed format. Required")
p <- optparse::add_option(p, c("-m", "--min-block-size"), default = 5, help = "The minimal size (base-pairs) of blocks used to estimate the enrichment. Default: %(default)s")
p <- optparse::add_option(p, c("-M", "--max-block-size"), default = 100, help = "The maximum size (base-pairs) of blocks used to estimate the enrichment. Default: %(default)s")
p <- optparse::add_option(p, c("--min-ov-size"), default = 20, help = "The minimal overlapped size betwen two intervals. Default: %(default)s")
p <- optparse::add_option(p, c("-o", "--out"), default = "hb_enrichment.csv", help = "Output prefix. Default: %(default)s")
p <- optparse::add_option(p, c("-p", "--plot"), action = "store_true", help = "A flag to save enrichment plots.")

opts <- optparse::parse_args(p, convert_hyphens_to_underscores = TRUE)
if (is.null(opts$qtl_db)) stop("-s/--qtl-db is required!")
if (is.null(opts$blocks)) stop("-b/--blocks is required!")
if (is.null(opts$ann_db)) stop("-a/--ann-db is required!")

# The main function
estimate_enrichment(opts) %>% data.table::fwrite(opts$out)
