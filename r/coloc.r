#!/usr/bin/env Rscript

#' Colocalization analysis by R/coloc
NULL

options(stringsAsFactors = FALSE, data.table.verbose = FALSE, future.globals.maxSize = 10000 * 1024^2)
suppressPackageStartupMessages({
  library(ggrepel)
  library(magrittr)
  library(tidyverse)
  library(patchwork)
  library(org.Hs.eg.db)
})


#' Fetch data for coloc analysis from the harmonized dataset
#'
fetch_coloc_data <- function(hm_dat, what = "exposure", warn_minp = 5e-5) {
  pos_cols <- c("SNP", "pos.exposure")
  src_cols <- c("beta", "se", "eaf", "samplesize")

  cols <- c(paste(src_cols, what, sep = "."), pos_cols)
  names(cols) <- c("beta", "varbeta", "MAF", "N", "snp", "position")

  tab <- hm_dat %>%
    dplyr::select(dplyr::all_of(cols)) %>%
    dplyr::filter(dplyr::if_all(dplyr::everything(), ~ !is.na(.x))) %>%
    dplyr::mutate(type = "quant", varbeta = varbeta ^ 2) %>%
    (function(dat) {
      res <- NULL
      if (nrow(dat)) {
        beta <- dplyr::select(dat, "snp", "beta") %>% tibble::deframe()
        varbeta <- dplyr::select(dat, "snp", "varbeta") %>% tibble::deframe()
        maf <- dplyr::select(dat, "snp", "MAF") %>% tibble::deframe()

        res <- list(beta = beta, varbeta = varbeta, snp = dat$snp, position = dat$position,
                    type = dat$type[1], N = max(dat$N), MAF = maf)
      }
      res
    })

  if (!is.null(tab)) coloc::check_dataset(tab, warn.minp = warn_minp)
  return(tab)
}


#' Colocalization anlaysis
colocalization_test <- function(hm_dat, min_snps = 10, ...) {
  if (nrow(hm_dat) < min_snps) return(NA)

  qry <- fetch_coloc_data(hm_dat, "exposure")
  ref <- fetch_coloc_data(hm_dat, "outcome")

  ccsum <- c()
  if (TRUE && !is.null(qry) && !is.null(ref)) {
    sink("/dev/null")
    ccsum <- coloc::coloc.abf(qry, ref, ...)$summary
    sink()
  }

  data.frame(
    nsnps = ccsum["nsnps"], H0 = ccsum["PP.H0.abf"], H1 = ccsum["PP.H1.abf"],
    H2 = ccsum["PP.H2.abf"], H3 = ccsum["PP.H3.abf"], H4 = ccsum["PP.H4.abf"]
  )
}


#' Fetch genomic features for the given region
fetch_genomic_features <- function(chrom, start, stop, gffpath, add_chr = TRUE, highlight = NULL) {
  if (add_chr && !stringr::str_detect(chrom, "chr")) chrom <- paste0("chr", chrom)
  which <- GenomicRanges::GRanges(paste0(chrom, ":", start, "-", stop))
  tar_type <- c("exon", "start_codon", "stop_codon", "five_prime_UTR", "three_prime_UTR")

  gene_tab <- rtracklayer::import(gffpath, which = which)
  if (length(gene_tab))
    gene_tab[, c("gene_name", "type", "tag", "gene_type", "transcript_support_level")] %>%
      as.data.frame() %>%
      apply(1, function(e) {
        type <- e["type"]

        gene_type <- e["gene_type"] %in% c("protein_coding", "lncRNA")
        is_canonical <- stringr::str_detect(e["tag"], "Ensembl_canonical")
        is_tsl_1 <- e["transcript_support_level"] %in% c(1, "1")
        if (type == "gene") {
          keep <- gene_type
        } else if (type %in% tar_type) {
          keep <- gene_type && is_canonical
        } else {
          return(NULL)
        }

        if (keep)
          return(data.frame(chrom = e$seqnames, start = e$start, end = e$end, strand = e$strand, width = e$width, gene = e$gene_name, type = e$type))

        return(NULL)
      }) %>%
      Reduce(rbind, .) %>%
      (function(e) if (!is.null(e)) dplyr::mutate(e, highlight = gene %in% highlight))
  else
    NULL
}


#' Plot colocalization
#'
#' @description Function to plot Figure x
plot_coloc <- function(load_path, min_h4 = 0.5, min_h3 = 0.25, override = FALSE, expand = 0, shift = 0, save_to = "./", ...) {
  cc_path <- file.path(load_path, "colocalization_test.csv")
  hm_path <- file.path(load_path, "harmonized_data.csv")
  if (!(file.exists(hm_path) && file.exists(cc_path))) stop("Missing colocalization_test.csv or harmonized_data.csv!")

  kwargs <- list(...)
  if (is.null(kwargs$bfile) || is.null(kwargs$plink_bin)) stop("You have to specify the path to plink and bfiles.")
  if (is.null(kwargs$gffpath)) stop("You have to specify the path to genomic feature files.")

  hm_tab <- data.table::fread(hm_path)
  cc_tab <- data.table::fread(cc_path) %>% dplyr::filter(H4 >= min_h4 | H3 >= min_h3)
  cmb_tab <- dplyr::inner_join(cc_tab, hm_tab, by = c("exposure", "outcome", "id.outcome"))
  height_dict <- c("gene" = 0.3, "exon" = 0.75, "start_codon" = 1.1, "stop_codon" = 1.1, "five_prime_UTR" = 0.5, "three_prime_UTR" = 0.5)

  dplyr::mutate(cmb_tab, outcome = paste(outcome, id.outcome, sep = "_")) %>%
  dplyr::group_by(exposure, outcome) %>%
    dplyr::summarise(fig_path = {
      per_exposure <- dplyr::cur_group()$exposure
      per_outcome <- dplyr::cur_group()$outcome

      fig_path <- file.path(save_to, paste("coloc", per_exposure, per_outcome, "pdf", sep = "."))
      if (!file.exists(fig_path) || override) {
        per_pair <- dplyr::cur_data_all() %>% dplyr::select(exposure, outcome, SNP, chr.exposure, pos.exposure, pval.outcome, pval.exposure)
        topsnp_id <- dplyr::slice_min(per_pair, pval.exposure)$SNP[1]
        tsnp_ld_vec <- rep(0, each = length(per_pair$SNP)) 
        names(tsnp_ld_vec) <- per_pair$SNP

        tryCatch({
          ldmat <- ieugwasr::ld_matrix_local(per_pair$SNP, kwargs$bfile, kwargs$plink_bin, FALSE)^2
          tsnp_ld_vec <- ldmat[, colnames(ldmat) == topsnp_id]
          }, error = function(e) print(e))

        if (length(tsnp_ld_vec) > 0) {
          pair_tab <- tidyr::pivot_longer(per_pair, cols = c("pval.outcome", "pval.exposure")) %>%
            dplyr::mutate(ld2topsnp = tsnp_ld_vec[SNP], ld2topsnp = dplyr::if_else(is.na(ld2topsnp), 0, ld2topsnp), value = -log10(value),
                          TopSNP = SNP == topsnp_id, name = dplyr::if_else(name == "pval.exposure", exposure, outcome))
          gf_chr <- unique(per_pair$chr.exposure)[1]
          gf_start <- as.integer(max(min(per_pair$pos.exposure) + (shift - expand) * 1000, 0))
          gf_end <- as.integer(max(per_pair$pos.exposure) + (shift + expand) * 1000)
          x_lims <- c(gf_start - 10, gf_end + 10)

          gftab <- fetch_genomic_features(gf_chr, gf_start, gf_end, kwargs$gffpath, FALSE, per_exposure)
          if (!is.null(gftab)) {
            gftab <- dplyr::mutate(gftab, height = height_dict[type]) %>%
              dplyr::mutate(start = dplyr::if_else(start <= gf_start, as.integer(gf_start), as.integer(start)),
                            width = dplyr::if_else((start + width) >= gf_end, as.integer(gf_end - start), as.integer(width)))

            y_coord_dict <- dplyr::select(gftab, gene) %>% dplyr::distinct() %>% dplyr::mutate(y_coord = row_number()) %>% tibble::deframe()
            gftab <- dplyr::mutate(gftab, gene_pos = y_coord_dict[gene] * 1.5)

            gf_plot <- ggplot() +
              geom_rect(aes(xmin = start, xmax = end, ymin = gene_pos - height, ymax = gene_pos + height, fill = strand),
                        data = gftab, position = "identity", stat = "identity") +
              geom_text(aes(start, gene_pos, label = gene, color = highlight), data = dplyr::filter(gftab, type == "gene"), hjust = 1.1, size = 2) +
              scale_fill_manual(name = NULL, values = c("-" = "darkblue", "+" = "darkred")) +
              scale_color_manual(name = NULL, values = c("TRUE" = "red", "FALSE" = "black")) +
              lims(x = x_lims) +
              theme_classic() +
              labs(x = "Position (base-pair)", y = NULL) +
              theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_blank())
          }

          topsnp_info <- dplyr::filter(pair_tab, SNP == topsnp_id)
          cc_plot <- ggplot() +
            geom_point(aes(x = pos.exposure, y = value, color = ld2topsnp), pair_tab, alpha = 0.85) +
            geom_point(aes(x = pos.exposure, y = value), topsnp_info, size = 3, shape = 23, fill = "purple") +
            geom_text_repel(aes(x = pos.exposure, y = value, label = SNP), topsnp_info, min.segment.length = 0) +
            geom_hline(yintercept = -log10(5e-8), linetype = "dashed") +
            facet_wrap(~name, ncol = 1, scales = "free_y", strip.position = "right") +
            scale_colour_gradientn(colours = c("darkblue", "blue", "orange", "red", "darkred"), name = quote("LD (" ~ R^2 ~ ")")) +
            labs(x = NULL, y = quote(~ -Log[10] ~ "(p-value)")) +
            lims(y = c(0, NA), x = x_lims) +
            theme_bw() +
            theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

          if (!is.null(gftab)) { cmb_plot <- cc_plot / gf_plot + plot_layout(heights = c(4, 1)) }
          else { cmb_plot <- cc_plot }

          ggsave(fig_path, cmb_plot, width = 7, height = 5)
        }
      }
      NULL
    })
}


#
## Main steps
#
proj_dir <- "/home/zzhang/Documents/projects/wp_bcg_eqtl"
in_dir <- file.path(proj_dir, "outputs/pseudo_bulk")

# Two model used in the analysis.
mode_vec <- c("normal", "interaction")

# Main cell types.
cell_type_vec <- c("Monocytes", "CD4T", "CD8T", "NK", "B") # , "pDC", "mDC")

# All comparisons used in the analysis.
condition_vec <- c("T0_LPS.vs.T0_RPMI", "T3m_LPS.vs.T3m_RPMI", "T3m_RPMI.vs.T0_RPMI")#, "T3m_LPS.vs.T0_RPMI")

# Reference genotypes panel for local clumping
eur_1kg_geno <- file.path(proj_dir, "inputs/reference/genotypes/GRCh38/EUR")

# Reference genomic feature file
gff_file <- file.path(proj_dir, "inputs/reference/Gencode/gencode.v41.basic.annotation.Ec.transcript.autosome.gff.gz")

override <- TRUE
in_file <- "run_2_harmonized_data.csv"
out_file <- "run_2_colocalization_test.csv"
for (mode in mode_vec) {
  in_base_dir <- file.path(proj_dir, "outputs/pseudo_bulk/harmonization", mode)
  out_base_dir <- file.path(proj_dir, "outputs/pseudo_bulk/coloc", mode)

  for (cell_type in cell_type_vec) {
    if (mode == "normal") {
      hm_dat_path <- file.path(in_base_dir, cell_type, in_file)
      hm_dat <- data.table::fread(hm_dat_path)

      cc_save_to <- file.path(out_base_dir, cell_type)
      if (!dir.exists(cc_save_to)) dir.create(cc_save_to, recursive = TRUE)

      if (!file.exists(cc_save_to) || override) {
        cat("[I]: Estimating co-localization by coloc for", cell_type, "...\n")
        ccres <- dplyr::group_by(hm_dat, exposure, outcome, id.outcome) %>%
          dplyr::summarise(cc_res = colocalization_test(dplyr::cur_data_all())) %>%
          data.table::as.data.table() %>%
          dplyr::rename_with(dplyr::starts_with("cc_res."), .fn = ~ stringr::str_remove(.x, "cc_res.")) %>%
          dplyr::filter(dplyr::if_all(dplyr::starts_with("H"), .fns = ~ !is.na(.x)))

        data.table::fwrite(ccres, file.path(cc_save_to, out_file))
      }
    } else {
      for (condition in condition_vec) {
        hm_dat_path <- file.path(in_base_dir, cell_type, condition, in_file)
        hm_dat <- data.table::fread(hm_dat_path)

        cc_save_to <- file.path(out_base_dir, cell_type, condition)
        if (!dir.exists(cc_save_to)) dir.create(cc_save_to, recursive = TRUE)

        cat("[I]: Estimating co-localization by coloc for", cell_type, condition, "...\n")
        ccres <- dplyr::group_by(hm_dat, exposure, outcome, id.outcome) %>%
          # dplyr::filter(dplyr::if_all(.fns = ~ !is.na(.x))) %>%
          dplyr::summarise(cc_res = colocalization_test(cur_data_all())) %>%
          data.table::as.data.table() %>%
          dplyr::rename_with(dplyr::starts_with("cc_res."), .fn = ~ stringr::str_remove(.x, "cc_res.")) %>%
          dplyr::filter(dplyr::if_all(dplyr::starts_with("H"), .fns = ~ !is.na(.x)))

        data.table::fwrite(ccres, file.path(cc_save_to, out_file))
      }
    }
  }
}
