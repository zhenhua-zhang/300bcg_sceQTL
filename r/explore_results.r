#!/usr/bin/env Rscript
# Rscript to explore the sceQTL results.
# Used libraries:
options(stringsAsFactors = FALSE, future.globals.maxSize = 10000 * 1024^2)

library(magrittr)
library(tidyverse)
library(data.table)
library(ggupset)
library(org.Hs.eg.db)
library(ggrepel)
library(patchwork)

if (FALSE) {
  # GWAS analysis of regulatory and function enrichment with LD correction
  library(garfield) # To do function enrichment by Garfield
  library(fgsea) # Fast gene set enrichment analysis
  library(TwoSampleMR) # Two sample Mendelian randomization
  library(coloc) # Colocalization analysis
  library(gvis) # Genonmic track visualization package.
  library(biomaRt)
  library(ieugwasr)
}



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


#' Plot a dot plot to show the concordance between two cell types.
#'
#' @description Function to plot Figure ?
#' @param
draw_dotplot <- function(x, y, joinby, labels, save_to = "concordance_dot_plot.png") {
  cmb_qtltab <- dplyr::inner_join(x, y, by = joinby) %>%
    dplyr::mutate(
      Significant = dplyr::case_when(
        is.na(global_corrected_pValue.x) & is.na(global_corrected_pValue.y) ~ "None",
        !is.na(global_corrected_pValue.x) & is.na(global_corrected_pValue.y) ~ labels[1],
        is.na(global_corrected_pValue.x) & !is.na(global_corrected_pValue.y) ~ labels[2],
        !is.na(global_corrected_pValue.x) & !is.na(global_corrected_pValue.y) ~ "Both",
      ),
      Significant = factor(Significant, levels = c("Both", labels, "None"))
    ) %>%
    dplyr::filter(p_value.x < 0.05 & p_value.y < 0.05) %>%
    dplyr::arrange(dplyr::desc(Significant))

  nr_sign <- table(cmb_qtltab$Significant)
  new_legend_labels <- paste0(names(nr_sign), "(", nr_sign, ")")
  names(new_legend_labels) <- names(nr_sign)

  color_map <- c("red", "darkgreen", "darkblue", "gray")
  names(color_map) <- c("Both", labels, "None")

  x_lab <- paste("Beta in", labels[1])
  y_lab <- paste("Beta in", labels[2])
  dot_plot <- cmb_qtltab %>%
    ggplot(aes(x = beta.x, y = beta.y)) +
    geom_point(aes(color = Significant), alpha = 0.5, size = 0.5) +
    geom_hline(linetype = "dotted", yintercept = 0) +
    geom_vline(linetype = "dotted", xintercept = 0) +
    scale_color_manual(name = "Significant in:", values = color_map, labels = new_legend_labels) +
    labs(x = x_lab, y = y_lab) +
    theme_classic()

  n_chars <- nchar(c(x_lab, y_lab, "Significant in:")) %>% max()
  h <- 7
  w <- n_chars / 18 + h

  ggplot2::ggsave(save_to, plot = dot_plot, width = w, height = h)
}


#' Plot a upset plot to show number of shared QTL across all runs
#'
#' @description Function to plot Figure ?
#' @param
draw_upsetplot <- function(qtltabs, save_to = "./", groupby = "feature_id", width = 10, height = 10, override = TRUE) {
  gup_path <- file.path(save_to, "shared_QTL_upset_plot.png")
  if (!file.exists(gup_path) || override) {
    qtltabs <- qtltabs %>%
      lapply(function(x) dplyr::filter(x, !is.na(global_corrected_pValue))) %>%
      Reduce(rbind, .) %>%
      dplyr::distinct()

    if (groupby == "feature_id") {
      qtltab_upset <- qtltabs %>%
        dplyr::group_by(feature_id) %>%
        dplyr::summarise(comparisons = list(dplyr::cur_data()$cell_type))
    } else if (groupby == "QTL") {
      qtltab_upset <- qtltabs %>%
      dplyr::group_by(QTL) %>%
      dplyr::summarise(comparisons = list(paste(dplyr::cur_data()$cell_type, dplyr::cur_data()$condition, sep = "@")))
    }

    g_up <- qtltab_upset %>%
      ggplot2::ggplot(ggplot2::aes(x = comparisons)) +
      ggplot2::geom_bar() +
      ggplot2::labs(x = NULL, y = NULL) +
      ggupset::scale_x_upset(n_intersections = 50)

    ggplot2::ggsave(gup_path, plot = g_up, width = width, height = height)
  } else {
    cat("[W]: Found", gup_path, " Skipping ...\n")
  }
}


#' Functional enrichment by GARFIELD
garfield_test <- function(eqtl_tab, garf_dir, gwas_trait = "trait", snp_chr_col = "snp_chromosome", snp_pos_col = "snp_position",
                          snp_pval_col = "global_corrected_pValue", save_to = "./garfield_out") {
  pval_dir <- file.path(save_to, "pval", gwas_trait)
  if (!dir.exists(pval_dir)) dir.create(pval_dir, recursive = TRUE)

  output_dir <- file.path(save_to, "output")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  file.symlink(file.path(garf_dir, c("annotation", "tags", "maftssd")), save_to)

  tar_cols <- c("CHR" = snp_chr_col, "BP" = snp_pos_col, "pval" = snp_pval_col)
  chroms <- eqtl_tab %>%
    dplyr::select(dplyr::all_of(tar_cols)) %>%
    dplyr::filter(!is.na(pval)) %>%
    dplyr::arrange(BP) %>%
    dplyr::group_by(CHR) %>%
    dplyr::summarise(save_to = {
      per_chr <- dplyr::cur_group()
      if (!stringr::str_detect(per_chr, "chr")) {
        per_chr <- paste0("chr", per_chr)
      }

      data.table::fwrite(dplyr::cur_data(), file.path(pval_dir, per_chr), col.names = FALSE, sep = " ")
    }) %>%
    dplyr::select(CHR) %>%
    unlist() %>%
    as.character()

  out_pref <- file.path(save_to, "output", gwas_trait)
  garfield::garfield.run(out_pref, data.dir = save_to, trait = gwas_trait, run.option = "prep", chrs = chroms)

  prep_file <- paste0(out_pref, ".prep")
  garfield::garfield.run(out_pref, data.dir = save_to, trait = "", run.option = "perm", prep.file = prep_file, chrs = chroms)

  perm_file <- paste0(out_pref, ".perm")
  garfield::garfield.plot(perm_file, num_perm = 100000, plot_title = gwas_trait, output_prefix = out_pref, filter = 1)
}


#' Functional enrichment by GoShifter
goshifter_test <- function() {
  cat("GoShifter is implemented by Python and has a CLI. Check https://github.com/immunogenomics/goshifter for more information.\n")
}


#' Plot heatmap for functional enrichment by GoShifter
#'
#' @description Function to plot Figure x
plot_goshifter <- function(fpath) {
  enrich_out_list <- list.files(fpath, pattern = "*.enrich", recursive = TRUE, full.name = TRUE)

  p_val_perm <- enrich_out_list %>%
    lapply(function(e) {
      fread(e) %>%
        dplyr::arrange(enrichment) %>%
        dplyr::mutate(rank = 1:nrow(.), prob = 1 - rank / nrow(.)) %>%
        dplyr::filter(nperm == 0) %$%
        prob
    }) %>%
    unlist()

  cell_types <- enrich_out_list %>% dirname() %>% basename()
  mode <- enrich_out_list %>% dirname() %>% dirname() %>% basename()

  cmb_tab <- basename(enrich_out_list) %>%
    stringr::str_remove(".nperm2000.enrich") %>%
    stringr::str_split(pattern = "-", simplify = TRUE) %>%
    as.data.frame() %>%
    dplyr::rename("Comparison" = "V1", "Annotation" = "V2") %>%
    dplyr::mutate(Cell_type = cell_types, Mode = mode, P_val_perm = p_val_perm, Significant = p_val_perm < 0.05)

  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(ggplot2::aes(x = Cell_type, y = Comparison, fill = P_val_perm), cmb_tab) +
    ggplot2::geom_point(ggplot2::aes(x = Cell_type, y = Comparison, shape = Significant), (cmb_tab %>% dplyr::filter(Significant))) +
    ggplot2::scale_fill_gradient(low = "darkred", high = "white") +
    ggplot2::facet_wrap(~Annotation, nrow = 1) +
    ggplot2::theme_classic()

  save_to <- file.path(fpath, "GoShifter_annotation_enrichment.png")
  ggplot2::ggsave(save_to, plot = p, width = 12, height = 4)
}


#' Harmonize variants
harmonize_vars <- function(eqtl, gwas, alt_af_db, eqtl_pval_col = "p_value", eqtl_pval = 0.1,
                           clump = TRUE, clump_kb = 20, clump_r2 = 0.8, clump_bfile = NULL, do_proxy = TRUE,
                           plink_bin = "plink", bcftools_bin = "bcftools") {
  # Set up bcftools and plink path
  gwasvcf::set_bcftools(bcftools_bin)
  gwasvcf::set_plink(plink_bin)

  # Select and rename columns in the eQTL summary statistics.
  tar_col <- c(
    "SNP" = "snp_id", "beta" = "beta", "se" = "beta_se", "effect_allele" = "assessed_allele", "eaf" = "maf", "Phenotype" = "gene_name",
    "chr" = "snp_chromosome", "pos" = "snp_position", "samplesize" = "n_samples", "pval" = eqtl_pval_col
  )

  # Load the eQTL summary statistics as the exposure dataset.
  eqtl_tab <- dplyr::select(eqtl, dplyr::all_of(tar_col)) %>%
    dplyr::filter(pval < eqtl_pval, !is.na(pval)) %>%
    dplyr::mutate(eaf = alt_af_db[SNP]) %>%
    TwoSampleMR::format_data(type = "exposure")

  # Sample size
  gwas_path <- lapply(gwas, function(e) e[[1]]) %>% unlist()
  sample_size_dict <- lapply(gwas, function(e) e[[2]]) %>% unlist()

  # Fetch public available GWAS.
  tar_snps <- unique(eqtl_tab$SNP)
  do_proxy <- ifelse(do_proxy && !is.null(clump_bfile), "no", "yes")
  gwas_tab <- lapply(gwas_path, function(path, .rsid, .pval) {
    outcome_info <- path %>% basename() %>% stringr::str_remove_all(".vcf.gz") %>% stringr::str_split("_", n = 3, simplify = TRUE)
    gwasvcf::query_gwas(path, rsid = .rsid, proxies = do_proxy, bfile = proxy_ld_file, threads = 4) %>% 
      gwasvcf::query_gwas(pval = .pval) %>%
      gwasglue::gwasvcf_to_TwoSampleMR("outcome") %>%
      dplyr::filter(!is.na(effect_allele.outcome), !is.na(other_allele.outcome), nchar(effect_allele.outcome) == 1, nchar(other_allele.outcome) == 1) %>%
      dplyr::mutate(outcome = outcome_info[2], id.outcome = outcome_info[3])
  }, .rsid = tar_snps, .pval = 0.5) %>%
  Reduce(rbind, .)

  # Harmonized exposure and outcome dataset.
  suppressMessages({ # Here we suppress messages from the package to make the stdout less mess.
    hm_tab <- TwoSampleMR::harmonise_data(eqtl_tab, gwas_tab) %>%
      dplyr::mutate(samplesize.outcome = dplyr::if_else(is.na(samplesize.outcome), as.integer(sample_size_dict[id.outcome]), as.integer(samplesize.outcome)))
  })

  # Do clumping
  if (clump && !is.null(clump_bfile)) {
    clump_snps <- dplyr::select(hm_tab, SNP, pval.exposure, exposure, outcome, mr_keep) %>% dplyr::filter(mr_keep) %>% dplyr::group_by(exposure, outcome)
    clump_cols <- c("rsid" = "SNP", "pval" = "pval.exposure", "id" = "exposure")
    hm_tab <- dplyr::summarise(clump_snps, clumped = {
      tryCatch({
        dplyr::cur_data_all() %>%
          dplyr::select(dplyr::all_of(clump_cols)) %>%
          ieugwasr::ld_clump(plink_bin = plink_bin, bfile = clump_bfile, clump_kb = clump_kb, clump_r2 = clump_r2) %>%
          dplyr::select(dplyr::all_of(c("SNP" = "rsid"))) %>%
          dplyr::mutate(clump_keep = TRUE)
      }, error = function(e) data.frame(SNP = NULL, clump_keep = NULL))
    }) %>%
      data.table::as.data.table() %>%
      dplyr::rename_with(.cols = dplyr::starts_with("clumped."), .fn = ~ stringr::str_remove_all(.x, "clumped.")) %>%
      dplyr::right_join(hm_tab, by = c("outcome", "exposure", "SNP")) %>%
      dplyr::mutate(clump_keep = dplyr::if_else(is.na(clump_keep), FALSE, clump_keep))
  }

  hm_tab
}


#' Mendelian randomization analysis
mr_test <- function(hm_dat, min_exp_pval = 5e-8, min_snps = 5, prune = FALSE) {
  # Do prune
  if (prune) hm_dat <- TwoSampleMR::power_prune(hm_dat)

  # Filter out SNPs not suitable for MR test and keep significant exposure SNPs
  if ("clump_keep" %in% colnames(hm_dat)) hm_dat <- dplyr::filter(hm_dat, mr_keep, clump_keep, pval.exposure < min_exp_pval)
  else hm_dat <- dplyr::filter(hm_dat, mr_keep, pval.exposure < min_exp_pval)

  # A list to store results
  mr_res <- list(hm_dat = NULL, mr = NULL, mr_single_snp = NULL, mr_pleio = NULL, mr_hetero = NULL)
  if (nrow(hm_dat)) {
    suppressMessages({
      mr_res[["mr"]] <- TwoSampleMR::mr(hm_dat) # Basic MR analysis
      mr_res[["mr_single_snp"]] <- TwoSampleMR::mr_singlesnp(hm_dat) # Single SNP test
      mr_res[["mr_pleio"]] <- TwoSampleMR::mr_pleiotropy_test(hm_dat) # Horizontal pleiotropy
      mr_res[["mr_hetero"]] <- TwoSampleMR::mr_heterogeneity(hm_dat) # Check heterogeneity
    })
  }

  return(mr_res)
}


#' Plot MR for single SNPs
#'
#' @description Function to plot Figure x
plot_mrres <- function(mrpath, save_to = "./") {
  hm_tab <- data.table::fread(file.path(mrpath, "harmonized_data.csv"))
  mr_all_tab <- data.table::fread(file.path(mrpath, "mendelian_randomization_test.csv"))

  mr_persnp_tab <- data.table::fread(file.path(mrpath, "mendelian_randomization_test_persnp.csv"))
  mr_persnp_tab %>%
    dplyr::group_by(exposure, outcome) %>%
    dplyr::summarise(mr_figs = {
      per_exposure <- dplyr::cur_group()$exposure
      per_outcome <- dplyr::cur_group()$outcome

      per_mr_persnp <- dplyr::cur_data_all() %>% dplyr::filter(!duplicated(SNP)) %>% as.data.frame()

      all_mr_sig <- per_mr_persnp %>% dplyr::filter(stringr::str_detect(SNP, "^All ")) %>% dplyr::mutate(is_sig = p < 0.05) %>% dplyr::select(is_sig) %>% unlist() %>% sum()
      any_snp_sig <- per_mr_persnp %>% dplyr::filter(!stringr::str_detect(SNP, "^All ")) %>% dplyr::mutate(is_sig = p < 0.05) %>% dplyr::select(is_sig) %>% unlist() %>% sum()

      to_plot <- !is.na(all_mr_sig) && all_mr_sig == 2 && !is.na(any_snp_sig) && any_snp_sig > 0
      if (to_plot) {
        per_outcome_xx <- stringr::str_remove_all(per_outcome, "[|]+ id:") %>% stringr::str_replace_all("[ ]+", "_")

        # Forest plot
        fig_path <- file.path(save_to, paste("mr_forest", per_exposure, per_outcome_xx, "pdf", sep = "."))
        if (!file.exists(fig_path)) {
          p <- TwoSampleMR::mr_forest_plot(per_mr_persnp)
          ggsave(fig_path, plot = p[[1]], width = 7, height = 7)
        }

        # MR plot
        fig_path <- file.path(save_to, paste("mr_scatter", per_exposure, per_outcome_xx, "pdf", sep = "."))
        if (!file.exists(fig_path)) {
          per_mr_all <- mr_all_tab %>% dplyr::filter(exposure == per_exposure, outcome == per_outcome) %>% as.data.frame()
          per_hm_tab <- hm_tab %>% dplyr::filter(exposure == per_exposure, outcome == per_outcome) %>% as.data.frame()

          p <- TwoSampleMR::mr_scatter_plot(per_mr_all, per_hm_tab)
          ggsave(fig_path, plot = p[[1]], width = 7, height = 7)
        }
      }

      NULL
    })

  return(NULL)
}


#' Fetch data for coloc analysis from the harmonized dataset
#'
fetch_coloc_data <- function(hm_dat, what, warn_minp = 5e-5) {
  pos_cols <- c("SNP", "pos.exposure")
  src_cols <- c("beta", "se", "eaf", "samplesize")

  cols <- c(paste(src_cols, what, sep = "."), pos_cols)
  names(cols) <- c("beta", "varbeta", "MAF", "N", "snp", "position")

  tab <- hm_dat %>%
    dplyr::select(dplyr::all_of(cols)) %>%
    dplyr::filter(dplyr::across(dplyr::everything(), ~ !is.na(.))) %>%
    dplyr::mutate(type = "quant", varbeta = varbeta^2) %>%
    (function(dat) {
      beta <- dplyr::select(dat, "snp", "beta") %>% tibble::deframe()
      varbeta <- dplyr::select(dat, "snp", "varbeta") %>% tibble::deframe()
      maf <- dplyr::select(dat, "snp", "MAF") %>% tibble::deframe()

      list(
        beta = beta, varbeta = varbeta, snp = dat$snp, position = dat$position,
        type = dat$type[1], N = max(dat$N), MAF = maf
      )
    })
  coloc::check_dataset(tab, warn.minp = warn_minp)

  return(tab)
}


#' Colocalization anlaysis
colocalization_test <- function(hm_dat, min_snps = 10) {
  if (nrow(hm_dat) < min_snps) return(NA)

  qry <- fetch_coloc_data(hm_dat, "exposure")
  ref <- fetch_coloc_data(hm_dat, "outcome")

  if (TRUE) {
    sink("/dev/null")
    ccsum <- coloc::coloc.abf(qry, ref)$summary
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
  gene_tab[, c("gene_name", "type", "tag", "gene_type", "transcript_support_level")] %>%
    as.data.frame() %>%
    apply(1, function(e) {
      type <- e["type"]

      gene_type <- e["gene_type"] %in% c("protein_coding", "lncRNA")
      is_canonical <- e["tag"] %>% stringr::str_detect("Ensembl_canonical")
      is_tsl_1 <- e["transcript_support_level"] %in% c(1, "1")
      if (type == "gene") keep <- gene_type
      else if (type %in% tar_type) keep <- gene_type && is_canonical # && is_tsl_1
      else return(NULL)

      if (keep) data.frame(chrom = e$seqnames, start = e$start, end = e$end, strand = e$strand, width = e$width, gene = e$gene_name, type = e$type)
      else return(NULL)
    }) %>%
    Reduce(rbind, .) %>%
    (function(e) if (!is.null(e)) dplyr::mutate(e, highlight = gene %in% highlight))
}


#' Plot colocalization
#'
#' @description Function to plot Figure x
plot_coloc <- function(load_path, min_h4 = 0.5, min_h3 = 0.25, override = FALSE, expand = 0, shift = 0, save_to = "./", ...) {
  # load_path <- "../outputs/pseudo_bulk/outcomes/normal/normal_B"
  cc_path <- file.path(load_path, "colocalization_test.csv")
  hm_path <- file.path(load_path, "harmonized_data.csv")
  if (!(file.exists(hm_path) && file.exists(cc_path))) stop("Missing colocalization_test.csv or harmonized_data.csv!")

  kwargs <- list(...)
  if (is.null(kwargs$bfile) || is.null(kwargs$plink_bin)) stop("You have to specify the path to plink and bfiles.")
  if (is.null(kwargs$gffpath)) stop("You have to specify the path to genomic feature files.")

  hm_tab <- data.table::fread(hm_path)
  cc_tab <- data.table::fread(cc_path) %>% dplyr::filter(H4 >= min_h4 | H3 >= min_h3)
  cmb_tab <- dplyr::inner_join(cc_tab, hm_tab, by = c("exposure", "outcome"))
  height_dict <- c("gene" = 0.3, "exon" = 0.75, "start_codon" = 1.1, "stop_codon" = 1.1, "five_prime_UTR" = 0.5, "three_prime_UTR" = 0.5)

  cmb_tab %>%
    dplyr::group_by(exposure, outcome) %>%
    dplyr::summarise(fig_path = {
      per_exposure <- dplyr::cur_group()$exposure
      per_outcome <- dplyr::cur_group()$outcome

      per_outcome_xx <- stringr::str_remove_all(per_outcome, "[|]+ id:") %>% stringr::str_replace_all("[ ]+", "_")
      fig_path <- file.path(save_to, paste("coloc", per_exposure, per_outcome_xx, "pdf", sep = "."))
      if (!file.exists(fig_path) || override) {
        per_pair <- dplyr::cur_data_all() %>% dplyr::select(exposure, outcome, SNP, chr, pos, pval.outcome, pval.exposure)
        topsnp_id <- dplyr::slice_min(per_pair, pval.exposure)$SNP[1]
        ldmat <- ieugwasr::ld_matrix_local(per_pair$SNP, kwargs$bfile, kwargs$plink_bin, FALSE) ^ 2

        tsnp_ld_vec <- ldmat[, colnames(ldmat) == topsnp_id]
        if (length(tsnp_ld_vec) > 0) {
          pair_tab <- per_pair %>%
            tidyr::pivot_longer(cols = c("pval.outcome", "pval.exposure")) %>%
            dplyr::mutate(
              ld2topsnp = tsnp_ld_vec[SNP], ld2topsnp = dplyr::if_else(is.na(ld2topsnp), 0, ld2topsnp), value = -log10(value),
              TopSNP = SNP == topsnp_id, name = dplyr::if_else(name == "pval.exposure", exposure, outcome)
            )
          gf_chr <- unique(per_pair$chr)[1]
          gf_start <- as.integer(max(min(per_pair$pos) + (shift - expand) * 1000, 0))
          gf_end <- as.integer(max(per_pair$pos) + (shift + expand) * 1000)
          x_lims <- c(gf_start - 10, gf_end + 10)

          gftab <- fetch_genomic_features(gf_chr, gf_start, gf_end, kwargs$gffpath, TRUE, per_exposure)
          if (!is.null(gftab)) {
            gftab <- dplyr::mutate(gftab, height = height_dict[type]) %>%
              dplyr::mutate(
                start = dplyr::if_else(start <= gf_start, as.integer(gf_start), as.integer(start)),
                width = dplyr::if_else((start + width) >= gf_end, as.integer(gf_end - start), as.integer(width))
              )

            y_coord_dict <- dplyr::select(gftab, gene) %>% dplyr::distinct() %>% dplyr::mutate(y_coord = row_number()) %>% tibble::deframe()
            gftab <- dplyr::mutate(gftab, gene_pos = y_coord_dict[gene] * 1.5)

            gf_plot <- ggplot() +
              geom_rect(aes(xmin = start, xmax = end, ymin = gene_pos - height, ymax = gene_pos + height, fill = strand),
                        data = gftab, position = "identity", stat = "identity") +
              geom_text(aes(start, gene_pos, label = gene, color = highlight), data = dplyr::filter(gftab, type == "gene"), hjust = 1.1, size = 2) +
              scale_fill_manual(name = NULL, values = c("-" = "darkblue", "+" = "darkred")) +
              scale_color_manual(name = NULL, values = c("TRUE" = "red", "FALSE" = "black")) +
              lims(x = x_lims) + theme_classic() + labs(x = "Position (base-pair)", y = NULL) +
              theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_blank())
          }

          topsnp_info <- dplyr::filter(pair_tab, SNP == topsnp_id)
          cc_plot <- ggplot() +
            geom_point(aes(x = pos, y = value, color = ld2topsnp), pair_tab, alpha = 0.85) +
            geom_point(aes(x = pos, y = value), topsnp_info, size = 3, shape = 23, fill = "purple") +
            geom_text_repel(aes(x = pos, y = value, label = SNP), topsnp_info, min.segment.length = 0) +
            geom_hline(yintercept = -log10(5e-8), linetype = "dashed") +
            facet_wrap(~name, ncol = 1, scales = "free_y", strip.position = "right") +
            scale_colour_gradientn(colours = c("darkblue", "blue", "orange", "red", "darkred"), name = quote("LD (" ~ R^2 ~ ")")) +
            labs(x = NULL, y = quote(~-Log[10] ~ "(p-value)")) +
            lims(y = c(0, NA), x = x_lims) + theme_bw() + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

          if (!is.null(gftab)) cmb_plot <- cc_plot / gf_plot + plot_layout(heights = c(4, 1))
          else cmb_plot <- cc_plot

          ggsave(fig_path, cmb_plot, width = 7, height = 5)
        }
      }
      NULL
    })
}


#' Plot multi-omics data genomic tracks
draw_multracks <- function() {
  return(0)
}


#
## Main steps.
#
proj_dir <- "/home/zzhang/Documents/projects/wp_bcg_eqtl"
in_dir <- file.path(proj_dir, "outputs/pseudo_bulk")

# Two model used in the analysis.
mode_vec <- c("normal", "interaction")

# Main cell types.
cell_type_vec <- c("Monocytes", "CD4T", "CD8T", "NK", "B") # , "pDC", "mDC")

# All comparisons used in the analysis.
condition_vec <- c("T0_LPS.vs.T0_RPMI", "T3m_LPS.vs.T0_RPMI", "T3m_LPS.vs.T3m_RPMI", "T3m_RPMI.vs.T0_RPMI")

# Public GWAS list
pub_gwas_list <- list(
  # Cancers
  "ieu-b-4965" = c("2021_ColorectalCancer_ieu-b-4965.vcf.gz", NA),
  "ieu-b-4809" = c("2021_ProstateCancer_ieu-b-4809.vcf.gz", NA),
  "ieu-b-4874" = c("2021_BladderCancer_ieu-b-4874.vcf.gz", NA),
  "ieu-a-1082" = c("2013_ThyroidCancer_ieu-a-1082.vcf.gz", NA),
  "ieu-b-4963" = c("2017_OvarianCancer_ieu-b-4963.vcf.gz", NA),
  "ieu-a-966" = c("2014_LungCancer_ieu-a-966.vcf.gz", NA),

  # Autoimmune disease
  # "ebi-a-GCST010681" = c("2022_Type1Diabetes_ebi-a-GCST010681.vcf.gz", NA), No GWAS VCF available, TODO: format available SS into GWAS VCF format.
  "ukb-d-M13_RHEUMA" = c("2018_RheumatoidArthritis_ukb-d-M13_RHEUMA.vcf.gz", NA),
  "ieu-a-31" = c("2015_InflammatoryBowelDisease_ieu-a-31.vcf.gz", NA),
  "ukb-b-17670" = c("2019_MultipleSclerosis_ukb-b-17670.vcf.gz", NA),
  "ieu-a-996" = c("2014_AtopicDermatitis_ieu-a-996.vcf.gz", NA),
  "ieu-a-32" = c("2015_UlcerativeColitis_ieu-a-32.vcf.gz", NA),
  "ieu-a-30" = c("2015_CrohnsDisease_ieu-a-30.vcf.gz", NA),
  "ukb-a-100" = c("2017_Psoriasis_ukb-a-100.vcf.gz", NA),

  # Infectious disease
  "ebi-a-GCST010780" = c("2020_COVID19Release4_ebi-a-GCST010780.vcf.gz", NA),

  # Brain disorders
  "ieu-b-5067" = c("2022_AlzheimerDiseases_ieu-b-5067.vcf.gz", NA),
  "ieu-b-42" = c("2014_Schizophrenia_ieu-b-42.vcf.gz", NA),

  # Other genetic-related disease
  "ebi-a-GCST006867" = c("2018_Type2Diabetes_ebi-a-GCST006867.vcf.gz", NA),
  "ukb-d-J10_ASTHMA" = c("2018_Asthma_ukb-d-J10_ASTHMA.vcf.gz", NA),
  "ieu-a-7" = c("2013_CoronaryHeartDisease_ieu-a-7.vcf.gz", NA),
  "ukb-a-107" = c("2017_GoutDisease_ukb-a-107.vcf.gz", NA),

  # Others traits
  "ieu-a-89" = c("2014_Height_ieu-a-89.vcf.gz", NA),
  "ieu-b-40" = c("2018_BodyMassIndex_ieu-b-40.vcf.gz", NA),
  "ieu-b-111" = c("2020_Triglycerides_ieu-b-111.vcf.gz", 441016),
  "ieu-b-109" = c("2020_HDLCholesterol_ieu-b-109.vcf.gz", 403943),
  "ieu-b-110" = c("2020_LDLCholesterol_ieu-b-110.vcf.gz", 440546)
) %>%
lapply(function(e) c(file.path(proj_dir, "inputs/public_gwas", e[1]), e[2]))

# Reference genotypes panel for local clumping
eur_1kg_geno <- file.path(proj_dir, "inputs/reference/genotypes/EUR")

# Allele frequecies of effect allele used in sceQTL mapping.
afdb_path <- file.path(proj_dir, "inputs/genotypes/300BCG_sub40_imp_hg38_ref_allele_frequency.csv")
afdb <- data.table::fread(afdb_path) %>% tibble::deframe()

# GARFIELD database
# gfout_dir <- file.path(proj_dir, "outputs/pseudo_bulk/outcomes/garfield")
# gfdb_dir <- file.path(proj_dir, "outputs/garfield-data-grch38")

# List to store results.
qtltab_list <- list()
gfres_list <- list()
mrres_list <- list()
ccres_list <- list()

mode <- "interaction"
for (mode in mode_vec) {
  for (cell_type in cell_type_vec) {
    if (mode == "normal") {
      root_dir <- file.path(in_dir, mode, cell_type)
      run_id <- paste(mode, cell_type, sep = "_")

      if (file.exists(file.path(root_dir, "top_qtl_results_all_FDR.txt"))) {
        cat("[I]: Loading results from", root_dir, "...\n")
        qtltab_list[[run_id]] <- load_eqtl_tab(root_dir, 0.1, other_info = c("cell_type" = cell_type))
      } else {
        cat("[W]: Not top QTL results available for", run_id, "Skipping ...\n")
      }
    } else {
      for (condition in condition_vec) {
        root_dir <- file.path(in_dir, mode, cell_type, condition)
        run_id <- paste(mode, cell_type, condition, sep = "_")

        if (file.exists(file.path(root_dir, "top_qtl_results_all_FDR.txt"))) {
          cat("[I]: Loading results from", root_dir, "...\n")
          qtltab_list[[run_id]] <- load_eqtl_tab(root_dir, 0.1, other_info = c("cell_type" = cell_type, "condition" = condition))
        } else {
          cat("[W]: Not top QTL results available for", run_id, "Skipping ...\n")
        }
      }
    }
  }

  save_to <- file.path(proj_dir, "outputs/pseudo_bulk/outcomes", mode)

  # Number of shared eGenes
  cat("[I]: Drawing upset plot", mode, "eQTL ...\n")
  draw_upsetplot(qtltab_list, save_to, "QTL", 9, 5)

  # Concordance between two runs.
  .tmp <- combn(names(qtltab_list), 2) %>%
    t() %>%
    as.data.frame() %>%
    apply(1, function(e) {
      x_ct <- e[1]
      y_ct <- e[2]

      per_fig <- file.path(save_to, "concordance", paste0(y_ct, ".vs.", x_ct, ".png"))
      per_lab <- str_remove(e, paste0(mode, "_")) %>% str_replace("_", "@")

      if (!file.exists(per_fig)) {
        cat("[I]: Drawing dotplot comparing", x_ct, "to", y_ct, "...\n")
        draw_dotplot(qtltab_list[[x_ct]], qtltab_list[[y_ct]], labels = per_lab, joinby = c("snp_id", "ensembl_gene_id"), save_to = per_fig)
      }
    })

  # Function enrichment, COLOC and MR
  .tmp <- lapply(names(qtltab_list), function(pr_run_id) {
    tmp_tab <- qtltab_list[[pr_run_id]] %>%
      (function(dat) {
        top_egenes <- dplyr::filter(dat, !is.na(global_corrected_pValue)) %>% dplyr::select(gene_name) %>% unlist() %>% as.vector()
        dplyr::filter(dat, gene_name %in% top_egenes, p_value < 0.1, !stringr::str_starts(gene_name, "HLA"))
      })

    pr_save_to <- file.path(save_to, pr_run_id)
    if (!dir.exists(pr_save_to)) dir.create(pr_save_to, recursive = TRUE)

    hm_save_to <- file.path(pr_save_to, "harmonized_data.csv")
    if (!file.exists(hm_save_to)) {
      cat("[I]: Harmonizing variants ...\n")
      hm_dat <- harmonize_vars(tmp_tab, pub_gwas_list, afdb, clump_bfile = eur_1kg_geno, plink_bin = "plink-quiet")
      data.table::fwrite(hm_dat, hm_save_to, verbose = FALSE, showProgress = FALSE)
    } else {
      hm_dat <- data.table::fread(hm_save_to)
    }

    mr_save_to <- file.path(pr_save_to, "mendelian_randomization_test.csv")
    mr_ps_save_to <- file.path(pr_save_to, "mendelian_randomization_test_persnp.csv")
    if (!file.exists(mr_save_to) || !file.exists(mr_ps_save_to)) {
      cat("[I]: Estimating causality by TwoSampleMR ...\n")
      mrres_list[[pr_run_id]] <- mr_test(hm_dat, min_exp_pval = 5e-6)

      data.table::fwrite(mrres_list[[pr_run_id]]$mr, mr_save_to, verbose = FALSE, showProgress = FALSE)
      data.table::fwrite(mrres_list[[pr_run_id]]$mr_single_snp, mr_ps_save_to, verbose = FALSE, showProgress = FALSE)
    }

    cc_save_to <- file.path(pr_save_to, "colocalization_test.csv")
    if (!file.exists(cc_save_to)) {
      cat("[I]: Estimating co-localization by coloc ...\n")
      ccres_list[[pr_run_id]] <- hm_dat %>%
        dplyr::group_by(exposure, outcome) %>%
        dplyr::filter(dplyr::if_all(.fns = ~ !is.na(.x))) %>%
        dplyr::summarise(cc_res = colocalization_test(cur_data_all())) %>%
        data.table::as.data.table() %>%
        dplyr::rename_with(dplyr::starts_with("cc_res."), .fn = ~ stringr::str_remove(.x, "cc_res.")) %>%
        dplyr::filter(dplyr::if_all(dplyr::starts_with("H"), .fns = ~ !is.na(.x)))

      data.table::fwrite(ccres_list[[pr_run_id]], cc_save_to, verbose = FALSE, showProgress = FALSE)
    }

    NULL
  })

  # Gene enrichment analysis
  .tmp <- lapply(names(qtltab_list), function(pr_run_id) {
    tmp_tab <- qtltab_list[[pr_run_id]] %>%
      (function(dat) {
        top_egenes <- dplyr::filter(dat, !is.na(global_corrected_pValue)) %>% dplyr::select(gene_name) %>% unlist() %>% as.vector()
        dplyr::filter(dat, gene_name %in% top_egenes, p_value < 0.1, !stringr::str_starts(gene_name, "HLA"))
      })

    pr_save_to <- file.path(save_to, pr_run_id)
    if (!dir.exists(pr_save_to)) dir.create(pr_save_to, recursive = TRUE)

    gene_info <- tmp_tab %>%
      dplyr::filter(!is.na(global_corrected_pValue)) %>%
      dplyr::select(gene_name, beta) %>%
      tibble::deframe() %>%
      sort(decreasing = TRUE) %>%
      (function(vec) {
        gene_tab <- clusterProfiler::bitr(names(vec), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") %>%
          dplyr::mutate(Beta = vec[SYMBOL])
      })

    gsea_save_to <- file.path(pr_save_to, "egene_function_enrichment_GSEA.csv")
    if (!file.exists(gsea_save_to)) {
      cat("[I]: Performing gene set enrichment analysis ...\n")

      # Enrichment analysis GSEA
      gsea_res <- gene_info %>%
        dplyr::select(ENTREZID, Beta) %>%
        dplyr::mutate(ENTREZID = as.integer(ENTREZID)) %>%
        tibble::deframe() %>%
        (function(vec) {
          msigdbr::msigdbr(species = "Homo sapiens") %>%
            dplyr::select(gs_name, entrez_gene) %>%
            clusterProfiler::GSEA(vec, TERM2GENE = ., pvalueCutoff = 1.0, minGSSize = 20, maxGSSize = 5000)
        })
        
      # Gene set enrichment analysis for customized gene sets.
      # gene_list <- gene_info %>%
      #   dplyr::select(dplyr::all_of(c("SYMBOL", "Beta"))) %>%
      #   tibble::deframe()

      # bg_tab <- {
      #   bgfp <- list(
      #     sepsis_nm = "Sepsis/Monocytes/sepsis_nm_degs.txt",
      #     sepsis_valerie = "Sepsis/Monocytes/sepsis_valerie_degs.txt",
      #     covid19_berlin = "Covid19/Monocytes/berlin1_severe_mild_degs.txt",
      #     covid19_bohn = "Covid19/Monocytes/covid2_severe_mild_degs.txt",
      #     train_immunity = "train_immunity/300BCG_TrainedImmunity_DEGs.txt",
      #     bcg_effect_mono = "bcg_effect/300BCG_BCG_DEGs_Monocytes.txt",
      #     bcg_effect_cd8t = "bcg_effect/300BCG_BCG_DEGs_CD8T.txt",
      #     bcg_effect_cd4t = "bcg_effect/300BCG_BCG_DEGs_CD4T.txt"
      #   )

      #   lapply(names(bgfp), function(e) {
      #     data.table::fread(file.path(proj_dir, "inputs/DEGs", bgfp[[e]])) %>%
      #       dplyr::filter(abs(avg_log2FC) > 0.15) %>%
      #       dplyr::select(p_val, avg_log2FC, pct.1, pct.2, p_val_adj, gene) %>%
      #       dplyr::mutate(direction = dplyr::if_else(avg_log2FC < 0, "Down", "Up"), trait = paste(direction, e, sep = "@"))
      #   }) %>%
      #   Reduce(rbind, .)
      # }

      # cm_gsea <- clusterProfiler::GSEA(gene_list, TERM2GENE = bg_tab[, c("trait", "gene")], minGSSize = 15, maxGSSize = 2000, pvalueCutoff = 1.0)
      # data.table::fwrite(gsea_res@result, gsea_save_to, verbose = FALSE, showProgress = FALSE)
    }

    NULL
  })
}


for (mode in mode_vec) {
  for (cell_type in cell_type_vec) {
    if (mode == "normal") {
      run_id <- paste(mode, cell_type, sep = "_")
      root_dir <- file.path(in_dir, "outcomes", mode, run_id)

      coloc_save_to <- file.path(root_dir, "coloc_plots")
      if (!dir.exists(coloc_save_to)) dir.create(coloc_save_to)
      plot_coloc(root_dir, save_to = coloc_save_to)

      mr_save_to <- file.path(root_dir, "mr_plots")
      if (!dir.exists(mr_save_to)) dir.create(mr_save_to)
      plot_mrres(root_dir, save_to = mr_save_to)
    } else {
      for (condition in condition_vec) {
        run_id <- paste(mode, cell_type, condition, sep = "_")
        root_dir <- file.path(in_dir, "outcomes", mode, run_id)

        coloc_save_to <- file.path(root_dir, "coloc_plots")
        if (!dir.exists(coloc_save_to)) dir.create(coloc_save_to)
        plot_coloc(root_dir, save_to = coloc_save_to)

        mr_save_to <- file.path(root_dir, "mr_plots")
        if (!dir.exists(mr_save_to)) dir.create(mr_save_to)
        plot_mrres(root_dir, save_to = mr_save_to)
      }
    }
  }
}
