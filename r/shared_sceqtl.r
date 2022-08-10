#!/usr/bin/env Rscript
library(ashr)
library(mashr)
library(tidyverse)
library(data.table)

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


#' Prepare data form mashr
prep_mashr_data <- function(tab_list) {
  if (is.null(names(tab_list))) stop("The input tab_list should be a named list!")

  mashr_data <- Reduce(function(x, y) {
    if (is.vector(x) && is.character(x)) {
      x_tab <- dplyr::select(tab_list[[x]], beta, beta_se, QTL)

      nm_vec <- c("beta", "QTL")
      names(nm_vec) <- c(x, "QTL")
      x_bhat <- x_tab %>% dplyr::select(beta, QTL) %>% dplyr::rename(nm_vec)

      nm_vec <- c("beta_se", "QTL")
      names(nm_vec) <- c(x, "QTL")
      x_shat <- x_tab %>% dplyr::select(beta_se, QTL) %>% dplyr::rename(nm_vec)

      x <- list(Bhat = x_bhat, Shat = x_shat)
    }

    y_tab <- dplyr::select(tab_list[[y]], beta, beta_se, QTL)

    nm_vec <- c("beta", "QTL")
    names(nm_vec) <- c(y, "QTL")
    y_bhat <- y_tab %>% dplyr::select(beta, QTL) %>% dplyr::rename(nm_vec)

    nm_vec <- c("beta_se", "QTL")
    names(nm_vec) <- c(y, "QTL")
    y_shat <- y_tab %>% dplyr::select(beta_se, QTL) %>% dplyr::rename(nm_vec)

    x$Bhat <- dplyr::inner_join(x$Bhat, y_bhat, by = "QTL")
    x$Shat <- dplyr::inner_join(x$Shat, y_shat, by = "QTL")

    return(x)
  }, names(tab_list))

  mashr_data %>% (function(mat) {
       rownames(mat$Bhat) <- mat$Bhat$QTL
       mat$Bhat <- dplyr::select(mat$Bhat, -QTL) %>% as.matrix()

       rownames(mat$Shat) <- mat$Shat$QTL
       mat$Shat <- dplyr::select(mat$Shat, -QTL) %>% as.matrix()

       mat
    })
}



#' Perform mashr analysis
#'
#'@description The function performs mashr analysis. The sceQTL results are huge, therefore, here we using random
#' selected effects to estimate null distribution and using the top-eqtl to estimate the data-driven covariance matrices.
#' The function implementation refers to mashr Vignettes (https://stephenslab.github.io/mashr/articles/eQTL_outline.html)
exec_mashr <- function(mdata, sig_pval = 0.05, n_rand = 20000, n_rows = NULL) {
  if (!is.null(n_rows)) 
    mdata <- mash_set_data(mdata$Bhat[1:n_rows, ], mdata$Shat[1:n_rows, ])
  else
    mdata <- mash_set_data(mdata$Bhat, mdata$Shat)

  # Prepare random dataset and strong effect dataset.
  mash_pair <- mash_1by1(mdata)
  strong_subset <- get_significant_results(mash_pair, sig_pval)
  random_subset <- sample(1:nrow(mdata$Bhat), n_rand)

  # Estimate correlation structure.
  temp_data <- mash_set_data(mdata$Bhat[random_subset, ], mdata$Shat[random_subset, ])
  vhat <- estimate_null_correlation_simple(temp_data) # Accounting correlations among conditions/measurements.
  rm(temp_data)

  data_random <- mash_set_data(mdata$Bhat[random_subset, ], mdata$Shat[random_subset, ], V = vhat)
  data_strong <- mash_set_data(mdata$Bhat[strong_subset, ], mdata$Shat[strong_subset, ], V = vhat)

  # Calculate data-driven covariances.
  u_pca <- cov_pca(data_strong, 5)
  u_ed <- cov_ed(data_strong, u_pca)

  # Fit the model.
  u_c <- cov_canonical(data_random)
  model_random <- mash(data_random, Ulist = c(u_ed, u_c), outputlevel = 1)

  # Compute posterior summaries.
  model_strong <- mash(data_strong, g = get_fitted_g(model_random), fixg = TRUE)

  # Collect results.
  list(alter_model = model_strong,
       null_model = model_random,
       locale_false_sign_rate = get_lfsr(model_strong),
       posterior_mean = get_pm(model_strong),
       posterior_std_dev = get_psd(model_strong),
       pos_prob = get_pp(model_strong),
       neg_prob = get_np(model_strong),
       shared_signals_by_magnitude = get_pairwise_sharing(model_strong),
       shared_signals_by_sign = get_pairwise_sharing(model_strong, factor = 0)
  )
}


#' Visualization
plot_mashr <- function(mm, share_factor = c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8), override = FALSE, save_to = "./") {
  # Propotion of (significant) signals shared by each pair of conditions (?)
  lapply(share_factor, function(f) {
    plot_path <- file.path(save_to, paste0("share_factor_", f, ".pdf"))
    mat <- get_pairwise_sharing(mm, factor = f) %>%
      (function(mat) {
         new_nm <- colnames(mat) %>% stringr::str_remove_all("normal_|interaction_")
         colnames(mat) <- new_nm
         rownames(mat) <- new_nm
         mat
      })

    if (!file.exists(plot_path) || override) {
      pdf(plot_path)
      pheatmap::pheatmap(mat, main = "Proportion of shared signals")
      dev.off()
    }

    mat
  })
}


#' Main steps
proj_dir <- "~/Documents/projects/wp_bcg_eqtl"
mode <- "pseudo_bulk"
model_vec <- c("normal", "interaction")
celltype_vec <- c("Monocytes", "CD4T", "CD8T", "NK", "B")
override <- TRUE

qtl_tab_list <- list()
for (model in model_vec) {
  for (celltype in celltype_vec) {
    qtl_tab_path <- file.path(proj_dir, "outputs", mode, model, celltype)
    run_id <- paste(model, celltype, sep = "_")
    qtl_tab_list[[run_id]] <- load_eqtl_tab(qtl_tab_path, fdr = 0.1)
  }

  out_dir <- file.path(proj_dir, "outputs", mode, "outcomes", model, "mashr")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  mashr_res_path <- file.path(out_dir, "mashr_res.Rdata")
  if (file.exists(mashr_res_path)) {
    load(mashr_res_path)
  } else {
    mashr_data <- prep_mashr_data(qtl_tab_list)
    mashr_res <- exec_mashr(mashr_data, sig_pval = 0.01)
    save(mashr_res, file = mashr_res_path)
  }

  plot_mashr(mashr_res$alter_model, seq(1, 9) / 10, save_to = out_dir)
  break
}
