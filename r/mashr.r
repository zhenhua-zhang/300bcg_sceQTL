#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, datatable.verbose = FALSE, datatable.showProgress = FALSE)
suppressPackageStartupMessages({
  # MASHR analysis
  library(ashr)
  library(mashr)

  # Data manipulation
  library(magrittr)
  library(tidyverse)
  library(data.table)

  # Plotting
  library(ComplexHeatmap)
  library(circlize)
})


#' Load eQTL summary statistics
#'
#'@description A utility function to load summary statistics from the disk
load_eqtl_tab <- function(indir, fdr = 0.05, max_p_val = NULL, gene_bl = NULL) {
  # Gene black list by default, if the parameter is NULL
  if (is.null(gene_bl)) {
    gene_bl <- c(
      "RPLP0", "RPLP1", "RPLP2", "RPL3", "RPL3L", "RPL5", "RPL6", "RPL7", "RPL7A", "RPL7L1", "RPL8", "RPL9", "RPL10",
      "RPL10A", "RPL10L", "RPL11", "RPL12", "RPL13", "RPL13A", "RPL14", "RPL15", "RPL17", "RPL18", "RPL18A", "RPL19",
      "RPL21", "RPL22", "RPL23", "RPL23A", "RPL24", "RPL26", "RPL26L1", "RPL27", "RPL27A", "RPL28", "RPL29", "RPL30",
      "RPL31", "RPL32", "RPL34", "RPL35", "RPL35A", "RPL36", "RPL36A", "RPL36AL", "RPL37", "RPL37A", "RPL38", "RPL39",
      "RPL39L", "UBA52", "RPL41",
      "RPSA", "RPS2", "RPS3", "RPS3A", "RPS4X", "RPS4Y1", "RPS4Y2", "RPS5", "RPS6", "RPS7", "RPS8", "RPS9", "RPS10",
      "RPS11", "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A", "RPS16", "RPS17", "RPS18", "RPS19", "RPS20", "RPS21",
      "RPS23", "RPS24", "RPS25", "RPS26", "RPS27", "RPS27A", "RPS27L", "RPS28", "RPS29", "FAU",
      "MRPL1", "MRPL2", "MRPL3", "MRPL4", "MRPL9", "MRPL10", "MRPL11", "MRPL12", "MRPL13", "MRPL14", "MRPL15", "MRPL16",
      "MRPL17", "MRPL18", "MRPL19", "MRPL20", "MRPL21", "MRPL22", "MRPL23", "MRPL24", "MRPL27", "MRPL28", "MRPL30",
      "MRPL32", "MRPL33", "MRPL34", "MRPL35", "MRPL36", "MRPL37", "MRPL38", "MRPL39", "MRPL40", "MRPL41", "MRPL42",
      "MRPL43", "MRPL44", "MRPL45", "MRPL46", "MRPL47", "MRPL48", "MRPL49", "MRPL50", "MRPL51", "MRPL52", "MRPL53",
      "MRPL54", "MRPL55", "MRPL58", "MRPL57", "GADD45GIP1", "MRPS30", "MRPS18A", "MRPS2", "MRPS5", "MRPS6", "MRPS7",
      "MRPS9", "MRPS10", "MRPS11", "MRPS12", "MRPS14", "MRPS15", "MRPS16", "MRPS17", "MRPS18C", "MRPS21", "MRPS22",
      "MRPS23", "MRPS24", "MRPS25", "MRPS26", "MRPS27", "MRPS28", "DAP3", "MRPS31", "MRPS33", "MRPS34", "MRPS35",
      "CHCHD1", "AURKAIP1", "PTCD3", "MRPS18B",
      "HLA-A", "HLA-B", "HLA-C", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1", "HLA-DPA2", "HLA-DPA3",
      "HLA-DPB1", "HLA-DPB2", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DQB3", "HLA-DRA", "HLA-DRB1",
      "HLA-DRB2", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-DRB6", "HLA-DRB7", "HLA-DRB8", "HLA-DRB9", "HLA-E",
      "HLA-F", "HLA-G", "HLA-H", "HLA-J", "HLA-K", "HLA-L", "HLA-N", "HLA-P", "HLA-S", "HLA-T", "HLA-U", "HLA-V",
      "HLA-W", "HLA-X", "HLA-Y", "HLA-Z", "HLA-F-AS1"
    )
  }

  # All associations
  all_qtl_tab <- file.path(indir, "qtl_results_all.txt") %>%
    data.table::fread(showProgress = FALSE) %>%
    dplyr::filter(!feature_id %in% gene_bl)

  # eGene, defined using the lowest p-value per feature if the q-value is <= fdr
  egene_tab <- dplyr::group_by(all_qtl_tab, feature_id) %>%
    dplyr::arrange(empirical_feature_p_value) %>%
    dplyr::slice(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(q_value = qvalue::qvalue_truncp(empirical_feature_p_value)$qvalues) %>%
    dplyr::filter(q_value <= fdr) %>%
    dplyr::mutate(QTL = paste(feature_id, snp_id, sep = "-"))

  # eGene q-values
  egene_qval_vec <- dplyr::pull(egene_tab, q_value, QTL)

  # Max p-value to define a significant eVariant
  max_p_by_fdr <- dplyr::select(egene_tab, q_value, empirical_feature_p_value) %>%
    dplyr::slice_max(order_by = q_value) %>%
    dplyr::pull(empirical_feature_p_value) %>%
    head(1)

  # eVariant. A much looser definition can be applied by setting max_p_val
  if (is.null(max_p_val)) {
    evariant_tab <- dplyr::filter(all_qtl_tab, empirical_feature_p_value <= max_p_by_fdr)
  } else {
    cat("[W]: Variants with p-value <=", max_p_val, "are also retained as eVariants.\n")
    evariant_tab <- dplyr::filter(all_qtl_tab, empirical_feature_p_value <= max_p_by_fdr | p_value <= max_p_val)
  }

  dplyr::mutate(evariant_tab, QTL = paste(feature_id, snp_id, sep = "-"), q_value = egene_qval_vec[QTL])
}


#' Prepare data form mashr
prep_mashr_data <- function(tab_list) {
  if (is.null(names(tab_list))) stop("The input tab_list should be a named list!")

  mashr_data <- Reduce(function(x, y) {
    if (is.vector(x) && is.character(x)) {
      x_tab <- dplyr::select(tab_list[[x]], beta, beta_se, QTL)

      nm_vec <- c("beta", "QTL")
      names(nm_vec) <- c(x, "QTL")
      x_bhat <- x_tab %>% dplyr::select(beta, QTL) %>% dplyr::rename(dplyr::all_of(nm_vec))

      nm_vec <- c("beta_se", "QTL")
      names(nm_vec) <- c(x, "QTL")
      x_shat <- x_tab %>% dplyr::select(beta_se, QTL) %>% dplyr::rename(dplyr::all_of(nm_vec))

      x <- list(Bhat = x_bhat, Shat = x_shat)
    }

    y_tab <- dplyr::select(tab_list[[y]], beta, beta_se, QTL)

    nm_vec <- c("beta", "QTL")
    names(nm_vec) <- c(y, "QTL")
    y_bhat <- y_tab %>% dplyr::select(beta, QTL) %>% dplyr::rename(dplyr::all_of(nm_vec))

    nm_vec <- c("beta_se", "QTL")
    names(nm_vec) <- c(y, "QTL")
    y_shat <- y_tab %>% dplyr::select(beta_se, QTL) %>% dplyr::rename(dplyr::all_of(nm_vec))

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
exec_mashr <- function(mdata, sig_pval = 5e-2, n_pcs = 5, null_zth = 3, factor = 0.75, n_rand = NULL, n_rows = NULL) {
  if (!is.null(n_rows)) {
    mdata <- mash_set_data(mdata$Bhat[1:n_rows, ], mdata$Shat[1:n_rows, ])
  } else {
    mdata <- mash_set_data(mdata$Bhat, mdata$Shat)
  }

  model_random <- NULL
  if (is.null(n_rand)) {
    cat("[I]: Using the combination data-driven and canonical variace matrices.\n")
    data_strong <- mash_1by1(mdata) %>% get_significant_results(sig_pval) # Obtain significant results.
    u_pca <- cov_pca(mdata, n_pcs, subset = data_strong) # Obtain initial data-driven covariance matrices
    u_ed <- cov_ed(mdata, u_pca, subset = data_strong) # Apply Extreme Deconvolution
    u_c = cov_canonical(mdata) # Obtain canonical covariance matrices
    model_strong <- mash(mdata, c(u_ed, u_c)) # Compute posterior summaries.
  } else {
    cat("[I]: Using the random-strong strategy to estimate null distribution.\n")
    n_rand <- as.integer(min(n_rand, nrow(mdata$Bhat) * 0.5))

    # Prepare random dataset and strong effect dataset.
    strong_subset <- mash_1by1(mdata) %>% get_significant_results(sig_pval)
    random_subset <- sample(1:nrow(mdata$Bhat), n_rand)

    # Estimate correlation structure.
    temp_data <- mash_set_data(mdata$Bhat[random_subset, ], mdata$Shat[random_subset, ])
    vhat <- estimate_null_correlation_simple(temp_data, z_thresh = null_zth) # Accounting correlations among conditions/measurements.

    # Obtain random and strong data
    data_random <- mash_set_data(mdata$Bhat[random_subset, ], mdata$Shat[random_subset, ], V = vhat)
    data_strong <- mash_set_data(mdata$Bhat[strong_subset, ], mdata$Shat[strong_subset, ], V = vhat)

    # Calculate data-driven covariances.
    u_pca <- cov_pca(data_strong, n_pcs)
    u_ed <- cov_ed(data_strong, u_pca)

    # Fit the model.
    u_c <- cov_canonical(data_random)
    model_random <- mash(data_random, c(u_ed, u_c), outputlevel = 1)

    # Compute posterior summaries.
    model_strong <- mash(data_strong, g = get_fitted_g(model_random), fixg = TRUE)
  }

  # Collect results.
  list(alter_model = model_strong,
       null_model = model_random,
       locale_false_sign_rate = get_lfsr(model_strong),
       posterior_mean = get_pm(model_strong),
       posterior_std_dev = get_psd(model_strong),
       pos_prob = get_pp(model_strong),
       neg_prob = get_np(model_strong),
       shared_signals_by_magnitude = get_pairwise_sharing(model_strong, factor = factor),
       shared_signals_by_sign = get_pairwise_sharing(model_strong, factor = 0)
  )
}


#' Visualization
plot_mashr <- function(mm, share_factor = 3:8 / 10, override = FALSE, width = 5, height = 5, colrange = c(0, 1), save_to = "./", token = "none") {
  # Propotion of (significant) signals shared by each pair of conditions (?)
  col_dict <- c("B" = "#E64B35", "CD4T" = "#4DBBD5", "CD8T" = "#00A087", "NK" = "#3C5488", "Monocytes" = "#F39B7F")
  col_fun_prop <- colorRamp2(colrange, c("gray99", "#2166AC"))
  lapply(share_factor, function(f) {
    plot_path <- file.path(save_to, paste0(token, ".share_factor_", f, ".pdf"))
    mat <- get_pairwise_sharing(mm, factor = f) %>%
      (function(mat) {
         rownames(mat) <- colnames(mat)
         mat
      })

    if (!file.exists(plot_path) || override) {
      annotations <- rownames(mat) %>% stringr::str_split(pattern = "\\.", n = 2, simplify = TRUE)
      celltypes <- annotations[, 1]
      celltypes_col <- col_dict[celltypes]
      row_labels <- annotations[, ifelse(mode == "normal", 1, 2)]

      if (mode =="normal") {
        row_ann <- HeatmapAnnotation(
          Celltype = celltypes, which = "row",
          show_annotation_name = c("Celltype" = FALSE), col = list(Celltype = celltypes_col)
        )
      } else {
        char_vec <- stringr::str_split(row_labels, "_", simplify = TRUE)
        .time_vec <- char_vec[, 1]
        .stim_vec <- char_vec[, 2]
        row_ann <- HeatmapAnnotation(
          Celltype = celltypes, Stimulation = .stim_vec, Time = .time_vec, which = "row",
          show_annotation_name = c("Celltype" = FALSE, Stimulation = FALSE, Time = FALSE),
          col = list(Celltype = celltypes_col, Stimulation = c("LPS" = "darkblue", "RPMI" = "darkred"),
                     Time = c("T0" = "gray", "T3m" = "black"))
        )
      }
      col_ann <- HeatmapAnnotation(
        Celltype = celltypes, which = "column", show_annotation_name = c("Celltype" = FALSE),
        col = list(Celltype = celltypes_col), show_legend = c("Celltype" = FALSE)
      )
      pdf(plot_path, width = width, height = height)
      hm <- Heatmap(mat, col = col_fun_prop, name = "Shared signals",
                    left_annotation = row_ann,
                    bottom_annotation = col_ann,
                    row_labels = row_labels,
                    row_names_gp = gpar(fontsize = 9),
                    row_names_max_width = max_text_width(rownames(mat), gp = gpar(fontsize = 9)),
                    column_labels = row_labels,
                    column_names_rot = 45, column_names_side = "top", column_dend_side = "bottom",
                    column_names_gp = gpar(fontsize = 9),
                    column_names_max_height = max_text_width(colnames(mat), gp = gpar(fontsize = 9)))
      draw(hm)
      dev.off()
    }

    mat
  })
}


#
## Main steps.
#
proj_dir <- "/home/zzhang/Documents/projects/wp_bcg_eqtl"
in_dir <- file.path(proj_dir, "outputs/pseudo_bulk/summary_statistic")

# Two model used in the analysis.
mode_vec <- c("normal", "interaction")#, "per_condition")

# Main cell types.
celltype_vec <- c("Monocytes", "CD4T", "CD8T", "NK", "B") # , "pDC", "mDC")
ctmap <- c("mono" = "Monocytes", "cd8" = "CD8T", "cd4" = "CD4T", "nk" = "NK", "b" = "B")
ctmap_r <- c("Monocytes" = "mono", "CD8T" = "cd8", "CD4T" = "cd4", "NK" = "nk", "B" = "b")

# All comparisons used in the analysis.
condition_vec <- c("T0_LPS.vs.T0_RPMI", "T3m_RPMI.vs.T0_RPMI", "T3m_LPS.vs.T3m_RPMI") #, "T3m_LPS.vs.T0_RPMI")
condmap <- c("rpmi" = "T0_RPMI", "lps" = "T0_LPS", "bcg" = "T3m_RPMI", "bcg.lps" = "T3m_LPS")
condmap_r <- c("T0_RPMI" = "rpmi", "T0_LPS" = "lps", "T3m_RPMI" = "bcg", "T3m_LPS" = "bcg.lps")

# Header of summary statistics per condition by QTLtools.
eqtl_tab_header <- c(
  "phe_id", "phe_chr", "phe_from", "phe_to", "phe_strd", "n_var_in_cis", "dist_phe_var", "var_id",
  "var_chr", "var_from", "var_to", "nom_pval", "r_squared", "slope", "slope_se", "best_hit"
)
names(eqtl_tab_header) <- paste0("V", 1:16)

override <- TRUE
fdr <- 0.05
token <- paste0("fdr", fdr)
max_npv <- 0.1

plot_size <- data.frame(normal = c(width = 4.5, height = 3.25), per_condition = c(width = 7, height = 5.5))
for (mode in mode_vec) {
  save_to <- file.path(proj_dir, "outputs/pseudo_bulk/mashr", mode)
  if (!dir.exists(save_to)) dir.create(save_to, recursive = TRUE)

  mashr_res_path <- file.path(save_to, paste("mashr_res", token, "Rdata", sep = "."))
  if (!file.exists(mashr_res_path) || override) {
    qtltab_list <- list()
    for (per_celltype in celltype_vec) {
      if (mode == "normal") {
        root_dir <- file.path(in_dir, mode, per_celltype)
        .run_id <- per_celltype

        cat("[I]: Loading summary statistic from", root_dir, "\n")
        if (file.exists(file.path(root_dir, "top_qtl_results_all_FDR.txt")))
          qtltab_list[[.run_id]] <- load_eqtl_tab(root_dir, max_p_val = 5e-1)
        else
          cat("[W]: No top QTL results available for", .run_id, "Skipping ...\n")
      } else if (mode == "interaction") {
        for (per_cmp in condition_vec) {
          root_dir <- file.path(in_dir, mode, per_celltype, per_condition)
        }
        next
      } else if (mode == "per_condition") {
        new_cols <- c("p_value" = "nom_pval", "beta" = "slope", "beta_se" = "slope_se")
        root_dir <- file.path(in_dir, mode)
        for (per_condition in condmap) {
          .run_id <- paste(per_celltype, per_condition, sep = ".")
          .per_cond <- condmap_r[per_condition]
          .per_ct <- ctmap_r[per_celltype]

          per_inf <- file.path(root_dir, "nominal", paste0(.per_cond, "_", .per_ct, "_nominal.txt"))
          cat("[I]: Loading", per_inf, "\n")
          if (file.exists(per_inf))
            qtltab_list[[.run_id]] <- data.table::fread(per_inf, nThread = 4) %>%
              dplyr::rename_with(.fn = ~ eqtl_tab_header[.x]) %>%
              dplyr::filter(nom_pval < max_npv) %>%
              dplyr::mutate(QTL = paste0(phe_id, "-", var_id)) %>%
              dplyr::select(QTL, dplyr::all_of(new_cols))
          else
            cat("[W]: No summary statistics available! ", per_inf)
        }
      }
    }

    tar_qtl <- qtltab_list %>% lapply(function(tab) dplyr::filter(tab, !is.na(q_value))$QTL) %>% unlist() %>% unique()
    mashr_data <- qtltab_list %>% lapply(function(tab) dplyr::filter(tab, QTL %in% tar_qtl)) %>% prep_mashr_data()
    mashr_res <- exec_mashr(mashr_data)
    save(mashr_res, file = mashr_res_path)
  } else {
    load(mashr_res_path)
  }

  pw <- plot_size["width", mode]
  ph <- plot_size["height", mode]
  shared_signals <- plot_mashr(
    mashr_res$alter_model, seq(0, 9) / 10, override = TRUE, width = pw, height = ph,
    colrange = c(0, 1.1), save_to = save_to, token = mode
  )
}
