#!/usr/bin/env Rscript
# Author : Zhenhua Zhang
# E-mail : zhenhua.zhang217@gmail.com
# Created: Feb 15, 2022
# Updated: Feb 15, 2022

# TODO:
#   1. Add documents for each functions
#   2. A more robust implementation for the `load_data` function


# Strings should not be treated as factors; trace back when encounters errors.
options(
  error = traceback, stringsAsFactors = F, datatable.verbose = F,
  datatable.showProgress = F
)


# This is the only package we need to load at the beginning
library(magrittr)


# Parse CLI options
parse_cli_opt <- function(args = NULL) {
  parser <- optparse::OptionParser(prog = "scrnaseq")
  parser <- optparse::add_option(parser,
    opt_str = c("-o", "--out-dir"), metavar = "STR", default = "scRNAseq_out",
    help = "Output folder. Default: %default"
  )
  parser <- optparse::add_option(parser,
    opt_str = c("-t", "--data-type"), metavar = "STR", default = "10X",
    help = "The input data type (10X or RC). Default: %default"
  )
  parser <- optparse::add_option(parser,
    opt_str = c("-d", "--design-matrix"), metavar = "FILE",
    help = "A file include experimental design of case/control, batches etc."
  )
  parser <- optparse::add_option(parser,
    opt_str = c("-T", "--tmp-dir"), metavar = "DIR", default = NULL,
    help = "A temporary directory. Default: %default"
  )

  if (is.null(args)) {
    args <- commandArgs(trailingOnly = T)
  }

  opts <- optparse::parse_args(parser,
    args = args, positional_arguments = T, convert_hyphens_to_underscores = T
  )

  if (length(opts$args) < 1 && !opts$options$help) {
    stop("At least one positional argument is required.")
  }
  return(opts)
}


#' @title
#' Load multiple data according to a master file including design informations.
#'
#' @param design_mtrx The design matrix.
#' @param ext_id Column names used to add further information to cell barcodes.
#' @param dtype Data type. 10X or read counts.
#' @param proj_name The project name.
load_data <- function(design_mtrx, ext_id = c("groups", "batches", "pools"),
                      dtype = "10X", proj_name = "scRNAseq", min_cells = 5,
                      min_features = 200, calc_mt_rate = T, mt_pat = "^MT-",
                      merge = TRUE, merge_by = "groups", tmpdir = NULL) {
  base::stopifnot(dtype %in% c("10X", "RC"))

  obj_list <- base::apply(design_mtrx, 1, FUN = function(per_design) {
    rc_path <- per_design["rcpath"] # Read counts file / dir path
    names(rc_path) <- NULL

    if (dtype == "10X") {
      count_mtrx <- Seurat::Read10X(data.dir = rc_path)
    } else if (dtype == "RC") {
      rc_file <- base::file.path(rc_path, "read_counts.csv")
      # NOTE: Columns are cell barcode, rows are features (e.g. gene symbol)
      count_mtrx <- data.table::fread(rc_file, tmpdir = tmpdir)
    }

    # Meta information by design matrix
    design_items <- names(per_design)
    meta_mtrx <- per_design[!names(per_design) %in% c("rcpath")] %>%
      base::t() %>%
      base::as.data.frame() %>%
      dplyr::slice(rep(1, ncol(count_mtrx))) %>%
      purrr::set_names(design_items[!design_items %in% c("rcpath", "mtpath")])

    # Extra meta information
    mt_path <- ifelse("mtpath" %in% design_items, per_design["mtpath"], "")
    if (file.exists(mt_path)) {
      meta_mtrx <- data.table::fread(mt_path, tmpdir = tmpdir) %>%
        cbind(meta_mtrx)
    }

    # The rownames of meta-info matrix should be the col names of count matrix.
    base::rownames(meta_mtrx) <- base::colnames(count_mtrx)

    # Creat a Seurat Object from read counts
    obj <- Seurat::CreateSeuratObject(
      counts = count_mtrx, proj_name = proj_name, meta.data = meta_mtrx,
      min.cells = min_cells, min.features = min_features
    )

    if (calc_mt_rate) {
      obj[["pMT_RNA"]] <- Seurat::PercentageFeatureSet(obj, pattern = mt_pat)
    }

    return(obj)
  })

  # NOTE: The EXT_ID should be unique
  design_mtrx["EXT_ID"] <- apply(design_mtrx,
    MARGIN = 1, FUN = function(e) paste(e[ext_id], collapse = "_")
  )
  names(obj_list) <- design_mtrx[, "EXT_ID"]

  # Groups in the design matrix
  groups <- unique(design_mtrx[, merge_by])
  do_merge <- merge &&
    length(merge_by) > 0 &&
    length(groups) != nrow(design_mtrx) &&
    length(obj_list) > 1

  if (do_merge) {
    obj_names <- lapply(groups, function(g) {
      subset_names <- design_mtrx[design_mtrx[, merge_by] == g, "EXT_ID"]

      paste(subset_names, collapse = ".")
    })

    obj_list <- lapply(setNames(groups, obj_names), function(g, .obj_list) {
      subset_names <- design_mtrx[design_mtrx[, merge_by] == g, "EXT_ID"]
      print(subset_names)

      if (length(subset_names) > 1) {
        first_idx <- subset_names[1]
        rest_idx <- subset_names[2:length(subset_names)]

        subset_obj <- merge(.obj_list[[first_idx]], .obj_list[rest_idx],
          add.cell.ids = subset_names, project = proj_name
        )

        return(subset_obj)
      }

      return(.obj_list[[subset_names]])
    }, .obj_list = obj_list)
  }

  return(obj_list)
}


#
draw_qc_plots <- function(obj, split_by = c("batches", "pools"),
                          save_to = "./") {
  # Features
  features <- c("nCount_RNA", "nFeature_RNA", "pMT_RNA")

  # Violin plots show
  g_vln <- Seurat::VlnPlot(obj, features = features) & ggplot2::labs(x = NULL)

  # Scatter plots show
  g_sct_1 <- Seurat::FeatureScatter(obj, "nCount_RNA", "nFeature_RNA")
  g_sct_2 <- Seurat::FeatureScatter(obj, "nCount_RNA", "pMT_RNA")

  g_comb <- g_vln / (g_sct_1 + g_sct_2) &
    ggplot2::theme_classic() &
    ggplot2::theme(legend.position = "none")

  ggplot2::ggsave(save_to, plot = g_comb, width = 7, height = 7)
}


# TODO: add an option for TMP dir
eval_de_genes <- function(obj_list, save_to, norm_method = "LogNormalize",
                          n_count_lim = c(200, NA), n_feature_lim = c(20, NA),
                          p_mt_lim = c(0, 0.5)) {
  obj_names <- names(obj_list)

  lapply(setNames(obj_names, obj_names), function(per_rec, .obj_list) {
    obj_name <- per_rec
    obj_perse <- .obj_list[[obj_name]]

    .save_to <- paste0(file.path(save_to, obj_name), ".rds")
    obj <- obj_perse %>%
      (function(.obj) {
        kept_cells <- .obj@meta.data %>%
          (function(md) {
            lims <- quantile(seq.int(1, nrow(md), 1), probs = c(0.05, 0.95)) %>%
              as.integer() %>%
              (function(e) seq(e[1], e[2], by = 1))

            keep_idx <- md$nCount_RNA %>% order() %in% lims
            rownames(md)[keep_idx]
          })

        .obj[, kept_cells]
      }) %>%
      Seurat::NormalizeData(
        normalization.method = norm_method, verbose = FALSE
      ) %>%
      Seurat::FindVariableFeatures(
        selection.method = "vst", nfeatures = 2000, verbose = FALSE
      ) %>%
      (function(.obj, label_top_n = 10) { # Show highly variable features
        # p1 <- Seurat::VariableFeaturePlot(.obj)
        # top_n_lables <- head(Seurat::VariableFeatures(.obj), label_top_n)
        # p2 <- Seurat::LabelPoints(.obj, points = top_n_lables, repel = TRUE)
        # p1 + p2
        return(.obj)
      }) %>%
      Seurat::ScaleData(features = rownames(.), verbose = FALSE) %>%
      Seurat::RunPCA(
        features = Seurat::VariableFeatures(.), verbose = FALSE
      ) %>%
      (function(.obj) { # Show PCA, dot plots and heatmap
        # Seurat::PCAPlot(.obj)
        return(.obj)
      }) %>%
      Seurat::FindNeighbors(dims = 1:10, verbose = FALSE) %>%
      Seurat::FindClusters(resolution = 0.5, verbose = FALSE) %>%
      Seurat::RunUMAP(dims = 1:10, verbose = FALSE) %>%
      (function(.obj) { # Show UMAP
        # Seurat::DimPlot(.obj)
        return(.obj)
      }) %>%
      (function(.obj) {
        base::saveRDS(.obj, .save_to)
        return(.obj)
      })

    obj
  }, .obj_list = obj_list)
}


assign_cell_type <- function(obj_list, which_label = "label.fine",
                             cache_dir = ".") {
  # Cache Annotation hub
  ah_cdir <- file.path(cache_dir, "annotation_hub")
  dir.create(ah_cdir, showWarnings = FALSE, recursive = TRUE)

  # Cache for experiment hub
  eh_cdir <- file.path(cache_dir, "experiment_hub")
  dir.create(eh_cdir, showWarnings = FALSE, recursive = TRUE)

  # Setup environmental variables, then celldex will use them for cache.
  Sys.setenv(ANNOTATION_HUB_CACHE = ah_cdir, EXPERIMENT_HUB_CACHE = eh_cdir)

  lapply(obj_list, function(per_obj) {
    # Normalized data by Seurat, the slot should be data or scale.data
    test <- Seurat::GetAssayData(per_obj, slot = "data")
    for (db_name in c("HPCA", "BE", "DICE")) {
      if (db_name == "HPCA") {
        db <- celldex::HumanPrimaryCellAtlasData()
      } else if (db_name == "BE") {
        db <- celldex::BlueprintEncodeData()
      } else {
        db <- celldex::DatabaseImmuneCellExpressionData()
      }

      # NovershternHematopoieticData
      # NOTE: It is likely the best option for bone marrow samples.

      singler_labels <- SingleR::SingleR(
        test = test, ref = db, labels = db[[which_label]],
        clusters = per_obj@meta.data$seurat_clusters,
        assay.type.test = 1
      )

      ref_idents_name <- paste("SingleR", db_name, sep = "_")
      matched_tag <- match(per_obj$seurat_clusters, rownames(singler_labels))
      per_obj[[ref_idents_name]] <- singler_labels$labels[matched_tag]
    }

    per_obj
  })
}


calc_mean_exp <- function(obj_list, method = "mean") {
  lapply(obj_list, function(per_obj) {
    Seurat::AverageExpression(per_obj, verbose = FALSE)
  })
}


base_path <- "/vol/projects/CIIM/300BCG/300BCG_scRNA/scRNAdata"
design_mtrx <- data.frame(
  batches = c("B1", "B1", "B1", "B1"), # required
  pools = c("P1", "P2", "P3", "P4"), # required
  groups = c("G1", "G1", "G2", "G2"),
  rcpath = file.path(base_path, c(
    "Batch_1_Pool_1count", "Batch_1_Pool_2count",
    "Batch_1_Pool_3count", "Batch_1_Pool_2count"
  ), "outs/filtered_feature_bc_matrix") # required
)


raw_obj <- load_data(design_mtrx)
eval_obj <- eval_de_genes(raw_obj, ".")
ann_obj <- assign_cell_type(eval_obj, which_label = "label.main")
avg_exp <- calc_mean_exp(ann_obj)
avg_exp[[1]]$RNA %>% head()


raw_qc_plots <- draw_qc_plots(obj)
eval_qc_plots <- draw_qc_plots(eval_obj)

write_to_disk <- function(mean_exp_matrix, fpath = "mean_exp_matrix.csv") {
}


main <- function() {
  cliargs <- parse_cli_opt()

  args <- cliargs$args
  opts <- cliargs$options

  out_dir <- opts$out_dir
  design_mtrx <- opts$design_matrix

  mean_exp_matrix <- load_data(design_mtrx) %>%
    eval_de_genes() %>%
    calc_mean_exp() %>%
    write_to_disk(out_dir)
}

main()
