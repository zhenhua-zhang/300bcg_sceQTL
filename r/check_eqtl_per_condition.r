# Code to collect conditional eQTL.
options(stringsAsFactors = FALSE, datatable.verbose = FALSE, datatable.showProgress = FALSE)
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

proj_dir <- "~/Documents/projects/wp_bcg_eqtl"
work_dir <- file.path(proj_dir, "outputs/pseudo_bulk/summary_statistic/per_condition")

data_dir <- "/vol/projects/wli/projects/genetics/results/qtltools_dos_1mb/nominal"
celltype_map <- c("mono" = "Monocytes", "cd8" = "CD8T", "cd4" = "CD4T", "nk" = "NK", "b" = "B")
condition_map <- c("rpmi" = "T0_RPMI", "lps" = "T0_LPS", "bcg" = "T3m_RPMI", "bcg.lps" = "T3m_LPS")
eqtl_tab_header <- c(
  "phe_id", "phe_chr", "phe_from", "phe_to", "phe_strd", "n_var_in_cis", "dist_phe_var", "var_id",
  "var_chr", "var_from", "var_to", "nom_pval", "r_squared", "slope", "slope_se", "best_hit"
)
names(eqtl_tab_header) <- paste0("V", 1:16)

overwrite <- FALSE
max_nom_pval <- 5e-6
file_pool <- file.path(data_dir, list.files(data_dir))

save_to <- file.path(work_dir, "eqtl.per_condition.csv")
if (!file.exists(save_to) || overwrite) {
  eqtl_nom_tab <- lapply(file_pool, function(p_if) {
      p_if_nm_vec <- stringr::str_split(basename(p_if), pattern = "_", simplify = TRUE)
      p_condition <- condition_map[p_if_nm_vec[1]]
      p_celltype <- celltype_map[p_if_nm_vec[2]]

      if (is.na(p_celltype) || is.na(p_condition)) return(NULL)

      data.table::fread(p_if, nThread = 4) %>%
        dplyr::rename_with(.fn = ~ eqtl_tab_header[.x]) %>%
        dplyr::filter(nom_pval < max_nom_pval) %>%
        dplyr::mutate(cell_type = p_celltype, condition = p_condition)
    }) %>%
    Reduce(rbind, .)
  data.table::fwrite(eqtl_nom_tab, save_to)
} else {
  eqtl_nom_tab <- data.table::fread(save_to)
}

candiate_eqtl_tab <- eqtl_nom_tab %>%
  dplyr::select(phe_id, var_id, nom_pval, slope, slope_se, best_hit, cell_type, condition) %>%
  tidyr::pivot_wider(names_from = c(cell_type, condition), values_from = c(nom_pval, slope, slope_se), names_sep = ".") %>%
  dplyr::filter(dplyr::if_any(dplyr::contains("nom_pval"), ~ !is.na(.x) & .x < 5e-8)) %>%
  dplyr::select(phe_id, var_id, sort(dplyr::starts_with("nom_pval")))
candiate_eqtl_tab %>% data.table::fwrite(file.path(work_dir, "combined_eqtl.all.csv"))

candiate_eqtl_tab %>%
  dplyr::filter(
    dplyr::if_all(dplyr::ends_with("T0_RPMI"), ~ is.na(.x)),
    dplyr::if_any(dplyr::ends_with("T0_LPS"), ~ !is.na(.x))
  ) %>%
  dplyr::group_by(phe_id) %>%
  dplyr::slice_head(n = 1) %>%
  as.data.frame()
