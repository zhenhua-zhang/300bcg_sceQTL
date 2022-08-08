#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, error = traceback)

library(magrittr)

#' Load genotypes. Fetch postion and GT field
load_gntp <- function(fpath, tmp_dir = ".") {
  dis_cols <- c("ID", "QUAL", "FILTER", "INFO", "FORMAT")
  pos_cols <- c("CHROM", "POS", "REF", "ALT")

  gntp_dtfm <- data.table::fread(fpath,
    skip = "CHROM", showProgress = FALSE, verbose = FALSE, tmpdir = tmp_dir
  ) %>%
    dplyr::rename("CHROM" = "#CHROM") %>%
    dplyr::select(-dplyr::one_of(dis_cols)) %>%
    dplyr::mutate(dplyr::across(
      -dplyr::one_of(pos_cols),
      function(.x) stringr::str_extract(.x, "[01][/|][01]")
    ))

  return(gntp_dtfm)
}


#' Load SouporCell results, including called genotypes and barcodes to clusters
load_soupercell_res <- function(path, keep_singlet = TRUE) {
  dis_cols <- c("ID", "QUAL", "FILTER", "INFO", "FORMAT")
  pos_cols <- c("CHROM", "POS", "REF", "ALT")

  sc_list <- lapply(path, function(per_path) {
    gntp_dtfm <- file.path(per_path, "cluster_genotypes.vcf") %>%
      load_gntp() %>%
      dplyr::rename_with(-dplyr::one_of(pos_cols),
        .fn = function(e) paste0("cluster_", e)
      )

    clst_dtfm <- file.path(per_path, "clusters.tsv") %>%
      data.table::fread() %>%
      dplyr::select(barcode, status, assignment) %>%
      dplyr::mutate(assignment = paste0("cluster_", assignment))

    if (keep_singlet) {
      clst_dtfm <- clst_dtfm %>% dplyr::filter(status == "singlet")
    }

    list(genotype = gntp_dtfm, assignments = clst_dtfm)
  })

  return(sc_list)
}



#' Assign barcode to donors
assign_bc <- function(dmres, gntp) {
  dis_cols <- c("CHROM", "POS", "REF_c", "ALT_c", "REF_r", "ALT_r")

  assign_list_per_donor <- lapply(dmres, function(per_dmres) {
    comb_gntp_dtfm <- per_dmres$genotype %>%
      dplyr::inner_join(gntp,
        by = c("CHROM", "POS"), suffix = c("_c", "_r")
      ) %>%
      dplyr::mutate(dplyr::across(
        -dplyr::one_of(dis_cols),
        ~ dplyr::case_when(
          .x %in% c("0/0", "0|0") ~ 0,
          .x %in% c("0/1", "1/0", "0|1", "1|0") ~ 1,
          .x %in% c("1/1", "1|1") ~ 2
        )
      )) %>%
      dplyr::mutate(dplyr::across(
        -dplyr::one_of(dis_cols),
        ~ dplyr::case_when(
          REF_c == REF_r & ALT_c == ALT_r ~ .x,
          REF_c == ALT_r & ALT_c == REF_r ~ (.x - 1) * (-1) + 1,
          T ~ -999
        )
      )) %>%
      dplyr::filter(dplyr::across(
        -dplyr::one_of(dis_cols),
        ~ .x >= 0 # No -999 is allowed
      ))

    n_total_snp <- nrow(comb_gntp_dtfm)
    tmp_dtfm <- comb_gntp_dtfm %>%
      colnames() %>%
      purrr::discard(function(e) e %in% dis_cols) %>%
      combn(2) %>%
      t() %>%
      as.data.frame() %>%
      purrr::set_names(c("left", "right")) %>%
      dplyr::filter(
        (stringr::str_detect(left, "cluster_") &
          stringr::str_detect(right, "300BCG")) |
          (stringr::str_detect(left, "300BCG") &
            stringr::str_detect(right, "cluster_"))
      ) %>%
      apply(1, function(per_row, .gntp_dtfm) {
        ctable <- .gntp_dtfm %>%
          dplyr::select(one_of(per_row)) %>%
          table()

        p_value_chi <- ctable %>% chisq.test() %$% p.value
        n_align_snp <- diag(ctable) %>% sum()
        match_ratio <- n_align_snp / n_total_snp

        names(per_row) <- NULL
        return(data.frame(
          cluster = per_row[1], donor = per_row[2], htest = "chisq",
          p_value = p_value_chi, match_ratio = match_ratio,
          n_matched_snp = n_align_snp, n_total_snp = n_total_snp
        ))
      }, .gntp_dtfm = comb_gntp_dtfm) %>%
      Reduce(rbind, .)

    return(tmp_dtfm)
  })

  for (nm in names(assign_list_per_donor)) {
    assign_list_per_donor[[nm]]["batch_pool"] <- nm
  }

  assign_dtfm_per_donor <- Reduce(rbind, assign_list_per_donor)

  assign_dtfm_per_cell <- lapply(names(dmres), function(per_dmres) {
    gntp_map <- assign_dtfm_per_donor %>%
      dplyr::filter(batch_pool == per_dmres) %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_max(match_ratio) %>%
      dplyr::select(cluster, donor) %>%
      dplyr::distinct() %>%
      tibble::deframe()

    print(gntp_map)
    dmres[[per_dmres]][["assignments"]] %>%
      dplyr::mutate(
        assigned_donor = gntp_map[assignment],
        batch_pool = per_dmres
      )
  }) %>%
    Reduce(rbind, .)

  return(
    list(donor_ass = assign_dtfm_per_donor, cell_ass = assign_dtfm_per_cell)
  )
}


proj_dir <- "~/Documents/projects/wp_bcg_eqtl"
in_dir <- file.path(proj_dir, "inputs")
out_dir <- file.path(proj_dir, "outputs")
temp_dir <- file.path(proj_dir, "temps")

gntp_dtfm <- file.path(in_dir, "genotypes/300BCG_sub40_imp_hg38_ids.vcf.gz") %>%
  load_gntp(tmp_dir = temp_dir)

dm_list <- "/vol/projects/wli/projects/300bcg/results/souporcell" %>%
  list.dirs(recursive = FALSE) %>%
  purrr::discard(function(e) grepl("_", e)) %>%
  purrr::set_names(basename(.)) %>%
  load_soupercell_res()

alist <- assign_bc(dm_list, gntp_dtfm)

donor_ass <- alist$donor_ass
save_to <- file.path(out_dir, "demultiplexing/Souporcell_dm_perdonor.csv")
donor_ass %>% data.table::fwrite(save_to, row.names = FALSE, quote = FALSE)

cell_ass <- alist$cell_ass
save_to <- file.path(out_dir, "demultiplexing/Souporcell_dm_percell.csv")
cell_ass %>% data.table::fwrite(save_to, row.names = FALSE, quote = FALSE)
