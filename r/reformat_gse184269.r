#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(magrittr)
  library(tidyverse)
  library(data.table)
  library(rtracklayer)

  library(GenomicRanges)
})

projdir <- "~/Documents/projects/wp_bcg_eqtl"

new_header <- c(
  "title", "geo_accession", "status", "submission_date", "last_update_date", "type", "channel_count",
  "tissue", "organism_ch1", "donor_id", "age", "sex", "treatment_protocol_ch1",
  "growth_protocol_ch1", "molecule_ch1", "extract_protocol_ch1", "label_ch1", "label_protocol_ch1",
  "taxid_ch1", "hyb_protocol", "scan_protocol", "description", "data_processing_kit", "data_processing_qc",
  "platform_id", "contact_name", "contact_laboratory", "contact_department", "contact_institute",
  "contact_address", "contact_city", "contact_state", "contact_zip", "contact_country",
  "supplementary_file_Grn_idat", "supplementary_file_Red_idat", "data_row_count"
)
tar_col <- paste0("V", seq(1, length(new_header)))
names(tar_col) <- new_header


read_cmd <- "zgrep -E '^!Sample_' ./GSE184269/GSE184269_series_matrix.txt.gz | tr -d '!'"
meta_info <- data.table::fread(cmd = read_cmd, header = FALSE) %>%
  data.table::transpose() %>%
  as.data.frame() %>%
  dplyr::filter(!V1 %in% "Sample_title") %>%
  dplyr::select(dplyr::all_of(tar_col[c(8, 10, 11, 12, 35)])) %>%
  dplyr::mutate(
    sample_id = stringr::str_extract(supplementary_file_Grn_idat, pattern = "[0-9]{12,12}_R[0-9]{2,2}C[0-9]{2,2}"),
    donor_id = stringr::str_remove(donor_id, "donor id: "),
    age = as.integer(stringr::str_remove(age, "age: ")),
    sex = stringr::str_remove(sex, "Sex: "),
    tissue = stringr::str_replace(tissue, "\303\257", "i"),
  ) %>%
  dplyr::select(donor_id, age, sex, tissue, sample_id) %>%
  dplyr::arrange(donor_id, tissue)

meta_info %>% data.table::fwrite("./GSE184269/meta_infomation.csv")

#
chain <- file.path(projdir, "inputs/reference/UCSC_GB/hg19ToHg38.over.chain") %>% import.chain()
meth_tab <- file.path(projdir, "inputs/methylation/GSE184269/GSE184269_Matrix_processed_FINAL.annotated.grch37.csv.gz") %>% 
  data.table::fread()

cpg_grch38_meth_tab <- meth_tab %>%
  dplyr::select(-dplyr::starts_with("annot."), -dplyr::ends_with("Detection_Pval")) %>%
  makeGRangesFromDataFrame(TRUE) %>%
  (function(tab) {seqlevelsStyle(tab) <- "UCSC"; tab}) %>%
  liftOver(chain) %>%
  as.data.frame() %>%
  dplyr::select(-group, -group_name)

cpg_grch38_meth_tab %>% data.table::fwrite("GSE184269_Matrix_processed_FINAL.grch38.csv")
