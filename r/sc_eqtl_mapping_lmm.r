#!/usr/bin/env Rscript
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: Apr 12, 2022
# Updated: Apr 12, 2022

# Some simulation code to test the LMM for sceQTL mapping.

library(data.table)
library(multidplyr)
library(tidyverse)
library(lmerTest)
library(magrittr)
library(Seurat)
library(lme4)

set.seed(31415)

tmpdir <- "~/Documents/projects/wp_bcg_eqtl/temps"

n_patients <- 40
n_genes <- 10
n_snps <- 10 * n_genes

gene_ids <- paste0("G", 1:n_genes)
patient_ids <- paste0("P", 1:n_patients)

bound <- sample(800:1200, 2) %>% sort()
n_cells <- sample(bound[1]:bound[2], n_patients)

cvrt_tab <- data.frame(
  PatientID = patient_ids,
  NCells = n_cells,
  Age = sample(18:100, n_patients, replace = TRUE),
  Sex = sample(c("M", "F"), n_patients, replace = TRUE)
)

# NCells to PatientID map
ctp_map <- cvrt_tab %>%
  select(NCells, PatientID) %>%
  deframe()

# Expression matrix per cell
expr_means <- rnorm(n_genes)
phtp_tab <- expand.grid(n_cells, expr_means) %>%
  rename("n_cells" = "Var1", "mean" = "Var2") %>%
  group_by(n_cells) %>%
  summarise(expmat = {
    cdata <- cur_data()
    mn <- cdata %$% mean
    nc <- cdata %$% n_cells %>% unique()

    lapply(mn, function(m) {
      rnorm(nc, m, 10)
    }) %>%
      Reduce(cbind, .) %>%
      as.data.frame()
  }) %>%
  as.data.table() %>%
  mutate(PatientID = ctp_map[as.character(n_cells)]) %>%
  select(-n_cells) %>%
  set_names(c(gene_ids, "PatientID"))

# The number of SNPs should be aroud n_snps * n_genes
gntp_tab <- expand.grid(
  patient_ids,
  paste0("rs", str_pad(1:n_snps, as.integer(log10(n_snps)) + 1, pad = "0"))
) %>%
  rename("PatientID" = "Var1", "SNPs" = "Var2") %>%
  mutate(
    dosage = 1 + rnorm(nrow(.), 0, 0.5),
    dosage = case_when(dosage < 0 ~ 0.0, dosage > 2 ~ 2, T ~ dosage),
  ) %>%
  pivot_wider(id_cols = SNPs, names_from = PatientID, values_from = dosage) %>%
  mutate(GeneID = sample(gene_ids, n_snps, TRUE)) %>%
  relocate(GeneID, .before = SNPs) %>%
  as.data.frame()


# Write
gntp_tab %>%
  group_by(GeneID) %>%
  group_walk(
    ~ fwrite(.x, sprintf("%s/%s.csv", tmpdir, .y$GeneID)),
    .keep = T
  )

files <- list.files(tmpdir, pattern = "*.csv", full.names = TRUE)

cluster <- new_cluster(10)
cluster_library(cluster, c("tidyverse", "data.table", "lmerTest", "lme4"))
cluster_assign_partition(cluster, filename = files)
cluster_send(cluster, my_data <- vroom::vroom(filename))
cluster_copy(cluster, c("phtp_tab", "cvrt_tab", "tmpdir"))

reg_tab <- party_df(cluster, "my_data") %>%
  group_by(GeneID, SNPs) %>%
  summarise(est = {
    tmp <- cur_group()
    gene_id <- tmp$GeneID %>% as.character()
    snp_id <- tmp$SNPs %>% as.character()

    sub_gntp <- cur_data() %>%
      pivot_longer(everything(), names_to = "PatientID", values_to = snp_id)

    sub_phtp <- phtp_tab %>% select(one_of(gene_id, "PatientID"))

    fml <- formula(paste(gene_id, "~", snp_id, "+ (1 | PatientID) + Age + Sex"))
    suppressMessages(
      sub_gntp %>%
        inner_join(sub_phtp, by = "PatientID") %>%
        inner_join(cvrt_tab, by = "PatientID") %>%
        lmer(fml, data = .) %>%
        summary() %>%
        coef() %>%
        as.data.frame() %>%
        mutate(Parameter = rownames(.)) %>%
        set_names(c("Est", "StdErr", "df", "tvalue", "Pvalue", "Param")) %>%
        filter(str_detect(Param, "rs"))
    )
  }) %>%
  collect() %>%
  as.data.table() %>%
  rename_with(.fn = function(e) str_remove_all(e, "est.")) %>%
  fwrite(file.path(tmpdir, "sc_eqtl.csv"))

# vim: set tw=200:
