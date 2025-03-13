#!/usr/bin/env Rscript
# File: determine_sample_size.r
# Author: Zhenhua Zhang
# E-mail: zhenha.zhang217@gmail.com
# Created: Mar 12, 2025
# Updated:

suppressPackageStartupMessages({
  library(ggplot2)
  library(designr) # designr, simLMM
  library(lme4) # lmer
  library(afex) # aov_car
  library(magrittr)
  library(data.table) # fread
})

set.seed(4096)

empirical_genotype_beta <- fread("/home/zzhang/Documents/projects/resources/GTEx/V8/Whole_Blood.v10.eGenes.txt") %>%
  dplyr::pull(slope) %>% abs() %>% quantile(probs = c(0.25))

empirical_intercept <- fread("/home/zzhang/Documents/projects/resources/GTEx/V8/Whole_Blood.v10.normalized_expression.bed") %>%
  dplyr::select(-c("#chr", "start", "end", "gene_id")) %>%
  as.matrix() %>%
  max()

empirical_tab <- fread("/home/zzhang/Documents/projects/resources/GTEx/V8/Whole_Blood.v10.normalized_expression.bed") %>%
  dplyr::select(-c("#chr", "start", "end")) %>%
  dplyr::slice_sample(n = 1000) %>%
  tidyr::pivot_longer(-gene_id, names_to = "Donor", values_to = "expr") %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarize(
    intercept_sd = dplyr::pick(expr, Donor) %>% (function(x) {
      m <- lm(expr ~ 1 + Donor)
      coeffcients <- coef(m)
      data.frame(sd = sd(coeffcients[2:length(coeffcients)]), mean = mean(coeffcients[2:length(coeffcients)]))
    })) %>%
  tidyr::unnest(cols = c(intercept_sd))

empirical_donor_residual_sd <- empirical_tab %>% dplyr::pull(sd) %>% mean()
empirical_donor_beta <- empirical_tab %>% dplyr::pull(mean) %>% mean()

n_sample_vec <- seq(11, 150, 1) # minimal sample size is 10, maximum is 100, step is 1
n_sims <- length(n_sample_vec) # number of simulations.
n_repeat <- 1000 # number of repeats per simulation

fix      <- c(empirical_intercept, empirical_genotype_beta, empirical_donor_beta) # Empirical statistics, 
sd_Subj  <- 0
sd_Res   <- empirical_donor_residual_sd
vc_sd <- list(sd_Subj, sd_Res)

per_sim_seeds <- rexp(n_sims * n_repeat) %>% `*`(100000000) %>% as.integer() %>% matrix(nrow = n_sims, ncol = n_repeat)

cof_tab <- NULL
for (i in 1:n_sims) {
  design <- fixed.factor("Genotype", levels = c("AA", "AT", "TT")) +
    random.factor("Conditoin", instances = 4) +
    random.factor("Donor", instances = n_sample_vec[i])

  dat <- design.codes(design) %>% dplyr::group_by(Conditoin, Donor) %>% dplyr::slice_sample(n = 1) %>% dplyr::ungroup()
  contrasts(dat$Genotype) <- c(0, 1, 2)
  dat$Genotype_n <- model.matrix(~ Genotype, dat)[, 3]

  n_donor <- length(unique(dat$Donor))
  for (j in 1:n_repeat) {
    cat("[I]: Simulation", i, "Repeat", j, "...\r\b")
    # simulate data
    fml <- formula("~ Genotype + (1 | Donor)")
    set.seed(per_sim_seeds[i, j])
    dat$expr <- simLMM(fml, data = dat, Fixef = fix, VC_sd = vc_sd, CP = 0.3, empirical = FALSE, verbose = FALSE)

    # LMM
    ww <- ""
    suppressMessages(suppressWarnings(
      LMM <- withCallingHandlers({
        lmer(expr ~ Genotype + (1 | Donor), data = dat, REML=FALSE, control = lmerControl(calc.derivs=FALSE))
      }, warning = function(w) { ww <<- w$message })
    ))

    cof_tab <- summary(LMM) %>% coef() %>% as.data.frame() %>%
      tibble::rownames_to_column("Coef") %>%
      dplyr::filter(Coef == "Genotype2") %>%
      dplyr::mutate(nsj = n_donor, perm_idx = j, singular = isSingular(LMM), warning = ww, noWarning = ww == "") %>%
      dplyr::rename(c("SE" = "Std. Error", "t" = "t value", "p" = "Pr(>|t|)")) %>%
      dplyr::mutate(sign = dplyr::if_else(p < 0.05, 1.0, 0.0)) %>%
      dplyr::bind_rows(cof_tab, .)
  }

  m0 <- loess(sign ~ nsj, data = cof_tab)
  max_power <- cof_tab %>% dplyr::mutate(pred = predict(m0)) %>% dplyr::pull(pred) %>% max()
  cat("[I]: maximum power is", max_power, "at sample size", n_sample_vec[i], "\n")
}
