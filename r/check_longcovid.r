#!/usr/bin/env Rscript
library(tidyverse)
library(lme4)
library(lmerTest)

slfn5 <- read.delim("../../outputs/pseudo_bulk/replication/longcovid_from_Qiuyao/slfn5.txt")
slfn5$Genotype_binary <- ifelse(slfn5$Genotype %in% c("GG", "GA"), "Gx", "AA")
head(slfn5)

m0 <- lmer(AveExp_cd4 ~ Genotype * Timepoint + (1 | donor), slfn5, REML = FALSE)
summary(m0)

m1 <- lmer(AveExp_cd4 ~ Genotype + Timepoint + (1 | donor), slfn5, REML = FALSE)
summary(m1)

m2 <- lmer(AveExp_B ~ Genotype + (Genotype | donor), slfn5)
summary(m2)
