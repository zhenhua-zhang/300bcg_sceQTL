#!/usr/bin/env bash
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: 2023 Mar 05
# Updated: 2023 Apr 05


# Raw data of 300OB
# /vol/projects/CIIM/meta_cQTL/data/Raw/300OB

# Imputed genotypes by Javi
# /vol/projects/CIIM/300OB/QTL/out/imputation/local/chr2.dose.vcf.gz

projdir=~/Documents/projects/wp_bcg_eqtl
workdir=$projdir/outputs/bmi_qtl

# #
# ## Genotype QC, pre-imputation. Done by Javi
# #
# raw_bfile=$projdir/inputs/300OB/genotypes/pre_imputation/300OB
# 
# # Sample QC
# # 1. Sex check
# # No sex chromosome available.
# if [[ do_sex_check == true ]]; then
#   plink --bfile $raw_bfile --check-sex --out $workdir/genotypes/quality_control/300OB.check_sex
#   awk 'NR == 1 || $x == "" {print}' >| $workdir/genotypes/quality_control/failed.sex_check.txt
# fi
# 
# # Genotype missing rates, and heterozigosity
# # Missing rate between 3% and 7%
# plink --bfile $raw_bfile --missing --out $workdir/genotypes/quality_control/300OB.call_rates
# awk 'NR==1 || $6>=0.1 {print}' $workdir/genotypes/quality_control/300OB.call_rates.imiss \
#   >| $workdir/genotypes/quality_control/failed.individual_call_rates.txt 
# awk 'NR==1 || $x<=0.xx {print}' $workdir/genotypes/quality_control/300OB.call_rates.lmiss \
#   >| $workdir/genotypes/quality_control/failed.snp_call_rates.txt
# 
# plink --bfile $raw_bfile --het --out $workdir/genotypes/quality_control/300OB.heterozigosity
# awk 'BEGIN {total = 0} {total += $3 / $4} END {print total / NR}' $workdir/genotypes/quality_control/300OB.heterozigosity.het
# 
# 
# # Identifying relatives
# # high_ld=$projdir/inputs/reference/high-LD-regions.txt
# # plink --bfile $raw_bfile --exclude $high_ld --range --indep-pairwise 50 5 0.2 --out $workdir/genotypes/quality_control/300OB.identify_relatives
# # plink --bfile $raw_bfile --extract $workdir/genotypes/quality_control/300OB.identify_relatives.prune.in --genome --out $workdir/genotypes/quality_control/300OB.identify_relatives
# 
# #
# ## Imputation. Done by Javi
# #

#
## Genotype QC, post-imputation
#
chrommap=$projdir/inputs/reference/chrom_map.txt
refvars=/vol/projects/BIIM/resources/snpDB/GRCh38_b150/00-All.vcf.gz
imputed=/vol/projects/CIIM/300OB/QTL/out/imputation/local/chr2.dose.vcf.gz
imputed_chr=$projdir/outputs/300OB/genotypes/chr2.dose.rm_chr.vcf.gz
imputed_chr_filtered=$projdir/outputs/300OB/genotypes/chr2.dose.rm_chr.filtered.vcf.gz

# Rename chromosome name: chrN to N
cat <<EOF >| $chrommap
chr1 1
chr2 2
chr3 3
chr4 4
chr5 5
chr6 6
chr7 7
chr8 8
chr9 9
chr10 10
chr11 11
chr12 12
chr13 13
chr14 14
chr15 15
chr16 16
chr17 17
chr18 18
chr19 19
chr20 20
chr21 21
chr22 22
EOF

bcftools annotate --rename-chrs $chrommap --threads 10 -o $imputed_chr $imputed \
  && bcftools index -t $imputed_chr


# Annotate and filter the variants.
bcftools annotate -c ID -a $refvars --threads 5 $imputed_chr \
  | bcftools norm -d snps --threads 5  \
  | bcftools view -O z -v snps -i 'MAF>=0.05 && R2>=0.3 && %ID!="." && COUNT(GT=="het")>=3 && COUNT(GT=="RR")>=3 && COUNT(GT=="AA")>=3' --threads 4 -o $imputed_chr_filtered \
  && bcftools index -t $imputed_chr_filtered

#
## Phenotypes
#
idmap_file=$projdir/inputs/300OB/phenotypes/id.code.rank.tsv
phenotype_file=$projdir/inputs/300OB/phenotypes/Leukocytenumbers.csv
pheno_dir=$projdir/outputs/300OB/phenotypes
qtl_genotype_file=$projdir/outputs/300OB/genotypes/qtl_genotypes.chr2.vcf.gz
vcf_sample_name_map=$projdir/outputs/300OB/genotypes/new_sample_name_map.txt
phenotype=BMI
interaction_covar=monocytes_percent
interaction_covar=Monocytes_absolute_10e09_per_l

head -1 $phenotype_file | tr "," "\n" | cat -n

Rscript - <<EOF
library(magrittr)

idmap <- data.table::fread("$idmap_file") %>% dplyr::mutate(exact_VCF_header = paste0("0_", \`VCF-header\`))
phenotype_tab <- data.table::fread("$phenotype_file")
merged_tab <- dplyr::inner_join(idmap, phenotype_tab, by = c("lab-code" = "Cxxx")) %>%
  dplyr::rename(lab_code = \`lab-code\`) %>%
  dplyr::arrange(lab_code)

outdir <- "$pheno_dir"

# Phenotypes
## For QTLtools
out_file <- file.path(outdir, "qtl_phenotype.QTLtools.$phenotype.bed")
merged_tab %>%
  dplyr::select(lab_code, \`$phenotype\`) %>%
  dplyr::mutate(\`#Chr\` = "ChrN", start = 0, end = 1, pid = "0", gid = "0", strand = "-") %>%
  tidyr::pivot_wider(values_from = "$phenotype", names_from = "lab_code") %>%
  data.table::fwrite(out_file, sep = '\\t')

## For PLINK
out_file <- file.path(outdir, "qtl_phenotype.PLINK.$phenotype.txt")
merged_tab %>%
  dplyr::select(dplyr::all_of(c("PID" = "lab_code", "IID" = "lab_code", "$phenotype"))) %>%
  dplyr::mutate(dplyr::across(\`$phenotype\`, ~dplyr::if_else(is.na(.x), -9, .x))) %>%
  data.table::fwrite(out_file, sep = ' ')

# Covariates
# For QTLtools
out_file <- file.path(outdir, "qtl_covariate.QTLtools.txt")
merged_tab %>%
  dplyr::select(lab_code, sex, age) %>%
  dplyr::mutate(sex = dplyr::if_else(sex == "Male", 1, 2)) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::mutate(id = c("#id", "sex", "age")) %>%
  dplyr::relocate(id) %>%
  data.table::fwrite(out_file, col.names = FALSE, sep = '\\t')

# For PLINK
out_file <- file.path(outdir, "qtl_covariate.PLINK.txt") # Non-interaction
merged_tab %>%
  dplyr::select(dplyr::all_of(c("PID" = "lab_code", "IID" = "lab_code", "sex", "age"))) %>%
  dplyr::mutate(sex = dplyr::if_else(sex == "Male", 1, 2)) %>%
  data.table::fwrite(out_file, sep = ' ')

out_file <- file.path(outdir, "qtl_covariate.PLINK.$interaction_covar.txt") # Interaction
merged_tab %>%
  dplyr::select(dplyr::all_of(c("PID" = "lab_code", "IID" = "lab_code", "sex", "age", "$interaction_covar"))) %>%
  dplyr::mutate(sex = dplyr::if_else(sex == "Male", 1, 2), dplyr::across(\`$interaction_covar\`, ~dplyr::if_else(is.na(.x), -9, .x))) %>%
  data.table::fwrite(out_file, sep = ' ')

# Sample ID mapping file, used to reset the sample names in VCF file.
if (!file.exists("$vcf_sample_name_map"))
  merged_tab %>%
    dplyr::select(exact_VCF_header, lab_code) %>%
    data.table::fwrite("$vcf_sample_name_map", col.names = FALSE, sep = ' ')
EOF

bcftools reheader -s $vcf_sample_name_map $imputed_chr_filtered \
  | bcftools view -S <(cut -f2 -d' ' $vcf_sample_name_map) -O z -o $qtl_genotype_file --force-samples \
  && bcftools index -t $qtl_genotype_file

#
## QTL mapping
#
# QTLtools
# Covariants, male:1, female:2, using the `gwas` subcommand
pheno_dir=$projdir/outputs/300OB/phenotypes
qtl_dir=$projdir/outputs/300OB/qtl

bgzip $pheno_dir/qtl_phenotype.QTLtools.$phenotype.bed
tabix $pheno_dir/qtl_phenotype.QTLtools.$phenotype.bed.gz

# QTLtools GWAS mode original data
QTLtools gwas \
  --vcf $qtl_genotype_file \
  --bed $pheno_dir/qtl_phenotype.QTLtools.$phenotype.bed.gz \
  --cov $pheno_dir/qtl_covariate.QTLtools.$interaction_covar.txt \
  --out $qtl_dir/QTLtools/qtl.gwas.$phenotype.txt

# QTLtools GWAS mode normalized data
QTLtools gwas \
  --vcf $qtl_genotype_file \
  --bed $pheno_dir/qtl_phenotype.QTLtools.$phenotype.bed.gz \
  --cov $pheno_dir/qtl_covariate.QTLtools.$interaction_covar.txt \
  --out $qtl_dir/QTLtools/qtl.gwas.$phenotype.txt \
  --normal

# QTLtools trans mode normalized data, nominal mode
QTLtools trans \
  --vcf $qtl_genotype_file \
  --bed $pheno_dir/qtl_phenotype.QTLtools.$phenotype.bed.gz \
  --cov $pheno_dir/qtl_covariate.QTLtools.$interaction_covar.txt \
  --out $qtl_dir/QTLtools/qtl.trans.normal.nominal.$phenotype \
  --threshold 5e-2 \
  --normal \
  --nominal

# QTLtools trans mode normalized data, permute mode
QTLtools trans \
  --vcf $qtl_genotype_file \
  --bed $pheno_dir/qtl_phenotype.QTLtools.$phenotype.bed.gz \
  --cov $pheno_dir/qtl_covariate.QTLtools.$interaction_covar.txt \
  --out $qtl_dir/QTLtools/qtl.trans.normal.permute.$phenotype \
  --threshold 5e-2 \
  --normal \
  --permute


# PLINK association mode using linear model
# Model: BMI ~ G + age + gender
plink \
  --linear \
  --vcf $qtl_genotype_file \
  --pheno $pheno_dir/qtl_phenotype.PLINK.$phenotype.txt \
  --covar $pheno_dir/qtl_covariate.PLINK.txt \
  --out $qtl_dir/PLINK/qtl.linear.$phenotype \
  --allow-no-sex


# Model: BMI ~ G x Monocyte-related-traits + age + gender. PLINK
plink \
  --linear interaction \
  --vcf $qtl_genotype_file \
  --pheno $pheno_dir/qtl_phenotype.PLINK.$phenotype.txt \
  --covar $pheno_dir/qtl_covariate.PLINK.$interaction_covar.txt \
  --parameters 1-4, 7 \
  --out $qtl_dir/PLINK/qtl.linear.$phenotype.by_$interaction_covar \
  --allow-no-sex

#
## Plotting
#

# QQplot of QTLtools results.
Rscript ~/tools/QTLtools/script/plotTrans.R \
  $projdir/outputs/300OB/qtl/QTLtools/QQplot.$phenotype.pdf \
  $projdir/outputs/300OB/qtl/QTLtools/qtl.trans.normal.nominal.$phenotype.hits.txt.gz \
  $projdir/outputs/300OB/qtl/QTLtools/qtl.trans.normal.nominal.$phenotype.bins.txt.gz \
  $projdir/outputs/300OB/qtl/QTLtools/qtl.trans.normal.permute.$phenotype.hits.txt.gz \
  $projdir/outputs/300OB/qtl/QTLtools/qtl.trans.normal.permute.$phenotype.bins.txt.gz
