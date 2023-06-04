# Scripts or tools used in the analysis

## Introduction
The repository deposits all customized R/Python/BASH code used for the analysis by study *xxx et. al.* (DOI: xxx) for the purpose of reproduction.

## Table of tools/packages used in the analysis

| Tool name | Version | Comments |
| :-------: | :-----: | :------: |
| PLINK     |         |          |
| BCFtools  |         |          |
| TODO      |         |          |


## Statistical analysis
The statistical analyses were performed by in-house R/Python script or published tools/packages.
The codes used to performed each analysis are listed and briefly introduced as the following.

### Single-cell RNA-seq analysis and cell type identification
1. Demultiplexing by `souporcell`
2. Single-cell analysis by `Seurat`
3. Cell type identification by `?`

### Identification of co-expression modules

We used weighted gene co-expression network analysis (WGCNA) to detect co-expression modules for each transcriptomic profiles of four identified cell types.
For more details about the analysis check R code at `r/wgcna.r`

### Definition of pseudo-time windows
`r/pseudotime.r` **TODO**

### Identification of single-cell expression quantitative traits loci (sceQTL)
`snakemake/300bcg_limix_qtl.smk`

### Estimation of shared signals
`r/mashr.r`

### Construction of reference genotype panel
To speed up the calculation of harmonization step, we prepared a reference panel using genetic variants from 1000 genome project.
The generated reference panel was also used to estimate enrichement of the identified QTL SNPs in epigenomic annotations.
The detailed steps is recorded in `bash/prepare_reference.sh`

### Harmonization of genetic variants
`r/harmonization.r`

### Mendelian randomization analysis
`bash/smr.sh`

### Colocalization analysis
`r/coloc.r`

### Enrichement of identified QTL SNPs in epigenomic annotations
We obtained cell-type-specific ATAC-seq peaks from three cell types, namely monocytes, CD8T, and NK cells, sorted from PBMC samples at d0 ("V1") in quantification*.csv.gz files.
`bash/enrich.sh` and `r/enrich.r`.


## Publicly available data and reference

### Table of publicly available GWAS summary statistics

### Table of reference data

| Reference name | Version | Availability | Comments |
| :------------: | :-----: | :----------: | :------: |
| Human genome   | GRCh38  |              |          |
| Genetic variants | GRCh38/build xxx | xxx | |


## Misc

SingleR reference

| Database                            | Species | Comments                                                                                                  |
| :---------------------------------: | :-----: | :-------------------------------------------------------------------------------------------------------: |
| Monaco immnune data                 | Human   | RNA-Seq profiling on 29 immune cell types consituting PBMCs sorted from 4 Singaporean-Chinese individuals |
| Human Primary Cell Atlas            | Human   | Microarray datasets derived from human primary cells                                                      |
| Blueprint / ENCODE                  | Human   | Bulk RNA-seq data for pure stroma and immune cells generated by Blueprint                                 |
| Database of Immnune Cell Expression | Human   | Bulk RNA-seq samples of sorted human immune cell populations                                              |


## Database

meta_information.sqlite, meta information about SNPs, GWAS, features, QTL mapping parameters
snp.sqlite, chrom, postion
feature.sqlite
eqtl.sqlite, per cell type per condition
harmonization.sqlite, per GWAS per table
colocalization.sqlite
mendelian_randomization.sqlite

MD5 sum per sqtlite DB

SNPs associated with IMPDH1
rs2288553, rs11770116, rs2288548, rs2288549, rs4731448, rs2278293, rs2278294, rs2228075


## Example eQTL

### rs7543002-CD55

**Reference**:

1. Upregulation of CD55 complement regulator in distinct PBMC subpopulations of COVID-19 patients is associated with suppression of interferon responses
2. Costimulation via CD55 on Human CD4+ T Cells Mediated by CD97
3. Beyond the Role of CD55 as a Complement Component
4. CD55 Deficiency, Early-Onset Protein-Losing Enteropathy and Thrombosis
5. Pathogenic implications for autoimmune mechanisms derived by comparative eQTL analysis of CD4+ versus CD8+ T cells
6. DNA methylation signatures reveal that distinct combinations of transcription factors specify human immune cell epigenetic identity. Methylation data.


## TODO

### 2023 Jan 26
- [x] ADCY3 ~ BMI + co-variables in 38 300BCG individuals with single-cell RNA-seq, Spearman's rho/linear mixed model.
      NOTE: no significant results
- [x] BMI ~ Genotype + co-variables in whole 300BCG cohort.

### 2023 Feb 14
- [x] Convert COVID-19 GWAS summary statistic into GWAS VCF format. Clump the SS and decide the top SNPs.

### 2023 Feb 17
- [x] Co-expression between all genes and CD55
- [x] Cell type specific eQTL of COVID-19 associated SNPs.
- [x] *trans*-eQTL for all CD55-coexpressed genes, where the SNPs are significant (?) associated with CD55. Only one gene in CD8T passed 5e-5 threshold.

### 2023 Feb 19
- [x] Trans-eQTL of CD55, where the trans loci are 1M closed to the TF that are potentially regulating CD55 expression by GeneCards.
- [x] Correlation between trans-eQTL Z-scores and CD55 cis-eQTL Z-scores, co-localization-like plot.

### 2023 Feb 21
- [-] Peak-to-gene link by MHH50 data or Hi-C by Javierre el al [PMID: 27863249](https://pubmed.ncbi.nlm.nih.gov/27863249)
- [ ] Methylation changes along the gene body and surroundings by Roy et al [PMID: 34706222](https://pubmed.ncbi.nlm.nih.gov/34706222)
- [ ] A table to collect all information about CD55/KANSL1 associated SNPs
- [ ] Phe-WAS for rs7543002 using Finngen dataset


### Candidate examples

1. Human SLFN5 is a transcriptional co-repressor of STAT1-mediated interferon responses and promotes the malignant phenotype in glioblastoma
