# Scripts and tools used in the analysis

<!--[[toc]]-->
[[toc]]

## 1 Introduction
The current repository deposits all customized R/Python/BASH code used for the analysis from **Zhang et. al. 2023** (doi: *TBA*) for the purpose of reproduction.
In the study, we identified single-cell expression quantitative trait loci (sceQTL) using a linear mixed model.
Downstream analysis identified regulatory programms in immune traits and common diseases, such trained immunity and COVID-19.

During immune cell activation, dynamic shifts in gene expression orchestrate cellular properties essential for effective functions.
Deciphering the intricate regulatory mechanisms underlying these processes holds promise for elucidating how genetic variants contribute to immune-related disorders.
We mapped genetic effects on gene expression (expression quantitative trait loci, eQTLs) using single-cell transcriptomes from 152 samples obtained from 38 healthy individuals, covering baseline state and lipopolysaccharide treatment either pre- or post- Bacillus Calmette-Guerin (BCG) vaccination.
Interestingly, we uncovered a monocyte eQTL linked to the LCP1, shedding light on inter-individual variation in trained immunity.
Furthermore, we elucidated genetic and epigenetic regulatory networks of SLFN5 and its pivotal role in COVID-19 pathogenesis by incorporating disease-associated loci, chromatin accessibility, and transcription factor binding affinities, aligning with the established function of SLFN5 in restricting virus replication during influenza infection.
Our study provides a paradigm to decipher the genetic underpinnings of complex traits by integrating single-cell eQTL with multi-omics data from patients and public databases.

## 2 Table of tools/packages used in the analysis

| Tool name         | Version    | Comments                                                                                         |
| :---------------- | ---------: | :----------------------------------------------------------------------------------------------- |
| BCFtools          | v1.12      | Genotype manipulation and quanlity control.                                                      |
| CellRanger        | v3.1.0     | Generate FASTQ files by demultiplexing BCL files and generate read counts from FASTQ files.      |
| PLINK             | v1.90b6.26 | Clumping, genotype manipulation, and kinship estimation                                          |
| SMR               | v1.3.1     | Summary-based Mendelian randomization                                                            |
| Souporcell        | v2.0       | scRNA-seq sample demultiplexing.                                                                 |
| limix_qtl         | -          | eQTL identification using linear mixed model.                                                    |
| Python            | v3.9.6     | Programming language.                                                                            |
| Python/matplotlib | v3.8.2     | Visualization                                                                                    |
| Python/numpy      | v1.26.2    | Dependency of matplotlib                                                                         |
| Python/pandas     | v2.1.3     | Dependency of matplotlib                                                                         |
| Python/pysam      | v0.18.0    | Python wrapper of HTSlib to work on VCF files.                                                   |
| R                 | v4.2.0     | Programming language.                                                                            |
| R/BoicManager     | v1.30.15   | R package management.                                                                            |
| R/CS-CORE         | v1.0.1     | Co-expression analysis.                                                                          |
| R/ComplexHeatmap  | v2.12.1    | Visualization.                                                                                   |
| R/Gviz            | v1.40.1    | Visualization.                                                                                   |
| R/MashR           | v0.2.57    | R implementation of multivariate adaptive shrinkage (MASH).                                      |
| R/PEER            | v1.3       | Bayesian approaches to infer hidden determinants, their effects from gene expression profiles.   |
| R/Seurat          | v4.1.1     | scRNA-seq procession, analysis, and cell type identification.                                    |
| R/TwoSampleMR     | v0.5.6     | Two sample Mendelian randomization.                                                              |
| R/clusterProfiler | v4.8.3     | Visualization and function enrichement of genes.                                                 |
| R/coloc           | v5.1.0.1   | Colocalization analysis by COLOC.                                                                |
| R/tidyverse       | v1.3.2     | Tidyverse for data manipulation and visualization. GGplot2, dplyr, tidyr, and purrr.             |
| R/SingleR         | v0.3.1     | Single-cell cell type identification.                                                            |

## 3 Singularity container to reproduce the analysis
To replicate the whole analysis, we packed the packages/tools used in the analysis into a `Singularity` container.
The configure file of the container is available in `singularity/sceQTL_300bcg.def`.
A precompiled `Singularity` container is available on [Zenodo](https://doi.org/10.5281/zenodo.xxx).

## 4 Analysis
The statistical analyses were performed by in-house R/Python script or published tools/packages.
The codes used to performed each analysis are listed and briefly introduced as the following.

### 4.1 Single-cell RNA-seq analysis and cell type identification
  - Read counts generated by `cellranger` pipeline
  - Demultiplexing by `souporcell`
  - Single-cell analysis by `Seurat`
  - Cell type identification by `SingleR`

For more information, please check [Li et al. 2023](https://doi.org/10.1016/j.celrep.2023.112487)

### 4.2 Identification of single-cell expression quantitative traits loci (sceQTL)
The eQTL of each cell type or conditions were estimated using a linear mixed model which is implemented in [`limix_qtl`](https://github.com/single-cell-genetics/limix_qtl).
An author recommended pipeline was implemented in `snakemake/300bcg_limix_qtl.smk` and `snakemake/scTI_Mono.smk` for general eQTL (consistent and response eQTL) and trained-immunity eQTL, respectively.
The input to `limix_qtl` was generated using a in-house R script which is implemented in `r/prepare_inputs_limix_eqtl.r`, while the kinship matrix was caculated by `PLINK` which is available in `bash/prepare_kinship_matrix.sh`

### 4.3 Estimation of shared signals
We used MASH (MASHR) to estimate the shared eQTL signals between each pair of eQTL summary statistics.
A strong-week strategy was used based on the [authors' recommendations](https://stephenslab.github.io/mashr/articles/eQTL_outline.html).
These steps are available by `r/mashr.r`.

### 4.4 Construction of reference genotype panel
To speed up the calculation of harmonization step, we prepared a reference panel using genetic variants from 1000 genome project.
The generated reference panel was also used to estimate enrichement of the identified QTL SNPs in epigenomic annotations.
The detailed steps is recorded in `bash/prepare_reference.sh`

### 4.5 Harmonization of genetic variants
To evaluate the biological function of the identified sceQTL, we harmonized these eVariants with publicly available GWAS summary statistics.
The harmonization steps are available in `r/harmonization.r`.

### 4.6 Summary based Mendelian randomization (SMR) analysis
We hypothesized that traits or diseases are mediated by the gene expressions.
This mediation effects was estimated using [SMR & HEIDI methods](https://yanglab.westlake.edu.cn/software/smr).
The `SMR` tool was used to prioritize genes underlying GWAS hits for follow-up functional studies.
These steps were implemented in `bash/smr.sh`.

### 4.7 Colocalization analysis
We used [COLOC](https://github.com/coloc/coloc) to identify the colocalized genes between sceQTL and GWAS loci.
The COLOC hypothesis shared common genetic casusal variants in the given region by the genetic colocalisation analysis of two potentially related phenotypes.
The statisitcal steps were implemented in `r/coloc.r`.

### 4.8 Single-cell Co-expression analysis
To determine the co-expression profiles we used [CS-CORE](https://github.com/ChangSuBiostats/CS-CORE), which is powered by [WGCNA](https://cran.r-project.org/web/packages/WGCNA/index.html).
`CS-CORE` is a R package for cell-type-specific co-expression inference from single cell RNA-sequencing data.
The co-expression inference steps were implemented in `r/plot_main_figures.r`.

### 4.9 Replication of the identified eQTLGen
To validate the identified eQTL (specificlly consistent eQTL), we checked their effect size in publicly available dataset [eQTLGen](https://www.eqtlgen.org).
Theese steps were implemented in `r/replication.r`.

### 4.10 Visualization codes
The figures in the manuscripts were generated using in-house R/Python codes. These codes are available in 
1. `r/PieDonut.R`
2. `r/plot_figure2.R`
3. `r/plot_main_figures.r`
4. `r/trackplot.r`
5. `py3/plot_asoc_rc.py`

### 4.11 Misc
Any additional scripts are recorded in the previous sections.
1. `r/collect_top_eqtl.r`. A `R` script to collect top eQTLs using permutation-passed method.
2. `bash/gwas_to_vcf.sh`. A `bash` script to convert GWAS summary statisitics to GWAS VCF format.


<!--
### 4.9 Enrichement of identified QTL SNPs in epigenomic annotations
We obtained cell-type-specific ATAC-seq peaks from three cell types, namely monocytes, CD8T, and NK cells, sorted from PBMC samples at d0 ("V1") in quantification*.csv.gz files.
`bash/enrich.sh` and `r/enrich.r`.
# Other
atac_seq.r, check_eqtl_per_condition.r, check_longcovid.r, coexpression.r, convert_seurat_into_loom.r, reformat_gse184269.r, explore_coloc.r
-->


## 5 Publicly available reference data

### 5.1 Table of reference data
| Reference name   | Version                | Availability                                                                                       |
| :--------------- | :--------------------- | :------------------------------------------------------------------------------------------------- |
| ADASTRA          | V5.1.3                 | https://adastra.autosome.org/bill-cipher                                                           |
| DICE             | Build Feb 25, 2022     | https://dice-database.org                                                                          |
| GTEx             | Release V8             | https://www.gtexportal.org/home                                                                    |
| GeneHancer       | January 2019 (V2)      | https://www.genecards.org                                                                          |
| Genetic variants | GRCh38p7/build 151     | https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz                  |
| Genomic features | Release 41 (GRCh38p13) | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz |
| eQTLGen          | Release 2019-12-23     | https://www.eqtlgen.org/phase1.html                                                                |

<!-- | Human genome     | GRCh38              |              | -->

### 5.2 Table of publicly available GWAS summary statistics
| Trait Category | Trait                      | ID                                   | TotalVariants | StudyType   | PMID   | SampleSize | CaseSize | ControlSize |
| :------------  | :------------------------- | :----------------------------------- | ------------: | :---------- | -----: | ---------: | -------: | ----------: |
| Autoimmune     | AlzheimerDiseases          | ieu-b-5067                           | 12321875      | CaseControl |        | NA         | 487331   | 954         |
| Autoimmune     | AsthmaDT2                  | UKB-a-255                            | 10545186      | Continuous  |        | NA         | NA       | NA          |
| Autoimmune     | AtopicDermatitis           | ieu-a-996                            | 9965822       | Continuous  |        | NA         | NA       | NA          |
| Autoimmune     | CrohnsDisease              | ieu-a-30                             | 11112707      | Continuous  |        | NA         | NA       | NA          |
| Autoimmune     | GoutDisease                | UKB-a-107                            | 10545186      | Continuous  |        | NA         | NA       | NA          |
| Autoimmune     | InflammatoryBowelDisease   | ieu-a-31                             | 11534250      | Continuous  |        | NA         | NA       | NA          |
| Autoimmune     | MultipleSclerosis          | UKB-b:17670                          | 9851866       | CaseControl |        | NA         | 461254   | 1679        |
| Autoimmune     | Psoriasis                  | UKB-a-100                            | 10545186      | Continuous  |        | NA         | NA       | NA          |
| Autoimmune     | RheumatoidArthritis        | ukb-d-M13_RHEUMA                     | 10079898      | Continuous  |        | NA         | NA       | NA          |
| Autoimmune     | SystemicLupusErythematosus | EBI-a-GCST003156                     | 7071163       | CaseControl |        | NA         | 9066     | 5201        |
| Autoimmune     | Type2Diabetes              | EBI-a-GCST006867                     | 5030727       | CaseControl |        | NA         | 1178     | 61714       |
| Autoimmune     | UlcerativeColitis          | ieu-a-32                             | 11094790      | Continuous  |        | NA         | NA       | NA          |
| Cancer         | BladderCancer              | ieu-b-4874                           | 9904926       | CaseControl |        | NA         | 372016   | 1279        |
| Cancer         | ColorectalCancer           | ieu-b-4965                           | 11738639      | CaseControl |        | NA         | 372016   | 5657        |
| Cancer         | LungCancer                 | ieu-a-966                            | 8827545       | Continuous  |        | NA         | NA       | NA          |
| Cancer         | OvarianCancer              | ieu-b-4963                           | 9822229       | CaseControl |        | NA         | 198523   | 1218        |
| Cancer         | ProstateCancer             | ieu-b-4809                           | 12097504      | CaseControl |        | NA         | 173493   | 9132        |
| Cancer         | ThyroidCancer              | ieu-a-1082                           | 571227        | Continuous  |        | NA         | NA       | NA          |
| General        | BodyMassIndex              | ieu-b-40                             | 2336260       | Continuous  |        | NA         | NA       | NA          |
| General        | DiastolicBloodPressure     | ieu-b-39                             | 7160619       | Continuous  |        | NA         | NA       | NA          |
| General        | HDLCholesterol             | ieu-b-109                            | 12321875      | Continuous  |        | NA         | NA       | NA          |
| General        | Height                     | ieu-a-89                             | 2546513       | Continuous  |        | NA         | NA       | NA          |
| General        | LDLCholesterol             | ieu-b-110                            | 12321875      | Continuous  |        | NA         | NA       | NA          |
| General        | SystolicBloodPressure      | ieu-b-38                             | 7088083       | Continuous  |        | NA         | NA       | NA          |
| General        | Triglycerides              | ieu-b-111                            | 12321875      | Continuous  |        | NA         | NA       | NA          |
| General        | Urate                      | met-a-345                            | 2544432       | Continuous  |        | NA         | NA       | NA          |
| Infection      | COVID19Release7            | HGI_A2_ALL_eur_leave23andme_20220403 | 12195847      | Continuous  |        | NA         | NA       | NA          |
| Others         | CoronaryHeartDisease       | ieu-a-7                              | 8479693       | Continuous  |        | NA         | NA       | NA          |
| Others         | ObesityClass1              | ieu-a-90                             | 2377582       | Continuous  |        | NA         | NA       | NA          |
| Others         | ObesityClass2              | ieu-a-91                             | 2328688       | Continuous  |        | NA         | NA       | NA          |
| Others         | ObesityClass3              | ieu-a-92                             | 2248088       | Continuous  |        | NA         | NA       | NA          |
| Others         | Schizophrenia              | ieu-b-42                             | 14237637      | CaseControl |        | NA         | 43456    | 33640       |

<!-- ### Identification of co-expression modules
  We used weighted gene co-expression network analysis (WGCNA) to detect co-expression modules for each transcriptomic profiles of four identified cell types.
  For more details about the analysis check R code at `r/wgcna.r`
-->


### 5.3 SingleR reference
| Database                            | Species | Comments                                                                                                  |
| :---------------------------------- | :------ | :-------------------------------------------------------------------------------------------------------- |
| Monaco immnune data                 | Human   | RNA-Seq profiling on 29 immune cell types consituting PBMCs sorted from 4 Singaporean-Chinese individuals |
| Human Primary Cell Atlas            | Human   | Microarray datasets derived from human primary cells                                                      |
| Blueprint / ENCODE                  | Human   | Bulk RNA-seq data for pure stroma and immune cells generated by Blueprint                                 |
| Database of Immnune Cell Expression | Human   | Bulk RNA-seq samples of sorted human immune cell populations                                              |

## 6 Check list
Scripts used to generate single-cell eQTL from 300BCG cohort.
  - [x] bash/estimate_independent_gwas_loci.sh
  - [x] bash/prepare_kinship_matrix.sh
  - [x] bash/gwas_to_vcf.sh
  - [x] bash/smr.sh
  - [x] py3/fetch_snp_info.py
  - [x] py3/plot_asoc_rc.py
  - [x] r/PieDonut.R
  - [x] r/harmonization.r
  - [ ] r/mashr.r
  - [ ] r/plot_main_figures.r
  - [x] r/prepare_inputs_limix_eqtl.r
  - [ ] r/replication.r
  - [ ] r/trackplot.r
  - [x] snakemake/300bcg_limix_qtl.smk
  - [ ] snakemake/scTI_Mono.smk
  - [ ] singularity/sceQTL_300bcg.def
