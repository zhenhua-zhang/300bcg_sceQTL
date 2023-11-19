# Scripts and tools used in the analysis

[[toc]]

## 1 Introduction
The repository deposits all customized R/Python/BASH code used for the analysis by study **Zhang et. al. 2023** (doi: *TBA*) for the purpose of reproduction.
In the study, we identified single-cell expression quantitative trait loci (sceQTL) using a linear mixed model.
Downstream analysis identified regulatory programms in immune traits and common diseases, such trained immunity and COVID-19.

## 2 Table of tools/packages used in the analysis

| Tool name  | Version | Comments |
| :--------- | ------: | :------- |
| PLINK      |         | Clumping, genotype manipulation, and kinship estimation |
| BCFtools   |         | Genotype manipulation and quanlity control |
| Souporcell |         | scRNA-seq sample demultiplexing       |
| limix_qtl  |         | eQTL identification |
| Seurat     |         | scRNA-seq procession, analysis, and cell type identification |
| Python     | 3.9.6   | Programming language |
| R          | 4.2.0   | Programming language |

## 3 Singularity container to reproduce the analysis
To replicate the whole analysis, we packed the packages/tools used in the analysis into a `Singularity` container.
The configure file of the container is available in `singularity/sceQTL_300bcg.def`.
A precompiled `Singularity` container is available in on [Zenodo](https://doi.org/10.5281/zenodo.xxx).

## 4 Analysis
The statistical analyses were performed by in-house R/Python script or published tools/packages.
The codes used to performed each analysis are listed and briefly introduced as the following.

### 4.1 Single-cell RNA-seq analysis and cell type identification
  - Demultiplexing by `souporcell`
  - Single-cell analysis by `Seurat`
  - Cell type identification by `?`

### 4.2 Identification of single-cell expression quantitative traits loci (sceQTL)
The eQTL of each cell type or conditions were estimated using a linear mixed model which is implemented in [`limix_qtl`](https://github.com/single-cell-genetics/limix_qtl).
An author recommended pipeline was implemented in `snakemake/300bcg_limix_qtl.smk`.

### 4.3 Estimation of shared signals
We used MASH (MASHR) to estimate the shared eQTL signals between each pair of eQTL summary statistics.
These steps are available by `r/mashr.r`.

### 4.4 Construction of reference genotype panel
To speed up the calculation of harmonization step, we prepared a reference panel using genetic variants from 1000 genome project.
The generated reference panel was also used to estimate enrichement of the identified QTL SNPs in epigenomic annotations.
The detailed steps is recorded in `bash/prepare_reference.sh`

### 4.5 Harmonization of genetic variants
To evaluate the biological function of the identified sceQTL, we harmonized these eVariants with publicly available GWAS summary statistics.
The harmonization steps are available in `r/harmonization.r`.

### 4.6 Summary based Mendelian randomization (SMR) analysis
`bash/smr.sh`

### 4.7 Colocalization analysis
`r/coloc.r`

### 4.8 Enrichement of identified QTL SNPs in epigenomic annotations
We obtained cell-type-specific ATAC-seq peaks from three cell types, namely monocytes, CD8T, and NK cells, sorted from PBMC samples at d0 ("V1") in quantification*.csv.gz files.
`bash/enrich.sh` and `r/enrich.r`.

## 5 Publicly available reference data

### 5.1 Table of reference data
| Reference name   | Version          | Availability | Comments |
| :--------------: | :--------------: | :----------: | :------: |
| Human genome     | GRCh38           |              |          |
| Genetic variants | GRCh38/build 108 |              |          |

### 5.2 Table of publicly available GWAS summary statistics
| Trait Category | Trait                      | ID                                   | TotalVariants | StudyType   | SampleSize | CaseSize | ControlSize |
| :------------  | :------------------------- | :----------------------------------- | ------------: | ----------: | ---------: | -------: | ----------: |
| Autoimmune     | AlzheimerDiseases          | ieu-b-5067                           | 12321875      | CaseControl | NA         | 487331   | 954         |
| Autoimmune     | AsthmaDT2                  | UKB-a-255                            | 10545186      | Continuous  | NA         | NA       | NA          |
| Autoimmune     | AtopicDermatitis           | ieu-a-996                            | 9965822       | Continuous  | NA         | NA       | NA          |
| Autoimmune     | CrohnsDisease              | ieu-a-30                             | 11112707      | Continuous  | NA         | NA       | NA          |
| Autoimmune     | GoutDisease                | UKB-a-107                            | 10545186      | Continuous  | NA         | NA       | NA          |
| Autoimmune     | InflammatoryBowelDisease   | ieu-a-31                             | 11534250      | Continuous  | NA         | NA       | NA          |
| Autoimmune     | MultipleSclerosis          | UKB-b:17670                          | 9851866       | CaseControl | NA         | 461254   | 1679        |
| Autoimmune     | Psoriasis                  | UKB-a-100                            | 10545186      | Continuous  | NA         | NA       | NA          |
| Autoimmune     | RheumatoidArthritis        | ukb-d-M13_RHEUMA                     | 10079898      | Continuous  | NA         | NA       | NA          |
| Autoimmune     | SystemicLupusErythematosus | EBI-a-GCST003156                     | 7071163       | CaseControl | NA         | 9066     | 5201        |
| Autoimmune     | Type2Diabetes              | EBI-a-GCST006867                     | 5030727       | CaseControl | NA         | 1178     | 61714       |
| Autoimmune     | UlcerativeColitis          | ieu-a-32                             | 11094790      | Continuous  | NA         | NA       | NA          |
| Cancer         | BladderCancer              | ieu-b-4874                           | 9904926       | CaseControl | NA         | 372016   | 1279        |
| Cancer         | ColorectalCancer           | ieu-b-4965                           | 11738639      | CaseControl | NA         | 372016   | 5657        |
| Cancer         | LungCancer                 | ieu-a-966                            | 8827545       | Continuous  | NA         | NA       | NA          |
| Cancer         | OvarianCancer              | ieu-b-4963                           | 9822229       | CaseControl | NA         | 198523   | 1218        |
| Cancer         | ProstateCancer             | ieu-b-4809                           | 12097504      | CaseControl | NA         | 173493   | 9132        |
| Cancer         | ThyroidCancer              | ieu-a-1082                           | 571227        | Continuous  | NA         | NA       | NA          |
| General        | BodyMassIndex              | ieu-b-40                             | 2336260       | Continuous  | NA         | NA       | NA          |
| General        | DiastolicBloodPressure     | ieu-b-39                             | 7160619       | Continuous  | NA         | NA       | NA          |
| General        | HDLCholesterol             | ieu-b-109                            | 12321875      | Continuous  | NA         | NA       | NA          |
| General        | Height                     | ieu-a-89                             | 2546513       | Continuous  | NA         | NA       | NA          |
| General        | LDLCholesterol             | ieu-b-110                            | 12321875      | Continuous  | NA         | NA       | NA          |
| General        | SystolicBloodPressure      | ieu-b-38                             | 7088083       | Continuous  | NA         | NA       | NA          |
| General        | Triglycerides              | ieu-b-111                            | 12321875      | Continuous  | NA         | NA       | NA          |
| General        | Urate                      | met-a-345                            | 2544432       | Continuous  | NA         | NA       | NA          |
| Infection      | COVID19Release7            | HGI_A2_ALL_eur_leave23andme_20220403 | 12195847      | Continuous  | NA         | NA       | NA          |
| Infection      | COVID19Release7            | HGI_A2_ALL_leave_23andme_20220403    | 11732502      | Continuous  | NA         | NA       | NA          |
| Others         | CoronaryHeartDisease       | ieu-a-7                              | 8479693       | Continuous  | NA         | NA       | NA          |
| Others         | ObesityClass1              | ieu-a-90                             | 2377582       | Continuous  | NA         | NA       | NA          |
| Others         | ObesityClass2              | ieu-a-91                             | 2328688       | Continuous  | NA         | NA       | NA          |
| Others         | ObesityClass3              | ieu-a-92                             | 2248088       | Continuous  | NA         | NA       | NA          |
| Others         | Schizophrenia              | ieu-b-42                             | 14237637      | CaseControl | NA         | 43456    | 33640       |

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
  - [x] bash/smr.sh
  - [x] py3/fetch_snp_info.py
  - [x] py3/plot_asoc_rc.py
  - [x] r/PieDonut.R
  - [ ] r/enrich.r
  - [x] r/harmonization.r
  - [ ] r/mashr.r
  - [ ] r/plot_main_figures.r
  - [ ] r/prepare_inputs_limix_eqtl.r
  - [ ] r/replication.r
  - [ ] r/trackplot.r
  - [x] snakemake/300bcg_limix_qtl.smk
  - [ ] singularity/sceQTL_300bcg.def
