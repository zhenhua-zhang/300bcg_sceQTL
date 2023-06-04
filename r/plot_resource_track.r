#!/usr/bin/env Rscript
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gamil.com
# Created: Aug 01, 2021
# Updated: Aug 11, 2021

options(error=traceback, stringsAsFactors=F)

# To plot the track plots

# Packages
library(Gviz)
library(GenomicInteractions)
library(rtracklayer)
library(AnnotationHub)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(biomaRt)
library(GenomicRanges)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)

library(data.table)
library(tidyverse)
library(magrittr)


# Workingspace
wkdir <- "../../outputs/scATAC-seq/enrichment"
version <- "ver_beta"

# SNP 2
rsid <- "rs6800484"
snp_contig <- 3
snp_position <- 46089651
eqtl_gene <- "CCR1"
# eqtl_gene <- "FYCO1"
# 
# # SNP 3
# rsid <- "rs10420007"
# snp_contig <- 19
# snp_position <- 4724057
# eqtl_gene <- "DPP9"


# SNP information
rsid <- "rs7255545"
snp_contig <- 19
snp_position <- 4724722
eqtl_gene <- "DPP9"

flank_len <- 330000
range_start <- snp_position - flank_len
range_end <- snp_position + flank_len


# Ideogram track
ideog_track <- IdeogramTrack(genome="hg38", chromosome=as.character(snp_contig))


# Genomic coordnation track
gaxis_track <- GenomeAxisTrack()


# ATAC-seq depth
atac_dep_file <- str_glue("{wkdir}/{version}/integration_atac_dep_track.csv")
atac_dep_track <- DataTrack(name="ATAC-seq\nMonocytes", background.panel="gray95", cex.axis=0.5)
if (file.exists(atac_dep_file)) {
  # NOTE: 1. When the depth <= 0.01, it's ground to 0.
  atac_dep_tab <- fread(atac_dep_file, verbose=F, showProgress=F)
  atac_dep_track <- atac_dep_tab %>%
    dplyr::filter(contig==snp_contig
                  & between(position, range_start, range_end)
                  & celltype=="cMono") %>%
    dplyr::mutate(depth=case_when(depth<=0.01 ~ 0, T ~ depth)) %>%
    makeGRangesFromDataFrame(seqnames.field="contig", start.field="position",
                             end.field="position", keep.extra.columns=T) %>%
    DataTrack(type="horizon", name=str_glue("ATAC-seq\nMonocytes"), background.panel="white",
              col.horizon="darkred", fill.horizon="darkred")
} else {
  message("Empty ATAC-seq track")
}


# Gene to peak links track
g2p_link_file <- str_glue("{wkdir}/{version}/integration_g2p_link_track.csv")
if (file.exists(g2p_link_file)) {
  g2p_link_tab <- fread(g2p_link_file)
} else {
  g2p_link_tab <- fread("../../inputs/scATAC-seq/p2g.merge.txt") %>%
    tidyr::separate(Peak, into=c("peak_contig", "peak_start", "peak_end"), sep="[:-]") %>%
    tidyr::separate(GeneLoci, into=c("gene_contig", "gene_start"), sep=":") %>%
    dplyr::mutate(gene_contig=str_remove(pattern="chr", gene_contig),
                  gene_start=as.integer(gene_start),
                  gene_end=gene_start,
                  peak_contig=str_remove(pattern="chr", peak_contig),
                  peak_start=as.integer(peak_start),
                  peak_end=as.integer(peak_start)) %>%
    dplyr::relocate(gene_end, .after="gene_start") %>%
    dplyr::select(-dist) %>%
    dplyr::filter(Corr>=0.45 & !peak_contig%in%c("X")) %T>%
    (function(x) {fwrite(x, g2p_link_file, sep=",", quote=F, row.names=F)})
}
g2p_link_tab %<>% dplyr::filter(Gene==eqtl_gene)
g2p_link_title <- str_glue("Peak to gene\n{eqtl_gene}")
if (nrow(g2p_link_tab)!=0) {
  # Peak anchors
  anchor_peak.0 <- g2p_link_tab %>%
    dplyr::select(peak_contig:peak_end) %>%
    makeGRangesFromDataFrame(seqnames.field="peak_contig", start.field="peak_start", end.field="peak_end", keep.extra.columns=T)
  # Gene anchor
  anchor_gene <- g2p_link_tab %>%
    dplyr::select(gene_contig:Corr) %>%
    makeGRangesFromDataFrame(seqnames.field="gene_contig", start.field="gene_start", end.field="gene_end", keep.extra.columns=T)
  g2p_link_track <- GenomicInteractions(anchor_peak.0, anchor_gene, counts=g2p_link_tab$Corr) %>%
    InteractionTrack(name=g2p_link_title, chromosome=as.character(snp_contig))
  displayPars(g2p_link_track) <- list(plot.outside=T, col.outside="gray",
                                      col.anchors.fill="darkblue",
                                      col.anchors.line="darkblue",
                                      col.interactions="darkblue",
                                      background.panel="gray95",
                                      interaction.dimension="height",
                                      interaction.measure ="counts",
                                      anchor.height=0.05)
} else {
  g2p_link_track <- AnnotationTrack(name=g2p_link_title)
}


# Hi-C track (public, the coordnation should be lifted to GRCh38)
hic_link_file <- str_glue("{wkdir}/{version}/integration_hic_link_track.csv")
if (file.exists(hic_link_file)) {
  hic_link_tab <- fread(hic_link_file)
} else {
  ahub <- AnnotationHub() # snapshotDate(): 2020-10-27
  ahub_chain <- subset(ahub, rdataclass=="ChainFile" & species=="Homo sapiens")
  ahub_chain <- ahub_chain[ahub_chain$title=="hg19ToHg38.over.chain.gz"][[1]]

  # Format the Hi-C dataset
  hic_link_tab <- fread("../../inputs/Hi-C/2016_cell_Javierre-etal_lineage-specific/PCHiC_peak_matrix_cutoff5.tsv") %>%
    dplyr::select(-one_of("Mac0", "Mac1", "Mac2", "Neu", "MK", "EP", "Ery",
                          "FoeT", "aCD4", "naCD4", "nCD4", "nCD8", "nB",
                          "clusterPostProb", "clusterID", "dist")) %>%
    tidyr::pivot_longer(cols=Mon:tB, values_to="score", names_to="celltype") %>%
    dplyr::mutate(baitChr=paste0("chr", baitChr), oeChr=paste0("chr", oeChr))

  # Liftover from GRCh37(hg19) to GRCh38 (UCSC)
  hic_link_tab <- hic_link_tab %>%
    makeGRangesFromDataFrame(seqnames.field="baitChr", start.field="baitStart", end.field="baitEnd", keep.extra.columns=T) %>%
    liftOver(ahub_chain) %>%
    unlist() %>%
    as.data.frame() %>%
    dplyr::select(-c(width, strand)) %>%
    dplyr::rename("baitChr"="seqnames", "baitStart"="start", "baitEnd"="end") %>%
    makeGRangesFromDataFrame(seqnames.field="oeChr", start.field="oeStart", end.field="oeEnd", keep.extra.columns=T) %>%
    liftOver(ahub_chain) %>%
    unlist() %>%
    as.data.frame() %>%
    dplyr::select(-c(width, strand)) %>%
    dplyr::rename("oeChr"="seqnames", "oeStart"="start", "oeEnd"="end") %>%
    dplyr::distinct() %>%
    dplyr::select(baitChr:baitName, oeChr:oeEnd, oeID:score) %>%
    dplyr::filter(!(baitChr %in% c("chrY", "chrX", "chrMT")), !(oeChr %in% c("chrY", "chrX", "chrMT"))) %>%
    dplyr::mutate(baitChr=str_remove(baitChr, "chr"), oeChr=str_remove(oeChr, "chr")) %>%
    dplyr::mutate(celltype=case_when(celltype=="Mon"  ~ "Mono",
                                     celltype=="tCD4" ~ "CD4T",
                                     celltype=="tCD8" ~ "CD8T",
                                     celltype=="tB"   ~ "B")) %>%
    (function(x) { fwrite(x, hic_link_file, sep="\t", quote=F, row.names=F); x })
  # removeCache(ahub) # Remove cache to give more space.
}
# Only use data for Monocytes
hic_link_tab <- hic_link_tab %>% dplyr::filter(celltype %in% c("Mono"))
# Otehr anchors not overlpped with the target SNP.
ext_hic_link_tab <- hic_link_tab %>%
  dplyr::filter(baitChr==snp_contig
                & !((baitStart<=snp_position & snp_position<=baitEnd)
                   | (oeStart<=snp_position & snp_position<=oeEnd)))
# First anchor
anchor_bait.0 <- ext_hic_link_tab %>%
  dplyr::select(baitChr:baitName) %>%
  makeGRangesFromDataFrame(seqnames.field="baitChr", start.field="baitStart", end.field="baitEnd", keep.extra.columns=T)
# Second anchor
anchor_oe.0 <- ext_hic_link_tab %>%
  dplyr::select(oeChr:score) %>%
  makeGRangesFromDataFrame(seqnames.field="oeChr", start.field="oeStart", end.field="oeEnd", keep.extra.columns=T)
# The interaction link
hic_link_track.0 <- GenomicInteractions(anchor_bait.0, anchor_oe.0, counts=ext_hic_link_tab$score) %>%
  InteractionTrack(name="PCHiC\nMonocytes", chromosome=as.character(snp_contig))
# Display parameters
displayPars(hic_link_track.0) <- list(plot.outside=T,
                                      col.outside="gray",
                                      col.interactions="lightblue",
                                      col.anchors.line="gray",
                                      interaction.dimension="height",
                                      interaction.measure ="counts",
                                      anchor.height=0.1)
# Anchors overlpped with the target SNP.
snp_hic_link_tab <- hic_link_tab %>%
  dplyr::filter(baitChr==snp_contig
                & ((baitStart<=snp_position & snp_position<=baitEnd)
                   | (oeStart<=snp_position & snp_position<=oeEnd)))
# First anchor
anchor_bait.1 <- snp_hic_link_tab %>%
  dplyr::select(baitChr:baitName) %>%
  makeGRangesFromDataFrame(seqnames.field="baitChr", start.field="baitStart", end.field="baitEnd", keep.extra.columns=T)
# Second anchor
anchor_oe.1 <- snp_hic_link_tab %>%
  dplyr::select(oeChr:score) %>%
  makeGRangesFromDataFrame(seqnames.field="oeChr", start.field="oeStart", end.field="oeEnd", keep.extra.columns=T)
# The interaction link
hic_link_track.1 <- GenomicInteractions(anchor_bait.1, anchor_oe.1, counts=snp_hic_link_tab$score) %>%
  InteractionTrack(chromosome=as.character(snp_contig))
# Display parameters
displayPars(hic_link_track.1) <- list(plot.outside=T, col.outside="gray",
                                      col.anchors.fill="red",
                                      col.anchors.line="black",
                                      col.interactions="red",
                                      interaction.dimension="height",
                                      interaction.measure ="counts",
                                      anchor.height=0.15)
hic_link_track <- OverlayTrack(list(hic_link_track.0, hic_link_track.1))


# Whole blood eQTL track (GTEx or eQTLGen)
which_eqtl_set <- "eqtlgen"
wbl_eqtl_file <- str_glue("{wkdir}/{version}/integration_wbl_eqtl_track.csv")
if (file.exists(wbl_eqtl_file)) {
  wbl_eqtl_tab <- fread(wbl_eqtl_file)
} else if (which_eqtl_set == "gtex") {
  # Fetch a data.frame of gene features
  # Here we use GENCODE basic and Ensembl Canonical transcripts
  # NOTE: 1. Please note, only protein coding genes on autosomals were included.
  #       2. No gene symbol available in GTEx dataset, the symbols were fetch by biomaRt.
  # TODO: 1. No rsid available, perhaps it's better to fetch rsid as well.
  ensembl <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl")
  gene_pos_tab <- getBM(attributes=c("chromosome_name", "ensembl_gene_id", "external_gene_name"),
                        filters=c("chromosome_name", "biotype", "transcript_gencode_basic", "transcript_is_canonical"),
                        values=list(as.character(1:22), c("protein_coding"), T, T),
                        mart=ensembl)

  wbl_eqtl_tab <- fread("../../inputs/GTEx/Whole_Blood.v8.EUR.signif_pairs.txt.gz") %>%
    tidyr::separate(variant_id, into=c("contig", "position", "refAllele", "altAllele", "build")) %>%
    dplyr::filter(nchar(refAllele)==1 & nchar(altAllele)==1 & maf>=0.1) %>%
    dplyr::mutate(phenotype_id=mapply(function(x) str_split(x, "\\.", simplify=T)[1, 1], phenotype_id),
                  position=as.integer(position),
                  contig=mapply(function(x) as.integer(gsub("chr", "", x)), x=contig)) %>%
    dplyr::inner_join(gene_pos_tab, by=c("phenotype_id"="ensembl_gene_id")) %>%
    dplyr::arrange(contig, position) %>%
    dplyr::rename("p_value_eqtl"="pval_nominal") %>%
    dplyr::select(contig, position, p_value_eqtl, external_gene_name) %>%
    dplyr::distinct() %>%
    (function(x) {fwrite(x, wbl_eqtl_file, sep=",", quote=F, row.names=F); x})
} else if (which_eqtl_set == "eqtlgen") {
  # NOTE: 1. Bonferroni P-value were used as p_value_eqtl.
  #       2. SNPs with non rsid were removed
  #       3. eQTLGen genome build is GRCh37, they were converted o GRCh38 by rsid.
  snps_b38 <- SNPlocs.Hsapiens.dbSNP144.GRCh38
  wbl_eqtl_tab <- fread("../../inputs/eQTLGen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz") %>%
    dplyr::filter(grepl(pattern="^rs", SNP) & (!SNPChr %in% c("X"))) %>%
    dplyr::select(one_of("SNP", "GeneSymbol", "BonferroniP"))

  wbl_eqtl_tab <- wbl_eqtl_tab %$% SNP %>%
    snpsById(snps_b38, ids=., ifnotfound="drop") %>%
    as.data.frame() %>%
    dplyr::inner_join(wbl_eqtl_tab, by=c("RefSNP_id"="SNP")) %>%
    dplyr::select(-one_of(c("strand", "alleles_as_ambig"))) %>%
    dplyr::rename("contig"="seqnames", "position"="pos",
                  "p_value_eqtl"="BonferroniP",
                  "external_gene_name"="GeneSymbol") %>%
    dplyr::mutate(contig=as.integer(contig), position=as.integer(position)) %>%
    dplyr::distinct() %T>%
    (function(x) {fwrite(x, wbl_eqtl_file, sep=",", quote=F, row.names=F)})
}
# Define y lims
egene_pv <- wbl_eqtl_tab %>%
  dplyr::filter(contig==snp_contig
                & between(position, range_start, range_end)
                & external_gene_name==eqtl_gene) %>%
  dplyr::mutate(p_value_eqtl=-log10(p_value_eqtl)) %$%
  p_value_eqtl
y_max_eqtl <- egene_pv %>% max()
y_min_eqtl <- egene_pv %>% min()
# eQTL SNP but non ASoC SNP
eqtl_snp_track.0 <- wbl_eqtl_tab %>%
  dplyr::filter(contig==snp_contig 
                & between(position, range_start, range_end) 
                & position!=snp_position
                & external_gene_name==eqtl_gene) %>%
  dplyr::mutate(p_value_eqtl=-log10(p_value_eqtl)) %>%
  makeGRangesFromDataFrame(seqnames.field="contig", start.field="position",
                           end.field="position", keep.extra.columns=T) %>%
  DataTrack(name=str_glue("Whole blood eQTL SNP\n{eqtl_gene}"),
            background.panel="gray95", cex.axis=0.5, cex=0.6, col="black",
            ylim=c(y_min_eqtl, y_max_eqtl))
# eQTL SNP and ASoC SNP
eqtl_snp_track.1 <- wbl_eqtl_tab %>%
  dplyr::filter(contig==snp_contig
                & position==snp_position
                & external_gene_name==eqtl_gene) %>%
  # If the target SNP isn't a eQTL SNP, we use a fake p-value = 1.
  (function(x) { if(nrow(x)>0) x
                 else data.frame(contig=snp_contig, position=snp_position, p_value_eqtl=1)}) %>%
  dplyr::mutate(p_value_eqtl=-log10(p_value_eqtl)) %>%
  makeGRangesFromDataFrame(seqnames.field="contig", start.field="position",
                           end.field="position", keep.extra.columns=T) %>%
  DataTrack(cex.axis=0.5, cex=2, col="red", ylim=c(y_min_eqtl, y_max_eqtl))
eqtl_snp_track <- OverlayTrack(list(eqtl_snp_track.0, eqtl_snp_track.1))


# Covid19 GWAS SNP track
gwas_snp_file <- str_glue("{wkdir}/{version}/integration_gwas_snp_track.csv")
if (file.exists(gwas_snp_file)) {
  gwas_snp_tab <- fread(gwas_snp_file)
} else {
  # Very severe respiratory confirmed covid vs population
  # NOTE: 1. SNPs with MAF < 0.1 were filtered out.
  gwas_snp_tab <- fread("../../inputs/covid19_gwas/covid19hg.org/COVID19_HGI_A2_ALL_leave_23andme_20210607.txt.gz") %>%
    dplyr::rename("contig"="#CHR", "position"="POS", "p_value_gwas"="all_inv_var_meta_p") %>%
    dplyr::filter(all_meta_AF>=0.1, p_value_gwas<0.05) %>%
    dplyr::select(contig, position, p_value_gwas) %>%
    (function(x) {fwrite(x, gwas_snp_file, sep=",", quote=F, row.names=F); x})
}
y_max_gwas <- gwas_snp_tab %>%
  dplyr::filter(contig==snp_contig&between(position, range_start, range_end)) %>%
  dplyr::mutate(p_value_gwas=-log10(p_value_gwas)) %$%
  p_value_gwas %>%
  max()
gwas_snp_track.0 <- gwas_snp_tab %>%
  dplyr::filter(contig==snp_contig&between(position, range_start, range_end)) %>%
  dplyr::mutate(p_value_gwas=-log10(p_value_gwas)) %>%
  makeGRangesFromDataFrame(seqnames.field="contig", start.field="position",
                           end.field="position", keep.extra.columns=T) %>%
  DataTrack(name="COVID19 GWAS SNP", background.panel="white", cex.axis=0.5,
            cex=0.6, col="black", ylim=c(0, y_max_gwas))
gwas_snp_track.1 <- gwas_snp_tab %>%
  dplyr::filter(contig==snp_contig & position==snp_position) %>%
  dplyr::mutate(p_value_gwas=-log10(p_value_gwas)) %>%
  makeGRangesFromDataFrame(seqnames.field="contig", start.field="position",
                           end.field="position", keep.extra.columns=T) %>%
  DataTrack(background.panel="white", cex.axis=0.5, cex=2, col="red",
            ylim=c(0, y_max_gwas))
gwas_snp_track <- OverlayTrack(list(gwas_snp_track.0, gwas_snp_track.1))


# ASoC SNP track
asoc_snp_file <- str_glue("{wkdir}/{version}/integration_asoc_snp_track.csv")
if (file.exists(asoc_snp_file)) {
  asoc_snp_tab <- fread(asoc_snp_file)
} else {
  # NOTE: 1. The saved table were splited by semi-colon (;)
  asoc_snp_tab <- fread("../../outputs/scATAC-seq/summary/ver_beta/pseudo-bulk_per-celltype_asoc-snps_annotated.csv") %>%
    dplyr::filter(p_value_adj<0.05) %>%
    dplyr::select(contig, position, variantID, celltype, condition) %>%
    dplyr::mutate(feature=paste(variantID, celltype, condition, sep=", ")) %>%
    dplyr::select(contig, position, feature) %>%
    (function(x) {fwrite(x, asoc_snp_file, sep=";", quote=F, row.names=F); x})
}
asoc_snp_track <- asoc_snp_tab %>%
  dplyr::filter(contig==snp_contig
                & between(position, range_start, range_end)
                & grepl(rsid, feature)) %>%
  makeGRangesFromDataFrame(seqnames.field="contig", start.field="position",
                           end.field="position", keep.extra.columns=T) %>%
  AnnotationTrack(name="ASoC\nSNP", showId=T, feature=.$feature, fill="red",
                  col="gray95", group=.$feature, just.group="right", lwd=6,
                  background.panel="gray95", fontsize.group=25, min.width=8)


# Ensembl gene track
gene_pos_track <- BiomartGeneRegionTrack(genome="hg38", chromosome=snp_contig,
                                         start=range_start, end=range_end,
                                         filters=list(transcript_is_canonical=T, biotype="protein_coding"),
                                         name="Ensembl gene",
                                         transcriptAnnotation="symbol",
                                         background.panel="white", fontsize.group=20)


# Indicate the location of ASoC SNPs
data_track <- HighlightTrack(trackList=list(atac_dep_track, g2p_link_track, hic_link_track,
                                            eqtl_snp_track, gwas_snp_track),
                             start=snp_position, width=1, col="blue",
                             chromosome=snp_contig, inBackground=F)


# The plot
pdf(str_glue("{wkdir}/{version}/integration_xomics_{rsid}_{eqtl_gene}.pdf"), width=16, height=12, onefile=F)
plotTracks(list(ideog_track, gaxis_track, data_track, asoc_snp_track, gene_pos_track),
           sizes=c(1, 1.5, 3, 4, 4, 4.5, 4.5, 2, 2),
           from=range_start, to=range_end,
           chromosome=as.character(snp_contig),
           fontsize=16, background.title="lightgray",
           col.title="black", col.axis="black")
dev.off()

