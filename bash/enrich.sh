#!/bin/bash

#
## Estimate top eQTL's enrichment against genomic features annotations
#

proj_dir=~/Documents/projects/wp_bcg_eqtl

mode=normal
cell_type=Monocytes
condition=

in_dir=$proj_dir/outputs/pseudo_bulk/$mode/$cell_type
out_dir=$proj_dir/outputs/pseudo_bulk/outcomes/annotation_enrichment/$mode/$cell_type/$condition

ref_dir=$proj_dir/inputs/reference
refdb=$ref_dir/genotypes/GRCh38/EUR

# Step 0. Prepare a referecen panel.
# Check more in ./prepare_reference.sh

# Step 1. Haplotype blocks identified using variants identified in 1000 Genomes Project (phase 3, GRCh38)
hb_file=$proj_dir/inputs/reference/haplotype_blocks/haplotype_blocks_EUR
if [[ ! -e $hb_file.blocks ]]; then plink --bfile $refdb --threads 2 --blocks no-pheno-req --out $hb_file; fi


# Step 2. Prepare sets of annotations
wkdir=$proj_dir/outputs/pseudo_bulk/outcomes/annotation_enrichment
if [[ ! -d $wkdir/annotation ]]; then mkdir -p $wkdir/annotation; fi

## Step 2.1. 300BCG ATAC-seq and ensembl annotations
skip=true
for rtype in enhancer promoter{,_flanking_region} {CTCF,TF}_binding_site open_chromatin_region; do
  if [[ $skip == true ]]; then continue; fi
  ## Prepare genomic feature file containing genomic annotations.
  zcat $proj_dir/inputs/annotations/Ensembl/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.gff.gz \
    | awk -v rtype=$rtype 'BEGIN {print "Chrom\tStart\tStop\tMetainfo"} /^[0-9]{1,2}/ && $3==rtype {print $1"\t"$4"\t"$5"\t"$9}' $gffile \
    | sort -k1,1g -k2,2g -k3,3g >| $wkdir/annotation/$rtype.NULL.ensembl.txt

  ## Prepare a plain file containing ATAC-seq peaks
  zcat $proj_dir/inputs/atac_seq/peaks_filtered.csv.gz \
    | awk -F, -v rtype=$rtype 'BEGIN {print "Chrom\tStart\tStop\tMetainfo"} $13==rtype {gsub(/chr/, "", $2); print $2"\t"$3"\t"$4"\treg_feature="$11}' \
    | sort -k1,1g -k2,2g -k3,3g >| $wkdir/annotation/$rtype.PBMC.300BCG_ATACseq.txt

  ## Prepare a plain file containing ATAC-seq peaks in monocytes
  ## These ATAC-seq peaks were only from d0, check more information at TODO:
  zcat $proj_dir/inputs/atac_seq/peaks_filtered_monocyte.csv.gz \
    | awk -F, -v rtype=$rtype 'BEGIN {print "Chrom\tStart\tStop\tMetainfo"} $13==rtype {gsub(/chr/, "", $2); print $2"\t"$3"\t"$4"\treg_feature="$11}' \
    | sort -k1,1g -k2,2g -k3,3g >| $wkdir/annotation/$rtype.monocytes.300BCG_ATACseq.txt

  zcat $proj_dir/inputs/atac_seq/peaks_filtered_cd8t.csv.gz \
    | awk -F, -v rtype=$rtype 'BEGIN {print "Chrom\tStart\tStop\tMetainfo"} $13==rtype {gsub(/chr/, "", $2); print $2"\t"$3"\t"$4"\treg_feature="$11}' \
    | sort -k1,1g -k2,2g -k3,3g >| $wkdir/annotation/$rtype.cd8t.300BCG_ATACseq.txt

  zcat $proj_dir/inputs/atac_seq/peaks_filtered_nkcell.csv.gz \
    | awk -F, -v rtype=$rtype 'BEGIN {print "Chrom\tStart\tStop\tMetainfo"} $13==rtype {gsub(/chr/, "", $2); print $2"\t"$3"\t"$4"\treg_feature="$11}' \
    | sort -k1,1g -k2,2g -k3,3g >| $wkdir/annotation/$rtype.nkcell.300BCG_ATACseq.txt
done

## Step 2.2 Genomic annotations of epigenomic markers by Roadmap Epigenome project.
declare -A celltype=(
  # [E030]="Primary_neutrophils_from_peripheral_blood"
  # [E035]="Primary_hematopoietic_stem_cells"
  # [E036]="Primary_hematopoietic_stem_cells_short_term_culture"
  # [E050]="Primary_hematopoietic_stem_cells_G-CSF-mobilized_Female"
  # [E051]="Primary_hematopoietic_stem_cells_G-CSF-mobilized_Male"
  # [E115]="Dnd41_TCell_Leukemia_Cell_Line"
  # [E116]="GM12878_Lymphoblastoid_Cells"
  # [E123]="K562_Leukemia_Cells"
  [E029]="Monocytes Primary_monocytes_from_peripheral_blood"
  [E032]="Bcell Primary_B_cells_from_peripheral_blood"
  [E034]="Tcells Primary_T_cells_from_peripheral_blood"
  [E037]="Thm2 Primary_T_helper_memory_cells_from_peripheral_blood_2"
  [E038]="Thn Primary_T_helper_naive_cells_from_peripheral_blood"
  [E039]="Thn Primary_T_helper_naive_cells_from_peripheral_blood"
  [E040]="Thm1 Primary_T_helper_memory_cells_from_peripheral_blood_1"
  [E043]="Th Primary_T_helper_cells_from_peripheral_blood"
  [E044]="Treg Primary_T_regulatory_cells_from_peripheral_blood"
  [E045]="Tem Primary_T_cells_effector_memory_enriched_from_peripheral_blood"
  [E046]="Nkcell Primary_Natural_Killer_cells_from_peripheral_blood"
  [E047]="CD8Tn Primary_T_CD8+_naive_cells_from_peripheral_blood"
  [E048]="CD8Tm Primary_T_CD8+_memory_cells_from_peripheral_blood"
  [E062]="Mononuclear Primary_mononuclear_cells_from_peripheral_blood"
  [E124]="MonocyteCD14 Monocytes-CD14+_RO01746_Primary_Cells"
)


# 1. Excluded states: 5_Tx5 6_T 7_Tx3 8_TxWk 18_EnhAc 19_DNase 21_Het 25_Quies 24_ReprPC
# 2. Included states (promoter): 2_PromU 3_PromD1 4_PromD2 22_PromP 23_PromBiv
# 2. Included states (enhancer):  10_TxEnh5 11_TxEnh3 12_TxEnhW 13_EnhA1 14_EnhA2 15_EnhAF 16_EnhW1 17_EnhW2
# 3. Included states (trannscription/TSS): 1_TssA 9_TxReg
skip=true
for idx in ${!celltype[@]}; do
  if [[ $skip == true ]]; then continue; fi
  per_ct="${celltype[$idx]}"
  per_ct_ab=$(cut -f1 -d\  <<<"$per_ct")
  per_ct_fn=$(cut -f2 -d\  <<<"$per_ct")

  for rtype in 1_TssA 9_TxReg 10_TxEnh5 11_TxEnh3 12_TxEnhW 13_EnhA1 14_EnhA2 15_EnhAF 16_EnhW1 17_EnhW2 2_PromU 3_PromD1 4_PromD2 22_PromP 23_PromBiv; do
    zcat $proj_dir/inputs/reference/Roadmap/$idx"_25_imputed12marks_hg38lift_dense.bed.gz" \
      | awk -v rtype=$rtype 'BEGIN {print "Chrom\tStart\tStop\tMetainfo"} $4 ~ rtype {gsub(/chr/, "", $1); print $1"\t"$2"\t"$3"\t."}' \
      | sort -k1,1g -k2,2g -k3,3g >| $wkdir/annotation/$rtype.$per_ct_ab.roadmap.txt
  done
done

# Step 2.3. GWAS summary statistics
skip=false
for per_gwas in $proj_dir/inputs/public_gwas/*.vcf.gz; do
  if [[ $skip == true ]]; then continue; fi
  gwas_name=$(basename $per_gwas | cut -f2 -d_)
  echo -e "[I]: Working on $gwas_name"
  zcat $per_gwas \
    | awk 'BEGIN {print "Chrom\tStart\tStop\tMetainfo"} /#/ {next} {split($10, ff, ":"); if (ff[3] <= (-log(5e-8)/log(10))) {next} else {print $1"\t"$2"\t"$2"\t"$3";"ff[3]}}' \
    | sort -k1,1g -k2,2g -k3,3g >| $wkdir/annotation/GWAS.$gwas_name.ieuopengwas.txt
  echo $wkdir/annotation/GWAS.$gwas_name.ieuopengwas.txt
done


# Step 3. Estimate the enrichment of all SNPs in given genomic features.
min_bs=2
max_bs=500
for cell_type in Monocytes CD{4,8}T B NK; do
  Rscript $proj_dir/scripts/r/epianno_enrich.r \
    --qtl-db $proj_dir/outputs/pseudo_bulk/normal/$cell_type \
    --blocks $proj_dir/inputs/reference/haplotype_blocks/haplotype_blocks_EUR.blocks.det \
    --ann-db $wkdir/annotation \
    --min-block-size $min_bs \
    --max-block-size $max_bs \
    --out $wkdir/$cell_type.${min_bs}-${max_bs}Kbp.csv
done
cd $wkdir
