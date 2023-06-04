#!/usr/bin/env bash
# r/coloc.r
# r/harmonization.r
# r/plot_main_figures.r
# r/trackplot.r

proj_dir=~/Documents/projects/wp_bcg_eqtl
source $proj_dir/scripts/.env/bin/activate

#
## Split pseudo bulk expression into different conditions.
#
workdir=$proj_dir/inputs/pseudo_bulk
for per_cell in Monocytes CD4T CD8T NK B; do
  per_inf=$workdir/$per_cell

  for per_cond in T0_RPMI T0_LPS T3m_RPMI T3m_LPS; do
    mkdir -p $per_inf/$per_cond

    # Phenotypes
    awk -v PERCOND=$per_cond -f- <<'EOF' $per_inf/phenotypes.tsv >| $per_inf/$per_cond/phenotypes.tsv
NR == 1 {
  chosen_cols[1] = 1
  for (ii = 2; ii <= NF; ii++) {
    chosen_cols[ii] = $ii ~ PERCOND ? 1 : 0
  }
} {
  for (ii = 1; ii < NF; ii++) if (chosen_cols[ii]) printf $ii"\t"
  if (chosen_cols[NF]) { print $ii } else { print "" }
}
EOF
    # Covariants
    grep -e sample_id -e $per_cond $per_inf/covariates.tsv | cut -f1,4,5 >| $per_inf/$per_cond/covariates.tsv
    grep -e sample_id -e $per_cond $per_inf/covariates_wpeer.tsv | cut -f1,4,5 >| $per_inf/$per_cond/covariates_wpeer.tsv

    # Sample mapping
    cut -f 1 $per_inf/$per_cond/covariates.tsv \
      | grep -v sample_id \
      | xargs -I % -n 1 grep % $per_inf/sample_mapping.tsv \
      >| $per_inf/$per_cond/sample_mapping.tsv
  done
done

# Fetch allele frequency
vcffile=$proj_dir/inputs/genotypes/300BCG_sub40_imp_hg38_ids_clean.vcf.gz
bcftools query -Hf '%CHROM,%ID,%POS,%REF,%ALT,%INFO/AF,[,%DS]\n' $vcffile \
  | sed 's/\[[0-9]*\]//g; s/:DS//g; s/^# \+//g' \
  | awk -F, '/^#/ {print; next;} {print}' \
  | head -1000 \
  >| $proj_dir/temps/genotype_test.csv


# Check window of current runs
windowSize = '100000'

# ADCY3 chr2:24819169-24920237, minus strand, GRCh38, 187 SNPs
# ADCY3 gene body
bcftools view 300BCG_sub40_imp_hg38_ids_clean.vcf.gz chr2:24819169-24920237 | grep -cv '^#'

# ADCY3 upstream 
# 100,000, 37 SNPs
bcftools view 300BCG_sub40_imp_hg38_ids_clean.vcf.gz chr2:24920237-25020237 | grep -cv '^#'
# 1,000,000, 600 SNPs
bcftools view 300BCG_sub40_imp_hg38_ids_clean.vcf.gz chr2:24920237-25920237 | grep -cv '^#'

# ADCY3 downstream
# 100,000, 41 SNPs
bcftools view 300BCG_sub40_imp_hg38_ids_clean.vcf.gz chr2:24719169-24819169 | grep -cv '^#'
# 100,000, 504 SNPs
bcftools view 300BCG_sub40_imp_hg38_ids_clean.vcf.gz chr2:23819169-24819169 | grep -cv '^#'

# Estimate the LD
bcftools view 300BCG_sub40_imp_hg38_ids_clean.vcf.gz chr2:23819169-24819169 | grep -v '^#' | head -1 | cut -f1-5
bcftools view 300BCG_sub40_imp_hg38_ids_clean.vcf.gz chr2:23819169-24819169 | grep -v '^#' | tail -1 | cut -f1-5

# Plot boxplots
# rs11687089:ADCY3, # Good common eQTL across different conditions.
# rs11731570:EXOSC9, # Good ieQTL example, but less information about the gene
# rs11080327:SLFN5, # Also good ieQTL, but not cell-type specific.
# rs11733586:TIFA, # Good example, regarding LPS response
# rs4761234:LYZ, # Good example. it is the monocyte marker
#
# Others
# rs1475642:BORCS7 rs7562347:COMMD1 rs10794331:DEAF1 rs2161548:ERAP2 rs2221936:EXOSC9 rs1367313:GFM1
# rs115613985:GZMK rs344352:HAGH rs2075846:HOMEZ rs339078:IMPDH1 rs2310003:IRF2 rs67260737:NUMA1
# rs4758197:PPFIBP2 rs10946198:RNASET2 rs1662745:RUNDC1 rs9292067:SNX18 rs1011731:SUCO
# rs62057151:ARL17B rs8138678:SMDT1 rs9410320:SEMA4D rs10792269:MS4A7 rs2298746:XRRA1 rs6938061:MDGA1
# rs10735234:GSTM3 rs10843881:DDX11 rs11150882:C17orf97 rs11958387:CENPK rs12972593:ZNF100
# rs13148861:FAM184B rs199499:ARL17B rs2290911:SH3YL1 rs2559854:CHPT1 rs2942377:MORN3 rs6002625:SMDT1
# rs62063281:KANSL1 rs62229264:CBR1 rs2288004:CORO1A rs2135923:CD55 rs6815614:DDX60L rs2906999:DTX2
# rs2066932:DIP2A rs312932:FZR1 rs4846032:AGTRAP rs9303316:NT5C3B rs35188965:SLC12A7 rs2036764:HOXB2
# rs7115854:ARL14EP rs7543002:CD55 rs10859822:NDUFA12
#

in_dir=$proj_dir/inputs/pseudo_bulk/all
python3 $proj_dir/scripts/py3/fetch_snp_info.py \
  -x rs11731570:EXOSC9 rs11080327:SLFN5 rs11733586:TIFA rs4761234:LYZ rs7543002:CD55 rs11687089:ADCY3 rs2564978:CD55 rs17653193:KANSL1 rs12949100:ORMDL3 \
  -g $proj_dir/inputs/genotypes/300BCG_sub40_imp_hg38_ids_clean.vcf.gz \
  -p $in_dir/phenotypes.tsv \
  -c $in_dir/covariates_celltype.tsv \
  -m $in_dir/sample_mapping.tsv \
  -o $proj_dir/outputs/pseudo_bulk/example_eQTL


#
## Estimate top eQTL's enrichment against genomic features annotations
#
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

## Step 2.1. Prepare 300BCG ATAC-seq and ensembl annotations
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

## Step 2.2 Prepare Genomic annotations of epigenomic markers by Roadmap Epigenome project.
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

# Step 2.3. Prepare GWAS summary statistics
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
min_bs=0.1
max_bs=1000
for min_bs in 0.1 0.25 0.5 1; do
  for max_bs in 250 500 750 1000; do
    for cell_type in Monocytes CD{4,8}T B NK; do
      out_dir=$wkdir/${min_bs}-${max_bs}Kbp
      if [[ ! -d $out_dir ]]; then mkdir -p $out_dir; fi

      Rscript $proj_dir/scripts/r/enrich.r \
        --qtl-db $proj_dir/outputs/pseudo_bulk/normal/$cell_type \
        --blocks $proj_dir/inputs/reference/haplotype_blocks/haplotype_blocks_EUR.blocks.det \
        --ann-db $wkdir/annotation \
        --min-block-size $min_bs \
        --max-block-size $max_bs \
        --out $out_dir/$cell_type.csv
    done
  done
done
cd $wkdir



#
## Harmonization between eQTL and GWAS sumstat
#
awk -F, -f- <<'EOF' harmonized_data.csv run_2_harmonized_data.csv >| temp.csv
  NR == FNR {print}
  NR != FNR && FNR == 1 {next}
  NR != FNR && FNR != 1 {
    printf $29","$32","$21","$13",";
    for (ii=1; ii<34; ii++) if (ii!=13&&ii!=21&&ii!=29&&ii!=32) printf $ii",";
    print $34;
  }
EOF

mv harmonized_data.csv run_1_harmonized_data.csv && gzip run_{1,2}_harmonized_data.csv && mv temp.csv harmonized_data.csv

awk -F, -f- <<'EOF' harmonized_data.csv
  function file_exists(file) { return (getline _ < file) >= 0 ? 1 : 0; }

  NR == 1 {HEADER = $0; next}
  {
    PER_GWAS_FN = $3"_harmonized_data.csv";
    if (file_exists(PER_GWAS_FN)) {
      print $0 >> PER_GWAS_FN
    } else {
      print HEADER > PER_GWAS_FN
    }
  }
EOF


#
## Plot ASoC of rs2564978
#
asoc_snps_for_rcplot=$proj_dir/outputs/COVID_MHH50/heter-snps_metainfo.two_conditions.csv
bamdir=/vol/projects/BIIM/Covid_50MHH/ASoC/outputs/scATAC-seq/run_v1/optdir/bam_pct_asoc/unmerged
snpid=rs2564978
snpid=rs11080327

python $proj_dir/scripts/py3/plot_asoc_rc.py \
  -r $asoc_snps_for_rcplot \
  -i $snpid \
  -b $bamdir \
  -f 40 \
  --fig-width 7 \
  -d Plasmablast pDC ncMono B cMono NK \
  -o $proj_dir/outputs/pseudo_bulk/example_eQTL

#
##
#
awk 'NR==1 {print; next} $4~"KANSL1/" {print}' $proj_dir/inputs/annotations/GeneHancer/GeneHancer.KANSL1.chr17_36438746_55814845.grch38.csv \
  >| KANSL1.GeneHancer.csv

while read -r line; do
  enhancer_start=$(cut -f10 <<<$line)
  enhancer_end=$(cut -f11 <<<$line)
  awk -F, -vPOS_START=$enhancer_start -vPOS_END=$enhancer_end \
    'NR==1 && POS_START == "geneHancerStart" {print; next} POS_START <= $19 && $19 <= POS_END {print}' \
    KANSL1.harmonized.csv
done < KANSL1.GeneHancer.csv >| KANSL1.enhancer_SNPs.csv

tar_snp_file=$proj_dir/outputs/pseudo_bulk/example_eQTL/tar_snps.txt
bcftools query -H -i "ID=@$tar_snp_file" -f '%CHROM,%POS,%ID,%REF,%ALT,%MAF,[,%GT]\n' $vcffile \
  >| $proj_dir/outputs/pseudo_bulk/example_eQTL/tar_snps.genotypes.csv

awk -f- <<'EOF' $proj_dir/outputs/pseudo_bulk/summary_statistic/interaction/*/*/qtl_results_all.txt \
  >| $proj_dir/outputs/pseudo_bulk/example_eQTL/CD55-rs2564978-interaction.csv
NR==1 {print "celltype\tcondition\t"$0; next}
$1~/CD55/ && $2~/rs2564978/ {n = split(FILENAME, fnarray, "/"); print fnarray[n-2]"\t"fnarray[n-1]"\t"$0}
EOF

awk -f- <<'EOF' $proj_dir/outputs/pseudo_bulk/summary_statistic/normal/*/qtl_results_all_FDR0.05.txt \
  >| $proj_dir/outputs/pseudo_bulk/example_eQTL/CD55-rs2564978-main.csv
NR==1 {print "celltype\t"$0; next}
$1~/CD55/ && $2~/rs2564978/ {n = split(FILENAME, fnarray, "/"); print fnarray[n-1]"\t"$0}
EOF

awk -f- <<'EOF' $proj_dir/outputs/pseudo_bulk/summary_statistic/interaction/*/*/qtl_results_all.txt \
  >| $proj_dir/outputs/pseudo_bulk/example_eQTL/interaction_eqtl.5e-7.csv
NR==1 {print "celltype\tcondition\t"$0; next}
$3<5e-7 {n = split(FILENAME, fnarray, "/"); print fnarray[n-2]"\t"fnarray[n-1]"\t"$0}
EOF

#
## add head
#
header_line=exposure,id.exposure,outcome,id.outcome,SNP,effect_allele.exposure,other_allele.exposure,effect_allele.outcome,other_allele.outcome,beta.exposure,beta.outcome,eaf.exposure,eaf.outcome,remove,palindromic,ambiguous,chr.outcome,pos.outcome,se.outcome,pval.outcome,samplesize.outcome,ncase.outcome,ncontrol.outcome,mr_keep.outcome,pval_origin.outcome,se.exposure,chr.exposure,pos.exposure,samplesize.exposure,pval.exposure,mr_keep.exposure,pval_origin.exposure,action,mr_keep
for x in *.gz; do
  zcat $x | awk -v HEADER_LINE=$header_line -F, 'NR==1 && $1~/exposure/ {print; next} NR==1&&$1!~/exposure/ {print HEADER_LINE; print; next} {print}' \
    | gzip >| tmp.gz \
    && mv -f tmp.gz $x;
done

# Update the sample size for HDL/LDL/Triglycerides
if [[ true ]]; then
  declare -A sample_size
  sample_size["HDLCholesterol"]=403943
  sample_size["LDLCholesterol"]=440546
  sample_size["Triglycerides"]=441016
  for per_key in "${!sample_size[@]}"; do
    per_ss=${sample_size[$per_key]}
    for per_hmfile in $(find . -name ${per_key}_harmonized_data.csv.gz); do
      zcat $per_hmfile \
        | awk -v SS=$per_ss -F, 'NR==1{print; next} {for(i=1;i<21;i++) printf $i","; printf SS","; for(i=22;i<34;i++) printf $i","; print $34}' \
        >| ${per_hmfile/.gz/}
      mv ${per_hmfile} $per_hmfile.bk
      gzip -f ${per_hmfile/.gz/}
    done
  done
fi
