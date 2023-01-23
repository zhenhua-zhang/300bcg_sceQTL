#!/usr/bin/env bash
# SMR, summary-based mendelian randomization.

workdir=~/Documents/projects/wp_bcg_eqtl

celltype=Monocytes
model=normal

temp_dir=$workdir/outputs/pseudo_bulk/smr/tempdir
out_dir=$workdir/outputs/pseudo_bulk/smr/database

# Required files
feature_ids=$temp_dir/features_ids.txt
snp_db=$temp_dir/snp_information.txt
gffdb=$workdir/inputs/reference/Gencode/gencode.v41.basic.annotation.Ec.transcript.autosome.gff


#
## Prepare a gene list with strand information
#
awk -f- <<'EOF' $gffdb | grep -v ENSG >| $workdir/inputs/reference/SMR/glist-hg38.sorted.strand
{
  split($9, ATTR, ";");
  split(ATTR[6], TMPLIST, "[=.]"); GENEID = TMPLIST[2];
  PROBEID_DICT[GENEID] = $1"\t"GENEID"\t0\t"$5"\t"GENEID"\t"$7;
  print $1" "$4" "$5" "GENEID" "$7
}
EOF

#
## Prepare GWAS summary
#
echo -e "[I]: Preapre GWAS summary ..."
for gwas_sum_raw in $workdir/inputs/public_gwas/*.gz; do
  # gwas_sum_raw=$workdir/inputs/public_gwas/2013_ObesityClass1_ieu-a-90.vcf.gz
  gwas_sum_ma=$(basename $gwas_sum_raw | sed "s/vcf.gz/ma/g")
  if [[ -e $out_dir/public_gwas/$gwas_sum_ma ]]; then continue; fi

  echo -e "SNP\tA1\tA2\tfreq\tb\tse\tp\tN" >| $out_dir/public_gwas/$gwas_sum_ma
  bcftools query -f '%ID\t%ALT\t%REF\t%AF\t[%ES\t%SE\t%LP\t%SS]\n' $gwas_sum_raw \
    | awk '{OFS="\t"; if($7==0) $7=1; else if($7>=999) $7=0; else $7=10^-$7; print}' \
    >> $out_dir/public_gwas/$gwas_sum_ma 
done


#
## Prepare eQTL summary
#
# SNP information
# bcftools query -f "%ID\t%REF\t%ALT\t%INFO/AF\n" ../../../inputs/genotypes/300BCG_sub40_imp_hg38_ids_clean.vcf.gz >| snp_information.txt
for celltype in B CD4T CD8T Monocytes NK; do
  # 1.1 Format the eQTL summary into .esd format.
  echo -e "[I]: Format eQTL VCF into ESD format ..."
  eqtl_sum_raw=$workdir/outputs/pseudo_bulk/summary_statistic/$model/$celltype/qtl_results_all.txt
  awk -v TEMP_DIR=$temp_dir -v CELLTYPE=$celltype -f- <<'EOF' $snp_db $eqtl_sum_raw
NR == FNR { REF_DICT[$1] = $2; ALT_DICT[$1] = $3; FREQ_DICT[$1] = $4; next; }
FNR > 1 {
  OUT_FILE = TEMP_DIR"/"CELLTYPE"/"$1".esd"
  if (getline < OUT_FILE < 0) { print "Chr\tSNP\tBp\tA1\tA2\tFreq\tBeta\tse\tp" > OUT_FILE; }
  CHR = $17; ID = $2; POS = $18; BETA = $4; SE = $5; PVAL = $3;
  A1 = $19;
  A2 = ALT_DICT[ID] == A1 ? REF_DICT[ID] : ALT_DICT[ID];
  FREQ = ALT_DICT[ID] == A1 ? FREQ_DICT[ID] : 1 - FREQ_DICT[ID];
  print CHR"\t"ID"\t"POS"\t"A1"\t"A2"\t"FREQ"\t"BETA"\t"SE"\t"PVAL >> OUT_FILE;
  LAST_FEATURE = $1;
}
EOF

  # 1.2 Prepare a file including probe information and file paths of the eQTL summary data, check 1.1
  # Columns are chromosome, probe ID(can be the ID of an exon or a transcript for RNA-seq data),
  # genetic distance (can be any arbitary value), physical position, gene ID, gene orientation,
  # and PathOfEsd
  echo -e "[I]: Prepare a file including probe information and file path of the eQTL summary"
  esd_file_list=$temp_dir/$celltype/esd_file_list.txt
  awk -v TEMP_DIR=$temp_dir -v CELLTYPE=$celltype -f- <<'EOF' $gffdb $feature_ids | sort -k1,1g -k4,4g >| $esd_file_list
BEGIN { print "Chr\tProbeID\tGeneticDistance\tProbeBp\tGene\tOrientation\tPathOfEsd"; }
NR == FNR {
  split($9, ATTR, ";");
  split(ATTR[6], TMPLIST, "[=.]"); GENEID = TMPLIST[2];
  PROBEID_DICT[GENEID] = $1"\t"GENEID"\t0\t"$5"\t"GENEID"\t"$7;
}
NR != FNR {
  if (PROBEID_DICT[$1]) print PROBEID_DICT[$1]"\t"TEMP_DIR"/"CELLTYPE"/"$1".esd";
  else next;
}
EOF

  # 1.3 Generate BESD dataset.
  echo -e "[I]: Generate BESD database from eQTL summary ..."
  smr --make-besd \
    --eqtl-flist $temp_dir/$celltype/esd_file_list.txt \
    --out $out_dir/sceqtl/$model/$celltype
done

#
## SMR analysis
#
bfile=$workdir/inputs/reference/genotypes/GRCh38//EUR
# gwas_sum_ma=$workdir/outputs/pseudo_bulk/smr/database/public_gwas/2013_ObesityClass1_ieu-a-90.ma
# eqtl_summary=$workdir/outputs/pseudo_bulk/smr/database/sceqtl/normal/$celltype

for per_eqtl in B CD4T CD8T Monocytes NK; do
  for gwas_summary in $out_dir/public_gwas/*.ma; do
    per_gwas=$(basename $gwas_summary | sed 's/.ma$//g')
    eqtl_summary=$workdir/outputs/pseudo_bulk/smr/database/sceqtl/normal/$per_eqtl
    # 2.1 SMR, single SNP, cis-eQTL
    echo -e "[I]: Perform SMR (single-SNP mode) for cis-eQTL ..."
    smr_outfile=$workdir/outputs/pseudo_bulk/smr/results/SMR.cis_single.$per_gwas.$per_eqtl
    smr_logfile=$workdir/outputs/pseudo_bulk/smr/results/SMR.cis_single.$per_gwas.$per_eqtl.log
    [[ ! -e $smr_outfile.smr ]] && smr --smr \
      --bfile $bfile \
      --gwas-summary $gwas_summary \
      --beqtl-summary $eqtl_summary \
      --thread-num 4 \
      --peqtl-smr 5e-6 \
      --out $smr_outfile \
      2>&1 1>| $smr_logfile

    # 2.2 SMR, multi SNP, cis-eQTL
    echo -e "[I]: Perform SMR (multi-SNP mode) for cis-eQTL ..."
    smr_outfile=$workdir/outputs/pseudo_bulk/smr/results/SMR.cis_multi.$per_gwas.$per_eqtl
    smr_logfile=$workdir/outputs/pseudo_bulk/smr/results/SMR.cis_multi.$per_gwas.$per_eqtl.log
    [[ ! -e $smr_outfile.msmr ]] && smr --smr-multi \
      --bfile $bfile \
      --gwas-summary $gwas_summary \
      --beqtl-summary $eqtl_summary \
      --thread-num 4 \
      --peqtl-smr 5e-6 \
      --out $smr_outfile \
      2>&1 1>| $smr_logfile

    # 2.3 SMR, single SNP, trans-eQTL
    echo -e "[I]: Perform SMR (single-SNP mode) for trans-eQTL ..."
    trans_win=1000
    smr_outfile=$workdir/outputs/pseudo_bulk/smr/results/SMR.trans_single.win_$trans_win.$per_gwas.$per_eqtl
    smr_logfile=$workdir/outputs/pseudo_bulk/smr/results/SMR.trans_single.win_$trans_win.$per_gwas.$per_eqtl.log
    [[ ! -e $smr_outfile.smr ]] && smr --trans \
      --trans-wind $trans_win \
      --bfile $bfile \
      --gwas-summary $gwas_summary \
      --beqtl-summary $eqtl_summary \
      --thread-num 4 \
      --peqtl-smr 5e-6 \
      --out $smr_outfile \
      2>&1 1>| $smr_logfile
  done
done



#
## Plot SMR results
#
gene_list=$workdir/inputs/reference/SMR/glist-hg38.sorted.strand
cmp_pool=(
  # ADCY3:2018_BodyMassIndex_ieu-b-40:Monocytes
  # ADCY3:2013_ObesityClass1_ieu-a-90:Monocytes
  # ADCY3:2013_ObesityClass2_ieu-a-91:Monocytes
  # RNASET2:2015_InflammatoryBowelDisease_ieu-a-31:CD8T
  # RNASET2:2015_UlcerativeColitis_ieu-a-32:CD8T
  # RNASET2:2015_CrohnsDisease_ieu-a-30:CD8T
  # RNASET2:2015_InflammatoryBowelDisease_ieu-a-31:CD4T
  # RNASET2:2015_UlcerativeColitis_ieu-a-32:CD4T
  # RNASET2:2015_CrohnsDisease_ieu-a-30:CD4T
  EXOSC9:2014_Height_ieu-a-89:Monocytes
)

for per_cmp in ${cmp_pool[@]}; do
  per_probe=$(cut -f1 -d: <<<$per_cmp)
  per_gwas=$(cut -f2 -d: <<<$per_cmp)
  per_eqtl=$(cut -f3 -d: <<<$per_cmp)

  smr --plot \
    --bfile $bfile \
    --gwas-summary $out_dir/public_gwas/$per_gwas.ma \
    --beqtl-summary $out_dir/sceqtl/normal/$per_eqtl \
    --probe $per_probe \
    --probe-wind 800 \
    --gene-list $gene_list \
    --peqtl-smr 5e-6 \
    --out $per_gwas.$per_eqtl

  Rscript - <<EOF
source("~/Downloads/smr_plot/plot/plot_SMR.r")
SMRData <- ReadSMRData("plot/$per_gwas.$per_eqtl.$per_probe.txt")
pdf("$per_gwas.$per_eqtl.$per_probe.SNP_effect.pdf")
SMREffectPlot(data=SMRData)
dev.off()

pdf("$per_gwas.$per_eqtl.$per_probe.locus_plot.pdf")
SMRLocusPlot(data=SMRData, smr_thresh=5e-2, heidi_thresh=0.05, plotWindow=800, max_anno_probe=16)
dev.off()
EOF
done
