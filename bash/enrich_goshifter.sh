# Annotation enrichment
# CIS-BP database version: 8-1-2019 Build version 2.00
proj_dir=~/Documents/projects/wp_bcg_eqtl

wkdir=$proj_dir/outputs/pseudo_bulk/enrichment
gffile=$proj_dir/inputs/annotations/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.gff

LDMATDIR=$HOME/Documents/projects/wp_bcg_eqtl/inputs/annotations/goshifter_ld
GOSHIFTERSIMG=$HOME/Documents/projects/wp_bcg_eqtl/scripts/singularity/goshifter.simg

if [[ ! -d $wkdir ]]; then mkdir -p $wkdir; fi

# Estimate enrichment by GoShifter
# GoShifter BASEDIR ANNOTYPE MODEL CELLTYPE [CONDITION]
GoShifter() {
  if [[ $@ -lt 4 ]]; then echo Not enough parameters! && return -1; fi

  local base_dir=$1
  local annotype=$2
  local model=$3
  local celltype=$4

  local in_dir=$base_dir/$model/$celltype
  local out_dir=$base_dir/enrichment/snp/$model/$celltype
  local ann_file=$base_dir/enrichment/annotation/$annotype.bed.gz

  if [[ -n $5 ]]; then
    local condition=$5
    local in_dir=$base_dir/$model/$celltype/$condition
    local out_dir=$base_dir/enrichment/snp/$model/$celltype/$condition
  fi

  if [[ ! -d $out_dir ]]; then mkdir -p $out_dir; fi

  awk 'NR==1{print "SNP\tChrom\tBP"; next} {print $1"\tchr"$16"\t"$17}' \
    $in_dir/top_qtl_results_all_FDR0.05.txt >| $out_dir/snp_pos.txt

  singularity exec --env APPEND_PATH=$HOME/tools/bin $GOSHIFTERSIMG goshifter.py \
    --ld $LDMATDIR \
    --permute 2000 \
    --rsquared 0.6 \
    --annotation $ann_file \
    --out $out_dir/$annotype \
    --snpmap $out_dir/snp_pos.txt \
    2>&1 >| $out_dir/$annotype.log
}


#
## Main steps
#

# No Significant results were found in the genomic annotations from ATAC-seq and ensembl

# Prepare genomic annotations
if [[ ! -d $wkdir/annotation ]]; then mkdir -p $wkdir/annotation; fi

for rtype in enhancer promoter{,_flanking_region} {CTCF,TF}_binding_site open_chromatin_region; do
  ## Prepare genomic feature file containing genomic annotations.
  awk -v rtype=$rtype '/^[0-9]{1,2}/ && $3==rtype {print "chr"$1"\t"$4"\t"$5"\t"$9}' $gffile \
    | sort -k1,1 -k2,2g -k3,3g | bgzip >| $wkdir/annotation/$rtype.ensembl.bed.gz
  tabix $wkdir/annotation/$rtype.ensembl.bed.gz

  ## Prepare a BED file containing ATAC-seq peaks
  zcat $proj_dir/inputs/atac_seq/peaks_filtered.csv.gz \
    | awk -F , -v rtype=$rtype '$13==rtype {print $2"\t"$3"\t"$4"\treg_type="$13";reg_feature="$11}' \
    | sort -k1,1 -k2,2g -k3,3g | bgzip >| $wkdir/annotation/$rtype.300BCG_ATACseq.bed.gz
  tabix $wkdir/annotation/$rtype.300BCG_ATACseq.bed.gz

  ## Prepare a BED file containing ATAC-seq peaks in monocytes
  zcat $proj_dir/inputs/atac_seq/peaks_filtered_monocyte.csv.gz \
    | awk -F , -v rtype=$rtype '$13==rtype {print $2"\t"$3"\t"$4"\treg_type="$13";reg_feature="$11}' \
    | sort -k1,1 -k2,2g -k3,3g | bgzip >| $wkdir/annotation/$rtype.300BCG_ATACseq.monocytes.bed.gz
  tabix $wkdir/annotation/$rtype.300BCG_ATACseq.monocytes.bed.gz

  zcat $proj_dir/inputs/atac_seq/peaks_filtered_cd8t.csv.gz \
    | awk -F , -v rtype=$rtype '$13==rtype {print $2"\t"$3"\t"$4"\treg_type="$13";reg_feature="$11}' \
    | sort -k1,1 -k2,2g -k3,3g | bgzip >| $wkdir/annotation/$rtype.300BCG_ATACseq.cd8t.bed.gz
  tabix $wkdir/annotation/$rtype.300BCG_ATACseq.cd8t.bed.gz

  zcat $proj_dir/inputs/atac_seq/peaks_filtered_nkcell.csv.gz \
    | awk -F , -v rtype=$rtype '$13==rtype {print $2"\t"$3"\t"$4"\treg_type="$13";reg_feature="$11}' \
    | sort -k1,1 -k2,2g -k3,3g | bgzip >| $wkdir/annotation/$rtype.300BCG_ATACseq.nkcell.bed.gz
  tabix $wkdir/annotation/$rtype.300BCG_ATACseq.nkcell.bed.gz
done


# for annotype in $wkdir/annotation/*.{monocytes,nkcell,cd8t}.bed.gz; do
for annotype in $wkdir/annotation/*.nkcell.bed.gz; do
  annotype=$(basename ${annotype/.bed.gz/})
  for mode in normal interaction; do
    for celltype in Monocytes CD{4,8}T B NK; do
      echo $annotype-$mode-$celltype
      if [[ $mode == "normal" ]]; then
        GoShifter $proj_dir/outputs/pseudo_bulk $annotype $mode $celltype
      else
        for condition in {T0_LPS,T3m_LPS,T3m_RPMI}.vs.T0_RPMI T3m_LPS.vs.T3m_RPMI; do
          echo GoShifter $proj_dir/outputs/pseudo_bulk $annotype $mode $celltype $condition
        done
      fi
    done
  done
done


mode="normal"
for annotype in $wkdir/annotation/*.txt; do
  annotype=$(basename ${annotype/.txt/})
  for celltype in Monocytes CD{4,8}T NK B; do
    GoShifter $proj_dir/outputs/pseudo_bulk/enrichment $annotype $mode $celltype
  done
done
