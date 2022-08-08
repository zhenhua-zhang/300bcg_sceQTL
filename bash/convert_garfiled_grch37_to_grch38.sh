#!/bin/bash
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: 2022 Jun 29

proj_dir=~/Documents/projects/wp_bcg_eqtl
out_dir=$proj_dir/outputs/garfield-data-grch38

liftOver_exe=/vol/projects/CIIM/resources/tools/liftOver/liftOver
src_gfdata=/vol/projects/CIIM/resources/tools/garfield-data
chain_file=/vol/projects/CIIM/resources/tools/liftOver/hg19ToHg38.over.chain.gz

mkdir -p $out_dir/{annotation,maftssd,output,tags/r0{1,8},pval,temp}

for chr in {1..22} X; do
  anno_file=$src_gfdata/annotation/chr$chr
  maftss_file=$src_gfdata/maftssd/chr$chr
  tags_r01_file=$src_gfdata/tags/r01/chr$chr
  tags_r08_file=$src_gfdata/tags/r08/chr$chr

  echo -e [I]: ------ Log for chr$chr ------

  # Vairant positions
  echo -e [I]: Format variant annotation into BED ...
  awk -F' ' -v CHROM=chr$chr \
    '{print CHROM"\t"$1-1"\t"$1"\t"CHROM";"$1";"$2}' $anno_file \
    >| $out_dir/annotation/chr$chr.GRCh37.bed

  echo -e [I]: LiftOver variant annotations from GRCh37 to GRCh38 ...
  $liftOver_exe \
    $out_dir/annotation/chr$chr.GRCh37.bed \
    $chain_file \
    $out_dir/annotation/chr$chr.GRCh38.bed \
    $out_dir/annotation/chr$chr.GRCh38.unlifted.bed

  echo -e [I]: Dump variant annotations in GRCh38 ...
  awk -F";" '{split($1, arr, "\t"); print arr[3]" "$3}' \
    $out_dir/annotation/chr$chr.GRCh38.bed \
    | sort -k1,1g -t " " \
    >| $out_dir/annotation/chr$chr

  # Coordination map
  echo -e [I]: Create a file mapping GRCh37 position to GRCh38 ones ...
  coord_map=$out_dir/temp/chr$chr.coordination_map.txt
  awk -F";" \
    'BEGIN{print "V37,V38"} {split($1, arr, "\t"); print $2" "arr[3]}' \
    $out_dir/annotation/chr$chr.GRCh38.bed \
    >| $coord_map

  echo -e [I]: Convert variant position in tags files ...
  rm -f $out_dir/tags/r0{1,8}/chr$chr
  awk -F" " -v OUTDIR=$out_dir/tags -v CHROM=chr$chr \
    -f- <<'EOF' $coord_map $tags_r01_file $tags_r08_file
NR == FNR {coordmap[$1] = $2; next}
NR != FNR {
  VAR_POS = coordmap[$1]

  LD_VAR_POS = ""
  split($2, ld_var_pos_37, ",")
  for (ii in ld_var_pos_37) {
    per_tag_pos = coordmap[ld_var_pos_37[ii]]
    if (length(per_tag_pos)) {
      if (LD_VAR_POS == "") {
        LD_VAR_POS = per_tag_pos
      } else {
        LD_VAR_POS = LD_VAR_POS","per_tag_pos
      }
    }
  }

  n = split(FILENAME, tmparr, "/")
  if (length(VAR_POS) && length(LD_VAR_POS)) {
    print VAR_POS" "LD_VAR_POS >> OUTDIR"/"tmparr[n-1]"/"CHROM
  }
}
EOF
  sort -k1,1g -t" " $out_dir/tags/r01/chr$chr > $out_dir/tags/r01/chr$chr.sorted
  mv -f $out_dir/tags/r01/chr$chr.sorted $out_dir/tags/r01/chr$chr 

  sort -k1,1g -t" " $out_dir/tags/r08/chr$chr > $out_dir/tags/r08/chr$chr.sorted
  mv -f $out_dir/tags/r08/chr$chr.sorted $out_dir/tags/r08/chr$chr 

  # TSS positions
  echo -e [I]: Format variant MAF/TSS distance into BED ...
  awk -F" " -v CHROM=chr$chr \
    '{print CHROM"\t"$1 + $3 - 1"\t"$1 + $3"\t"$1";"$2";"$3}' $maftss_file \
    >| $out_dir/maftssd/chr$chr.GRCh37.bed

  echo -e [I]: LiftOver variant MAF/TSS distance from GRCh37 to GRCh38 ...
  $liftOver_exe \
    $out_dir/maftssd/chr$chr.GRCh37.bed \
    $chain_file \
    $out_dir/maftssd/chr$chr.GRCh38.bed \
    $out_dir/maftssd/chr$chr.GRCh38.unlifted.bed

  echo -e [I]: Convert variant position in MAF/TSS distance files ...
  awk -F" " \
    -f- <<'EOF' $coord_map $out_dir/maftssd/chr$chr.GRCh38.bed \
    | grep -v '#' \
    | sort -k1,1g -t" "  \
    >| $out_dir/maftssd/chr$chr
NR == FNR {coord_map[$1] = $2; next}
NR != FNR {
  split($0, arr_1, "\t")
  split(arr_1[4], arr_2, ";")

  VAR_POS = coord_map[arr_2[1]]
  VAR_MAF = arr_2[2]
  VAR_TSSD = arr_1[3] - VAR_POS

  if (length(VAR_POS) && length(VAR_TSSD)) {
    print VAR_POS" "VAR_MAF" "VAR_TSSD
  } else {
    print "#"
  }
}
EOF

  echo -e [I]: Clean up ...
  rm -f $out_dir/{annotation,maftssd}/chr$chr.GRCh3{7,8}*
  rm -f $out_dir/temp/chr$chr.coordination_map.txt

  echo -e "[I]: ------ Done for chr$chr ------\n"
done
