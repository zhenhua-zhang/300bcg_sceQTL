#!/usr/bin/env bash
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: 2022 May 15
# Updated: 2022 May 15

proj_dir=~/Documents/projects/wp_bcg_eqtl

echo "[I]: Filter genotypes to make sure every genotype have at least three samples."
bcftools view \
  -O z \
  -v snps \
  -i '%ID!="." && COUNT(GT=="het")>2 && COUNT(GT=="RR")>2 && COUNT(GT=="AA")>2' \
  --threads 4 \
  -o $proj_dir/inputs/genotypes/300BCG_sub40_imp_hg38_ids_clean.vcf.gz \
  $proj_dir/inputs/genotypes/300BCG_sub40_imp_hg38_ids.vcf.gz


echo "[I]: Convert VCF into PLINK fam/bim/bed format..."
plink --vcf $proj_dir/inputs/genotypes/300BCG_sub40_imp_hg38_ids_clean.vcf.gz \
  --make-bed \
  --out $proj_dir/inputs/genotypes/300BCG_sub40_imp_hg38_ids_clean

echo "[I]: Filter genotypes..."
plink --bfile $proj_dir/inputs/genotypes/300BCG_sub40_imp_hg38_ids_clean \
  --maf 5e-2 \
  --hwe 1e-6 \
  --hwe-all \
  --indep-pairwise 250 50 0.2 \
  --make-bed \
  --out $proj_dir/inputs/genotypes/300BCG_sub40_imp_hg38_clean

echo "[I]: Calculate kinship matrix..."
plink --bfile $proj_dir/inputs/genotypes/300BCG_sub40_imp_hg38_clean \
  --make-rel square \
  --out $proj_dir/inputs/genotypes/300BCG_sub40_imp_hg38_clean

echo "[I]: Merge .rel.id and .rel matrix..."
awk -f- <<'EOF' $proj_dir/inputs/genotypes/300BCG_sub40_imp_hg38_clean.rel{.id,} \
  >| $proj_dir/inputs/kinships/300BCG_sub40_imp_hg38_clean.kinship
NR==FNR {map[NR]=$1; next}
NR>FNR {
  if (header_printed == 0) {
    printf "sample_id"
    for (i in map) {printf "\t"map[i]}
    printf "\n"
    header_printed = 1
  }
  print map[FNR]"\t"$0
}
EOF

echo "[I]: Clean up..."
rm -f \
  $proj_dir/inputs/genotypes/300BCG_sub40_imp_hg38_ids{,_clean}.{bed,bim,fam,nosex,log} \
  $proj_dir/inputs/genotypes/300BCG_sub40_imp_hg38_clean.{nosex,log,prune.{out,in},rel{,.id}}

echo "Done!"
