#!/usr/bin/env
#Author: Zhenhua Zhang
#E-mail: zhenhua.zhang3@helmholtz-hzi.de
#Created: 2023 Apr 28
#Updated: 2023 Apr 28

proj_dir=~/Documents/projects/wp_bcg_eqtl
temp_dir=$proj_dir/temps

reference_panel=/vol/projects/BIIM/resources/1000G/EUR/GRCh38/EUR
nr_independent_loci=$proj_dir/outputs/public_gwas/number_of_indenpendent_loci_per_gwas.csv

for per_gwas in $proj_dir/inputs/public_gwas/*.gz; do
  gwas_name=$(basename $per_gwas | sed 's/.vcf.gz//g')
  temp_file=$temp_dir/public_gwas_${gwas_name}.txt
  awk -f- <<'EOF' <(zcat $per_gwas) >| $temp_file
$0~/^##/ {next}
$0~/^#CHROM/ {print "SNP P"; next}
{
  n = split($9, name_vec, ":");
  split($10, value_vec, ":");
  for (ii in name_vec) {
    FORMAT_MAP[name_vec[ii]] = value_vec[ii]
  }

  p_value = 10 ** (-FORMAT_MAP["LP"])
  print $3" "p_value
}
EOF

  plink --bfile $reference_panel \
    --clump $temp_file \
    --clump-p1 5e-8 \
    --clump-p2 5e-2 \
    --clump-r2 0.1 \
    --clump-kb 500 \
    --out ${temp_file/.txt/}

  n_independent_loci=$(grep -c . ${temp_file/.txt/.clumped})

  if [[ ! -e $nr_independent_loci ]]; then
    echo -e "GWAS,nr_independent_loci" >| $nr_independent_loci
  fi

  gwas_name=$(basename $per_gwas | sed 's/.vcf.gz//g')
  echo -e "$gwas_name,$n_independent_loci" >> $nr_independent_loci

  rm -f ${temp_file/.txt/}*
done
