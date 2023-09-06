#!/usr/bin/env
# File: estimate_independent_gwas_loci.sh
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang3@helmholtz-hzi.de
# Created: Apr 28, 2023
# Updated: Aug 10, 2023

proj_dir=~/Documents/projects/wp_bcg_eqtl
temp_dir=$proj_dir/temps

reference_panel=/vol/projects/BIIM/resources/1000G/EUR/GRCh38/EUR
indep_loci_file=$proj_dir/outputs/public_gwas/all_of_indenpendent_loci_per_gwas.txt

for per_gwas in $proj_dir/inputs/public_gwas/*.gz; do
  gwas_name=$(basename $per_gwas | sed 's/.vcf.gz//g')
  temp_file=$temp_dir/public_gwas_${gwas_name}.txt
  awk -f- <<'EOF' <(zcat $per_gwas) >| $temp_file
$0 ~ /^##/ {next}
$0 ~ /^#CHROM/ {print "SNP P"; next}
{
  n = split($9, name_vec, ":");
  split($10, value_vec, ":");
  for (ii in name_vec) { FORMAT_MAP[name_vec[ii]] = value_vec[ii] }
  p_value = 10 ** (-FORMAT_MAP["LP"]);
  print $3" "p_value
}
EOF

  plink --bfile $reference_panel --clump $temp_file --clump-p1 5e-5 --clump-p2 1e-2 --clump-r2 0.1 --clump-kb 500 --out ${temp_file/.txt/} 2> /dev/null

  gwas_name=$(basename $per_gwas | sed 's/.vcf.gz//g')
  if [[ ! -e $indep_loci_file ]]; then echo -e "GWAS CHR F SNP BP P TOTAL NSIG S05 S01 S001 S0001 SP2" >| $indep_loci_file; fi
  awk -v GWAS=$gwas_name 'NR==1 || $0 ~ /^$/ {next} {gsub(/[ ]+/, " "); print GWAS""$0}' ${temp_file/.txt/.clumped} >> $indep_loci_file

  rm -f ${temp_file/.txt/}*
done
