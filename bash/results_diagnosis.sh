# Create a table of zero ratio per gene
proj_dir=~/Documents/projects/wp_bcg_eqtl

awk -f- <<'EOF' $proj_dir/inputs/phenotypes/B_T0_LPS.tsv \
  | sort -k2,2g -k3,3g \
  >| $proj_dir/outputs/limix_qtl/diagnosis/zero_ratio_per_gene_B_T0_LPS.tsv
BEGIN {
  print "FeatureId\tZerosNr\tZerosRatio"
} {
  nr_zeros = 0
  for (i=2; i<=NF; i++) {
    nr_zeros += $i == 0
  }

  ra_zeros = nr_zeros/(NF - 1.0)
  if (ra_zeros != 1 && $1 != "feature_id")
    print $1"\t"nr_zeros"\t"ra_zeros
}
EOF
