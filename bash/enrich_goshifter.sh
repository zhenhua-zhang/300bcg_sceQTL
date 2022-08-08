# Annotation enrichment
# CIS-BP database version: 8-1-2019 Build version 2.00
proj_dir=~/Documents/projects/wp_bcg_eqtl

wkdir=$proj_dir/outputs/pseudo_bulk/archived/20220622
reg_file=$proj_dir/inputs/annotations/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.gff
for regtype in enhancer promoter{,_flanking_region} {CTCF,TF}_binding_site open_chromatin_region; do
  awk -v tar_reg=$regtype \
    -f- <<'EOF' $reg_file >| $proj_dir/inputs/annotations/$regtype.bed
    /^[0-9]{1,2}/ && $3 == tar_reg {print "chr"$1"\t"$4"\t"$5"\t"$9}
EOF

  for mode in normal interaction; do
    for cell_type in Monocytes CD{4,8}T B NK; do
      for condition in T0_RPMI.vs.{T0_LPS,T3m_RPMI} T3m_RPMI.vs.T3m_LPS; do
        awk 'NR==1 {print "SNP\tChrom\tBP"; next} {print $1"\tchr"$16"\t"$17}' \
          $wkdir/$mode/$cell_type/$condition/top_qtl_results_all_FDR0.05.txt \
          >| $wkdir/$mode/$cell_type/$condition-snp_pos.txt

        singularity exec goshifter.simg goshifter.py \
          --snpmap $wkdir/$mode/$cell_type/$condition-snp_pos.txt \
          --annotation $proj_dir/inputs/annotations/$regtype.bed.gz \
          --permute 2000 \
          --ld $proj_dir/inputs/annotations/goshifter_ld \
          --out $wkdir/$mode/$cell_type/$condition-$regtype
      done
    done
  done
done
