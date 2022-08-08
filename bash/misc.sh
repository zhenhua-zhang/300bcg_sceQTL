#!/bin/bash
proj_dir=~/Documents/projects/wp_bcg_eqtl

source $proj_dir/scripts/.env/bin/activate

mode=normal
cell_type=B
condition=T0_LPS.vs.T0_RPMI

for cell_type in Monocytes CD4T CD8T B NK; do
  for condition in T0_LPS.vs.T0_RPMI T3m_LPS.vs.T0_RPMI T3m_RPMI.vs.T0_RPMI T3m_LPS.vs.T3m_RPMI; do
    in_dir=$proj_dir/inputs/pseudo_bulk/$cell_type/$condition
    out_dir=$proj_dir/outputs/pseudo_bulk/outcomes/$mode/${mode}_${cell_type}_${condition}/box_plots

    case $condition in
      T0_LPS.vs.T0_RPMI | T3m_LPS.vs.T3m_RPMI)
        group_by=stim; cond_map="0:RPMI 1:LPS" ;;
      T3m_LPS.vs.T0_RPMI | T3m_RPMI.vs.T0_RPMI)
        group_by=time; cond_map="0:T0 1:T3m" ;;
    esac

    mkdir -p $out_dir
    python $proj_dir/scripts/py3/plot_snp_effects.py \
      -x "IRF2:rs2310003" "ADCY3:rs11687089" \
      -g $proj_dir/inputs/genotypes/300BCG_sub40_imp_hg38_ids_clean.vcf.gz \
      -p $in_dir/phenotypes.tsv \
      -c $in_dir/covariates.tsv \
      -m $in_dir/sample_mapping.tsv \
      --cond-col $group_by \
      --cond-map $cond_map \
      --fig-size 3 7 \
      -o $out_dir

    echo $out_dir
  done
done
