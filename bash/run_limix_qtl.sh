#!/usr/bin/env bash
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: 2022 Apr 22
# Updated: 2023 Feb 08

# A function to set comparison pairs.
check_term() {
  case $1 in
    T3m_LPS.vs.T0_LPS | T3m_RPMI.vs.T0_RPMI | T3m_LPS.vs.T0_RPMI | time)
      echo time ;;
    T3m_LPS.vs.T3m_RPMI | T0_LPS.vs.T0_RPMI | stim)
      echo stim ;;
    *)
      echo "." ;;
  esac
}


# Debug mode.
DEBUG=false
if [[ $DEBUG == true ]]; then njobs=1; else njobs=200; fi


# Working path.
proj_dir=~/Documents/projects/wp_bcg_eqtl


#
## Create chunks file
#
if [[ ! -e $proj_dir/inputs/annotations/chunks_file.txt ]]; then
  annfile=$proj_dir/inputs/annotations/annotations_hg38.tsv
  chunkfile=$proj_dir/inputs/annotations/chunks_file.txt

  Rscript - $annfile $chunkfile <<'EOF'
library(magrittr)

options <- commandArgs(trailingOnly = TRUE)
if (length(options) != 2) stop("Require two positional options: in-file and out-file!")

anntab <- data.table::fread(options[1])
max_pos <- max(anntab$start, anntab$end) + 1000
anntab %>%
  dplyr::mutate(tss = dplyr::if_else(feature_strand == "+", start, end)) %>%
  dplyr::arrange(chromosome, tss) %>%
  dplyr::group_by(chromosome) %>%
  dplyr::mutate(group_id = dplyr::row_number() %/% 15) %>%
  dplyr::group_by(chromosome, group_id) %>%
  dplyr::summarize(start = min(tss) - 1000) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(chromosome) %>%
  dplyr::mutate(start = start, end = dplyr::lead(start, default = max_pos) - 1, chunk = paste0(chromosome, ":", start, "-", end)) %>%
  dplyr::ungroup() %>%
  dplyr::select(chunk) %>%
  data.table::fwrite(options[2], col.names = FALSE)

cat("Done\n")
EOF
fi


# Prepare inputs for each pairs
if [[ prepare_inputs == true ]]; then
  for cell_type in Monocytes CD4T CD8T NK B; do
    for per_cmp in {T0_LPS,T3m_LPS,T3m_RPMI}.vs.T0_RPMI T3m_LPS.vs.T3m_RPMI; do
      case_cond=$(cut -f1 -d. <<<$per_cmp)
      ctrl_cond=$(cut -f3 -d. <<<$per_cmp)

      per_celltype=$proj_dir/inputs/pseudo_bulk/$cell_type
      mkdir -p $per_celltype/$per_cmp

      # Expression matrix
      awk -v CC=$case_cond -v CT=$ctrl_cond \
        -f- <<'EOF' $per_celltype/phenotypes.tsv >| $per_celltype/$per_cmp/phenotypes.tsv
NR == 1 {
  split($0, header, "\t")
  for (i in header) {
    if (i == 1) {
      printf header[i]
      tar_col[i] = i
    } else if (header[i] ~ CC || header[i] ~ CT) {
      printf "\t"header[i]
      tar_col[i] = i
    }
  }
  printf "\n"
}
NR >= 2 {
  for (i in tar_col) {
    if (i == 1) { printf $tar_col[i] } else { printf "\t"$tar_col[i] }
  }
  printf "\n"
}
EOF

    # General covariates
    grep -e sample_id -e $case_cond -e $ctrl_cond \
      $per_celltype/covariates.tsv >| $per_celltype/$per_cmp/covariates.tsv

    # covariates with PEERs
    grep -e sample_id -e $case_cond -e $ctrl_cond \
      $per_celltype/covariates_wpeer.tsv >| $per_celltype/$per_cmp/covariates_wpeer.tsv

    # Sample to individuals mapping
    grep -e $case_cond -e $ctrl_cond \
      $per_celltype/sample_mapping.tsv >| $per_celltype/$per_cmp/sample_mapping.tsv
    done
  done
fi


# Use a conda environment
eval "$(conda shell.bash hook)"
conda deactivate
conda activate limix_qtl

#
## Run the pipeline by Snakemake, pseudo-bulk
#
# Submit jobs to do the regression analysis.
mode=pseudo_bulk
smfile=$proj_dir/scripts/snakemake/300bcg_limix_qtl.smk
smconf=$proj_dir/scripts/snakemake/configs
# for model in normal interaction per_condition; do
for model in per_condition; do
  # for per_celltype in Monocytes CD4T CD8T NK B; do
  for per_celltype in CD8T CD4T NK; do
    out_dir=$proj_dir/outputs/$mode/summary_statistic/$model/$per_celltype

    # Save the log files.
    if [[ ! -d $out_dir/logs ]]; then mkdir -p $out_dir/logs; fi

    # Unlock the control folder of Snakemake.
    snakemake --unlock -d $out_dir -s $smfile -C run_mode=. cell_type=. condition=. eval_model=. inter_term=.

    if [[ $model == normal ]]; then
      # Include all samples to estimate shared genetic effects
      snakemake -nrc 1 -d $out_dir -s $smfile --profile $smconf \
        -C run_mode=$mode cell_type=$per_celltype condition=. eval_model=$model use_peer=true \
        1>&2 >| $out_dir/logs/${model}_smdry.log

      # Run
      snakemake -j $njobs -d $out_dir -s $smfile --profile $smconf \
        -C run_mode=$mode cell_type=$per_celltype condition=. eval_model=$model use_peer=true
    elif [[ $model == interaction ]]; then
      # Runs pair-wise comparison for the "interaction" model. {T0_LPS,T3m_LPS,T3m_RPMI}.vs.T0_RPMI T3m_LPS.vs.T3m_RPMI
      for per_cmp in {T0_LPS,T3m_RPMI}.vs.T0_RPMI T3m_LPS.vs.T3m_RPMI; do
        icterm=$(check_term $per_cmp)

        # Snakemake dry-run log.
        snakemake -nrc 1 -d $out_dir -s $smfile --profile $smconf \
          -C run_mode=$mode cell_type=$per_celltype condition=$per_cmp eval_model=$model inter_term=$icterm use_peer=true \
          1>&2 >| $out_dir/logs/${model}_${per_cmp}_smdry.log

        # Run
        snakemake -j $njobs -d $out_dir -s $smfile --profile $smconf \
          -C run_mode=$mode cell_type=$per_celltype condition=$per_cmp eval_model=$model inter_term=$icterm use_peer=true
      done
    elif [[ $model == per_condition ]]; then
      # Runs per condition.
      for per_cond in T0_LPS T0_RPMI T3m_LPS T3m_RPMI; do
        snakemake -nrc 1 -d $out_dir -s $smfile --profile $smconf \
          -C run_mode=$mode cell_type=$per_celltype condition=$per_cond eval_model=$model inter_term=. use_peer=true \
          1>&2 >| $out_dir/logs/${model}_${per_cond}_smdry.log

        # Run
        snakemake -j $njobs -d $out_dir -s $smfile --profile $smconf \
          -C run_mode=$mode cell_type=$per_celltype condition=$per_cond eval_model=$model inter_term=. use_peer=true
        break
      done
    fi
  done
  break
done
