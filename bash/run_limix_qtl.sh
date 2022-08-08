#!/usr/bin/env bash
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: 2022 Apr 22
# Updated: 2022 Jul 27


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

if [[ $DEBUG == true ]]; then
  njobs=10
else
  njobs=200
fi


# Working path.
proj_dir=~/Documents/projects/wp_bcg_eqtl


# Prepare inputs for each pairs
if [[ prepare_inputs == true ]]; then
  for celltype in Monocytes CD4T CD8T NK B; do
    for cmpair in {T0_LPS,T3m_LPS,T3m_RPMI}.vs.T0_RPMI T3m_LPS.vs.T3m_RPMI; do
      case_cond=$(cut -f1 -d. <<<$cmpair)
      ctrl_cond=$(cut -f3 -d. <<<$cmpair)

      ctype=$proj_dir/inputs/pseudo_bulk/$celltype
      mkdir -p $ctype/$cmpair

      # Expression matrix
      awk -v CC=$case_cond -v CT=$ctrl_cond \
        -f- <<'EOF' $ctype/phenotypes.tsv >| $ctype/$cmpair/phenotypes.tsv
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
      $ctype/covariates.tsv >| $ctype/$cmpair/covariates.tsv

    # covariates with PEERs
    grep -e sample_id -e $case_cond -e $ctrl_cond \
      $ctype/covariates_wpeer.tsv >| $ctype/$cmpair/covariates_wpeer.tsv

    # Sample to individuals mapping
    grep -e $case_cond -e $ctrl_cond \
      $ctype/sample_mapping.tsv >| $ctype/$cmpair/sample_mapping.tsv
    done
  done
fi


# Use a conda environment
eval "$(conda shell.bash hook)"
conda deactivate
conda activate limix_qtl


# Run the pipeline by Snakemake
mode=pseudo_bulk
tmstmp=$(date +%Y%m%d%H%M%S)
smfile=$proj_dir/scripts/snakemake/300bcg_limix_qtl.smk

# Submit jobs to do the regression analysis.
# for model in normal interaction; do
for model in normal; do
  for ctype in Monocytes CD4T CD8T NK B; do
    out_dir=$proj_dir/outputs/$mode/$model/$ctype

    # Save the log files.
    if [[ ! -d $out_dir/logs ]]; then mkdir -p $out_dir/logs; fi

    # Unlock the control folder of Snakemake.
    snakemake --unlock -s $smfile \
      -C runMode=. cellType=. compPair=. evalModel=. interTerm=.

    if [[ $model == normal ]]; then
      # Include all samples to estimate shared genetic effects
      dr_log=$out_dir/logs/${tmstmp}_${model}_smdry.log
      snakemake -r -n -c 1 -s $smfile \
        -C runMode=$mode cellType=$ctype compPair=. evalModel=$model usePEER=true \
        1>&2 >| $dr_log

      # Run
      errf=$out_dir/logs/%j-${tmstmp}_${model}_sbatch.err
      outf=$out_dir/logs/%j-${tmstmp}_${model}_sbatch.out
      snakemake -k -j $njobs -w 60 -s $smfile \
        -C runMode=$mode cellType=$ctype compPair=. evalModel=$model usePEER=true \
        --cluster 'sbatch -t 1:59:0 -p cpu -o '$outf' -e '$errf' -J '$model-$ctype' -N 1 --cpus-per-task 1 --ntasks 1 --mem 4G'
    else
      # Runs pair-wise comparison for the "interaction" model.
      for cmppair in {T0_LPS,T3m_LPS,T3m_RPMI}.vs.T0_RPMI T3m_LPS.vs.T3m_RPMI; do
        icterm=$(check_term $cmppair)

        # Snakemake dry-run log.
        dr_log=$out_dir/logs/${tmstmp}_${model}_${cmppair}_smdry.log
        snakemake -r -n -c 1 -s $smfile \
          -C runMode=$mode cellType=$ctype compPair=$cmppair evalModel=$model interTerm=$icterm usePEER=true \
          1>&2 >| $dr_log

        # Run
        errf=$out_dir/logs/%j-${tmstmp}_${model}_${cmppair}_sbatch.err
        outf=$out_dir/logs/%j-${tmstmp}_${model}_${cmppair}_sbatch.out
        snakemake -k -j $njobs -w 60 -s $smfile \
          -C runMode=$mode cellType=$ctype compPair=$cmppair evalModel=$model interTerm=$icterm usePEER=true \
          --cluster 'sbatch -t 1:59:0 -p cpu -o '$outf' -e '$errf' -J '$model-$ctype-$cmppair' -N 1 --cpus-per-task 1 --ntasks 1 --mem 4G'
        # break
      done
    fi
    # break
  done
  # break
done
