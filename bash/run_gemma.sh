#!/bin/bash
#SBATCH --mem 4G
#SBATCH --array 1-5
#SBATCH --time 5:0:0
#SBATCH --partition cpu
#SBATCH --cpus-per-task 4
#SBATCH --output %A_%a-run_gemma.log

proj_dir=~/Documents/projects/wp_bcg_eqtl

source $proj_dir/scripts/.env/bin/activate

celltype=Monocytes
ia_term=time

# Prepare inputs gene expression matrix
python3 $proj_dir/scripts/py3/prepare_gemma_input.py \
  -g $proj_dir/inputs/genotypes/300BCG_sub40_imp_hg38_ids_clean.vcf.gz \
  -p $proj_dir/inputs/pseudo_bulk/$celltype/phenotypes.tsv \
  -c $proj_dir/inputs/pseudo_bulk/$celltype/covariates_wpeer.tsv \
  -m $proj_dir/inputs/pseudo_bulk/$celltype/sample_mapping.tsv \
  -f $proj_dir/inputs/reference/Gencode/gencode.v41.basic.annotation.Ec.transcript.autosome.gff.gz \
  -I $ia_term \
  -o $proj_dir/temps/$celltype

# Combined all genotypes and annotations into one
cat $proj_dir/temps/$celltype/block_*/*_genotype.txt >| $proj_dir/temps/genotype.txt
cat $proj_dir/temps/$celltype/block_*/*_genotype_annotation.txt >| $proj_dir/temps/genotype_annotation.txt

# Relatedness matrix from all used genotypes
gemma -gk 2 \
  -g $proj_dir/temps/Gemma/genotype.txt \
  -a $proj_dir/temps/Gemma/genotype_annotation.txt \
  -p $proj_dir/temps/Gemma/blocks/block_00000/PNRC2_phenotype.txt \
  -o kinship \
  -outdir $proj_dir/temps/Gemma


for per_input in $proj_dir/temps/Gemma/blocks/block_*/*_genotype.txt; do
  per_block=$(dirname $per_input)
  per_gene=$(basename $per_input | cut -f1 -d_)

  # Run GEMMA of linear mixed model
  gemma -lmm 3 \
    -g $per_block/${per_gene}_genotype.txt \
    -a $per_block/${per_gene}_genotype_annotation.txt \
    -p $per_block/${per_gene}_phenotype.txt \
    -k $proj_dir/temps/Gemma/kinship.abs.cXX.txt \
    -c $proj_dir/temps/Gemma/covariate.txt \
    -o ${per_gene}_QTL_common \
    -outdir $per_block

  gemma -lmm 1 \
    -g $per_block/${per_gene}_genotype.txt \
    -a $per_block/${per_gene}_genotype_annotation.txt \
    -p $per_block/${per_gene}_phenotype.txt \
    -k $proj_dir/temps/Gemma/kinship.abs.cXX.txt \
    -c $proj_dir/temps/Gemma/covariate.txt \
    -o ${per_gene}_QTL_interaction \
    -gxe $proj_dir/temps/Gemma/interaction_covariates-$ia_term.txt \
    -outdir $per_block

  break

  if [[ $count -le 8 ]]; then
    count=$(( $count + 1 )) && continue
  else
    count=0 && wait
  fi
done
