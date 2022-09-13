#!/bin/bash

proj_dir=~/Documents/projects/wp_bcg_eqtl
pct=Monocytes

source $proj_dir/scripts/.env/bin/activate

# Prepare inputs gene expression matrix
python3 $proj_dir/scripts/py3/prepare_gemma_input.py \
  -g abc \
  -p abc \
  -f /vol/projects/zzhang/projects/wp_bcg_eqtl/inputs/reference/Gencode/gencode.v41.basic.annotation.gff3.gz \
  -c abc \
  -o $proj_dir/temps/mean_genotype.geno

#
vcfpath=$proj_dir/inputs/genotypes/300BCG_sub40_imp_hg38_ids_clean.vcf.gz
plink --vcf $vcfpath --recode bimbam --out $proj_dir/outputs/genotypes/BIMBAM/300BCG_sub40_imp_hg38_clean

## Prepare genotypes


# Run GEMMA of linear mixed model
gemma \
  -g $snppath \
  -p $exppath \
  -a $annpath \
  -k $kinpath \
  -gxe $gbepath \
  -lmm \
  -notsnp \
  -o $outdir

gemma -gk \
  -g mouse_hs1940.geno.txt.gz \
  -p mouse_hs1940.pheno.txt \
  -a mouse_hs1940.anno.txt \
  -o mouse_hs1940

gemma -lmm \
  -g mouse_hs1940.geno.txt.gz \
  -p mouse_hs1940.pheno.txt \
  -a mouse_hs1940.anno.txt \
  -k output/mouse_hs1940.cXX.txt \
  -n 1 \
  -gxe mouse_hs1940.gxe.txt \
  -o mouse_hs1940_CD8MCH_lmm
