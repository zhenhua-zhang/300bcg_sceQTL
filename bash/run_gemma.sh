#!/bin/bash

proj_dir=~/Documents/projects/wp_bcg_eqtl
pct=Monocytes

## Prepare genotypes

python3 prepare_gemma_input.py \
  -g $vcfpath \
  -p $exppath \
  -w $winsize \
  -o $


# Run GEMMA of linear mixed model
gemma \
  -g $snppath \
  -p $exppath \
  -a $annpath \
  -k $kinpath \
  -lmm \
  -gxe $gbepath \
  -o

less mouse_hs1940.pheno.txt

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
