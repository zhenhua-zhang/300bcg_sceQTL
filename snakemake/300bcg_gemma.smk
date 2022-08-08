#
## Parameters
#

import os

projdir = config["projdir"]
sigimg = None


rule all:
  top_qtls


rule run_post_processing:
  input:
    ""
  output:
    ""
  shell:
    """
    singularity exec {sigimg}
    """


# Run gemma under linear mixed model.
rule run_gemma_lmm:
  input:
    geno = "{projdir}/inputs/{chunk}.geno.txt.gz",
    anno = "{projdir}/inputs/{chunk}.annot.txt",
    pheno = "{projdir}/inputs/{chunk}.pheno.txt",
    covar = "{projdir}/inputs/{chunk}.covar.txt"
    kinship = "{projdir}/outputs/{chunk}.cXX.txt",
  output:
    "{projdir}/outputs/{chunk}.assoc.txt"


# Create a genotype kinship matrix.
rule run_gemma_gk:
  input:
    geno_pc = "{projdir}/inputs/{chunk}.geno.txt.gz",
    anno_pc = "{projdir}/inputs/{chunk}.annot.txt",
    pheno_pc = "{projdir}/inputs/{chunk}.pheno.txt"
  output:
    kinship = "{projdir}/outputs/{chunk}.cXX.txt"
  params:
    binddir = ""
  shell:
    """
    singularity -b {binddir} exec {sigimg} \
        gemma -gk \
        -g {input.geno_pc} \
        -a {input.anno_pc} \
        -p {input.pheno_pc} \
        -o {output.kinship}
    """


# Prepare gemma input files.
rule prepare_gemma_input:
  input:
    ttl_geno_file = "",
    ttl_pheno_file = "",
    ttl_covar_file = "",
  output:
    ""
  params:
    outdir = f"{projdir}"
  shell:
    """
    singularity -b {binddir} exec {sigimg} \
        python3 prepare_gemma_input.py \
        -g {input.} \
        -p {input.} \
        -c {input.} \
        -o {params.}
        {chunk}
    """


# Create a chunk file
rule prepare_chunks:
  input:
    gtf_file = ""
  output:
    chunks_file = ""
  shell:
    """
    """
