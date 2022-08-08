"""
Snakefile for genotyping

Author: Zhenhua Zhang
E-mail: zhenhua.zhang217@gmail.com
Created: 2022 May 10
Updated: 2022 May 17

Adapted from: Marc Jan Bonder
Affiliation: EMBL-EBI
Date: Tuesday 10 February 2019
"""

import os

from subprocess import run

import pandas as pd

shell.prefix("set -euo pipefail;")


def _multi_arg_start(flag, files):
    flag += " "
    return " ".join(flag + f for f in files)


def _multi_arg_end(flag, files):
    flag = " "+flag
    return " ".join(f + flag for f in files)


def _multi_arg_both_ends(flag1, flag2, files):
    flag1 += " "
    flag2 = " "+flag2
    return " ".join(flag1 + f + flag2 for f in files)


def flatenChunk(chunk):
    return chunk.replace(":", "_").replace("-", "_")


def extendChunk(chunk):
    relChunk = chunk.pop()
    chunkSplitted = relChunk.split("_")
    return chunkSplitted[0]+":"+chunkSplitted[1]+"-"+chunkSplitted[2]


# Variables
projDir = os.environ.get("HOME") + '/Documents/projects/wp_bcg_eqtl'


# Inputs
# These input files are the same for each run.
inputFolder = f"{projDir}/inputs"
chunkFile = f"{inputFolder}/chunks_file.txt"
annotationFile = f"{inputFolder}/annotations/annotations_hg38.tsv"
genotypeFile = f"{inputFolder}/genotypes/300BCG_sub40_imp_hg38_clean"
kinshipFile = f"{inputFolder}/kinships/300BCG_sub40_imp_hg38_clean.kinship"


# XXX: using the --config to set these variables
runMode = config["runMode"]
cellType = config["cellType"]
evalModel = config["evalModel"]

if config["compPair"] == ".":
    compPair = ""
    inRootDir = f"{inputFolder}/{runMode}/{cellType}"
    outputFolder = f"{projDir}/outputs/{runMode}/{evalModel}/{cellType}".rstrip("/")
else:
    compPair = config["compPair"]
    inRootDir = f"{inputFolder}/{runMode}/{cellType}/{compPair}"
    outputFolder = f"{projDir}/outputs/{runMode}/{evalModel}/{cellType}/{compPair}".rstrip("/")

# These input files are different for each run.
phenotypeFile = f"{inRootDir}/phenotypes.tsv"
sampleMappingFile =  f"{inRootDir}/sample_mapping.tsv"

if "usePEER" in config and config["usePEER"] == "true":
  covariateFile = f"{inRootDir}/covariates_wpeer.tsv"
else:
  covariateFile = f"{inRootDir}/covariates.tsv"


# Tools
limixSimg = f"{projDir}/scripts/tools/limix20220510.simg"

limixPath = f"{projDir}/scripts/tools/limix_qtl"
postProcess = f"{limixPath}/Limix_QTL/post_processing/minimal_postprocess.py"
multTestCorrect = f"{limixPath}/Limix_QTL/post_processing/multTestCorrect_PipeLine.R"

intermParam = ""
if evalModel == "interaction":
    qtlAnalysis = f"{limixPath}/Limix_QTL/run_interaction_QTL_analysis.py"
    intermParam = "-it " + config["interTerm"]
elif evalModel == "hurdle":
    qtlAnalysis = f"{limixPath}/Limix_QTL/run_QTL_analysis_hurdle.py"
else:
    qtlAnalysis = f"{limixPath}/Limix_QTL/run_QTL_analysis.py"


# Parameters
numberOfPermutations = '1000'
minorAlleleFrequency = '0.1'
relatednessScore = '0.8'
windowSize = '100000'
FDR = '0.05'


# Outputs
finalQTLRun = f"{outputFolder}/qtl_results_all.txt"
topQTL = f"{outputFolder}/top_qtl_results_all.txt"
correctedQTL = f"{outputFolder}/top_qtl_results_all_FDR{FDR}.txt"

if not os.path.exists(outputFolder): # Make output folder if it doesn't exists
  os.makedirs(outputFolder, mode=0o755, exist_ok=True)


with open(chunkFile, 'r') as f:
    chunks = [x.strip() for x in f.readlines()]


qtlOutput = []
for chunk in chunks:
    #print(chunk)
    processedChunk = flatenChunk(chunk)
    #print(processedChunk)
    processedChunk = expand(outputFolder + '/{chunk}.finished', chunk=processedChunk)
    #print(processedChunk)
    qtlOutput.append(processedChunk)


## flatten these lists
qtlOutput = [filename for elem in qtlOutput for filename in elem]


#finalQtlRun = [filename for elem in finalQtlRun for filename in elem]

rule all:
    input:
        qtlOutput,finalQTLRun,topQTL,correctedQTL

rule run_qtl_mapping:
    input:
        af = annotationFile,
        pf = phenotypeFile,
        cf = covariateFile,
        rf = kinshipFile,
        smf = sampleMappingFile
    output:
        # outputFolder + '{chunk}.finished'
        outputFolder.rstrip("/") + "/{chunk}.finished"
    params:
        img = limixSimg,
        cmd = qtlAnalysis,
        gen = genotypeFile,
        od = outputFolder,
        np = numberOfPermutations,
        maf = minorAlleleFrequency,
        rs = relatednessScore,
        w = windowSize,
        interms = intermParam
    run:
        chunkFull = extendChunk({wildcards.chunk})
        shell(
            " singularity exec {params.img} python {params.cmd}"
            " -w {params.w} "
            " -gm gaussnorm "
            " -af {input.af} "
            " -pf {input.pf} "
            " -cf {input.cf} "
            " -rf {input.rf} "
            " -smf {input.smf} "
            " -gr {chunkFull} "
            " -od {params.od} "
            " -pg {params.gen} "
            " -np {params.np} "
            " -rs {params.rs} "
            " -maf {params.maf} "
            " {params.interms} "
#           " -rc " # Regressing out covariants will throw a error of non-matched matrix shape (152, 5) vs (6,) at run_interaction_QTL_analysis.py:584
            " -c "
            " && touch {output}")

rule aggregate_qtl_results:
    input:
        IF = outputFolder,
        OF = outputFolder,
        finalFiles = qtlOutput
    output:
        finalQTLRun
    params:
        img = limixSimg,
        cmd = postProcess
    run:
        shell(
            " singularity exec {params.img} python {params.cmd}"
            " -id {input.IF} "
            " -od {input.OF} "
            " -sfo ")

rule top_feature:
    input:
        IF = outputFolder,
        OF = outputFolder,
        finalFiles = qtlOutput
    output:
        topQTL
    params:
        img = limixSimg,
        cmd = postProcess
    run:
        shell(
            " singularity exec {params.img} python {params.cmd}"
            " -id {input.IF} "
            " -od {input.OF} "
            " -tfb "
            " -sfo ")

rule multi_test_correction:
    input:
        IF = outputFolder,
        OF = outputFolder,
        finalFile = topQTL
    output:
        correctedQTL
    params:
        img = limixSimg,
        cmd = multTestCorrect,
        fdr = FDR
    run:
        shell(
            "singularity exec {params.img} Rscript {params.cmd} "
            " {input.IF} "
            " {input.OF} "
            " {params.fdr} ")
