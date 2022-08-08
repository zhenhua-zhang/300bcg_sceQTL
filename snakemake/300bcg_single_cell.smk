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
import re
import glob

from subprocess import run
from os.path import join

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
# XXX: using the --cofig to set these variables
cell_type = config["cell"]
time_point = config["time"]
stimulation = config["stim"]

runBarcode = "_".join([cell_type, time_point, stimulation])

inputFile = runBarcode + ".tsv"
inputFolder = join(projDir, "inputs")
chunkFile = join(inputFolder, "chunks_file.txt")
annotationFile = join(inputFolder, "annotations/annotations_hg38.tsv")
genotypeFile = join(inputFolder, "genotypes/300BCG_sub40_imp_hg38_clean")
phenotypeFile = join(inputFolder, "phenotypes", inputFile)
covariateFile = join(inputFolder, "covariates", inputFile)
kinshipFile = join(inputFolder, "kinships/300BCG_sub40_imp_hg38_clean.kinship")
sampleMappingFile = join(inputFolder, "sample_mapping", inputFile)


# Tools
limixSimg = join(projDir, "scripts/tools/limix20220510.simg")

limixPath = join(projDir, "scripts/tools/limix_qtl")
qtlAnalysis = join(limixPath, "Limix_QTL/run_QTL_analysis.py")
postProcess = join(limixPath, "Limix_QTL/post_processing/minimal_postprocess.py")
multTestCorrect = join(limixPath, "Limix_QTL/post_processing/multTestCorrect_PipeLine.R")


# Parameters
numberOfPermutations = '1000'
minorAlleleFrequency = '0.1'
relatednessScore = '0.8'
windowSize = '100000'
FDR = '0.05'


# Outputs
outputFolder = join(projDir, "outputs/limix_qtl", runBarcode)
finalQTLRun = join(outputFolder, "qtl_results_all.txt")
topQTL = join(outputFolder, "top_qtl_results_all.txt")
correctedQTL = join(outputFolder, "top_qtl_results_all_FDR" + FDR + ".txt")

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
        join(outputFolder, '{chunk}.finished')
    params:
        img = limixSimg,
        cmd = qtlAnalysis,
        gen = genotypeFile,
        od = outputFolder,
        np = numberOfPermutations,
        maf = minorAlleleFrequency,
        rs = relatednessScore,
        w = windowSize,
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
            " -rc "
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
