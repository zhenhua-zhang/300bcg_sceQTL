
import glob
import os
from subprocess import run
import pandas as pd
import re
from os.path import join

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

#Variables
chunkFile = './scTI_chunks2.txt'
genotypeFile = './genotype'
annotationFile = './anno.nochr.txt'
phenotypeFile = '.scMonoTI.txt'
covariateFile = './scCov.txt'
kinshipFile = './scKinship.txt'
sampleMappingFile = './scSampleMapping.txt'

numberOfPermutations = '1000'
minorAlleleFrequency = '0.1'
windowSize = '1000000'
FDR = '0.1'
outputFolder = './scTIMono/'
finalQTLRun = './scTIMono/qtl_results_all.txt'
topQTL = './top_qtl_results_all.txt'
correctedQTL = "".join(['./top_qtl_results_all_FDR',FDR,'.txt'])



with open(chunkFile,'r') as f:
    chunks = [x.strip() for x in f.readlines()]

qtlOutput = []
for chunk in chunks:
    #print(chunk)
    processedChunk = flatenChunk(chunk)
    #print(processedChunk)
    processedChunk=expand(outputFolder+'{chunk}.finished',chunk=processedChunk )
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
        outputFolder + '{chunk}.finished'
    params:
        gen=genotypeFile,
        od = outputFolder,
        np = numberOfPermutations,
        maf = minorAlleleFrequency,
        w = windowSize,
    run:
        chunkFull = extendChunk({wildcards.chunk})
        shell(
	    " ./Python-3.9.8/bin/python3.9 ./run_QTL_analysis.py "
            " --plink {params.gen} "
            " -af {input.af} "
            " -pf {input.pf} "
            " -cf {input.cf} "
            " -od {params.od} "
            " -rf {input.rf} "
            " --sample_mapping_file {input.smf} "
            " -gr {chunkFull} "
            " -np {params.np} "
            " -maf {params.maf} "
            " -c -gm gaussnorm "
            " -w {params.w} "
            " -rs 0.8 "
            " -rc ")
        shell("touch {output}")

rule aggregate_qtl_results:
    input:
        IF = outputFolder,
        OF = outputFolder,
        finalFiles = qtlOutput
    output:
        finalQTLRun
    run:
        shell(
 	    " ./Python-3.9.8/bin/python3.9 ./minimal_postprocess.py "
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
    run:
        shell(
	    " ./Python-3.9.8/bin/python3.9 ./minimal_postprocess.py "
            "-id {input.IF} "
            "-od {input.OF} "
            "-tfb "
            "-sfo ")

rule multi_test_correction:
    input:
        IF = outputFolder,
        OF = outputFolder,
        finalFile = topQTL
    output:
        correctedQTL
    params:
        fdr = FDR
    run:
        shell(
	    "./R-4.0.2/bin/Rscript ./multTestCorrect_PipeLine.R " 
            "{input.IF} "
            "{input.OF} "
            "{params.fdr} ")
