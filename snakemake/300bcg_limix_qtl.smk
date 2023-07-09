"""
Snakefile for eQTL mapping

Author: Zhenhua Zhang
E-mail: zhenhua.zhang217@gmail.com
Created: 2022 May 10
Updated: 2023 Feb 08

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


def dynamic_slurm_time_limit(wc):
    long_long_long_time_chunk = [ "6_31300000_32300000" ]
    long_long_time_chunk = [
        "6_30300000_31300000", "6_32300000_33300000", "7_2200000_3210000",
        "8_143000000_144000000", "11_0_1130000", "11_118000000_119000000", "11_62400000_63400000",
        "11_65400000_66400000", "12_6070000_7070000", "17_1150000_2160000", "17_4180000_5200000",
        "17_7220000_8230000", "17_75000000_76000000", "17_81100000_82100000", "19_12200000_13200000",
        "19_17200000_18200000", "19_18200000_19300000", "19_35400000_36400000", "19_36400000_37400000",
        "19_4140000_5150000", "19_48500000_49500000", "19_51500000_52500000", "19_52500000_53500000",
        "19_56500000_57600000", "19_57600000_58600000", "19_7160000_8170000",
    ]
    long_time_chunk = [
        "1_161000000_162000000", "2_27000000_28000000", "7_99700000_101000000", "8_144000000_145000000",
        "9_136000000_137000000", "9_137000000_138000000", "14_23400000_24400000", "16_0_1020000",
        "16_1020000_2030000", "16_2030000_3030000", "16_4030000_5030000", "17_44600000_45700000",
        "17_82100000_83200000", "19_10200000_11200000", "19_1110000_2120000", "19_1110000_2120000",
        "19_2120000_3130000", "19_3130000_4140000", "19_43400000_44500000", "19_45500000_46500000",
        "19_46500000_47500000", "19_48600_1110000", "19_49500000_50500000", "19_6150000_7160000",
        "22_38400000_39400000", "22_49700000_50800000",
    ]

    is_long_long_long_job = wc.chunk in long_long_long_time_chunk
    if is_long_long_long_job:
        return "47:59:59"

    is_long_long_job = wc.chunk in long_long_time_chunk
    if is_long_long_job:
        return '23:59:59'

    is_long_job = wc.chunk in long_time_chunk
    if is_long_job:
        return '11:59:59'

    return '4:59:59'


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
run_mode = config["run_mode"]
cell_type = config["cell_type"]
eval_model = config["eval_model"]

if config["condition"] == ".":
    condition = ""
    inRootDir = f"{inputFolder}/{run_mode}/{cell_type}".rstrip("/")
    outputFolder = f"{projDir}/outputs/{run_mode}/summary_statistic/{eval_model}/{cell_type}".rstrip("/")
else:
    condition = config["condition"]
    inRootDir = f"{inputFolder}/{run_mode}/{cell_type}/{condition}".rstrip("/")
    outputFolder = f"{projDir}/outputs/{run_mode}/summary_statistic/{eval_model}/{cell_type}/{condition}".rstrip("/")

# These input files are different for each run.
phenotypeFile = f"{inRootDir}/phenotypes.tsv"
sampleMappingFile =  f"{inRootDir}/sample_mapping.tsv"

if "use_peer" in config and config["use_peer"] == "true":
  covariateFile = f"{inRootDir}/covariates_wpeer.tsv"
else:
  covariateFile = f"{inRootDir}/covariates.tsv"


# Tools
limixSimg = f"{projDir}/scripts/tools/limix20220510.simg"
limixPath = f"{projDir}/scripts/tools/limix_qtl"
postProcess = f"{limixPath}/Limix_QTL/post_processing/minimal_postprocess.py"
multTestCorrect = f"{limixPath}/Limix_QTL/post_processing/multTestCorrect_PipeLine.R"

intermParam = ""
if eval_model == "interaction":
    qtlAnalysis = f"{limixPath}/Limix_QTL/run_interaction_QTL_analysis.py"
    intermParam = "-it " + config["inter_term"]
elif eval_model == "hurdle":
    qtlAnalysis = f"{limixPath}/Limix_QTL/run_QTL_analysis_hurdle.py"
else:
    qtlAnalysis = f"{limixPath}/Limix_QTL/run_QTL_analysis.py"


# Parameters
# numberOfPermutations = '10000'
# minorAlleleFrequency = '0.05'
numberOfPermutations = '1000'
minorAlleleFrequency = '0.1'
relatednessScore = '0.8'
windowSize = '1000000'
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
    input: qtlOutput, finalQTLRun, topQTL, correctedQTL

rule run_qtl_mapping:
    input:
        af = annotationFile, pf = phenotypeFile, cf = covariateFile, rf = kinshipFile, smf = sampleMappingFile
    output:
        # outputFolder + '{chunk}.finished'
        outputFolder.rstrip("/") + "/{chunk}.finished"
    params:
        cmd = qtlAnalysis,
        gen = genotypeFile,
        od = outputFolder,
        np = numberOfPermutations,
        maf = minorAlleleFrequency,
        rs = relatednessScore,
        w = windowSize,
        interms = intermParam
    resources:
      time = dynamic_slurm_time_limit
    run:
        chunkFull = extendChunk({wildcards.chunk})
        shell(
            " singularity exec {limixSimg} python {params.cmd}"
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
    input: IF = outputFolder, OF = outputFolder, finalFiles = qtlOutput
    output: finalQTLRun
    params: cmd = postProcess
    shell: "singularity exec {limixSimg} python {params.cmd} -id {input.IF} -od {input.OF} -sfo"

rule top_feature:
    input: IF = outputFolder, OF = outputFolder, finalFiles = qtlOutput
    output: topQTL
    params: cmd = postProcess
    shell: "singularity exec {limixSimg} python {params.cmd} -id {input.IF} -od {input.OF} -tfb -sfo"

rule multi_test_correction:
    input: IF = outputFolder, OF = outputFolder, finalFile = topQTL
    output: correctedQTL
    params: cmd = multTestCorrect, fdr = FDR
    shell: "singularity exec {limixSimg} Rscript {params.cmd} {input.IF} {input.OF} {params.fdr}"
