"""
Snakefile for  genotyping

Author: Marc Jan Bonder
Affiliation: EMBL-EBI
Date: Tuesday 10 February 2019
#Run: snakemake -s snakemake/limix_scqtl_snakemake_totmean_cellcount_kinship.smk --jobs 1000 --latency-wait 30 --cluster-config cluster.json --cluster 'bsub -q {cluster.queue} -n {cluster.n} -R "rusage[mem={cluster.memory}]" -M {cluster.memory} -o /hps/nobackup/hipsci/scratch/ComparingQtlMapping/SingleCell/PseudoBulk/totmean/Run_Output_PCA20_88_log_TPM_scater_libsize_206_cellcount_kinship_chr2/mapping.o -e /hps/nobackup/hipsci/scratch/ComparingQtlMapping/SingleCell/PseudoBulk/totmean/Run_Output_PCA20_88_log_TPM_scater_libsize_206_cellcount_kinship_chr2/mapping.e' --keep-going

"""

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
chunkFile = '/hps/nobackup/hipsci/scratch/trans_eqtls/QTL_Mapping/splitFiles/Ensembl_75_Limix_Annotation_FC_Gene_step100.txt'
genotypeFile = '/hps/nobackup/hipsci/scratch/genotypes/imputed/REL-2018-01/Full_Plink/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.20180102.genotypes.norm.renamed.recode.vcf.gz' #Can be .bed or .bgen, no extenion needed
annotationFile = '/hps/nobackup/hipsci/scratch/ComparingQtlMapping/Annotation/Ensembl_75_Limix_Annotation_FC_Gene.txt'
phenotypeFile = '/hps/nobackup/hipsci/scratch/ComparingQtlMapping/SingleCell/PseudoBulk/totmean/Run_Output_PCA20_88_log_TPM_scater_libsize_206_AllTissues_MJ/pseudobulk_platebasedmean_endoProjectMerged.exp.tsv.gz'
covariateFile = '/hps/nobackup/hipsci/scratch/ComparingQtlMapping/SingleCell/PseudoBulk/totmean/Run_Output_PCA20_88_log_TPM_scater_libsize_206_AllTissues_MJ/endoProjectMerged.PC20.tsv'
kinshipFile = '/hps/nobackup/hipsci/scratch/genotypes/imputed/REL-2018-01/Full_Filtered_Plink-f/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.20170327.genotypes.norm.renamed.recode.filtered.rel'
sampleMappingFile = '/hps/nobackup/hipsci/scratch/ComparingQtlMapping/SingleCell/PseudoBulk/totmean/Run_Output_PCA20_88_log_TPM_scater_libsize_206_AllTissues_MJ/GTE_EndoDiff.txt' #Not needed
numberOfPermutations = '1000'
minorAlleleFrequency = '0.1'
windowSize = '100000'
FDR = '0.05'
outputFolder = '/hps/nobackup/hipsci/scratch/ComparingQtlMapping/SingleCell/PseudoBulk/totmean/Run_Output_PCA20_88_log_TPM_scater_libsize_206_AllTissues_MJ/'

finalQTLRun = '/hps/nobackup/hipsci/scratch/ComparingQtlMapping/SingleCell/PseudoBulk/totmean/Run_Output_PCA20_88_log_TPM_scater_libsize_206_AllTissues_MJ/qtl_results_all.txt'
topQTL = '/hps/nobackup/hipsci/scratch/ComparingQtlMapping/SingleCell/PseudoBulk/totmean/Run_Output_PCA20_88_log_TPM_scater_libsize_206_AllTissues_MJ/top_qtl_results_all.txt'
correctedQTL = "".join(['/hps/nobackup/hipsci/scratch/ComparingQtlMapping/SingleCell/PseudoBulk/totmean/Run_Output_PCA20_88_log_TPM_scater_libsize_206_AllTissues_MJ/top_qtl_results_all_FDR',FDR,'.txt'])

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
            " singularity exec /hps/nobackup/stegle/users/galvari/containers/limix206_qtl.simg python /hps/nobackup/stegle/users/galvari/src/hipsci_pipeline/limix_QTL_pipeline/run_QTL_analysis.py "
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
            " singularity exec /hps/nobackup/stegle/users/galvari/containers/limix206_qtl.simg python /hps/nobackup/stegle/users/galvari/src/hipsci_pipeline/post-processing_QTL/minimal_postprocess.py "
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
            " singularity exec /hps/nobackup/stegle/users/galvari/containers/limix206_qtl.simg python /hps/nobackup/stegle/users/galvari/src/hipsci_pipeline/post-processing_QTL/minimal_postprocess.py "
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
            "singularity exec /hps/nobackup/stegle/users/galvari/containers/limix206_qtl.simg Rscript /hps/nobackup/stegle/users/galvari/src/hipsci_pipeline/post-processing_QTL/multTestCorrect_PipeLine.R "
            "{input.IF} "
            "{input.OF} "
            "{params.fdr} ")