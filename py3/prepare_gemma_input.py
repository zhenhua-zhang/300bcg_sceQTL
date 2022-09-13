#!/usr/bin/env python3
# Because phenotype is required, so the order of the samples in the rest files is the same as phenotype's

from argparse import ArgumentParser
from collections import OrderedDict

import pysam as ps
import pandas as pd

GFF_HEADER = ["Chrom", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes"]

def filter_record(rec):
    return all([
        "transcript_support_level=1" in rec.Attributes,
        "Ensembl_canonical" in rec.Attributes,
        "not_best_in_genome_evidence" not in rec.Attributes,
    ])


def select_attr(rec, kept_attr=["gene_name", "gene_type", "gene_id", "tag"]):
    attr_dict = dict([x.split("=") for x in rec.Attributes.split(";")])
    for pka in kept_attr:
        rec[pka] = attr_dict.get(pka, None)

    return rec.loc[["Chrom", "Start", "End"] + kept_attr]


def parse_genomic_features(path, out_pref=None):
    global GFF_HEADER
    gftab = pd.read_csv(path, sep="\t", comment="#", header=None, names=GFF_HEADER).query("Type == 'transcript'")
    gftab = gftab.loc[gftab.apply(filter_record, axis=1)].apply(select_attr, axis=1)

    return gftab


def prep_genotype(path, p2gmap=None, out_pref = "./genotype"):
    variter = ps.VariantFile(path).fetch()
    with open(out_pref.strip("/") + ".csv", mode="w") as outfh:
        for prec in variter:
            snp_id = prec.id

            # Dealing with minor and major alleles
            ref_allele = prec.ref
            alt_allele, *_ = prec.alts if prec.alts else (None, None)
            alt_freq, *_ = prec.info.get("AF")

            if alt_freq and alt_freq > 0.5:
                major, minor = alt_allele, ref_allele
            else:
                major, minor = ref_allele, alt_allele

            if p2gmap:
                smpds = [prec.samples.get(ps, {"DS": None}).get("DS") for ps in p2gmap.values()]
            else:
                smpds = [pfmt.get("DS", None) for _, pfmt in prec.samples.items()]

            # mean genotype format (for more, check GEMMA manual at Mean Genotype File).
            mgrec = ",".join([str(x) for x in [snp_id, minor, major] + smpds]).strip("\n") + "\n"

            outfh.write(mgrec)


def prep_phenotype(path, out_pref="./phenotype"):
    pheno_tab = pd.read_csv(path)


def prep_covariate(path, p2gmap=None):
    cov_tab = pd.read_csv(path)


def main(opts):
    # p2gmap = prep_phenotype(opts.phenotype, opts.out_pref)
    parse_genomic_features(opts.genomic_feature)

    if opts.p2g_map is not None:
        p2gmap = pd.read_csv(opts.p2g_map, sep = "\t")
        p2gmap = OrderedDict(zip(p2gmap.sample_id, p2gmap.genotype_id))

    if opts.genotype is not None:
        prep_genotype(opts.genotype, p2gmap, opts.out_pref)


if __name__ == "__main__":
    par = ArgumentParser(description="Prepare GEMMA input files.")

    par.add_argument("-g", "--genotype", required=True, metavar="FILE", help="Genotypes in VCF/BCF format. Required")
    par.add_argument("-p", "--phenotype", required=True, metavar="FILE",
                     help="Phenotype in tab-delimited format, with rows representing phenotypes and columns representing samples. "
                          "A columns named feature_id is required to indicate the phenotype name. Required")
    par.add_argument("-f", "--genomic-feature", required=True, metavar="FILE",
                     help="Genomic features file including target position information. Required")
    par.add_argument("-c", "--covariate", default=None, metavar="FILE",
                     help="Covariate in tab-delimited format, with columns as covariate names and orws as samples. "
                          "A columns named sample_id is required. Default: %(default)s")
    par.add_argument("-m", "--p2g-map", default=None, metavar="FILE",
                     help="Phenotype to genotype map. If it is missing the script supposes the names are matched across all files. Default: %(default)s")
    par.add_argument("-F", "--flank-size", default=2.5e5, metavar="INT",
                     help="The window-size used to fetch SNPs for each feature. Default: %(default)s")
    par.add_argument("-o", "--out-pref", default="./gemma_input", metavar="STR",
                     help="The output prefix, path can be included. Default: %(default)s")

    opts = par.parse_args()

    main(opts)
