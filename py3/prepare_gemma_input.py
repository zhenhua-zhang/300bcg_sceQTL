#!/usr/bin/env python3
# Because phenotype is required, so the order of the samples in the rest files is the same as phenotype's

# TODO: 1-based or 0-based

import os
from argparse import ArgumentParser

import pysam as ps
import pandas as pd


class GFFTree:
    def __init__(self, in_path):
        self._in_path = in_path
        self._gff_records = []

    def fetch(self, chrom=None, begin=None, end=None, bin_size=1024):
        gff_file = ps.TabixFile(self._in_path, parser=ps.asGFF3())
        subgfr = gff_file.fetch(chrom, begin, end)

        record_tmp = []
        for idx, per_rec in enumerate(subgfr):
            attr_dict = dict([x.split("=") for x in per_rec.attributes.split(";")])
            record_tmp.append((per_rec.contig, per_rec.start, per_rec.end, attr_dict["gene_name"]))

            if (idx + 1) % bin_size == 0:
                self._gff_records.append(record_tmp)
                record_tmp = []
        if record_tmp: self._gff_records.append(record_tmp)
        gff_file.close()

    @property
    def records(self):
        return self._gff_records


class Variant:
    def __init__(self, in_path):
        self._in_path = in_path
        self._var_records = []
        self._var_annos = []

    def fetch(self, chrom, begin, end, samples=None):
        var_file = ps.VariantFile(self._in_path)

        subvars = var_file.fetch(chrom, begin, end)
        for idx, prec in enumerate(subvars):

            # Dealing with minor and major alleles
            ref_allele = prec.ref
            alt_allele, *_ = prec.alts if prec.alts else (None, None)
            alt_freq, *_ = prec.info.get("AF")

            if alt_freq and alt_freq > 0.5:
                major, minor = alt_allele, ref_allele
            else:
                major, minor = ref_allele, alt_allele

            if samples:
                smpds = [prec.samples.get(ps, {"DS": None}).get("DS") for ps in samples]
            else:
                smpds = [pfmt.get("DS", None) for _, pfmt in prec.samples.items()]

            if idx == 0:
                if samples:
                    colname = ["Chrom", "Minor", "Major"] + samples
                else:
                    colname = ["Chrom", "Minor", "Major"] + [str(x) for x in prec.samples.keys()]

                self._var_records.append(colname)

            per_record = [prec.id, minor, major] + smpds
            self._var_records.append(per_record)
            self._var_annos.append([prec.id, prec.pos, prec.chrom])

        var_file.close()

    def write_to(self, out_path):
        var_path = f"{out_path}_genotype.txt"
        with open(var_path, "w") as outhand:
            for per_rec in self._var_records:
                line = "\t".join([str(x) for x in per_rec]).strip("\n") + "\n"
                outhand.write(line)

        ann_path = f"{out_path}_genotype_annotation.txt"
        with open(ann_path, "w") as outhand:
            for per_rec in self._var_annos:
                line = ",".join([str(x) for x in per_rec]).strip("\n") + "\n"
                outhand.write(line)

    @property
    def records(self):
        return self._var_records


def fread(fpath, n_test=20, **kwargs):
    delims = {"\t": [], ",": [], ";": [], "|": [], " ": []}
    with open(fpath) as inhand:
        for idx, line in enumerate(inhand):
            if idx == n_test: break
            delims["\t"].append(line.count("\t"))
            delims[","].append(line.count(","))
            delims[";"].append(line.count(";"))
            delims["|"].append(line.count("|"))
            delims[" "].append(line.count(" "))

        obs = [(d, max(c)) for d, c in delims.items() if min(c) == max(c) and max(c) != 0]

    if len(obs) != 1: raise ValueError("Ambugious field delimiter, you have to specify it manully.")
    if "sep" not in kwargs: kwargs["sep"] = obs[0][0]

    return pd.read_csv(fpath, **kwargs)


def main(opts):
    out_dir = opts.out_dir
    if not os.path.exists(out_dir): os.makedirs(out_dir, exist_ok=True)

    # Expression matrix. Gene per row, sample per column.
    # The first column should be the gene name or any unique identifier indicating the gene.
    expr_tab = fread(opts.phenotype, index_col="feature_id")
    avail_genes = expr_tab.index.to_list()
    avail_samples = expr_tab.columns.to_list()

    sample_orders = f"{out_dir}/sample_orders.txt"
    with open(sample_orders, "w") as outhand:
        outhand.write("\n".join(avail_samples) + "\n")

    # Covariates. Sample/individual per row, trait per column.
    # The first column should be the sample id or any unique identifier indicating the sample/individual.
    if opts.covariate:
        covar_tab = fread(opts.covariate, index_col="sample_id")
        covar_save_to = f"{out_dir}/covariate.txt"
        (covar_tab
         .assign(intercept=1)
         .loc[avail_samples, ["intercept"] + covar_tab.columns.to_list()]
         .to_csv(covar_save_to, sep="\t", index=False, header=False))

        covar_order = f"{out_dir}/covariate_orders.txt"
        with open(covar_order, "w") as outhand:
            outhand.write("\n".join(["intercept"] + covar_tab.columns.to_list()) + "\n")

        if opts.ia_var:
            itvar_save_to = f"{out_dir}/interaction_covariates-{opts.ia_var}.txt"
            covar_tab.loc[avail_samples, opts.ia_var].to_csv(itvar_save_to, index=False, header=False)

    # Phenotype to genotype map.
    # Two columns, first is the genotype id in the VCF file, the second is the sample/individual id in phenotype and covariate files.
    if opts.p2g_map:
        p2g_tab = fread(opts.p2g_map, header=None, index_col=1)
        p2g_map = p2g_tab.to_dict()[0]
        tar_geno_ids = [p2g_map[pid] for pid in avail_samples]
    else:
        p2g_map = None
        tar_geno_ids = None

    # Genomic feature tree including gene coordiantions.
    # The features were assigned bins to ensure equal size of chuncks regarding number of genes.
    all_gfr = GFFTree(opts.genomic_feature)
    all_gfr.fetch()

    for idx, per_block in enumerate(all_gfr.records):
        dir_name = f"{out_dir}/block_{idx:>05}"
        if not os.path.exists(dir_name): os.makedirs(dir_name, exist_ok=True)

        for chrom, start, stop, gene_name in per_block:
            if gene_name not in avail_genes: continue

            chrom = f"chr{chrom}"
            start = max(start - opts.flank_size, 0)
            stop = stop + opts.flank_size

            all_vars = Variant(opts.genotype)
            all_vars.fetch(chrom, start, stop, tar_geno_ids)

            genotype_save_to = f"{dir_name}/{gene_name}"
            all_vars.write_to(genotype_save_to)

            phenotype_save_to = f"{dir_name}/{gene_name}_phenotype.txt"
            expr_tab.loc[gene_name, avail_samples].to_csv(phenotype_save_to, sep="\t", index=False, header=False)


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
                     help="Phenotype to genotype map. If it is missing the script supposes the names are matched across all files."
                          " Default: %(default)s")
    par.add_argument("-I", "--ia-var", default=None, metavar="STR",
                     help="The covariate used for the interaction. Default: %(default)s")
    par.add_argument("-F", "--flank-size", default=2.5e5, metavar="INT",
                     help="The window-size used to fetch SNPs for each feature. Default: %(default)s")
    par.add_argument("-o", "--out-dir", default="./gemma_inputs", metavar="STR",
                     help="The output dir. Default: %(default)s")

    opts = par.parse_args()

    main(opts)
