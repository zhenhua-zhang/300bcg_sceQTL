#!/usr/bin/env python3
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: 2022 May 23
# Updated: 2022 May 23

from argparse import ArgumentParser

import pysam as ps
import pandas as pd

def get_opts():
    par = ArgumentParser()
    par.add_argument("-x", "--pg-pairs", nargs="*", required=True,
                     help="Phenotype-genotype pairs. E.g., rs115613985:GZMK."
                     " Multiple pairs should be splitted by space. Required")
    par.add_argument("-g", "--genotype-file", required=True,
                     help="File includes variants in VCF format. Required")
    par.add_argument("-p", "--phenotype-file", required=True,
                     help="The file includes phenotypes for each sample. The"
                     " program supposes the rows are phenotypes and columns are"
                     " sample id, the first column are name or ids of"
                     " phenotypes with name 'feature_id'. Required")
    par.add_argument("-c", "--covariate-file", default=None,
                     help="The file include covariates. The file is supposed to"
                     " have a header with sample id ('sample_id') as the first"
                     " column. Each column contains one covariate while each"
                     " row contains covariates per sample. Default: None")
    par.add_argument("-m", "--sample-map-file", default=None,
                     help="The map between genotype and phenotype IDs. The file"
                     " is supposed to be a non-header two-column text file. The"
                     " first column contains genotype ID which must exist in"
                     " the genotype file and the second column contains sample"
                     " ID which must exist in the phenotype file. If no file"
                     " given (default), the map will be inferred from phenotype"
                     " file. Default: None")
    par.add_argument("-o", "--out-dir", default="./",
                     help="Output folder. Default: %(default)s")

    return par.parse_args()


def fetch_var_by_rsid(gntp_file: str, snp_ids: tuple):
    var_info = {}

    vcffile = ps.VariantFile(gntp_file)
    for line in vcffile.fetch():
        per_snp_id = line.id
        if per_snp_id not in snp_ids:
            continue
        else:
            snp_order = line.alleles
            for smp_idx, per_smp in line.samples.items():
                gntp = "".join(sorted(per_smp.alleles, key=lambda x: snp_order.index(x)))

                if per_snp_id in var_info:
                    var_info[per_snp_id].update({smp_idx: gntp})
                else:
                    var_info[per_snp_id] = {smp_idx: gntp}

        if len(var_info) == len(snp_ids):
            break

    if not vcffile.is_closed:
        vcffile.close()

    return var_info


def main():
    opt = get_opts()
    pg_pairs = opt.pg_pairs
    gntp_file = opt.genotype_file
    phtp_file = opt.phenotype_file
    cvrt_file = opt.covariate_file
    smap_file = opt.sample_map_file
    out_dir = opt.out_dir

    pg_pairs = [pp.split(":") for pp in pg_pairs]
    snp_ids = tuple([pp[0] for pp in pg_pairs])
    tar_features = list(set([pp[1] for pp in pg_pairs]))

    cvrt_mat = pd.read_table(cvrt_file, header=0, index_col="sample_id")
    phtp_mat = pd.read_table(phtp_file, header=0, index_col="feature_id")
    var_info = fetch_var_by_rsid(gntp_file, snp_ids)
    if smap_file:
        kept_samples = cvrt_mat.index.intersection(phtp_mat.columns)
        map_mat = pd.read_table(smap_file, header=None, index_col=False)
        smap = map_mat.set_index(1).filter(kept_samples, axis=0).loc[:, 0]
    else:
        smap = phtp_mat.columns.to_series()
    gntp_mat = pd.DataFrame(var_info).loc[smap, :].set_index(smap.index)
    phtp_mat = phtp_mat.loc[tar_features, :].T
    info_mat = pd.concat([phtp_mat, gntp_mat, cvrt_mat], axis=1)

    for per_snp, per_feature in pg_pairs:
        tar_cols = [per_snp, per_feature] + cvrt_mat.columns.to_list()
        save_to = f"{out_dir}/{per_feature}-{per_snp}.QTL_info.csv"
        info_mat.loc[:, tar_cols].to_csv(save_to)


if __name__ == "__main__":
    main()
