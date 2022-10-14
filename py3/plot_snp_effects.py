#!/usr/bin/env python3
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: 2022 May 23
# Updated: 2022 May 23

from argparse import ArgumentParser
from collections import OrderedDict

import pysam as ps
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def get_opts():
    par = ArgumentParser()
    par.add_argument("-x", "--pg-pairs", nargs="*", required=True,
                     help="Phenotype-genotype pairs. E.g., IRF7:rs7113204."
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
    par.add_argument("-G", "--group-by", default=None,
                     help="Covariates that are used to group the feature level."
                     " It is useful when the condition plays potential role. "
                     " Default: None")
    par.add_argument("-m", "--sample-map-file", default=None,
                     help="The map between genotype and phenotype IDs. The file"
                     " is supposed to be a non-header two-column text file. The"
                     " first column contains genotype ID which must exist in"
                     " the genotype file and the second column contains sample"
                     " ID which must exist in the phenotype file. If no file"
                     " given (default), the map will be inferred from phenotype"
                     " file. Default: None")
    par.add_argument("--cond-col", default=None,
                     help="Condition column. Default: %(default)s")
    par.add_argument("--cond-map", nargs="*", default=None,
                     help="Rename the condition to given name and sort it in"
                     " the presented order. The syntax is OLD-NAME:NEW-NAME,"
                     " where multiple ones should be splitted by space.")
    par.add_argument("--drop-outlier", action="store_true",
                     help="Whether drop outliers.")
    par.add_argument("--fig-size", nargs=2, default=[7, 7], type=float,
                     help="The figure size in format w h. Default: %(default)s")
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
                gntp = "".join(
                    sorted(per_smp.alleles, key=lambda x: snp_order.index(x)))

                if per_snp_id in var_info:
                    var_info[per_snp_id].update({smp_idx: gntp})
                else:
                    var_info[per_snp_id] = {smp_idx: gntp}

        if len(var_info) == len(snp_ids):
            break

    if not vcffile.is_closed:
        vcffile.close()

    return var_info


def draw_snpeff(
    mat, phtp_col, gntp_col, cond_col, cond_order=None, drop_outlier=True,
    save_to="scatterplot.png", fig_size=(7, 7)
):
    hue_order = mat[gntp_col].drop_duplicates().sort_values().to_list()
    mat.sort_values(gntp_col, inplace=True)

    do_dodge = True
    if cond_col is None:
        do_dodge = False
        cond_col = gntp_col

    fig, axe = plt.subplots()
    axe = sns.boxplot(
        x=cond_col, y=phtp_col, order=cond_order, hue=gntp_col,
        hue_order=hue_order, data=mat, ax=axe, dodge=do_dodge,
        showfliers=(not drop_outlier)
    )
    y_lims = axe.get_ylim()

    axe = sns.stripplot(
        x=cond_col, y=phtp_col, order=cond_order, hue=gntp_col,
        hue_order=hue_order, data=mat, ax=axe, dodge=do_dodge,
        linewidth=1, edgecolor="gray"
    )

    if drop_outlier:
        axe.set_ylim(y_lims[0] * 0.9, y_lims[1] * 1.1)

    hdl, lbl = axe.get_legend_handles_labels()
    plt.legend(hdl[0:len(hue_order)], lbl[0:len(hue_order)], title=gntp_col)

    axe.set_ylabel("Level of " + phtp_col)
    axe.set_xlabel("Conditions")
    axe.spines["top"].set_visible(False)
    axe.spines["right"].set_visible(False)

    fig.set_figwidth(fig_size[0])
    fig.set_figheight(fig_size[1])
    fig.set_tight_layout(True)

    fig.savefig(save_to)


def main():
    opt = get_opts()
    pg_pairs = opt.pg_pairs
    gntp_file = opt.genotype_file
    phtp_file = opt.phenotype_file
    cvrt_file = opt.covariate_file
    smap_file = opt.sample_map_file
    drop_outlier = opt.drop_outlier
    cond_col = opt.cond_col
    cond_map = opt.cond_map
    out_dir = opt.out_dir
    fig_size = opt.fig_size

    pg_pairs = [pp.split(":") for pp in pg_pairs]
    snp_ids = tuple([pp[1] for pp in pg_pairs])
    tar_features = [pp[0] for pp in pg_pairs]

    if cond_map:
        cond_map = OrderedDict([pp.split(":") for pp in cond_map])

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

    cond_order = None
    if cond_map and cond_col:
        cond_order = list(cond_map.values())
        info_mat[cond_col] = (info_mat[cond_col]
                              .apply(lambda x, cm: cm[str(x)], cm=cond_map))

    for per_feature, per_snp in pg_pairs:
        save_to = f"{out_dir}/{per_feature}-{per_snp}.jpg"
        draw_snpeff(
            info_mat, per_feature, per_snp, cond_col, cond_order, drop_outlier,
            save_to, fig_size
        )


if __name__ == "__main__":
    main()
