#!/usr/bin/env python3
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: May 21, 2021
# Updated: Jun 01, 2021

import logging
import argparse
import itertools

from collections import OrderedDict

import matplotlib.pyplot as plt
from pysam import AlignmentFile

# Plasmablast were removed as nearly no ASoC reads can be retrieved.
BLACKLIST_CELLTYPE = ["Plasmablast"]


class FieldMissingError(Exception):
    def __init__(self, msg):
        self.msg = msg
        super(FieldMissingError, self).__init__(msg)

    def __str__(self):
        return f"Missing required fields: {self.msg}!"


class HeaderMissingError(Exception):
    def __init__(self, msg=0):
        self.msg = msg
        super(HeaderMissingError, self).__init__(msg)

    def __str__(self):
        return f"Missing header line! Or you skipped it by {self.msg} lines."


class MultipleSNPIdError(Exception):
    def __init__(self, msg=0):
        self.msg = msg
        super(MultipleSNPIdError, self).__init__(msg)

    def __str__(self):
        return f"Multiple SNP ID given, which is not supported yet! {self.msg}"


class MissingSNPIdError(Exception):
    def __init__(self, msg=0):
        self.msg = msg
        super(MissingSNPIdError, self).__init__(msg)

    def __str__(self):
        return f"Multiple SNP ID given, which is not supported yet! {self.msg}"


def check_file_header(header_fields):
    """Check the header of a file.
    """
    for required_fields in ["contig", "position", "variantID", "refAllele",
                            "altAllele"]:
        if required_fields not in header_fields:
            raise FieldMissingError(required_fields)


def check_rcpool(rcpool):
    """Check whether the read counts pool is good.
    """
    snp_rsid, celltype, condition = [], [], []
    for _snp_rsid, _celltype, _condition, _, _ in rcpool.keys():
        snp_rsid.append(_snp_rsid)
        celltype.append(_celltype)
        condition.append(_condition)

    snp_rsid = list(set(snp_rsid))
    celltype = list(set(celltype))
    condition = list(set(condition))

    n_snp_rsid = len(snp_rsid)
    logging.info(f"Current SNP rsids: {', '.join(snp_rsid)}")

    if n_snp_rsid > 1:
        raise MultipleSNPIdError(", ".join(snp_rsid))
    elif n_snp_rsid < 1:
        raise MissingSNPIdError("Missing SNP rsid!")

    return sorted(celltype), sorted(condition)


def parse_region_file(fpath, sep="\t", skip=0, header=1, add_chr=True, discard_celltype=BLACKLIST_CELLTYPE, bampath="./"):
    """Parse position file."""
    snp_pos_map = {}
    if skip >= header:
        raise HeaderMissingError(skip)

    with open(fpath, "r") as fhandle:
        header_fields = None
        for idx, line in enumerate(fhandle):
            if idx < skip:
                continue

            if idx == header - 1:
                header_fields = line.strip().split(sep)
                check_file_header(header_fields)
                continue

            snp_info = line.strip().split(sep)
            chrom, pos, snpid, ref, alt, donor_id, celltype, condition, *_ = snp_info

            if celltype in discard_celltype:
                continue

            chrom = "chr" + chrom if add_chr else chrom
            bamfile = f'{bampath}/{donor_id}-{celltype}_asoc.bam'
            rec_key = (snpid, celltype, condition, chrom, ref, alt, int(pos))
            if rec_key not in snp_pos_map:
                snp_pos_map[rec_key] = [bamfile]
            else:
                snp_pos_map[rec_key].append(bamfile)

    return snp_pos_map


def fetch_read_counts(snpid, region_tab, flank_len=200, threads=4, discard_zero=True):
    readcounts = {}
    for key, htsfile_pool in region_tab.items():
        if key in readcounts:
            raise KeyError(f"Duplicated key: {key}")

        if snpid == key[0]:
            _, celltype, condition, chrom, ref, alt, pos = key
            pile_start, pile_end = pos - flank_len, pos + flank_len
            asoc_rc = OrderedDict()
            asoc_rc[ref] = {}
            asoc_rc[alt] = {}
            flank_rc = {}

            for htsfile in htsfile_pool:
                mode = "rb" if htsfile.endswith(".bam") else "r"
                with AlignmentFile(htsfile, mode, threads=threads) as samhd:
                    for per_col in samhd.pileup(chrom, pile_start, pile_end):
                        cur_pos = per_col.pos + 1
                        if cur_pos in flank_rc:
                            flank_rc[cur_pos] += per_col.nsegments
                        else:
                            flank_rc[cur_pos] = per_col.nsegments

                        if cur_pos == pos:
                            for per_seg in per_col.pileups:
                                base = per_seg.alignment.query_sequence[per_seg.query_position]
                                segpos = per_seg.alignment.get_reference_positions()
                                if base.upper() not in [ref, alt]: continue
                                for per_pos in segpos:
                                    if per_pos in asoc_rc[base]:
                                        asoc_rc[base][per_pos] += 1
                                    else:
                                        asoc_rc[base][per_pos] = 1

            rc_key = (snpid, celltype, condition, chrom, pos)
            readcounts[rc_key] = {"asoc_rc": asoc_rc, "flank_rc": flank_rc}

    return readcounts


# def plot_rc(rcpool, figheight=5, figwidth=6, cond_order=["Convalescent", "Hospitalized"]):
def plot_rc(rcpool, figheight=5, figwidth=6, cond_order=["Post", "Mild", "Severe"]):
    """Plot ASoC read counts."""
    celltype, condition = check_rcpool(rcpool)
    condition.sort(key=lambda x: cond_order.index(x))

    n_celltype, n_condition = len(celltype), len(condition)
    if n_celltype == n_condition == 1:
        logging.info("Only one panel will be plotted!")

    axes_idx = list(itertools.product(range(n_condition), range(n_celltype)))

    snpid = None
    lines, labels = None, None
    fig, axes = plt.subplots(ncols=n_celltype, nrows=n_condition, sharex=True, sharey=True, squeeze=False, subplot_kw={"frame_on": True})
    for (snpid, cur_cell, cur_cond, cur_chrom, cur_pos), val in rcpool.items():
        celltype_idx = celltype.index(cur_cell)
        condition_idx = condition.index(cur_cond)
        cur_axe = axes[condition_idx][celltype_idx]

        asoc_rc, flank_rc = val["asoc_rc"], val["flank_rc"]
        flank_reads_pos = sorted(list(flank_rc.keys()))
        flank_reads_counts = [flank_rc[x] for x in flank_reads_pos]
        flank_reads_counts_ = [-x for x in flank_reads_counts]
        cur_axe.fill_between(flank_reads_pos, flank_reads_counts_, flank_reads_counts, lw=0, alpha=0.2, color="#00BA38", label="Gross read counts")

        ref, alt = list(asoc_rc.keys())
        color_dict = {ref: "#00BFC4", alt: "#C77CFF"}
        report_dict = OrderedDict()

        report_dict["chrom"] = cur_chrom
        report_dict["pos"] = cur_pos
        report_dict["snpid"] = snpid

        asoc_reads_pos = sorted(list(set(list(asoc_rc[ref].keys()) + list(asoc_rc[alt].keys()))))
        for pp in asoc_reads_pos:
            if pp not in asoc_rc[ref]: asoc_rc[ref][pp] = 0
            if pp not in asoc_rc[alt]: asoc_rc[alt][pp] = 0

        ref_counts = [asoc_rc[ref][x] for x in asoc_reads_pos]
        alt_counts = [-asoc_rc[alt][x] for x in asoc_reads_pos]
        cur_axe.fill_between(asoc_reads_pos, 0, ref_counts, lw=0, color=color_dict[ref], label=ref + " allele counts")
        cur_axe.fill_between(asoc_reads_pos, 0, alt_counts, lw=0, color=color_dict[alt], label=alt + " allele counts")

        report_dict["cond"] = cur_cond
        report_dict["celltype"] = cur_cell
        report_dict["refcounts"] = asoc_rc[ref][cur_pos]
        report_dict["altcounts"] = asoc_rc[alt][cur_pos]
        logging.info(",".join([str(x) for x in report_dict.values()]))

        cur_axe.scatter(y=0, x=cur_pos, s=2, c="k", marker="|", label=snpid)

        if lines is None and labels is None:
            lines, labels = cur_axe.get_legend_handles_labels()

    for xi, yi in axes_idx:
        axes[xi][yi].set_xticks([])
        axes[xi][yi].tick_params(length=0.01, labelsize=5)

        cur_cond, cur_cell = condition[xi], celltype[yi]
        if yi == 0:
            axes[xi][yi].set_ylabel(cur_cond)

        if xi == len(condition) - 1:
            axes[xi][yi].set_xlabel(cur_cell)

        for spine_pos in ["top", "left", "right", "bottom"]:
            axes[xi][yi].spines[spine_pos].set_linewidth(0.0)

    fig.set_tight_layout(True)
    fig.subplots_adjust(wspace=-0.01, hspace=-0.01)
    fig.set_figheight(figheight)
    fig.set_figwidth(figwidth)
    fig.legend(lines, labels, ncol=1, fontsize="xx-small")

    return fig


def main():
    """The main entry of the script."""
    # CLI arguments
    parser = argparse.ArgumentParser(description="Plot ASoC read counts.")
    parser.add_argument("-r", "--region-file", required=True, metavar="FILE", help="A tab-delimited file including SNP positions.")
    parser.add_argument("-i", "--snp-rsids", required=True, metavar="RSID", nargs="+", help="SNP rsid to be ploted.")
    parser.add_argument("-b", "--bam-path", default="./", metavar="DIR", help="The path to BAM files. Default: %(default)s")
    parser.add_argument("-f", "--flank-len", default=200, type=int, metavar="INTEGER", help="The length of flank sequence in up/down-stream. of the ASoC SNP." " Default: %(default)s")
    parser.add_argument("-p", "--threads", default=1, type=int, metavar="INTEGER", help="Number of threads used for (de)compression. Default: %(default)s")
    parser.add_argument("-D", "--discard-zero", action="store_true", help="Discard cell type in which no reads counted.")
    parser.add_argument("-d", "--discard-celltype", nargs="*", default=None, metavar="STRING", help="Cell types will be discarded if region files contains cell type. Default: %(default)s")
    parser.add_argument("-a", "--add-chr", action="store_true", help="Add 'chr' the chromosome.")
    parser.add_argument("--seqtech", choices=["se", "pe"], default="pe", metavar="[se, pe]", help="Sequencing technology. Not in usage yet.")
    parser.add_argument("--fig-width", type=float, default=6, metavar="FLOAT", help="The width of the figure. Default: %(default)s")
    parser.add_argument("--fig-height", type=float, default=5, metavar="FLOAT", help="The height of the figure. Default: %(default)s")
    parser.add_argument("-o", "--outdir", default="./", metavar="PATH", help="The output directory. Default: %(default)s")

    options = parser.parse_args()
    region_file = options.region_file
    flank_len = options.flank_len
    bam_path = options.bam_path
    threads = options.threads
    outdir = options.outdir
    snp_rsids = options.snp_rsids
    add_chr = options.add_chr
    discard_zero = options.discard_zero
    discard_celltype = options.discard_celltype
    fig_width = options.fig_width
    fig_height = options.fig_height

    logging.basicConfig(format="{levelname: ^8}| {asctime} | {message}", style="{", datefmt="%Y%m%d,%H:%M:%S", level=logging.INFO)

    if discard_celltype: BLACKLIST_CELLTYPE.extend(discard_celltype)
    logging.info(f"Celltype in blacklist: {', '.join(list(set(BLACKLIST_CELLTYPE)))}")

    region_tab = parse_region_file(region_file, sep=",", add_chr=add_chr, bampath=bam_path)
    for snpid in snp_rsids:
        rcpool = fetch_read_counts(snpid, region_tab, flank_len, threads, discard_zero)
        rcplot = plot_rc(rcpool, figheight=fig_height, figwidth=fig_width)

        rcplot.savefig(f"{outdir}/{snpid}_asoc_read_depth.pdf")


if __name__ == "__main__":
    main()
