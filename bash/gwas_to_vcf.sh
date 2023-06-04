#!/usr/bin/bash
set -Eeu -o pipefail

echo -e "Convert COVID-19 HGI summary statistic into GWAS VCF format."
#
## NOTE:
#  1. Reference genome version GRCH38

work_dir=~/Documents/projects/wp_bcg_eqtl/temps/GWAS2VCF
tool_dir=~/tools/gwas2vcf

# Install the package
# NOTE: the pysam 0.15.2 does not work under Python3.9.6, so using the latest 0.20.0
git clone https://github.com/MRCIEU/gwas2vcf.git $tool_dir

python3 -m venv $tool_dir/.env
source $tool_dir/.env/bin/activate
pip install -U pip
pip install -r $tool_dir/requirements.txt
pip install 'git+https://github.com/bioinformed/vgraph@v1.4.0#egg=vgraph'

# GWAS summary statistics, COVID-19, HGI, Very severe respiratory confirmed covid vs. population, GRCh38
wget -cP $work_dir https://storage.googleapis.com/covid19-hg-public/freeze_7/results/20220403/pop_spec/sumstats/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv.gz
gunzip $work_dir/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv.gz

# Reference genome sequence, GRCh38/hg38/b38
wget -cP $work_dir https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
wget -cP $work_dir https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai
wget -cP $work_dir https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict

# Reference genomic variants, dbSNP b153, GRCh38/hg38/b38
wget -cP $work_dir http://fileserve.mrcieu.ac.uk/dbsnp/dbsnp.v153.hg38.vcf.gz
wget -cP $work_dir http://fileserve.mrcieu.ac.uk/dbsnp/dbsnp.v153.hg38.vcf.gz.tbi

# Parameters
cat <<EOF >| $work_dir/parameters.json
{
  "chr_col": 0,
  "pos_col": 1,
  "oa_col": 2,
  "ea_col": 3,
  "snp_col": 14,
  "beta_col": 6,
  "se_col": 7,
  "ncontrol_col": 10,
  "pval_col": 8,
  "eaf_col": 13,
  "delimiter": "\\t",
  "header": true,
  "build": "GRCh38"
}
EOF


# Prepare an alias file for chromosomes


parameters=$work_dir/parameters.json
summary_statistic=$work_dir/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv
reference_genome=$work_dir/Homo_sapiens_assembly38.fasta
dbsnp_variants=$work_dir/dbsnp.v153.hg38.vcf.gz
vcf_out_file=$work_dir/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.vcf.gz
out_id="COVID19_HGI_A2_ALL_eur_leave23andme_20220403"

python $tool_dir/main.py \
  --id ${out_id} \
  --json ${parameters} \
  --ref ${reference_genome} \
  --dbsnp ${dbsnp_variants} \
  --data ${summary_statistic} \
  --out ${vcf_out_file} \
  --alias $tool_dir/alias-hg38.txt


#
## Clump GWAS summary
#
sumstat_file=/vol/projects/zzhang/projects/wp_bcg_eqtl/temps/GWAS2VCF/COVID19_HGI_A2_ALL_eur_leave23andme_20220403.tsv
refgeno_file=/vol/projects/BIIM/resources/1000G/GRCh38/EUR
plink \
  --clump $sumstat_file \
  --bfile $refgeno_file \
  --clump-p1 5e-8 \
  --clump-p2 5e-2 \
  --clump-r2 0.5 \
  --clump-kb 250 \
  --clump-snp-field rsid \
  --clump-field all_inv_var_meta_p \
  --out COVID19_HGI_A2_ALL_eur_leave23andme_20220403.independent_snp


# Clean up
# echo "Clean up..." && rm -rf $work_dir && echo "Clean up successful!" || echo "Clean up failed!"
