#!/bin/bash
#SBATCH --mem 2G
#SBATCH --time 5:0:0
#SBATCH --qos normal
#SBATCH --partition cpu
#SBATCH --cpus-per-task 8
#SBATCH --output /home/zzhang/Documents/projects/wp_bcg_eqtl/inputs/reference/genotypes/GRCh38/%j-%u-prepare_reference.log
#SBATCH --job-name prepare-reference-genotyeps
#SBATCH --mail-user zzhang

#
## Prepare reference, including variants
#

# NOTE:
# 1. The reference panel were built on common SNPs (MAF >= 0.05) identified in 1000 Genomes Projects (release 20190312 biallelic SNV and INDEL)
# 2. The variants were annotated (i.e., rs-id assignment) against to common variants (MAF >= 0.05, GRCh38p7 build 151) from NCBI.
# 3. Only variants on autosome were included.
# 4. The script download all required files by wget and check MD5 signature.
# 5. The sample information of 1000 Genomes projects should be downloaded manually and placed in the same path of the current script file.
# 6. The script was designed for a job submission by slurm but can also be run directly in a proper shell, e.g., BASH/SH
# 7. For more information please contact: Zhenhua Zhang <zhenhua.zhang217@gmail.com>

n_cpus=8
subpop=EUR
proj_dir=~/Documents/projects/wp_bcg_eqtl
wk_dir=$proj_dir/inputs/reference/genotypes/GRCh38

# Download common reference genetic variants
# If we do not have the reference genetic variants we download them from NCBI's FTP site.
base_url=https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF
if [[ ! -e $wk_dir/00-common_all.vcf.gz ]]; then
  wget -qcP $wk_dir $base_url/00-common_all.vcf.gz
  wget -qcP $wk_dir $base_url/00-common_all.vcf.gz.tbi
  wget -qcP $wk_dir $base_url/00-common_all.vcf.gz.md5

  md5_local=$(md5sum $wk_dir/00-common_all.vcf.gz | cut -f1 -d\ )
  md5_remote=$(cut -f1 -d\  $wk_dir/00-common_all.vcf.gz.md5)
  if [[ $md5_local != $md5_remote ]]; then echo -e "[E]: Local and remote are not the same!"; fi
fi

# Download sample information
echo -e "[W]: Go to https://www.internationalgenome.org/data-portal/sample to download the sample information"
awk -v pop=$subpop -F$'\t' '$6==pop{print $1}' $wk_dir/igsr_samples.tsv >| $wk_dir/igsr_${subpop}_samples.tsv

# Download 1000 genome genotypes
base_url=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL

md5_summary_file=$wk_dir/20190312_biallelic_SNV_and_INDEL_MANIFEST.txt
if [[ ! -e $md5_summary_file ]]; then wget -qcP $wk_dir $base_url/20190312_biallelic_SNV_and_INDEL_MANIFEST.txt; fi

for chr_nr in {1..22}; do
  in_pcvcf=ALL.chr$chr_nr.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz

  echo -e "[I]: Downloading genotypes for chromosome "$chr_nr
  if [[ ! -e $wk_dir/$in_pcvcf ]]; then wget -qcP $wk_dir $base_url/$in_pcvcf; fi
  md5_local=$(md5sum $wk_dir/$in_pcvcf | cut -f1 -d\ )
  md5_remote=$(grep -E $in_pcvcf'[^.]' $md5_summary_file | cut -f3 -d$'\t')
  [[ $md5_local != $md5_remote ]] && echo -e "[E]: md5sum does not match: $in_pcvcf, $md5_local (local) vs $md5_remote (remote)"

  if [[ ! -e $wk_dir/$in_pcvcf.tbi ]]; then wget -qcP $wk_dir $base_url/$in_pcvcf.tbi; fi
  md5_local=$(md5sum $wk_dir/$in_pcvcf.tbi | cut -f1 -d\ )
  md5_remote=$(grep -E $in_pcvcf.tbi $md5_summary_file | cut -f3 -d$'\t')
  [[ $md5_local != $md5_remote ]] && echo -e "[E]: md5sum does not match: $in_pcvcf, $md5_local (local) vs $md5_remote (remote)"

  # Only keep SNPs with MAF >= 0.05
  echo -e "[I]: Filtering SNPs for chromosome "$chr_nr
  out_pcvcf=${in_pcvcf/ALL/$subpop}
  if [[ ! -e $out_pcvcf ]]; then
    ~/tools/bin/bcftools view -v snps -q 0.05:minor -S $wk_dir/igsr_${subpop}_samples.tsv --threads $n_cpus --force-samples -o $wk_dir/$out_pcvcf $wk_dir/$in_pcvcf
    ~/tools/bin/bcftools index -f -t --threads $n_cpus $wk_dir/$out_pcvcf
  fi

  # Annotate the SNPs using the common SNPs from NCBI
  echo -e "[I]: Annotating SNPs for chromosome "$chr_nr
  in_pcvcf=$out_pcvcf
  out_pcvcf=${in_pcvcf/.vcf.gz/.anno.vcf.gz}
  if [[ ! -e $out_pcvcf ]]; then
    ~/tools/bin/bcftools annotate -c ID -a $wk_dir/00-common_all.vcf.gz --threads $n_cpus $wk_dir/$in_pcvcf | bcftools view -k -o $wk_dir/$out_pcvcf
    ~/tools/bin/bcftools index -f -t --threads $n_cpus $wk_dir/$out_pcvcf
  fi
done

# Concate all processed genotypes.
echo -e "[I]: Concating SNPs for chromosome: $(seq 1 22)"
~/tools/bin/bcftools concat -a -D -o $wk_dir/$subpop.shapeit2_integrated_snv_v2a_27022019.GRCh38.phased.anno.vcf.gz $subpop.chr*.anno.vcf.gz

# Covnert the vcf file into bed/bim/fam file
echo -e "[I]: Converting VCF into bim/fam/bed format ..."
~/tools/bin/plink --vcf $wk_dir/$subpop.shapeit2_integrated_snv_v2a_27022019.GRCh38.phased.anno.vcf.gz --make-bed --out $wk_dir/$subpop

echo -e "Job done at "$(date)
echo -e To clean up temporary files do:"\n"rm *.vcf.gz *.vcf.gz.tbi *.md5 *.txt igsr_${subpop}_samples.tsv 20190312_biallelic_SNV_and_INDEL_MANIFEST.txt
