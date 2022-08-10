#!/usr/bin/env bash
# A script to download publicly available GWAS summary statistics (SS) from IEU OpenGWAS project.
# The SS is in a variant of VCF format, GWAS VCF, for more information, please refer https://github.com/MRCIEU/gwas-vcf-specification


# Selected publicly available GWAS summary statistics which are deposited by IEU OpenGWAS project.
declare -A gwas_dict=(
  # Trait, year, abbrev., sample_size
  # Autoimmune disease
  # ["ebi-a-GCST005536"]="Type1Diabetes 2015 T1D"
  ["ieu-a-31"]="InflammatoryBowelDisease 2015 IBD"
  ["ieu-a-30"]="CrohnsDisease 2015 CrD"
  ["ieu-a-32"]="UlcerativeColitis 2015 UC"
  ["ieu-a-996"]="AtopicDermatitis 2014 AD"
  ["ukb-a-100"]="Psoriasis 2017 Ps"
  ["ukb-b-17670"]="MultipleSclerosis 2019 MS"
  ["ukb-d-M13_RHEUMA"]="RheumatoidArthritis 2018 RA"

  # Immune disorders
  ["ukb-d-J10_ASTHMA"]="Asthma 2018 Asthma"

  # Infectious disease
  ["ebi-a-GCST010780"]="COVID19Release4 2020 COVID"

  # Cancers
  ["ieu-a-966"]="LungCancer 2014 LC"
  ["ieu-b-4809"]="ProstateCancer 2021 PC"
  ["ieu-a-1082"]="ThyroidCancer 2013 TC"
  ["ieu-b-4963"]="OvarianCancer 2017 OC"
  ["ieu-b-4965"]="ColorectalCancer 2021 CC"
  ["ieu-b-4874"]="BladderCancer 2021 BC"

  # Brain disorders
  ["ieu-b-5067"]="AlzheimerDiseases 2022 AlD"
  ["ieu-b-42"]="Schizophrenia 2014 Sch"

  # Other genetic-related disease
  ["ukb-a-107"]="GoutDisease 2017 GD"
  ["ebi-a-GCST006867"]="Type2Diabetes 2018 T2D"
  ["ieu-a-7"]="CoronaryHeartDisease 2013 ChD"

  # Others traits
  ["ieu-a-89"]="Height 2014 Height"
  ["ieu-b-40"]="BodyMassIndex 2018 BMI"
  ["ieu-b-109"]="HDLCholesterol 2020 HDL"
  ["ieu-b-110"]="LDLCholesterol 2020 LDL"
  ["ieu-b-111"]="Triglycerides 2020 Tri"
)

# The base URL used to download the summary statistics
base_url="https://gwas.mrcieu.ac.uk/files"

out_dir=~/Documents/projects/wp_bcg_eqtl/inputs/public_gwas
if [[ -d $out_dir ]]; then mkdir -p $out_dir; fi

for ieu_id in ${!gwas_dict[*]}; do
  trait=$(cut -f1 -d" " <<<${gwas_dict[$ieu_id]})
  year=$(cut -f2 -d" " <<<${gwas_dict[$ieu_id]})
  abbrev=$(cut -f3 -d" " <<<${gwas_dict[$ieu_id]})
  save_to=$out_dir/${year}_${trait}_${ieu_id}

  if [[ -f $save_to.vcf.gz ]]; then continue; fi
  wget -c -O $save_to.vcf.gz $base_url/$ieu_id/$ieu_id.vcf.gz
  wget -c -O $save_to.vcf.gz.tbi $base_url/$ieu_id/$ieu_id.vcf.gz.tbi
done

for x in *.gz; do
  (printf "%s: " $x; zcat $x | grep -v "^#" -m 1 | cut -f8-) | grep -v AF=
done
