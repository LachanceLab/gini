# 2 - download_clump_panUKB.sh

# Bash script that extracts panUKB GWAS summary statistics for each trait, uses
# PLINK to LD prune for independent GWAS hits, and saves the file using the
# Rscript in filter_GWAS_independent.R.

### Directories and paths ###
pwd -P
# sets directory to plink 1.9
dir_plink="~/plink1_9"
# sets path to 1000 Genomes binary files (or another prefered LD panel, preferably
# composed of just Europeans). Leave as just the prefix before the chromosome
# number and the file extensions
dir_1kG='/storage/coda1/p-jlachance6/0/shared/1kG_phase3_GRCh37/'
loc_bfile=${dir_1kG}'1kG_phase3_GRCh37_ALL'
# sets path to text file containing FID and IID of individuals to use in LD panel
loc_LD_IIDs='/storage/coda1/p-jlachance6/0/shared/1kG_phase3_GRCh37/EUR_IIDs.txt'

# sets directory where panUKB GWAS summary statistics will be saved
dir_generated_data='../generated_data/'
dir_sf=${dir_generated_data}'panUKB_sf/'
mkdir -p $dir_sf
# sets path to the traits_list.txt file
loc_traits_list='../input_data/traits_list.txt'


module load r/4.2.1-tidy

# Reads traits_list and extracts phenotype name and AWS summary file link
IFS=$'\n'
phenotypes=()
while read -r line; do phenotypes+=("$line"); done <<< "$(awk -F'\t' '(NR>1 && $37 == "TRUE") {print $1}' $loc_traits_list)"
aws_links=()
while read -r line; do aws_links+=("$line"); done <<< "$(awk -F'\t' '(NR>1 && $37 == "TRUE") {print $31}' $loc_traits_list)"

# Loops through each trait
length=${#phenotypes[@]}
for i in $(seq 0 $((length-1)))
do
phenotype=${phenotypes[$i]}
echo ========== [${i} ${phenotype}] ==========

# downloads large panUKB GWAS summary file (zipped, about 2.5GB)
echo [${i} ${phenotype}] Downloading panUKB summary file
aws_link=${aws_links[$i]}
wget -qO ${dir_sf}${phenotype}_raw_sf.tsv.bgz ${aws_link}

# unzips panUKB GWAS summary file
echo [${i} ${phenotype}] Unzipping summary file
zcat ${dir_sf}${phenotype}_raw_sf.tsv.bgz > ${dir_sf}${phenotype}_raw_sf.tsv
rm -v ${dir_sf}${phenotype}_raw_sf.tsv.bgz

# transforms summary file to allow PLINK to properly interact with it
echo [${i} ${phenotype}] Transforming summary file for PLINK
index=$(awk -F'\t' '{for (i=1; i<=NF; i++) if ($i == "neglog10_pval_meta") print i; exit}' ${dir_sf}${phenotype}_raw_sf.tsv)
# if index is empty (no meta GWAS was done), it instead looks for pval_EUR
if [ -z "$index" ]; then index=$(awk -F'\t' '{for (i=1; i<=NF; i++) if ($i == "neglog10_pval_EUR") print i; exit}' ${dir_sf}${phenotype}_raw_sf.tsv); fi

awk -F'\t' -v col="$index" '{if (NR == 1) {print "SNP\tpval\t"$0} else if ($col != "NA"){print $1":"$2"_"$3"_"$4"\t"10^-$col"\t"$0}}' ${dir_sf}${phenotype}_raw_sf.tsv > ${dir_sf}${phenotype}_sf.tsv

# runs --clump in PLINK to LD prune GWAS summary results
echo [${i} ${phenotype}] Running --clump in PLINK

pval_cutoff=0.0001
r2_cutoff=0.2

~/plink1_9/plink \
--bfile ${loc_bfile} \
--allow-extra-chr \
--keep ${loc_LD_IIDs} \
--clump ${dir_sf}${phenotype}_sf.tsv \
--clump-snp-field SNP \
--clump-field pval \
--clump-p1 ${pval_cutoff} \
--clump-r2 ${r2_cutoff} \
--out ${dir_sf}${phenotype}_sf

# runs Rscript to filter full summary statistics to just independent SNPs
echo [${i} ${phenotype}] Filtering GWAS summary file for independent SNPs
Rscript helper_functions/filter_GWAS_independent.R ${phenotype}

# removes all trait-specific files except for filtered GWAS summary statistics
rm -v ${dir_sf}${phenotype}_raw_sf.*
rm -v ${dir_sf}${phenotype}_sf.*

done

echo DONE WITH ALL TRAITS


####
# Make list of all unique SNPs in the above summary files
loc_out_SNPs=${dir_generated_data}"all_sf_SNPs.txt"

# Loop through each file in the directory
for file in "${dir_sf}"*.txt; do
  echo ${file}
  # Extract the chromosome and position columns and append to output file
  awk 'BEGIN {FS="\t"; OFS=" "} NR>1 {print $3, $4, $4, 1}' "${file}" >> "${loc_out_SNPs}"
done

# Remove duplicate rows from the output file
sort "${loc_out_SNPs}" | uniq > "${loc_out_SNPs}.tmp"
mv "${loc_out_SNPs}.tmp" "${loc_out_SNPs}"

####
# Groups 1kG samples by continent
loc_out_IDs=${dir_1kG}'1kG_continent_IDs.txt'
    
####
# Gets continent-specific allele frequencies within 1kG for significant SNPs
~/plink1_9/plink \
--bfile ${loc_bfile} \
--allow-extra-chr \
--extract range ${loc_out_SNPs} \
--keep ${loc_out_IDs} \
--within ${loc_out_IDs} \
--freq \
--out ${dir_generated_data}AFs_1kG_sfSNPs


