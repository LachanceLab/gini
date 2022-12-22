# 2 - create_ancestry_AFs.sh

# Loops through each chromosome and calculates allele frequencies, stratified by ancestry

### Directories and paths ###

# sets working directory
cd ./
# sets directory to plink
dir_plink="~/plink1_9"
# sets path to UKB binary files. Leave as just the prefix before the chromosome
# number and the file extensions
loc_bfile_prefix="/directory/UKB_QC_chr"
# sets directory of generated_data from previous scripts
dir_generated_data="../generated_data"


### Code ###

for i in {1..22};
do

${dir_plink}/plink \
--bfile ${loc_bfile_prefix}${i} \
--keep ${dir_generated_data}/pop_ALLrc_IIDs.txt \
--within ${dir_generated_data}/pop_ALLrc_IIDs.txt \
--freq \
--out ${dir_generated_data}/allele_frequencies/pop_ALL_AFs_chr$i

echo Done with chromosome ${i}

done