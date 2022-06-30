# 2 - create_ancestry_AFs.sh

# Loops through each chromosome and calculates allele frequencies, stratified by ancestry

for i in {1..22};
do

~/plink1_9/plink \
--bfile /directory/UKB_QC_chr$i \
--keep ../generated_data/pop_ALLrc_IIDs.txt \
--within ../generated_data/pop_ALLrc_IIDs.txt \
--freq \
--out ../generated_data/allele_frequencies/pop_ALL_AFs_chr$i

done