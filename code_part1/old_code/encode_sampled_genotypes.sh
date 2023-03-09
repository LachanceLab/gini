# 6 - encode_sampled_genotypes.sh

# This script runs a PLINK command that encodes the genotypes of the sampled
# indviduals for the top 100 bin SNPs so that R can read this data easily

# A --merge-list file is needed where each line is the path to the bfile for a
# chromosome. The path needs to be a prefix, meaning the file extension is
# automatically appended by PLINK
# https://www.cog-genomics.org/plink/1.9/data#merge_list

### Direcories and paths ###

# sets working directory
cd /project/code/

# sets directory to plink
dir_plink="~/plink1_9"
# sets path to UKB binary files. Leave as just the prefix before the chromosome
# number and the file extensions
loc_bfile_prefix="/directory/UKB_QC_chr"
# sets path to list of IIDs sampled to be used in calculation, generated in assign_ancestries_to_IIDs.R
loc_keep="../generated_data/pop_sampled_IIDs.txt"
# set path to list of rsIDs of SNPs in top 100 bins, generated in create_table.R
loc_extract="../generated_data/top100bin_SNPs_rsIDs.txt"
# sets directory where encoded genotype files will be outputted to
dir_out="../generated_data/encoded_genotypes/"

### Code ###

mkdir ${dir_out}

for i in {1..22};
do

${dir_plink}/plink \
--bfile ${loc_bfile_prefix}${i} \
--keep ${loc_keep} \
--extract ${loc_extract} \
--recode A \
--recode-allele ${loc_extract} \
--out ${dir_out}encoded_genotype_chr${i}

echo Done with chromosome ${i}

done

