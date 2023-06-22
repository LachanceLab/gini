# 4 - calculate_divergence.sh

# Computes PGS for all traits on a subset of the UKB population using PLINK2

loc_bfile_prefix="/storage/scratch1/1/ncarvalho6/03-06_UKB/UKB_bfiles/UKB_BED_AFTER_QC_removed182IDs_chr"
# sets path to list of IIDs sampled to be used in calculation, generated in assign_ancestries_to_IIDs.R
loc_keep="../generated_data/pop_sampled_IIDs.txt"
# sets directory where encoded genotype files will be outputted to
dir_out="../generated_data/polygenic_scores/"
mkdir -p $dir_out

loc_betas_csv="../input_data/PGS-effects-PLR.csv"
loc_out_prefix=${dir_out}"ALL_traits-PGS_chr"

loc_betas_txt="${loc_betas_csv::-3}txt"
awk -F',' 'BEGIN{OFS="\t";} {$1=$1;print $0}' ${loc_betas_csv} > ${loc_betas_txt}
NF=$(head -n 1 ${loc_betas_txt} | awk -F'\t' '{if (NR=1) print NF}')
score_cols=$(seq 6 1 $NF)


for i in {1..22}; do

echo Chromosome ${i}

~/plink2/plink2 \
--bfile ${loc_bfile_prefix}${i} \
--keep ${loc_keep} \
--score ${loc_betas_txt} 1 5 header-read ignore-dup-ids cols=maybefid,phenos,scoresums \
--score-col-nums ${score_cols} \
--out ${loc_out_prefix}${i}

done