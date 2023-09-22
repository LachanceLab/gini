# 3 compute_genetic_liability.sh

loc_bfile_prefix="/storage/scratch1/1/ncarvalho6/03-06_UKB/UKB_bfiles/UKB_BED_AFTER_QC_removed182IDs_chr"
dir_simulations="../generated_data/simulations/"

dir_out=${dir_simulations}genetic_liabilities/
mkdir -p $dir_out
loc_out_prefix=${dir_out}genetic_liabilities_chr
loc_betas=${dir_simulations}all_true_betas.txt

NF=$(head -n 1 ${loc_betas} | awk -F'\t' '{if (NR=1) print NF}')
score_cols=$(seq 12 1 $NF)


for i in {1..22}; do

echo Chromosome ${i}

~/plink2/plink2 \
--bfile ${loc_bfile_prefix}${i} \
--score ${loc_betas} 5 4 header-read ignore-dup-ids cols=maybefid,phenos,scoresums \
--score-col-nums ${score_cols} \
--out ${loc_out_prefix}${i}

done
