# 6 create_GWAS_scripts.sh
dir_code3="/storage/coda1/p-jlachance6/0/shared/gini/gini/code_part3/"
loc_bfile_prefix="/storage/scratch1/1/ncarvalho6/03-06_UKB/UKB_bfiles/UKB_BED_AFTER_QC_removed182IDs_chr"
dir_simulations="../generated_data/simulations/"
dir_GWAS_ranges=${dir_simulations}GWAS_ranges/

dir_GWAS_scripts="GWAS_scripts/"
mkdir -p ${dir_GWAS_scripts}

dir_out=${dir_simulations}GWAS_results/
mkdir -p ${dir_out}

# gets number of GWAS to do
loc_simphenos=${dir_simulations}simphenos_tbl.txt
N_GWAS=$(($(wc -l < ${loc_simphenos}) - 1))

loc_pheno=${dir_simulations}pheno_values.txt
# 
# for i in $(seq 1 $N_GWAS); do
# # only does trial 1 for now
# if [ $(( (i + 4) % 5 )) -eq 0 ]; then
# echo $i

loc_SNPs=${dir_simulations}GWAS_ranges.txt
loc_out=${dir_out}sim_GWAS_chr

for j in {1..22}; do

echo '#!/bin/bash
#SBATCH -Jsim_GWAS_chr'${j}'                    # Job name
#SBATCH --account=gts-jlachance6-biocluster                 # charge account
#SBATCH -N1 --ntasks-per-node=4                 # Number of nodes and cores per node required
#SBATCH --mem-per-cpu=32G                        # Memory per core
#SBATCH -t24:00:00                                    # Duration of the job (Ex: 15 mins)
#SBATCH -qinferno                               # QOS Name
#SBATCH -o'${dir_GWAS_scripts}'sim_GWAS_chr'${j}'.out                         # Combined output and error messages file
cd '${dir_GWAS_scripts}'                            # Change to working directory


~/plink2/plink2 \
--bfile '${loc_bfile_prefix}${j}' \
--pheno '${dir_code3}${loc_pheno}' \
--linear allow-no-covars \
--extract range '${dir_code3}${loc_SNPs}' \
--bed-border-kb 1000 \
--out '${dir_code3}${loc_out}${j}'

' > ${dir_GWAS_scripts}sim_GWAS_chr${j}.sh
sbatch ${dir_GWAS_scripts}sim_GWAS_chr${j}.sh
echo sim_GWAS_chr${j}.sh

done

# --pheno-col-nums $((2 + i)) \
# fi
# 
# done