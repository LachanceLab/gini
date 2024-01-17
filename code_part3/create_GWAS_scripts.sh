# 5 create_GWAS_scripts.sh
dir_code3="/storage/coda1/p-jlachance6/0/shared/gini/gini/code_part3/"
loc_bfile_prefix="/storage/scratch1/1/ncarvalho6/03-06_UKB/UKB_bfiles/UKB_BED_AFTER_QC_removed182IDs_chr"
dir_simulations="../generated_data/simulations/"

dir_GWAS_scripts=${dir_code3}"GWAS_scripts/"
mkdir -p ${dir_GWAS_scripts}

dir_out=${dir_simulations}GWAS_results_PC/
mkdir -p ${dir_out}

loc_pheno=${dir_simulations}pheno_values.txt
loc_covar=${dir_simulations}covar_PCs.txt

loc_SNPs=${dir_simulations}GWAS_ranges.txt
loc_out=${dir_out}sim_GWAS_chr

for j in {1..22}; do

hrs=$((70 - (j-1)*2))

echo '#!/bin/bash
#SBATCH -Jchr'${j}'_sim_GWAS                    # Job name
#SBATCH --account=gts-jlachance6-biocluster                 # charge account
#SBATCH -N1 --ntasks-per-node=4                 # Number of nodes and cores per node required
#SBATCH --mem-per-cpu=6G                        # Memory per core
#SBATCH -t'${hrs}':00:00                                    # Duration of the job (Ex: 15 mins)
#SBATCH -qinferno                               # QOS Name
#SBATCH -o'${dir_GWAS_scripts}'sim_GWAS_chr'${j}'.out                         # Combined output and error messages file
cd '${dir_GWAS_scripts}'                            # Change to working directory


~/plink2/plink2 \
--bfile '${loc_bfile_prefix}${j}' \
--keep '${dir_code3}${dir_simulations}'pop_GWAS2.txt \
--pheno '${dir_code3}${loc_pheno}' \
--covar '${dir_code3}${loc_covar}' \
--linear allow-no-covars hide-covar cols=+a1freq \
--out '${dir_code3}${loc_out}${j}'

# 

' > ${dir_GWAS_scripts}sim_GWAS_chr${j}.sh
sbatch ${dir_GWAS_scripts}sim_GWAS_chr${j}.sh
echo sim_GWAS_chr${j}.sh

done