# 8 LD_clump_simGWAS.sh
dir_code3="/storage/coda1/p-jlachance6/0/shared/gini/gini/code_part3/"
dir_sims="../generated_data/simulations/"

dir_1kG='/storage/coda1/p-jlachance6/0/shared/1kG_phase3_GRCh37/'
loc_bfile=${dir_1kG}'1kG_phase3_GRCh37_ALL'
# sets path to text file containing FID and IID of individuals to use in LD panel
loc_LD_IIDs=${dir_1kG}'EUR_IIDs.txt'

dir_results=${dir_sims}GWAS_results/
mkdir -p ${dir_results}clumped/

dir_scripts=${dir_code3}clump_scripts/
mkdir -p ${dir_scripts}

loc_simphenos=${dir_sims}simphenos_tbl.txt
mapfile -t -s 1 -O 1 traitnames < <(awk -F'\t' '{print $5}' $loc_simphenos)


# LD clump

pval_cutoff=0.00001
r2_cutoff=0.2

for traitname in "${traitnames[@]}"; do

echo '#!/bin/bash
#SBATCH -Jclump_'${traitname}'                    # Job name
#SBATCH --account=gts-jlachance6-biocluster                 # charge account
#SBATCH -N1 --ntasks-per-node=1                 # Number of nodes and cores per node required
#SBATCH --mem-per-cpu=16G                        # Memory per core
#SBATCH -t0:05:00                                    # Duration of the job (Ex: 15 mins)
#SBATCH -qinferno                               # QOS Name
#SBATCH -o'${dir_scripts}'clump_'${traitname}'.out                         # Combined output and error messages file
cd '${dir_scripts}'                            # Change to working directory


echo '${traitname}'

~/plink1_9/plink \
--bfile '${loc_bfile}' \
--allow-extra-chr \
--keep '${loc_LD_IIDs}' \
--clump '${dir_code3}${dir_results}combined/sim_GWAS_ALL.ph_${traitname}.glm.linear' \
--clump-snp-field ID \
--clump-field P \
--clump-p1 '${pval_cutoff}' \
--clump-r2 '${r2_cutoff}' \
--out '${dir_code3}${dir_results}clumped/sim_GWAS_ALL.ph_${traitname}'

' > ${dir_scripts}clump_${traitname}.sh
#sbatch ${dir_scripts}clump_${traitname}.sh
echo clump_${traitname}.sh

done