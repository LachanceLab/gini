# 7 combine_GWAS_results.sh

dir_code3="/storage/coda1/p-jlachance6/0/shared/gini/gini/code_part3/"
#loc_bfile_prefix="/storage/scratch1/1/ncarvalho6/03-06_UKB/UKB_bfiles/UKB_BED_AFTER_QC_removed182IDs_chr"
dir_sims=${dir_code3}"../generated_data/simulations/"

dir_results=${dir_sims}GWAS_results_PC/
mkdir -p ${dir_results}combined/
loc_simphenos=${dir_sims}simphenos_tbl.txt

dir_GWAS_scripts="combine_scripts/"
mkdir -p ${dir_GWAS_scripts}

# 1. Read the 5th column from the tab-separated file in $loc_simphenos into an array 'traitnames'
mapfile -t -s 1 -O 1 traitnames < <(awk -F'\t' '{print $5}' $loc_simphenos)

# 2. Loop through each traitname in $traitnames
for traitname in "${traitnames[@]}"; do

cat <<- EOF > ${dir_GWAS_scripts}combine_${traitname}.sh
#!/bin/bash
#SBATCH -Jcombine_${traitname}                    # Job name
#SBATCH --account=gts-jlachance6-biocluster                 # charge account
#SBATCH -N1 --ntasks-per-node=1                 # Number of nodes and cores per node required
#SBATCH --mem-per-cpu=2G                        # Memory per core
#SBATCH -t0:05:00                                    # Duration of the job (Ex: 15 mins)
#SBATCH -qinferno                               # QOS Name
#SBATCH -o${dir_GWAS_scripts}combine_${traitname}.out                         # Combined output and error messages file
cd ${dir_GWAS_scripts}                            # Change to working directory

echo ${traitname}
# Initialize an empty file for the combined results for this traitname
> "${dir_results}combined/sim_GWAS_ALL.ph_${traitname}.glm.linear"

# Loop through each chromosome i : 1-22
for i in \$(seq 1 22); do
  
  # If i is 1, copy the whole file including the header
  if [ \$i -eq 1 ]; then
    awk -F'\t' 'NR==1 {print \$0 "\tvarid"}' "${dir_results}sim_GWAS_chr\${i}.ph_${traitname}.glm.linear" > "${dir_results}combined/sim_GWAS_ALL.ph_${traitname}.glm.linear"
  fi
  # Skip the header and append to the combined file
  awk -F'\t' 'NR>1 && \$8=="ADD" {print \$0 "\t" \$1 ":" \$2 "_" \$4 "_" \$5}' "${dir_results}sim_GWAS_chr\${i}.ph_${traitname}.glm.linear" >> "${dir_results}combined/sim_GWAS_ALL.ph_${traitname}.glm.linear"
done
  
EOF

sbatch ${dir_GWAS_scripts}combine_${traitname}.sh
echo ${dir_GWAS_scripts}combine_${traitname}.sh
done