#PBS -N download_clump_panUKB
#PBS -l nodes=1:ppn=1
#PBS -l pmem=16gb
#PBS -l walltime=48:00:00
#PBS -A GT-jlachance6
#PBS -q inferno
#PBS -j oe
#PBS -o download_clump_panUKB.out

cd $PBS_O_WORKDIR

module load r/4.2.1-tidy

dir_wd='/storage/home/hcoda1/1/ncarvalho6/gini_shared/UKBB-Genetic-Architecture/code_part1/'
dir_sf='/storage/home/hcoda1/1/ncarvalho6/gini_shared/UKBB-Genetic-Architecture/generated_data/panUKB_sf/'
loc_traits_list='/storage/home/hcoda1/1/ncarvalho6/gini_shared/UKBB-Genetic-Architecture/code_part1/traits_list.txt'
loc_bfile='/storage/home/hcoda1/1/ncarvalho6/scratch/12-19_1000G_phase3/phase3kG'

# Set the field separator to a newline character
IFS=$'\n'
# Read the output of the awk command into an array
phenotypes=()
while read -r line; do phenotypes+=("$line"); done <<< "$(awk -F'\t' '(NR>1 && $21 == "TRUE") {print $1}' $loc_traits_list)"
aws_links=()
while read -r line; do aws_links+=("$line"); done <<< "$(awk -F'\t' '(NR>1 && $21 == "TRUE") {print $15}' $loc_traits_list)"

pval_cutoff=0.0001
length=${#phenotypes[@]}
#mkdir ${dir_sf}
for i in $(seq 0 4)
#for i in $(seq 0 $((length-1)))
do
phenotype=${phenotypes[$i]}
echo ========== [${i} ${phenotype}] ==========

echo [${i} ${phenotype}] Downloading panUKB summary file
aws_link=${aws_links[$i]}
wget -qO ${dir_sf}${phenotype}_raw_sf.tsv.bgz ${aws_link}

echo [${i} ${phenotype}] Unzipping summary file
zcat ${dir_sf}${phenotype}_raw_sf.tsv.bgz > ${dir_sf}${phenotype}_raw_sf.tsv
rm -v ${dir_sf}${phenotype}_raw_sf.tsv.bgz

echo [${i} ${phenotype}] Transforming summary file for PLINK
index=$(awk -F'\t' '{for (i=1; i<=NF; i++) if ($i == "pval_meta") print i; exit}' ${dir_sf}${phenotype}_raw_sf.tsv)
# if index is empty (no meta GWAS was done), it instead looks for pval_EUR
if [ -z "$index" ]; then index=$(awk -F'\t' '{for (i=1; i<=NF; i++) if ($i == "pval_EUR") print i; exit}' ${dir_sf}${phenotype}_raw_sf.tsv); fi
#echo $index

awk -F'\t' -v col="$index" '{if (NR == 1) {print "SNP\tpval\t"$0} else if ($col != "NA"){print $1":"$2"\t"exp($col)"\t"$0}}' ${dir_sf}${phenotype}_raw_sf.tsv > ${dir_sf}${phenotype}_sf.tsv

echo [${i} ${phenotype}] Running --clump in PLINK

~/plink1_9/plink \
--bfile ${loc_bfile} \
--allow-extra-chr \
--keep /storage/home/hcoda1/1/ncarvalho6/scratch/12-19_1000G_phase3/EUR_IDs.txt \
--clump ${dir_sf}${phenotype}_sf.tsv \
--clump-snp-field SNP \
--clump-field pval \
--clump-p1 ${pval_cutoff} \
--clump-r2 0.2 \
--out ${dir_sf}${phenotype}_sf

echo [${i} ${phenotype}] Filtering GWAS summary file for independent SNPs

Rscript ${dir_wd}filter_GWAS_independent.R ${phenotype}

rm -v ${dir_sf}${phenotype}_raw_sf.*
rm -v ${dir_sf}${phenotype}_sf.*

done

echo DONE WITH ALL TRAITS

