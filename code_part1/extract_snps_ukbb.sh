#! /usr/bin/bash

#Example input: ./extract_snps_ukbb.sh -t 4 -p <plink_path> -u <ukbb_path> -s <snp_dir> -k <keep_iids> -o <output_directory>

threads=1
#add --keep command

function HELP {
	echo "Use -t flag for number of threads"
	echo "Use -p flag for path to plink"
	echo "Use -u flag for path to UKBB data"
	echo "Use -s flag for directory of SNP.txt files"
	echo "Use -k flag to indicate which iids to keep during extraction"
	echo "Use -o flag for output directory"
	exit 2
}

while getopts "t:p:u:s:k:o:v" option; do 
	case $option in
		t) threads=$OPTARG;;
		p) plink_path=$OPTARG;;
		u) ukbb_path=$OPTARG;;
		s) snp_dir=$OPTARG;;
		k) keep_iids=$OPTARG;;
		o) out_dir=$OPTARG;;
		v)set -x;;
		\?) HELP;;
	esac
done

mkdir $out_dir

for FILE in $snp_dir/*;
do 
	extension="${FILE##*/}"
	extension=${extension:0:-13}
	cd $out_dir
	mkdir $extension
	echo "####################Running ${extension}####################"
	for i in $(seq 1 22); do 
		${plink_path} --bfile ${ukbb_path}/UKB_BED_AFTER_QC_removed182IDs_chr${i} --keep ${keep_iids} --recode A --recode-allele $FILE --extract $FILE --out ${out_dir}/${extension}/genotype_matrix_chr${i} --threads $threads
	done
done