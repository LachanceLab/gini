# Goes to input_data folder
dir_repo='/storage/home/hcoda1/1/ncarvalho6/gini_shared/UKBB-Genetic-Architecture/'
dir_1kG=${dir_repo}'input_data/1kG_phase3_GRCh37/'
cd ${dir_1kG}

# Downloads variant manifest file from PanUKB
wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/full_variant_qc_metrics.txt.bgz

# Decompresses variant manifest file
zcat full_variant_qc_metrics.txt.bgz > full_variant_qc_metrics.txt
rm full_variant_qc_metrics.txt.bgz
awk 'NR>1 {print $6}' full_variant_qc_metrics.txt > panUKB_SNPs_varIDs.txt

# Downloads each 'integrated variant call set' file from 1kG phase3 (GRCh37) data
# https://www.internationalgenome.org/data-portal/data-collection/phase-3
# This loop takes over an hour to run
for i in {9..22}; do
  echo =======[ CHROMOSOME ${i} ]=======
  
  # Downloads 1kG vcf file and unzips it
  url="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
  file_gz="ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
  file=${file_gz:0:-3}
  wget ${url}
  gunzip -f ${file_gz}
  
  # Uses plink2 to convert to bed, also filters to panUKBB SNPs
  ~/plink2/plink2 \
  --vcf ${file} \
  --extract panUKB_SNPs_varIDs.txt \
  --allow-extra-chr \
  --set-missing-var-ids @:#_\$r_\$a \
  --new-id-max-allele-len 23 missing \
  --max-alleles 2 \
  --keep-allele-order \
  --make-bed \
  --out 1kG_phase3_GRCh37_chr${i}
  
  # Removes original vcf file to save space
  rm ${file}
  
  # Adds bfile prefix to list used later for merging by PLINK
  echo 1kG_phase3_GRCh37_chr${i} >> list_1kG_phase3_GRCh37.txt
done

# Merges bfiles into one, should be about 25M SNPs
~/plink1_9/plink \
--merge-list list_1kG_phase3_GRCh37.txt \
--make-bed \
--out 1kG_phase3_GRCh37_ALL

# Downloads 1kG sample file
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
# Extracts just the FID, IID, and continent code
awk 'BEGIN {OFS=" "} NR>1 {print "0", $1, $3}' integrated_call_samples_v3.20130502.ALL.panel > 1kG_continent_IDs.txt
# Extracts the FID and IID of EUR samples only
awk '$3 == "EUR" {print $1, $2}' 1kG_continent_IDs.txt > EUR_IIDs.txt