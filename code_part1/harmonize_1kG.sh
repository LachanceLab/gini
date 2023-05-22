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
for i in {1..22}; do
  url="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
  file_gz="ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
  file=${file_gz:0:-3}
  wget ${url}
  gunzip -f ${file_gz}
  
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
done


