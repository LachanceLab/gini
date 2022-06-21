# 3 - append_betas_AFs.R

# Appends the allele frequencies calculated by PLINK onto Prive et al.'s effect
# weights (a.k.a. betas) for each SNP for each trait and then generates a file
# with the appended data for each trait

### Packages and directories ----

library(tidyverse)
library(data.table)

# sets location to allele frequencies generated in previous PLINK script
dir_AFs <- "../generated_data/allele_frequencies/"
# sets location to Prive et al.'s list of SNPs and their betas, from:
# https://figshare.com/articles/dataset/Effect_sizes_for_215_polygenic_scores/14074760/2?file=31619351
loc_betas <- "../prive_data/PGS-effects-PLR.csv"
# sets directory where all the betas+appended AFs will be saved to
dir_out <- "../generated_data/betas_and_AFs/"


### Code ----

# reads the 'PGS-effects-PLR.csv' file and extract the traits (codes)
betas_wide <- as_tibble(fread(loc_betas))
codes <- colnames(betas_wide)[6:length(colnames(betas_wide))]

# creates tibble with correctly formatted column names
betas <- tibble(
  rsid=as.character(),
  chr=as.numeric(),
  pos=as.numeric(),
  a0=as.character(),
  a1=as.character(),
  A1=as.character(),
  A2=as.character()
)
for (code in codes) {
  betas[,as.character(code)] <- as.numeric()
}
col_ancestries <- colnames(AFs[[1]])[5:length(colnames(AFs[[1]]))]
for (col_ancestry in col_ancestries) {
  betas[,col_ancestry] <- as.numeric()
}

# defines function for cases where the reference alleles are mismatched between
# PLINK's allele frequencies and Prive et al.'s betas
flip_AF <- function(AF) {return(1-AF)}

# appends ancestry-specific allele frequencies to list of all SNPs
AFs <- list()
for (chr_i in 1:22) {
  loc_AFs_chr <- paste0(dir_AFs,"pop_ALL_AFs_chr",chr_i,".frq.strat")
  
  AFs_chr <- as_tibble(fread(loc_AFs_chr, drop=c("CHR"))) %>%
    select(SNP, A1, A2, CLST, MAF) %>%
    pivot_wider(names_from=CLST,values_from=MAF,names_prefix="VarFreq_")
  
  betas_wide_chr <- betas_wide %>% filter(chr == chr_i) %>%
    left_join(AFs_chr, by=c("rsid"="SNP"))
  
  betas <- betas %>%
    add_row(betas_wide_chr %>% filter(a0==A1)) %>%
    add_row(betas_wide_chr %>%
              filter(a0!=A1) %>%
              mutate_at(col_ancestries,flip_AF))
  
  print(paste("Loaded allele frequencies for chromosome",chr_i))
}

# goes through each trait and filters SNPs to just those that are significantly
# correlated for that trait and then creates a .txt file containing the SNPs,
# chromosome and positions, alleles, effect weights (betas), and allele 
# frequencies for each trait
for (code in codes) {
  betas_trait <- betas[which(betas[code] != 0),] %>%
    select(chrom=chr,rsID=rsid,chr_position=pos,A1=a0,A2=a1,
           effect_weight=!!as.name(code), starts_with("VarFreq_"))
  
  loc_out <- paste0(dir_out,code,"-betasAFs.txt")
  write.table(betas_trait,loc_out,sep="\t",quote=FALSE,row.names=FALSE)
  print(paste("Created betas+AFs file for",code))
}