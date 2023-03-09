# 7 - calculate_divergence.R

# This script calculates the PGS for each sampled individual and calculates the
# level of divergence in PGS using the ANOVA F-stat. Warning: this script can
# take a while; you may want to split the for loop into multiple parallel jobs

### Libraries and directories ####
library(tidyverse)
library(data.table)

# sets working directory
setwd("./")

# sets location of traits_table
loc_table <- "../generated_data/traits_table.txt"
# sets directory for generated_data in general, the encoded genotypes made by
# PLINK, and the summary files with appended allele frequencies
dir_generated_data <- "../generated_data/"
dir_encoded_geno <- "../generated_data/encoded_genotypes/"
dir_betasAFs <- "../generated_data/betas_and_AFs/"


## Code ####

# loads traits table and extract prive_codes (traits)
traits_table <- as_tibble(fread(loc_table))
codes <- traits_table$prive_code

# loads the SNPs in the top 100 bins of the traits
loc_file <- paste0(dir_generated_data,"top100bin_SNPs.txt")
top100bin_SNPs <- as_tibble(fread(loc_file))

# loads the list of IIDs randomly sampled from each population
loc_file <- paste0(dir_generated_data,"pop_sampled_IIDs.txt")
pop_sampled <- as_tibble(fread(loc_file, fill=TRUE)) %>% select(-V4) %>%
  rename(FID=V1, IID=V2, ancestry=V3)

f_stats <- c()
p_values <- c()
# loops through each trait, calculates a PGS for the sampled individuals,
# and computes the ANOVA F-stat and p-value
for (code in codes) {
  # extracts the SNPs that are in the top 100 bins for this trait
  SNPs_sig <- (top100bin_SNPs %>% filter(prive_code == code))$rsID
  
  # loads the summary file of this trait
  loc_sf <- paste0(dir_betasAFs,code,"-betasAFs.txt")
  sf <- as_tibble(fread(loc_sf)) %>%
    select(-starts_with("varFreq")) %>%
    filter(rsID %in% SNPs_sig) %>%
    mutate(rsID_A2 = paste0(rsID,"_",A2))
  
  # loads the encoded genotypes for the selected significant SNPs
  geno <- tibble()
  # loops through each chromosome
  for (chr in 1:22) {
    # filters to just significant SNPs and genotype information for chromosome
    cols <- (sf %>% filter(chrom==chr))$rsID_A2
    loc_geno <- paste0(dir_encoded_geno,"encoded_genotype_chr",chr,".raw")
    geno_chr <- as_tibble(fread(loc_geno, select=c("IID",cols)))
    
    if (nrow(geno) == 0) {geno <- geno_chr
    } else {geno <- geno %>% left_join(geno_chr, by="IID")}
  }
  
  dps <- c()
  # loops through each individual in the sample and computes their PGS
  for (i in 1:nrow(geno)) {
    geno[is.na(geno)] <- 0
    geno_vec <- geno[i,-1] %>% unlist(use.names=FALSE)
    dot_product <- (sf[["effect_weight"]] %*% geno_vec)[1]
    dps[i] <- dot_product
    if (i %% 1000 == 0) {print(paste(code,":: Calculated PRS for",i,"individuals"))}
  }
  # appends the PGS to the table containing the IIDs and ancestries of the
  # sampled individuals
  PRS <- geno %>% select(IID) %>% mutate(PRS = dps)
  pop_sampled <- pop_sampled %>% left_join(PRS, by="IID") %>%
    filter(ancestry != "Ashkenazi")
  
  # computes an ANOVA and extracts F-stat and p-value
  aov_model <- aov(PRS ~ ancestry, data=pop_sampled)
  f_stats <- c(f_stats, summary(aov_model)[[1]][1,4])
  p_values <- c(p_values, summary(aov_model)[[1]][1,5])
  
  # changes name of column containing the PRS for this trait to the prive_code
  colnames(pop_sampled)[which(colnames(pop_sampled)=="PRS")] <- code
}
# appends the f-stat and p-value to the traits_table and saves to system
traits_table$f_stat <- f_stats
traits_table$p_value_f <- p_values
write.table(traits_table,loc_table,sep="\t",quote=FALSE,row.names=FALSE)

# saves the PGS for all sampled individuals and all traits to a file
loc_out <- paste0(dir_generated_data,"pop_sampled_PRSs.txt")
write.table(pop_sampled,loc_out,row.names=FALSE,quote=FALSE)