# Libraries and paths ####
library(tidyverse)
library(data.table)

dir_LD <- "../input_data/ldscores/"
dir_sims <- "../generated_data/simulations/"

# Code ####
# loops through each population
pops <- c("EAS","AMR","CSA","AFR","EUR")
for (i in 1:length(pops)) {
  pop <- pops[i]
  loc_LD <- paste0(dir_LD,"UKBB.",pop,".ldscores_FULL_adj.txt")
  
  # reads pop-specific LD dara
  LD_pop <- as_tibble(fread(loc_LD))
  LD_pop <- LD_pop %>% select(CHR, BP, A0, A1, rsid, varid, AF)
  colnames(LD_pop)[ncol(LD_pop)] <- paste0("AF_",pop)
  
  # filters running SNP list to SNPs found in LD_pop
  if (i==1) {LD_all <- LD_pop
  } else {LD_all <- LD_all %>% inner_join(LD_pop, by=c("CHR","BP","A0","A1","rsid","varid"))}
}
# about 8.74M SNPs
LD_all <- LD_all %>% 
  filter(AF_EUR > 0.001) %>%
  arrange(CHR, BP, A0, A1)

# writes to system
loc_out <- paste0(dir_sims,"LD_all_SNPs.txt")
fwrite(LD_all, loc_out, sep="\t")
