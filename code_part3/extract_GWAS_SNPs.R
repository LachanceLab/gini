# 5 extract_GWAS_SNPs.R

# Libraries and paths ####
library(tidyverse)
library(data.table)

dir_sims <- "../generated_data/simulations/"
dir_out <- paste0(dir_sims, "GWAS_ranges/")
dir.create(dir_out, showWarnings = FALSE)

# Code ####

# loads true betas
true_betas <- as_tibble(fread(paste0(dir_sims, "all_true_betas.txt")))
# loads simmed phenotype table
simphenos_tbl <- as_tibble(fread(paste0(dir_sims,"simphenos_tbl.txt")))

# for (i in 1:nrow(simphenos_tbl)) {
#   print(i)
#   trait_tb <- true_betas[true_betas[[11+i]] != 0,]
#   # formats for plink --extract range
#   trait_range <- trait_tb %>% mutate(one = 1, BP2 = BP) %>%
#     select(CHR, BP, BP2, one, varid, A0, A1)
#   
#   loc_out <- paste0(dir_out, "GWAS_range_",i,".txt")
#   fwrite(trait_range, loc_out, sep=" ", col.names=FALSE)
# }

GWAS_ranges <- true_betas %>% mutate(one = 1, BP2 = BP) %>%
  select(CHR, BP, BP2, one, varid, A0, A1)

loc_out <- paste0(dir_sims, "GWAS_ranges.txt")
fwrite(GWAS_ranges, loc_out, sep=" ", col.names=FALSE)