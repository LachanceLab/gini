library(tidyverse)
library(data.table)
source("../code_part1/helper_functions/helper_functions.R")

setting <- '_p8_r01'
dir_combined <- "../generated_data/simulations/GWAS_results_PC/combined/"
loc_simphenos <- "../generated_data/simulations/simphenos_tbl.txt"

simphenos_tbl <- as_tibble(fread(loc_simphenos))


settings <- c('_p5_r2','_p8_r01')
for (setting in settings) {
  dir_clumped <- paste0("../generated_data/simulations/GWAS_results_PC/clumped",setting,'/')
  
  simphenos_tbl$gini_obs <- as.numeric(NA)
  simphenos_tbl$gini_obs_nopad <- as.numeric(NA)
  
  # traitname. <- "h1_Mc5000_t1"
  
  for (i in 1:nrow(simphenos_tbl)) {
    traitname. <- simphenos_tbl$traitname[i]
    print(traitname.)
    
    loc_clumped <- paste0(dir_clumped, "sim_GWAS_ALL.ph_",traitname.,".clumped")
    clumped <- as_tibble(fread(loc_clumped))
    
    if (nrow(clumped) > 25000) {
      print(paste("Skipping",traitname.))
      next
    }
    
    loc_combined <- paste0(dir_combined, "sim_GWAS_ALL.ph_",traitname.,".glm.linear")
    combined <- as_tibble(fread(loc_combined))
    
    
    indep_sf <- combined[combined$varid %in% clumped$SNP,] %>%
      mutate(gvc = 2 * BETA^2 * A1_FREQ * (1 - A1_FREQ))
    
    gvc_list <- pad_zeros(indep_sf$gvc, 500)
    gini <- get_gini(gvc_list)
    
    simphenos_tbl$gini_obs[i] <- gini
    simphenos_tbl$gini_obs_nopad[i] <- get_gini(indep_sf$gvc)
  }
  
  colnames(simphenos_tbl)[(-1:0)+ncol(simphenos_tbl)] <- paste0(colnames(simphenos_tbl)[(-1:0)+ncol(simphenos_tbl)],setting)
}