library(tidyverse)
library(data.table)

dir_combined <- '/storage/coda1/p-jlachance6/0/shared/gini/gini/generated_data/simulations/GWAS_results_PC/combined/'

loc_simphenos <- '/storage/coda1/p-jlachance6/0/shared/gini/gini/generated_data/simulations/simphenos_tbl.txt'
simphenos <- as_tibble(fread(loc_simphenos))
simphenos$lambda <- as.numeric(NA)


for (i in 1:nrow(simphenos)) {

  if (simphenos$trial[i] != 1) {next}
  
  traitname <- simphenos$traitname[i]
  
  loc_data <- paste0(dir_combined,'sim_GWAS_ALL.ph_',traitname,'.glm.linear')
  
  data.full <- as_tibble(fread(loc_data, fill=TRUE))
  
  data.full <- data.full %>%
    mutate(T_STAT = as.numeric(T_STAT),
           chisq = T_STAT^2) %>%
    filter(!is.na(T_STAT))
  
  lambda <- median(data.full$chisq)/qchisq(0.5,1)
  print(paste(i,traitname,round(lambda,2)))
  
  simphenos$lambda[i] <- lambda
}

loc_simphenos <- '/storage/coda1/p-jlachance6/0/shared/gini/gini/generated_data/simulations/simphenos_tbl.txt'
fwrite(simphenos, loc_simphenos, sep='\t')
