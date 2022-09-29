### Packages and directories ----

library(tidyverse)
library(data.table)
source("../code_part1/helper_functions.R")

# sets working directory
setwd("./")

# directory where summary files with betas+AFs for each trait are stored
dir_summary_files <- "../generated_data/betas_and_AFs/"
# directory where previous generated data is kept
dir_gen_data <- "../generated_data/"
# sets the location of the traits table
loc_table <- "../generated_data/traits_table.txt"


### Code ----

# reads traits table
traits_table <- as_tibble(fread(loc_table)) %>%
  mutate(summed_gvc100 = as.numeric(NA),
         summed_gvc100_raw = as.numeric(NA))

# gets list of top 100 bin SNPs for each trait
loc_top100bins <- paste0(dir_gen_data,"top100bin_SNPs.txt")
top100bins <- as_tibble(fread(loc_top100bins))

cols_keep <- c("rsID","effect_weight","VarFreq_United")
codes <- traits_table$prive_code


for (i in 1:length(codes)) {
  code <- codes[i]
  loc_betas_AFs <- paste0(dir_summary_files,"/",code,"-betasAFs.txt")
  betas_AFs <- as_tibble(fread(loc_betas_AFs, select=cols_keep))
  
  trait_rsIDs <- (top100bins %>% filter(prive_code == code))$rsID
  
  trait_SNPs <- betas_AFs %>% filter(rsID %in% trait_rsIDs) %>%
    get_h2("effect_weight","VarFreq_United")
  
  summed_gvc100 <- sum(trait_SNPs$h2)
  summed_gvc100_raw <- summed_gvc100
  slice <- traits_table %>% filter(prive_code==code)
  if (slice$trait_type[[1]]=="binary") {
    summed_gvc100 <- summed_gvc100_raw * slice$liab_coef[[1]]
  }
  ldpred2_h2 <- (slice)$ldpred2_h2[[1]]
  traits_table$summed_gvc100[i] <- summed_gvc100
  traits_table$summed_gvc100_raw[i] <- summed_gvc100_raw
  print(paste(code,summed_gvc100, summed_gvc100_raw, ldpred2_h2))
}

## Saving the traits_table
#write.table(traits_table,loc_table,sep="\t",quote=FALSE,row.names=FALSE)


######

traits_table <- as_tibble(fread(loc_table))

ggplot(traits_table %>%
         filter(
           (prevalence > 0.01 | trait_type=="quantitative"),
           summed_gvc100_raw < 1
         ),
       aes(x=ldpred2_h2, y=summed_gvc100_raw)) +
  geom_point(aes(color=trait_type)) +
  facet_wrap(~ trait_type)
