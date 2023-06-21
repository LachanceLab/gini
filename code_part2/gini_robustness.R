# gini_robustness.R
# Calculates gini using different parameters to check robustness

### Libraries and directories ####
library(tidyverse)
library(data.table)
source("../code_part1/helper_functions/helper_functions.R")

# sets working directory
setwd("./")

# directory where panUKB GWAS summary statistics were saved to
dir_sf <- "../generated_data/panUKB_sf/"
# sets location to traits_table.txt
loc_table <- "../generated_data/traits_table.txt"

# reads traits_table and extract quantitative trait codes
traits_table <- as_tibble(fread(loc_table))
qcodes <- (traits_table %>% filter(GWAS_trait_type=="quantitative",
                                   PGS_trait_type=="quantitative"))$prive_code

# parameter options
thresholds <- c(50, 100, 250, 500, 750, 1000)
# makes empty table for collecting data
robustness_tbl <- tibble(
  prive_code = as.character(),
  pop = as.character(),
  threshold = as.numeric(),
  gini = as.numeric(),
  n_sig_SNPs = as.numeric(),
  n_sig_SNPs_pop = as.numeric(),
  sum_gvc_top = as.numeric(),
  sum_gvc_all = as.numeric()
)

ii <- 0
# mega-loop
for (code in qcodes) {
  # reads summary file
  #slice <- traits_table %>% filter(prive_code==code)
  filename <- paste0(code,"_sf_indep.txt")
  loc_summary_file <- paste0(dir_sf,filename)
  sf <- as_tibble(fread(loc_summary_file))
  # gets number of SNPs
  n_sig_SNPs <- nrow(sf)
  
  # loops through available gvc columns
  cols_gvc <- colnames(sf %>% select(starts_with("gvc_")))
  pops2use <- sapply(cols_gvc, function(x) substring(x, 5, nchar(x))) %>% unname()
  for (pop in pops2use) {
    col_gvc <- paste0("gvc_",pop)
    sf_pop <- sf %>%
      select(gvc_pop = !!enquo(col_gvc)) %>%
      arrange(-gvc_pop) %>%
      drop_na()
    # gets number of SNPs and sum_gvc (for this population)
    n_sig_SNPs_pop <- nrow(sf_pop)
    sum_gvc_all <- sum(sf_pop$gvc_pop)
    
    for (threshold in thresholds) {
      # gets top n SNPs, padding if necessary
      gvc_list <- (sf_pop %>% filter(row_number() <= threshold))$gvc_pop %>%
        pad_zeros(threshold)
      
      # calculates gini
      gini <- get_gini(gvc_list)
      # gets the sum of gvc among top SNPs
      sum_gvc_top <- sum(gvc_list)
      
      # records data
      robustness_tbl <- robustness_tbl %>% add_row(
        prive_code = code,
        pop = pop,
        threshold = threshold,
        gini = gini,
        n_sig_SNPs = n_sig_SNPs,
        n_sig_SNPs_pop = n_sig_SNPs_pop,
        sum_gvc_top = sum_gvc_top,
        sum_gvc_all = sum_gvc_all
      )
      
      print(paste(ii,code, pop, threshold, round(gini,3), sep="    "))
      ii <- ii + 1
    }
  }
}

# saves results
loc_out <- paste0("../generated_data/gini_robustness_data.txt")
fwrite(robustness_tbl, loc_out, sep="\t")
