# generate_multiple_ginis.R

# An early version of this script was used to generate multiple ginis using
# different settings in the gini calculation in order to compare their results.
# This data was used to select an ideal combination of settings

### Libraries and directories ####
library(tidyverse)
library(data.table)
source("../code_part1/helper_functions.R")

# sets working directory
setwd("./")

# sets location to traits_table.txt
loc_table <- "../generated_data/traits_table.txt"

# location to a file we generated that vastly speeds up the process of binning
# can be obtained from our github under ~/generated_data/
loc_chr_max_bps <- "../code_part1/chr_max_bps.txt"

# sets directory of betas with appended allele frequencies
dir_betas <- "../generated_data/betas_and_AFs/"

# sets directory of generated_data
dir_generated_data <- "../generated_data/"


### Code ###
traits_table <- as_tibble(fread(loc_table))

# loads a file that contains the max base pair position for each chromosome
chr_max_bps <- as_tibble(fread(loc_chr_max_bps))

# will loop through each of these
codes <- traits_table$prive_code
ancestries <- sapply(str_split(colnames(traits_table)[substring(colnames(traits_table),1,5)=="gini_"],"gini_"),"[",2)
bin_sizes <- c(1, 1E4, 1E5, 1E6)
bin_summary_methods <- c("sum","max")
thresholds <- c(10, 25, 50, 100, 150, 250, 500, 1000, 5000, 20000, 50000, 200000)

# used to estimate the progress bar in the console output
max_i <- length(codes) * length(ancestries) * length(bin_sizes) * length(bin_summary_methods) * (length(thresholds)+1)

# Helper function
record_gini <- function(ginis, i, prive_code, n_snps, ancestry, bin_size, n_significant_bins, bin_summary_method, threshold, threshold_zero_padding, gini) {
  ginis <- ginis %>% add_row(
    prive_code = prive_code,
    n_snps = n_snps,
    ancestry = ancestry,
    bin_size = bin_size,
    n_significant_bins = n_significant_bins,
    bin_summary_method = bin_summary_method,
    threshold = threshold,
    threshold_zero_padding = threshold_zero_padding,
    gini = gini
  )
  pct_done <- (i / max_i) * 100
  print(paste(prive_code, bin_size, ancestry, bin_summary_method, threshold, threshold_zero_padding,
              ":: About",round(pct_done,3),"% done."))
  
  return(ginis)
}

ginis <- tibble(
  prive_code = as.character(),
  n_snps = as.numeric(),
  ancestry = as.character(),
  bin_size = as.numeric(),
  n_significant_bins = as.numeric(),
  bin_summary_method = as.character(),
  threshold = as.numeric(),
  threshold_zero_padding = as.logical(),
  gini = as.numeric()
)

i <- 0
# The Hypernested Loop
for (code in codes) {
  
  loc_betas <- paste0(dir_betas,code,"-betasAFs.txt")
  if (!file.exists(loc_betas)) {next}
  
  print(paste("Reading summary file for",code))
  sf <- as_tibble(fread(loc_betas))
  n_snps <- nrow(sf)
  
  for (bin_size in bin_sizes) {
    
    if (bin_size > 1) {sf_binID <- bin_snps(sf, bin_size)}
    else {sf_binID <- sf %>% mutate(bin_ID = row_number())}
    
    for (ancestry in ancestries) {
      
      col_AF <- paste0("VarFreq_", ancestry)
      sf_h2 <- sf_binID %>%
        select(effect_weight, !!as.name(col_AF), bin_ID) %>%
        drop_na() %>%
        get_h2("effect_weight", col_AF)
      
      for (bin_summary_method in bin_summary_methods) {
        
        temp <- sf_h2
        if (bin_size > 1) {sf_binned <- get_data_binned(sf_h2, bin_summary_method)}
        else {sf_binned <- sf_h2}
        n_significant_bins <- nrow(sf_binned)
        
        threshold_zero_padding <- FALSE
        
        for (threshold in thresholds) {
          
          if (threshold < n_significant_bins) {
            sf_filtered <- sf_binned %>% filter(rank <= threshold)
            h2_list <- sf_filtered$h2
            gini <- get_gini(h2_list)
            i <- i + 1
            ginis <- record_gini(ginis, i, code, n_snps, ancestry, bin_size, n_significant_bins, bin_summary_method, threshold, threshold_zero_padding, gini)
          }
          
          if (threshold > n_significant_bins & !threshold_zero_padding) {
            threshold_zero_padding <- TRUE
            h2_list <- c(rep(0, (threshold - n_significant_bins)), sf_binned$h2)
            gini2 <- get_gini(h2_list)
            
            h2_list <- sf_binned$h2
            gini1 <- get_gini(h2_list)
            i <- i + 1
            ginis <- record_gini(ginis, i, code, n_snps, ancestry, bin_size, n_significant_bins, bin_summary_method, n_significant_bins, FALSE, gini1)
            i <- i + 1
            ginis <- record_gini(ginis, i, code, n_snps, ancestry, bin_size, n_significant_bins, bin_summary_method, threshold, TRUE, gini2)
          } else if (threshold > n_significant_bins & threshold_zero_padding) {
            h2_list <- c(rep(0, (threshold - n_significant_bins)), sf_binned$h2)
            gini <- get_gini(h2_list)
            i <- i + 1
            ginis <- record_gini(ginis, i, code, n_snps, ancestry, bin_size, n_significant_bins, bin_summary_method, threshold, threshold_zero_padding, gini)
          }
        }
      }
    }
  }
}
loc_out <- paste0(dir_generated_data,"ginis_extra.txt")
write.table(ginis,loc_out,sep="\t",row.names=FALSE, quote=FALSE)