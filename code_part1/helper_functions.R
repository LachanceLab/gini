# helper_functions.R

# This script just contains functions used by multiple scripts. There is no need
# to run this script as the other files source directly to this one.

### Libraries and directories ####
library(tidyverse)

## Functions ###

# function that creates a dictionary used for giving SNP bins IDs
create_bin_dict = function(bin_size) {
  chr_bins <- chr_max_bps %>%
    mutate(
      bins_per = ceiling(max_bps / bin_size),
      bins = cumsum(bins_per)
    ) %>% select(chr, bins)
  
  bin_dict = list("1" = 0)
  for (i in 2:22) {bin_dict[as.character(i)] <-  chr_bins$bins[i-1]}
  
  return(bin_dict)
}
# function that assigns each SNP a binID before the actual binning function
bin_snps = function(data, bin_size) {
  bin_dict <- create_bin_dict(bin_size)
  data_snp_bins <- data[0,] %>%
    mutate(bin_ID = as.numeric())
  for (chr in 1:22) {
    chr_bins_floor <- bin_dict[[chr]]
    chr_data <- data %>% filter(chrom == chr) %>% mutate(
      bin_ID = ceiling(chr_position / bin_size) + chr_bins_floor
    )
    data_snp_bins <- data_snp_bins %>% add_row(chr_data)
  }
  return(data_snp_bins)
}
# function that groups SNPs together by bin, combining their SNP h^2
get_data_binned = function(data_snp_bins, method="sum") {
  data_binned <- data_snp_bins %>%
    group_by(bin_ID)
  if (method == "sum") {
    data_binned <- data_binned %>% summarise(h2 = sum(h2),bin_n_snps = n())
  } else if (method == "max") {
    data_binned <- data_binned %>% summarise(h2 = max(h2),bin_n_snps = n())
  }
  data_binned <- data_binned %>% arrange(h2) %>%
    mutate(
      rank = nrow(data_binned) - row_number() + 1,
      rank_percentile = rank / max(rank)
    )
  return(data_binned)
}
# function that calculates and appends individual SNP gvc using betas and allele frequencies
get_gvc <- function(data_AF, col_beta, col_AF) {
  pop_data <- data_AF %>%
    mutate(
      gvc = 2 * data_AF[[col_beta]]**2 * data_AF[[col_AF]] * (1 - data_AF[[col_AF]]) 
    ) %>%
    arrange(gvc) %>%
    mutate(
      rank = nrow(data_AF) - row_number() + 1,
      rank_percentile = rank / max(rank)
    )
  return(pop_data)
}
# function that computes the gini of a list of values in ascending order (can't be all zeros either)
get_gini <- function(list) {
  # Adapted from: https://github.com/oliviaguest/gini
  list <- sort(list)
  n <- length(list)
  numerator <- 0
  for (i in 1:n) {numerator <- numerator + (2*i - n - 1)*list[i]}
  denominator <- n * sum(list)
  
  G <- numerator/denominator
  return(G)
}
