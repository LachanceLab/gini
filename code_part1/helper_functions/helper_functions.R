# helper_functions.R

# This script just contains functions used by multiple scripts. There is no need
# to run this script as the other files source directly to this one.

# Libraries and directories ####
library(tidyverse)

# Functions ####

# function that calculates and appends individual SNP gvc using betas and allele frequencies
# it is likely in your interest to remove NA beta/AF values before calling function
get_gvc <- function(data_AF, col_beta, col_AF) {
  pop_data <- data_AF %>%
    mutate(
      gvc = 2 * data_AF[[col_beta]]**2 * data_AF[[col_AF]] * (1 - data_AF[[col_AF]]) 
    )
  return(pop_data$gvc) #only returns gvc column
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
# function that pads a list with extra zeros in the beginning to meet some threshold
# number of elements N, or returns a list of the N highest values in the list
pad_zeros <- function(list, threshold) {
  list <- sort(list)
  if (length(list) >= threshold) {
    list_out <- list[(length(list)-threshold+1):length(list)]
  } else {
    list_out <- c(rep(0, threshold - length(list)), list)
  }
  return(list_out)
}
