# 4 - create_table.R

# This script creates the traits table from which most proceeding analysis take
# their data from. This includes joining all of Prive et al.'s trait info into
# one, and calculating gini, portability, and PRS divergence for each trait

### Libraries and directories ####
library(tidyverse)
library(data.table)

# sets working directory
setwd("./")

# sets directory where necessary prive files are stored. The following files
# need to be in this directory with the same file names listed here:

# 'phenotype-description.csv' download from Prive et al's github:
# https://github.com/privefl/UKBB-PGS/blob/main/phenotype-description.xlsx
# NOTE: MUST CONVERT TO READABLE FORMAT FOR R (.csv, .tsv, .txt, etc.)

# 'phenotype-info.csv' download from Prive et al.'s github:
# https://github.com/privefl/UKBB-PGS/blob/main/phenotype-info.csv

# 'pred-cor-PLR.csv' download from Prive et al.'s figshare:
# https://figshare.com/articles/dataset/Effect_sizes_for_215_polygenic_scores/14074760/2?file=31619357
dir_prive_data <- "../prive_data/"

# location to a file we generated that vastly speeds up the process of binning
# can be obtained from our github under ~/generated_data/
loc_chr_max_bps <- "../generated_data/chr_max_bps.txt"
# directory where summary files with betas+AFs for each trait are stored
dir_summary_files <- "../generated_data/betas_and_AFs/"
# location where the final traits_table should be saved to
loc_traits_table <- "../generated_data/traits_table.txt"

### Code ----

# sets location to files used later in the script
loc_phenotype_description <- paste0(dir_prive_data,"phenotype-description.csv")
loc_phenotype_info <- paste0(dir_prive_data,"phenotype-info.csv")
loc_pcor <-paste0(dir_prive_data,"pred-cor-PLR.csv")

## Joining trait descriptions ####

# not all traits listed in Prive et al.'s tables had betas computed, so we are
# going to restrict the tables to just the taits that do have them
filenames <- dir(dir_summary_files)
codes <- str_replace(filenames,"-betasAFs.txt","")

# reads prive's phenotype description and info files
prive_info <- as_tibble(fread(loc_phenotype_info))
prive_description <- as_tibble(fread(loc_phenotype_description))
traits_table <- prive_description %>%
  filter(phenotype %in% codes) %>%
  left_join(prive_info, by=c("phenotype"="pheno")) %>%
  rename("prive_code"="phenotype") %>%
  mutate(trait_type = ifelse(is.na(N),"binary","quantitative"))

# consolidates trait groups into fewer categories
groups_consolidated <- list(
  "psychological" = c("psychiatric disorders"),
  "diseases" = c("circulatory system","dermatologic","digestive",
                 "endocrine/metabolic","genitourinary","hematopoietic",
                 "musculoskeletal","neoplasms","neurological",
                 "psychiatric disorders","respiratory","sense organs",         
                 "symptoms"),
  "biological measures" = c("biological measures"),
  "lifestyle/environment" = c("lifestyle and environment"),
  "physical measures" = c("injuries & poisonings","physical measures","sex-specific factors")
)
for (i in 1:nrow(traits_table)) {
  for (group_consolidated in names(groups_consolidated)) {
    if (traits_table$group[i] %in% groups_consolidated[[group_consolidated]]) {
      traits_table$group_consolidated[i] <- group_consolidated
      break
    }
  }
}
psychological_codes <- c("fluid_intelligence")
traits_table[traits_table$prive_code %in% psychological_codes,"group_consolidated"] <- "psychological"

## Joining Prive et al.'s partial correlation values ####

pcors <- as_tibble(fread(loc_pcor))
pcors$pop <- gsub('United Kingdom', 'United', pcors$pop)
pcors <- pcors %>%
  rename(Nsize = N) %>%
  arrange(pop) %>%
  pivot_wider(
    names_from = pop,
    names_glue = "{.value}_{pop}",
    values_from = c(Nsize,pcor,inf,sup)
  )
traits_table <- traits_table %>% left_join(pcors, by=c("prive_code"="pheno"))

## Generating gini for each ####

# loads a file that contains the max base pair position for each chromosome
chr_max_bps <- as_tibble(fread(loc_chr_max_bps))
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
# function that assigns each SNP a binID before binning function
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
# function that calculates and appends individual SNP h2 using betas and
# allele frequencies
get_h2 <- function(data_AF, col_beta, col_AF) {
  pop_data <- data_AF %>%
    dplyr::mutate(
      h2 = 2 * data_AF[[col_beta]]**2 * data_AF[[col_AF]] * (1 - data_AF[[col_AF]]) 
    ) %>%
    arrange(h2) %>%
    dplyr::mutate(
      rank = nrow(data_AF) - row_number() + 1,
      rank_percentile = rank / max(rank)
    )
  return(pop_data)
}
# function that computes the gini of a list of values in ascending order
get_gini <- function(list) {
  # Adapted from: https://github.com/oliviaguest/gini
  n <- length(list)
  numerator <- 0
  for (i in 1:n) {numerator <- numerator + (2*i - n - 1)*list[i]}
  denominator <- n * sum(list)
  
  G <- numerator/denominator
  return(G)
}

# expands traits_table for gini calculation
pop_centers <- read.csv(
  "https://raw.githubusercontent.com/privefl/UKBB-PGS/main/pop_centers.csv",
  stringsAsFactors = FALSE)
ancestries <- sort(pop_centers$Ancestry)
ancestries[9] <- "United" # in order to match other data
for (ancestry in ancestries) {
  col_gini <- paste0("gini_",ancestry)
  traits_table[col_gini] <- as.numeric(NA)
}

# settings used for calculating gini
bin_size <- 100000              # bp size of bins
bin_summary_method <- "sum"     # how heritability is combined for a bin
threshold_zero_padding <- TRUE  # whether traits with less than 100 significant
                                # bins are padded with 0 heritability bins
threshold <- 100                # top # of bins (by heritability) to include

# loops through each trait and calculates gini for each
for (i in 1:nrow(traits_table)) {
  code <- traits_table$prive_code[i]
  filename <- paste0(code,"-betasAFs.txt")
  loc_summary_file <- paste0(dir_summary_files,filename)
  
  print(paste("Calculating gini for",code))
  summary_file <- as_tibble(fread(loc_summary_file))
  n_snps <- nrow(summary_file)
  
  if (bin_size > 1) {summary_file <- bin_snps(summary_file, bin_size)}
  else {summary_file <- summary_file %>% mutate(bin_ID = row_number())}
  
  for (ancestry in ancestries) {
    col_AF <- paste0("VarFreq_",ancestry)
    col_gini <- paste0("gini_",ancestry)
    pop_data <- summary_file %>%
      select(effect_weight,!!as.name(col_AF), bin_ID) %>%
      drop_na()
    pop_data[[col_AF]] <- as.numeric(pop_data[[col_AF]])
    
    pop_data <- get_h2(pop_data, "effect_weight",col_AF)
    
    if (bin_size > 1) {pop_data <- get_data_binned(pop_data, bin_summary_method)}
    
    n_significant_bins <- nrow(pop_data)
    
    if (threshold < n_significant_bins) {
      h2_list <- (pop_data %>% filter(rank <= threshold))$h2
    } else {
      h2_list <- c(rep(0, (threshold - n_significant_bins)), pop_data$h2)
    }
    pop_gini <- get_gini(h2_list)
    
    traits_table[i,col_gini] <- pop_gini
  }
}

## Calculating portability indices ####

# obtains mean PC distance between ancestries
prive_PC <- pop_centers %>% select(PC1:PC16)
prive_dist_to_UK <- as.matrix(dist(prive_PC))[,1]
distances <- tibble(
  population = pop_centers$Ancestry,
  prive_dist_to_UK = prive_dist_to_UK
) %>% arrange(population)
distances$population <- str_replace(distances$population,"United Kingdom","United")


portability_indices <- c()
portability_index_SEs <- c()
portability_index_Ps <- c()
codes <- traits_table$prive_code
# loops through trait and calculates the slope of the line of best fit for
# for relative predictive performance against PC distance from UK
for (code in codes) {
  temp <- traits_table %>% filter(prive_code == code) %>%
    select(starts_with("pcor_")) %>%
    pivot_longer(
      cols = starts_with("pcor_"),
      names_prefix = "pcor_",
      names_to = "population",
      values_to = "pcor"
    ) %>% mutate(
      relative_pcor = pcor / (traits_table %>% filter(prive_code == code))$pcor_United[1]
    ) %>% left_join(distances, by="population")
  
  lin_model <- lm((temp$relative_pcor - 1) ~ 0 + temp$prive_dist_to_UK)
  portability_index <- summary(lin_model)$coefficients[1,1]
  portability_index_SE <- summary(lin_model)$coefficients[1,2]
  portability_index_P <- summary(lin_model)$coefficients[1,4]
  
  portability_indices <- append(portability_indices,portability_index)
  portability_index_SEs <- append(portability_index_SEs,portability_index_SE)
  portability_index_Ps <- append(portability_index_Ps,portability_index_P)
  
  print(paste(code,portability_index))
}
# appends portability index statistics to traits_table
traits_table$portability_index <- portability_indices
traits_table$portability_index_SE <- portability_index_SEs
traits_table$portability_index_P <- portability_index_Ps

## Calculating F_statistic (PRS divergence) ####

# insert code here


## Saving the traits_table
write.table(traits_table,loc_traits_table,sep="\t",quote=FALSE,row.names=FALSE)
