# 5 - create_table.R

# This script creates the traits table from which most proceeding analyses take
# their data from. This includes joining all of Prive et al.'s trait info into
# one, and calculating gini, portability, and PRS divergence for each trait

### Libraries and directories ####
library(tidyverse)
library(data.table)
source("../code_part1/helper_functions.R")

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

# location to the genetic recombination map in hg19 format, created by liftover.py
loc_map <- "../other_data/aau1043_datas3_hg19"
# location to a file we generated that vastly speeds up the process of binning
# can be obtained from our github under ~/generated_data/
loc_chr_max_bps <- "../code_part1/chr_max_bps.txt"
# directory where summary files with betas+AFs for each trait are stored
dir_summary_files <- "../generated_data/betas_and_AFs/"
# directory where the traits_table and other intermediate output files will be
# saved to
dir_out <- "../generated_data/"

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
traits_table$group_consolidated <- as.character(NA)
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

## Generating gini for each trait ####

# loads a file that contains the max base pair position for each chromosome
chr_max_bps <- as_tibble(fread(loc_chr_max_bps))

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
# creates empty column for recombination rate (in units of cM per Mb)
traits_table$cMperMb <- as.numeric(NA)

# settings used for calculating gini
bin_size <- 100000              # bp size of bins
bin_summary_method <- "sum"     # how heritability is combined for a bin
threshold_zero_padding <- TRUE  # whether traits with less than 100 significant
                                # bins are padded with 0 heritability bins
threshold <- 100                # top # of bins (by heritability) to include
average_window <- bin_size      # size of averaging window for recombination rate

# loads recombination map
rec_map <- as_tibble(fread(loc_map)) %>% filter(Chr != "chrX")
rec_map$Chr <- as.numeric(substring(rec_map$Chr,4,5))

# loops through each trait and calculates its gini, its recombination rate,
# and keeps track of the SNPs in the top 100 bins in each trait for the UK
top100bin_SNPs <- tibble(
  chrom = as.numeric(),
  rsID = as.character(),
  chr_position = as.numeric(),
  A1 = as.character(),
  A2 = as.character(),
  prive_code = as.character()
)
for (i in 1:nrow(traits_table)) {
  code <- traits_table$prive_code[i]
  filename <- paste0(code,"-betasAFs.txt")
  loc_summary_file <- paste0(dir_summary_files,filename)
  
  summary_file <- as_tibble(fread(loc_summary_file)) %>%
    bin_snps(bin_size)
  n_snps <- nrow(summary_file)
  
  for (ancestry in ancestries) {
    col_AF <- paste0("VarFreq_",ancestry)
    col_gini <- paste0("gini_",ancestry)
    pop_data <- summary_file %>% drop_na() %>%
      get_h2("effect_weight", col_AF)
    
    pop_data_binned <- get_data_binned(pop_data, bin_summary_method)
    
    n_significant_bins <- nrow(pop_data_binned)
    
    if (threshold < n_significant_bins) {
      h2_list <- (pop_data_binned %>% filter(rank <= threshold))$h2
    } else {
      h2_list <- c(rep(0, (threshold - n_significant_bins)), pop_data_binned$h2)
    }
    pop_gini <- get_gini(h2_list)
    
    traits_table[i,col_gini] <- pop_gini
    
    if (ancestry == "United") {
      # extracts significant SNPs
      pop_data_binned <- pop_data_binned %>%
        arrange(-h2) %>%
        filter(row_number() <= threshold) %>%
        arrange(bin_ID)
      
      pop_data_sig <- pop_data %>%
        filter(bin_ID %in% pop_data_binned$bin_ID) %>%
        select(chrom, rsID, chr_position, A1, A2) %>%
        mutate(prive_code = code)
      top100bin_SNPs <- top100bin_SNPs %>% add_row(pop_data_sig)
      
      # calculates recombination rate
      pop_data <- pop_data %>%
        filter(bin_ID %in% pop_data_binned$bin_ID) %>%
        arrange(bin_ID)
      
      recombination_vector <- c()
      for (chr_i in 1:22) {
        rec_map_chr <- rec_map %>% filter(Chr == chr_i)
        pop_data_chr <- pop_data %>% filter(chrom == chr_i)
        if (nrow(pop_data_chr) == 0) {next}
        for (j in 1:nrow(pop_data_chr)) {
          start_region <- (pop_data_chr$chr_position[j]) - average_window/2
          end_region <- (pop_data_chr$chr_position[j]) + average_window/2 - 1
          rec_map_bins <- rec_map_chr %>% filter(Begin <= end_region,
                                                 End > start_region)
          
          if (nrow(rec_map_bins) == 0) {SNP_cMperMb <- NA}
          else if (max(rec_map_bins$cMperMb) == 0) {SNP_cMperMb <- 0}
          else {
            rec_map_bins$Begin[1] <- start_region
            rec_map_bins$End[nrow(rec_map_bins)] <- end_region
            rec_map_bins$bp_range <- (rec_map_bins$End - rec_map_bins$Begin)
            rec_map_bins$bp_range <- rec_map_bins$bp_range / sum(rec_map_bins$bp_range)
            
            SNP_cMperMb <- sum(rec_map_bins$bp_range * rec_map_bins$cMperMb) / sum(rec_map_bins$bp_range)
            
            
          }
          recombination_vector <- c(recombination_vector,SNP_cMperMb)
        }
      }
      
      pop_data$cMperMb <- recombination_vector
      trait_cMperMb <- sum(pop_data$h2 * pop_data$cMperMb, na.rm=TRUE) / sum(pop_data$h2)
      traits_table$cMperMb[i] <- trait_cMperMb
    }
  }
  print(paste0("Trait ", code,
               " :: Gini_United = ", round(traits_table$gini_United[i],3),
               " :: Rec Rate = ", round(traits_table$cMperMb[i],3)))
}
loc_out <- paste0(dir_out,"top100bin_SNPs.txt")
write.table(top100bin_SNPs, loc_out, row.names = FALSE, quote = FALSE, sep="\t")
# saves a version of this file with just rsIDs for PLINK to use later
loc_out <- paste0(dir_out,"top100bin_SNPs_rsIDs.txt")
write.table(top100bin_SNPs %>% select(rsID, A2) %>% distinct(),
            loc_out, row.names = FALSE, col.names = FALSE, quote = FALSE, sep="\t")

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

# THIS IS DONE IN THE NEXT TWO SCRIPTS:
# encode_sampled_genotypes.sh
# calculate_divergence.R


## Saving the traits_table
loc_out <- paste0(dir_out,"traits_table.txt")
write.table(traits_table,loc_out,sep="\t",quote=FALSE,row.names=FALSE)
