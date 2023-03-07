# 5 - create_table.R

# This script creates the traits table from which most proceeding analyses take
# their data from. This includes joining all of Prive et al.'s trait info into
# one, and calculating Gini, portability, and PRS divergence for each trait

### Libraries and directories ####
library(tidyverse)
library(data.table)
source("../code_part1/helper_functions/helper_functions.R")
source("../code_part1/helper_functions/winners_curse_functions.R")

# sets working directory
setwd("./")

# sets directory where necessary input files are stored.
# Read input_data/README.md for more information on necessary files
dir_input_data <- "../input_data/"
# directory where summary files with betas+AFs for each trait are stored
# (might be unused in the future?)
#dir_summary_files <- "../generated_data/betas_and_AFs/"
# directory where panUKB GWAS summary statistics were saved to
dir_sf <- "../generated_data/panUKB_sf/"
# directory where the traits_table and other intermediate output files will be
# saved to
dir_out <- "../generated_data/"

# location to the genetic recombination map in hg19 format, created by liftover.py
loc_map <- paste0(dir_input_data,"aau1043_datas3_hg19")
# location to a file we generated that speeds up the process of binning
# can be obtained from our github under ~/code_part1/
#loc_chr_max_bps <- paste0(dir_input_data,"chr_max_bps.txt")



### Code ----

# sets location to files used later in the script
loc_phenotype_description <- paste0(dir_input_data,"phenotype-description.csv")
loc_phenotype_info <- paste0(dir_input_data,"phenotype-info.csv")
loc_pcor <- paste0(dir_input_data,"pred-cor-PLR.csv")
loc_short_labels <-  paste0(dir_input_data,"traits_list.txt")

## Joining trait descriptions ####

# not all traits listed in Prive et al.'s tables had betas/PGS computed, so we
# are going to restrict the tables to just the traits that do have them
pcors <- as_tibble(fread(loc_pcor))
codes <- pcors$pheno %>% unique
#filenames <- dir(dir_summary_files)
#codes <- str_replace(filenames,"-betasAFs.txt","")

# coef_to_liab function from package bigsnpr
coef_to_liab <- function (K_pop, K_gwas = 0.5) {
  z <- stats::dnorm(stats::qnorm(min(K_pop, 1 - K_pop)))
  (K_pop * (1 - K_pop)/z)^2/(K_gwas * (1 - K_gwas))
}

# reads prive's phenotype description and info files
prive_info <- as_tibble(fread(loc_phenotype_info))
prive_description <- as_tibble(fread(loc_phenotype_description))
short_labels <- as_tibble(fread(loc_short_labels))
# defines table containing all info related to each trait
traits_table <- prive_description %>%
  filter(phenotype %in% codes) %>%
  left_join(short_labels %>% select(-group, -description), by="phenotype") %>%
  left_join(prive_info, by=c("phenotype"="pheno")) %>%
  dplyr::rename("prive_code"="phenotype") %>%
  mutate(PGS_trait_type = ifelse(is.na(N),"binary","quantitative"),
         GWAS_trait_type = ifelse(!is.na(n_controls_EUR),"binary","quantitative"),
         prevalence = ifelse(is.na(N),N_case / (N_case+N_control),as.numeric(NA))) %>%
  rowwise() %>% mutate(
         N_total = sum(N,N_case,N_control,na.rm=TRUE),
         liab_coef = coef_to_liab(prevalence, prevalence))

# filter traits to just those with prevalence > 1% (binary only) and respective
# panUKB GWAS data (should have 163 traits remaining)
traits_table <- traits_table %>%
  filter((prevalence > 0.01) | (PGS_trait_type=="quantitative")) %>%
  filter(!is.na(wget))

# consolidates trait groups into fewer categories
groups_consolidated <- list(
  "diseases" = c("circulatory system","dermatologic","digestive",
                 "endocrine/metabolic","genitourinary","hematopoietic",
                 "musculoskeletal","neoplasms","neurological",
                 "respiratory","sense organs", "symptoms"),
  "biological measures" = c("biological measures"),
  "lifestyle/psychological" = c("lifestyle and environment","psychiatric disorders"),
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
# manually adds trait to lifestyle/psychological category
psychological_codes <- c("fluid_intelligence")
traits_table[traits_table$prive_code %in% psychological_codes,"group_consolidated"] <- "lifestyle/psychological"

## Joining Prive et al.'s partial correlation values ####
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

## Generating gini and recombination rate for each trait ####
# (in early stages of the project, gvc)

# loads a file that contains the max base pair position for each chromosome
#chr_max_bps <- as_tibble(fread(loc_chr_max_bps))

# expands traits_table for gini calculation, obtains list of ancestries
pop_centers <- read.csv(
  "https://raw.githubusercontent.com/privefl/UKBB-PGS/main/pop_centers.csv",
  stringsAsFactors = FALSE)
#ancestries <- sort(pop_centers$Ancestry)
#ancestries[9] <- "United" # in order to match other data
# for (ancestry in ancestries) {
#   col_gini <- paste0("gini_",ancestry)
#   traits_table[col_gini] <- as.numeric(NA)
# }
# creates empty column for recombination rate (in units of cM per Mb)
traits_table$cMperMb <- as.numeric(NA)

# settings used for calculating gini
threshold_zero_padding <- TRUE  # whether traits with less than 100 significant
                                # bins are padded with 0 gvc bins
threshold <- 100                # top # of bins (by gvc) to include
average_window <- 100000        # size of averaging window for recombination rate
pval_cutoff <- 1E-5             # pval cutoff to use for Winner's Curse correction

# loads recombination map, filters out X chromosome, converts chr to number
rec_map <- as_tibble(fread(loc_map)) #%>% filter(Chr != "chrX")
#rec_map$Chr <- as.numeric(substring(rec_map$Chr,4,5))
rec_map$Chr <- substring(rec_map$Chr,4,5)

# loops through each trait and calculates its gini, its recombination rate,
# and keeps track of the SNPs in the top 100 bins in each trait for the UK
top_indep_SNPs <- tibble(
  chrom = as.character(),
  chr_position = as.numeric(),
  SNP = as.character(),
  A1 = as.character(),
  A2 = as.character(),
  prive_code = as.character()
)
pop_ginis <- tibble(
  prive_code = as.character(),
  beta_pop = as.character(),
  af_pop = as.character(),
  gini = as.numeric()
)
# loops through each trait
for (i in 1:nrow(traits_table)) {
  
  # extracts trait-related information
  code <- traits_table$prive_code[i]
  slice <- traits_table %>% filter(prive_code==code)
  
  # reads summary file for trait
  filename <- paste0(code,"_sf_indep.txt")
  loc_summary_file <- paste0(dir_sf,filename)
  sf <- sf_raw <- as_tibble(fread(loc_summary_file)) %>%
    group_by(SNP) %>% filter(pval == min(pval)) %>% ungroup() %>%
    filter(chr != "X") %>%
    mutate(chr = as.numeric(chrom))
  n_sig_SNPs <- nrow(sf)
  
  # filters to SNPs with known allele frequencies
  sf <- sf %>% filter(if_all(starts_with("af_"), ~!is.na(.)),
                      pval < pval_cutoff)
  n_sig_SNPs_allpops <- nrow(sf)
  
  # gets list of populations GWASs were done on
  beta_pops <- substring(colnames(sf %>% select(starts_with("beta_"))),6)
  af_pops <- substring(colnames(sf %>% select(starts_with("af_cases_"))),10)
  trait_type <- "binary"
  if (length(af_pops) == 0) {
    trait_type <- "quantitative"
    af_pops <- substring(colnames(sf %>% select(starts_with("af_"))),4)
  }
  af_pops <- af_pops[af_pops != "meta_hq"]
  
  # defines columns to use in gvc calculations
  if ("meta" %in% beta_pops) {
    col_beta <- "beta_meta"
    col_se <- "se_meta"
  } else {
    col_beta <- "beta_EUR"
    col_se <- "se_EUR"
  }
  
  traits_table[i, "n_sig_SNPs"] <- n_sig_SNPs
  traits_table[i, "n_sig_SNPs_allpops"] <- n_sig_SNPs_allpops
  # gets list of populations GWASs were done on
  GWAS_pops <- substring(colnames(sf %>% select(starts_with("beta_"))),6)
  
  # defines columns to use in gvc calculations
  col_AF <- "af_EUR"
  if ("meta" %in% GWAS_pops) {
    col_beta <- "beta_meta"
    col_se <- "se_meta"
  } else {
    col_beta <- "beta_EUR"
    col_se <- "se_EUR"
  }
  
  ##
  for (pop in af_pops) {
    col_AF <- paste0("af_",pop)
    
    n_pop <- switch(pop,
                    "meta" = "full_cohort_both_sexes",
                    "meta_hq" = "hq_cohort_both_sexes",
                    pop)
    
    n_cases <- slice[1,paste0("n_cases_",n_pop)][[1]]
    if (trait_type=="binary") {
      if (pop == "meta") {
        n_controls <- rowSums(slice %>% select(starts_with("n_controls_")), na.rm=TRUE)
      } else {
        n_controls <- slice[1,paste0("n_controls_",n_pop)][[1]]
      }
    } else {
      n_controls <- 0
    }
    
    #n_controls <- rowSums(slice %>% select(starts_with("n_controls_")), na.rm=TRUE)
    n_total <- n_cases + n_controls
    
    if (!(col_AF %in% colnames(sf))) {
      af_cases_pop <- paste0("af_cases_",pop)
      af_controls_pop <- paste0("af_controls_",pop)
      sf[,col_AF] <- (n_cases * sf[,af_cases_pop] + n_controls * sf[,af_controls_pop]) / n_total
    }
    
    if (pop == "meta" | all(af_pops == c("EUR"))) {
      sf_WC <- sf %>% mutate(discovery.n = n_total) %>%
        select(discovery.beta = !!enquo(col_beta),
               discovery.se = !!enquo(col_se),
               discovery.n, discovery.freq = !!enquo(col_AF))
      sf_WC <- as_tibble(correct_winners_curse(as.data.frame(sf_WC), pval_cutoff))
      sf[,col_beta] <- sf_WC$debiased.beta.mle
    }
    
    
    # gets gvc for each SNP
    sf <- get_gvc(sf, col_beta, col_AF)
    if (pop == "meta" | all(af_pops == c("EUR"))) {
      sum_gvc_all <- sum(sf$gvc)
      sf_raw <- sf_raw%>%
        select(-starts_with("WC_beta."),
               -starts_with("AF."),
               -starts_with("gvc.")) %>%
        left_join(sf %>% select(SNP, WC_beta = !!enquo(col_beta), AF = !!enquo(col_AF), gvc),
                  by="SNP")
      write.table(sf_raw, loc_summary_file, sep="\t", quote=FALSE, row.names = FALSE)
    }
    
    gvc_list <- pad_zeros(sf$gvc, threshold)
    
    if (pop == "meta" | all(af_pops == c("EUR"))) {sum_gvc_top <- sum(gvc_list)}
    
    # calculates gini
    gini <- get_gini(gvc_list)
    
    print(paste(code, pop, gini))
    pop_ginis <- pop_ginis %>% add_row(
      prive_code = code,
      beta_pop = substring(col_beta,6),
      af_pop = pop,
      gini = gini
    )
    if (substring(col_beta,6) != pop) {next}
    
    traits_table[i,"gini_panUKB"] <- gini
    traits_table[i,"sum_gvc_all"] <- sum_gvc_all
    traits_table[i,"sum_gvc_top"] <- sum_gvc_top
    
    # extracts the top significant SNPs
    sf_top <- sf %>% filter(rank <= threshold) %>%
      rename(chrom = chr, chr_position = pos, A1 = ref, A2 = alt) %>%
      select(chrom, chr_position, SNP, A1, A2, gvc)
    
    sf_top_out <- sf_top %>%
      select(chrom, chr_position, SNP, A1, A2) %>%
      mutate(prive_code = code, chrom = as.character(chrom))
    top_indep_SNPs <- top_indep_SNPs %>%
      add_row(sf_top_out)
    
    
    ##########
    ##########
    ##########
    
    # calculates recombination rate
    recombination_vector <- c()
    for (chr_i in c(1:22,"X")) {
      # extracts recombination map and significant SNPs for chromosome
      rec_map_chr <- rec_map %>% filter(Chr == chr_i)
      sf_top_chr <- sf_top %>% filter(chrom == chr_i)
      # skips to next chromosome if no significant SNPs found in chromosome
      if (nrow(sf_top_chr) == 0) {next}
      # loops through each significant SNP
      for (j in 1:nrow(sf_top_chr)) {
        # sets window size around SNP to average out recombination in
        start_region <- (sf_top_chr$chr_position[j]) - average_window/2
        end_region <- (sf_top_chr$chr_position[j]) + average_window/2 - 1
        rec_map_bins <- rec_map_chr %>% filter(Begin <= end_region,
                                               End > start_region)
        # manually sets recombination rate to 0 or NA if no recombination data
        # found in window surrounding SNP
        if (nrow(rec_map_bins) == 0) {SNP_cMperMb <- NA}
        else if (max(rec_map_bins$cMperMb) == 0) {SNP_cMperMb <- 0}
        else {
          # calculates proportion of window size occupied by each region in
          # the recombination rate map
          rec_map_bins$Begin[1] <- start_region
          rec_map_bins$End[nrow(rec_map_bins)] <- end_region
          rec_map_bins$bp_range <- (rec_map_bins$End - rec_map_bins$Begin)
          rec_map_bins$bp_range <- rec_map_bins$bp_range / sum(rec_map_bins$bp_range)
          
          # calculates recombination rate as an arithmetic mean weighted by
          # proportion of window size occupied by recombination region
          SNP_cMperMb <- sum(rec_map_bins$bp_range * rec_map_bins$cMperMb) / sum(rec_map_bins$bp_range)
        }
        recombination_vector <- c(recombination_vector,SNP_cMperMb)
      }
    }
    
    # sets trait recombination rate as arithmetic mean of significant SNPs'
    # recombination rate weighted by gvc
    sf_top$cMperMb <- recombination_vector
    trait_cMperMb <- sum(sf_top$gvc * sf_top$cMperMb, na.rm=TRUE) / sum(sf_top$gvc)
    traits_table$cMperMb[i] <- trait_cMperMb
    # prints progress to console
    print(paste0("Trait ", code,
                 " :: Gini_panUKB = ", round(traits_table$gini_panUKB[i],3),
                 " :: Rec Rate = ", round(traits_table$cMperMb[i],3)))
    
  }
}
# binds the population-specific Ginis to traits_table
pop_ginis2 <- pop_ginis %>%
  pivot_wider(
    names_from = af_pop,
    names_prefix = "gini_",
    values_from = gini
  )
traits_table <- traits_table %>% left_join(pop_ginis2, by = "prive_code")
# saves list of significant SNPs
top_indep_SNPs <- top_indep_SNPs %>% distinct()
loc_out <- paste0(dir_out,"top_indep_SNPs.txt")
write.table(top_indep_SNPs, loc_out, row.names = FALSE, quote = FALSE, sep="\t")
# saves a version of this file with just rsIDs for PLINK to use later
loc_out <- paste0(dir_out,"top_indep_SNPs_rsIDs.txt")
write.table(top_indep_SNPs %>% select(SNP, A2) %>% distinct(),
            loc_out, row.names = FALSE, col.names = FALSE, quote = FALSE, sep="\t")

## Calculating portability indices ####

# obtains mean PC distance between ancestries, similar to Prive et al.
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
  
  # extracts best fit line slope, standard error, and p-value
  lin_model <- lm((temp$relative_pcor - 1) ~ 0 + temp$prive_dist_to_UK)
  portability_index <- summary(lin_model)$coefficients[1,1]
  portability_index_SE <- summary(lin_model)$coefficients[1,2]
  portability_index_P <- summary(lin_model)$coefficients[1,4]
  
  portability_indices <- append(portability_indices,portability_index)
  portability_index_SEs <- append(portability_index_SEs,portability_index_SE)
  portability_index_Ps <- append(portability_index_Ps,portability_index_P)
  
  print(paste(code,portability_index))
}
# sets a hard cap of 0 for portability
portability_indices[portability_indices > 0] <- 0
# appends portability index statistics to traits_table
traits_table$portability_index <- portability_indices
traits_table$portability_index_SE <- portability_index_SEs
traits_table$portability_index_P <- portability_index_Ps

## Calculating F_statistic (PGS divergence) ####

# THIS IS DONE IN THE NEXT TWO SCRIPTS:
# encode_sampled_genotypes.sh
# calculate_divergence.R

## temporary solution: import old divergence statistics:
traits_table2 <- as_tibble(fread("../generated_data/traits_table_ASHG.txt"))
traits_table2 <- traits_table2 %>% select(prive_code, f_stat) %>%
  mutate(log_F = log10(f_stat))
traits_table <- traits_table %>% left_join(traits_table2, by="prive_code")

## Saving the traits_table to system
loc_out <- paste0(dir_out,"traits_table.txt")
write.table(traits_table,loc_out,sep="\t",quote=FALSE,row.names=FALSE)
