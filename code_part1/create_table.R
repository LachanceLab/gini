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
# directory where panUKB GWAS summary statistics were saved to
dir_sf <- "../generated_data/panUKB_sf/"
# directory where the traits_table and other intermediate output files will be
# saved to
dir_out <- "../generated_data/"

# location to the genetic recombination map in hg19 format, created by liftover.py
loc_map <- paste0(dir_input_data,"aau1043_datas3_hg19")

# Set to TRUE if you want to calculate rec_rate (old analysis)
calculate_rec_rate <- FALSE
# Set to TRUE if 1kG AFs have already been appended (skips loading step)
already_appended_1kG <- TRUE

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
  rowwise() %>% mutate(N_total = sum(N,N_case,N_control,na.rm=TRUE))

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

# expands traits_table for gini calculation, obtains list of ancestries
pop_centers <- read.csv(
  "https://raw.githubusercontent.com/privefl/UKBB-PGS/main/pop_centers.csv",
  stringsAsFactors = FALSE)
traits_table$cMperMb <- as.numeric(NA)

# settings used for calculating gini
threshold <- 500                # top # of SNPs (by gvc) to include
average_window <- 100000        # size of averaging window for recombination rate
pval_cutoff <- 1E-5             # pval cutoff to use for Winner's Curse correction
all_pops <- c("meta_hq","meta","EUR","AFR","AMR","CSA","EAS")

if (calculate_rec_rate) {
  # loads recombination map, filters out X chromosome, converts chr to number
  rec_map <- as_tibble(fread(loc_map)) #%>% filter(Chr != "chrX")
  #rec_map$Chr <- as.numeric(substring(rec_map$Chr,4,5))
  rec_map$Chr <- substring(rec_map$Chr,4,5)
}

# loops through each trait and calculates its gini, its recombination rate,
# and keeps track of the SNPs in the top 100 bins in each trait for the UK
top_indep_SNPs <- tibble(
  chrom = as.character(),
  chr_position = as.numeric(),
  SNP = as.character(),
  ref = as.character(),
  alt = as.character(),
  prive_code = as.character(),
  pop = as.character(),
  gvc = as.numeric()
)
pop_ginis <- tibble(
  prive_code = as.character(),
  beta_pop = as.character(),
  af_pop = as.character(),
  gini = as.numeric()
)
# loads up 1kG allele frequencies
if (!already_appended_1kG) {
  loc_1kG <- paste0(dir_out,"AFs_1kG_sfSNPs.frq.strat")
  AFs_1kG <- as_tibble(fread(loc_1kG)) %>% select(-MAC, -NCHROBS) %>%
    pivot_wider(
      names_from = "CLST",
      names_prefix = "AF_1kG_",
      values_from = "MAF",
    ) %>%
    separate(SNP, c(NA,"BP"),sep=":", remove=FALSE) %>%
    separate(BP, c("BP", NA, NA), sep="_", remove=TRUE) %>%
    mutate(BP = as.numeric(BP)) %>% rename(AF_1kG_CSA = AF_1kG_SAS)
}

# loops through each trait
for (i in 1:nrow(traits_table)) {
  
  # extracts trait-related information
  code <- traits_table$prive_code[i]
  slice <- traits_table %>% filter(prive_code==code)
  
  # reads summary file for trait
  filename <- paste0(code,"_sf_indep.txt")
  loc_summary_file <- paste0(dir_sf,filename)
  sf <- as_tibble(fread(loc_summary_file)) %>%
    group_by(SNP) %>% filter(pval == min(pval)) %>% ungroup() %>%
    filter(chr != "X") %>% mutate(chr = as.numeric(chr)) %>% unique() %>%
    select(-ends_with("MID")) # removed since not found in 1kG
  
  # appends 1kG AFs to summary files if not already present
  if (!("AF_1kG_EUR" %in% colnames(sf))) {
    sf <- sf %>%
      left_join(AFs_1kG %>% select(CHR, BP, A1, A2, starts_with("AF_1kG_")),
                by=c("chr"="CHR","pos"="BP","ref"="A2","alt"="A1"))
  }
  n_sig_SNPs <- nrow(sf)
  traits_table[i, "n_sig_SNPs"] <- n_sig_SNPs
  
  # gets list of populations GWASs were done on
  beta_pops <- substring(colnames(sf %>% select(starts_with("beta_"), -contains("meta2use"))),6)
  af_pops <- substring(colnames(sf %>% select(starts_with("af_cases_"), -contains("meta2use"))),10)
  trait_type <- "binary"
  if (length(af_pops) == 0) {
    trait_type <- "quantitative"
    af_pops <- substring(colnames(sf %>% select(starts_with("af_", ignore.case = FALSE), -contains("meta2use"))),4)
  }
  af_pops <- beta_pops # temporary fix
  # defines columns to use in gvc calculations later
  if ("meta_hq" %in% beta_pops) { meta2use <- "meta_hq"
  } else if ("meta" %in% beta_pops) { meta2use <- "meta"
  } else { meta2use <- "EUR" }
  col_beta <- paste0("beta_",meta2use)
  col_se <- paste0("se_",meta2use)
  
  # for quantitative traits, makes backup of AFs before 1kG imputation
  if (trait_type=="quantitative") {
    for (pop in all_pops) {
      col_AF_panUKB <- paste0("AF_panUKB_",pop)
      col_AF_1kG <- paste0("AF_1kG_",pop)
      col_AF <- paste0("af_",pop)
      # makes backup of AFs before 1kG imputation
      if (pop %in% af_pops) { sf[,col_AF_panUKB] <- sf[,col_AF]
      } else if (pop!="meta_hq" & pop!="meta") {sf[,col_AF] <- as.numeric(NA)} 
      # imputes missing AF_panUKB with AF_1kG
      if (pop!="meta_hq" & pop!="meta") {
        sf[is.na(sf[[col_AF]]),col_AF] <- sf[is.na(sf[[col_AF]]),col_AF_1kG]
      }
    }
  } else if (trait_type=="binary") {
    for (pop in all_pops) {
      col_AF_1kG <- paste0("AF_1kG_",pop)
      col_AF <- paste0("af_",pop)
      af_cases_pop <- paste0("af_cases_",pop)
      af_controls_pop <- paste0("af_controls_",pop)
      if (pop %in% af_pops) {
        # uses meta cases/controls ratio for meta_hq
        n_pop <- switch(pop,
                        "meta" = "full_cohort_both_sexes",
                        "meta_hq" = "full_cohort_both_sexes", #"hq_cohort_both_sexes",
                        pop)
        # extracts number of cases
        n_cases <- slice[1,paste0("n_cases_",n_pop)][[1]]
        # extracts number of controls, for meta/meta_hq, sums all population controls
        if ((pop == "meta_hq") | (pop == "meta")) {
          n_controls <- rowSums(slice %>% select(starts_with("n_controls_")), na.rm=TRUE)
        } else { n_controls <- slice[1,paste0("n_controls_",n_pop)][[1]] }
        # gets total number of cases/controls
        n_total <- n_cases + n_controls
        if (pop == meta2use) {n_total_meta <- n_total}
        
        # adds pop-specific case+control AFs based on weighted average AFs
        sf[,col_AF] <- (n_cases * sf[,af_cases_pop] + n_controls * sf[,af_controls_pop]) / n_total
      } else if (pop!="meta_hq" & pop!="meta") {sf[,col_AF] <- as.numeric(NA)} 
      # imputes missing AF_panUKB with AF_1kG
      if (pop!="meta_hq" & pop!="meta") {
        sf[is.na(sf[[col_AF]]),col_AF] <- sf[is.na(sf[[col_AF]]),col_AF_1kG]
      }
    }
  }
  # filters to SNPs with AF-data for every pop and under p-val cutoff
  col_afs <- paste0("af_",af_pops)
  sf <- sf %>% filter(if_all(all_of(col_afs), ~!is.na(.)), pval < pval_cutoff)
  
  # for the chosen meta2use column, uses Winner's Curse correction for its beta
  col_AF <- "af_meta2use"
  sf[,col_AF] <- sf[,paste0("af_",meta2use)]
  sf_WC <- sf %>% mutate(discovery.n = n_total_meta) %>%
    select(discovery.beta = !!enquo(col_beta),
           discovery.se = !!enquo(col_se),
           discovery.n, discovery.freq = !!enquo(col_AF))
  sf_WC <- as_tibble(correct_winners_curse(as.data.frame(sf_WC), pval_cutoff))
  col_beta <- "beta_meta2use_WC"
  sf[,col_beta] <- sf_WC$debiased.beta.mle
  # gets gvc for each SNP
  sf[,"gvc_meta2use"] <- get_gvc(sf, col_beta, col_AF)
  # for meta/meta_hq: gets the sum of gvc and joins columns of WC-corrected meta
  # beta, meta AF, and subsequent gvc
  sum_gvc_all <- sum(sf$gvc_meta2use)
  
  # loops through meta + all populations
  for (pop in c(af_pops, all_pops[-c(1:2)]) %>% unique()) {
    col_AF <- paste0("af_",pop)
    col_gvc <- paste0("gvc_",pop)
    # gets gvc for each SNP
    sf[,col_gvc] <- get_gvc(sf, col_beta, col_AF)
    gvc_list <- pad_zeros(sf[[col_gvc]], threshold)
    # calculates gini
    gini <- get_gini(gvc_list)
    # gets the sum of gvc among top SNPs
    if (pop == meta2use) {
      sum_gvc_top <- sum(gvc_list)
      traits_table[i,"gini_panUKB"] <- gini
    }
    # adds data to table
    print(paste(code, pop, round(gini,4)))
    pop_ginis <- pop_ginis %>% add_row(
      prive_code = code,
      beta_pop = substring(col_beta,6),
      af_pop = pop,
      gini = gini
    )
    # extracts the top significant SNPs and adds to growing list
    sf_top <- sf %>% arrange(desc(!!sym(col_gvc))) %>%
      filter(row_number() <= threshold)
    sf_top_out <- sf_top %>%
      rename(chrom = chr, chr_position = pos) %>%
      mutate(prive_code = code, chrom = as.character(chrom), pop=pop) %>%
      select(chrom, chr_position, SNP, ref, alt, prive_code, pop, gvc=!!enquo(col_gvc))
    top_indep_SNPs <- top_indep_SNPs %>% add_row(sf_top_out)
  }
  # saves sf back to system (WARNING: overrides file)
  #fwrite(sf, loc_summary_file, sep="\t")
  
  ##############################################################################
  # calculates recombination rate #
  if (calculate_rec_rate) {
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
  }
  ##############################################################################
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
loc_out <- paste0(dir_out,"top_indep_SNPs.txt")
fwrite(top_indep_SNPs, loc_out, sep="\t")

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
# calculate_divergence.sh
# calculate_divergence.R

## Saving the traits_table to system
loc_out <- paste0(dir_out,"traits_table.txt")
fwrite(traits_table,loc_out,sep="\t")
