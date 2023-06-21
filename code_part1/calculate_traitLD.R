# 8 - calculate_traitLD.R

# calculates a weighted average LD score for SNPs associated with a trait, as
# well as metrics of population-differences in LD scores


## libraries and directories

library(tidyverse)
library(data.table)

dir_generated_data <- "../generated_data/"
dir_sf <- paste0(dir_generated_data,"panUKB_sf/")
dir_input_data <- "../input_data/"
dir_ldscores <- paste0(dir_input_data, "ldscores/")

pops <- c("EUR","AFR","AMR","CSA","EAS")
# set to TRUE if code has already been run once, meaning that the LD scores
# have already been adjusted for AF (saves time)
ldscores_already_ajdusted <- TRUE
# set to TRUE if code has already been run once, meaning that the top independent
# SNPs file already contains the ld-scores for each SNP
top_SNPs_already_appended <- FALSE

## Functions

# function for adjsting LD scores by AF, as described in Gazal et al. 2017:
# https://www.nature.com/articles/ng.3954
adjust_LDscores <- function(ldscores, col_AF="AF", col_ld_score = "ld_score") {
  
  # puts SNPs into 0.05 AF bins, transforms LD score distribution within each
  # bin into standard normal distribution
  ldscores <- ldscores %>%
    mutate(AF_bin = floor(20*!!as.name(col_AF)) /20) %>%
    group_by(AF_bin) %>%
    mutate(ld_score_pct = rank(!!as.name(col_ld_score)) / n()) %>%
    ungroup() %>%
    mutate(ld_score_adj = qnorm(ld_score_pct)) %>%
    filter(!is.infinite(ld_score_adj))
  # sets LD score for MAF < 0.05 SNPs to 0 (skipping this Gazal et al. step)
  #ldscores[ldscores[,col_AF][[1]] < 0.05,"ld_score_adj"] <- 0
  #ldscores[ldscores[,col_AF][[1]] > 0.95,"ld_score_adj"] <- 0
  
  return(ldscores)
}

## CODE

# reads trait_table and top_indep_SNPs from create_table.R
loc_traits_table <- paste0(dir_generated_data,"traits_table.txt")
traits_table <- as_tibble(fread(loc_traits_table))
loc_top_indep_SNPs <- paste0(dir_generated_data,"top_indep_SNPs.txt")
top_indep_SNPs <- as_tibble(fread(loc_top_indep_SNPs))

# appends ld scores to top_indep_SNPs (if not done already)
if (!top_SNPs_already_appended) {
  for (pop in pops) {
    # reads LD scores file
    loc_ldscores <- paste0(dir_input_data, "ldscores/UKBB.",pop,".ldscores_FULL_adj.txt")
    print(paste0(pop, ": Reading ldscores"))
    ldscores <- as_tibble(fread(loc_ldscores))
    # adds adjusted LD scores column if not done already
    if (!ldscores_already_ajdusted) {
      print(paste0(pop, ": AF-adjusting ldscores"))
      ldscores <- ldscores %>% adjust_LDscores()
      # saves table with adjusted LD scores to system
      print(paste0(pop, ": Writing adjusted ldscores to system"))
      fwrite(ldscores, loc_ldscores, sep="\t")
    }
    # joins LD scores to top_indep_SNPs
    print(paste0(pop, ": Joining ldscores to top_indep_SNPs"))
    top_indep_SNPs <- top_indep_SNPs %>%
      left_join(ldscores %>% select(CHR,BP, A0, A1, ld_score, ld_score_adj), 
                by=c("chrom"="CHR","chr_position"="BP", "ref"="A0","alt"="A1"))
    # corrects column names
    colnames(top_indep_SNPs) <- c(colnames(top_indep_SNPs)[1:(ncol(top_indep_SNPs)-2)],
                                  paste0("ld_score_",pop),paste0("ld_score_adj_",pop))
  }
  # writes LDscore-appended top_indep_SNPs back to system
  fwrite(top_indep_SNPs, loc_top_indep_SNPs, sep="\t")
}
# sees how many top SNPs for each pop is missing LD data within each pop
for (pop in pops) {
  col_ld <- paste0("ld_score_",pop)
  print(pop)
  top_indep_SNPs %>% group_by(pop) %>%
    summarize(prop_NA = sum(is.na(!!as.name(col_ld))) / n()) %>% print()
}
# Loops through each trait to compute LD-related metrics
for (i in 1:nrow(traits_table)) {
  # reads summary file
  code <- traits_table$prive_code[i]
  filename <- paste0(code,"_sf_indep.txt")
  loc_summary_file <- paste0(dir_sf,filename)
  sf <- as_tibble(fread(loc_summary_file))
  # filters top_indep_SNPs to just those for this trait
  top_trait <- top_indep_SNPs %>% filter(prive_code==code)
  # loops through each population
  for (pop in pops) {
    col_AF <- paste0("af_", pop)
    col_gvc <- paste0("gvc_", pop)
    # filters to only top SNPs for this population
    top_trait_pop <- top_trait %>% filter(pop == .GlobalEnv$pop)
    # appends LD scores to top SNPs summary file data
    sf_top_pop <- sf %>%
      filter(SNP %in% top_trait_pop$SNP) %>%
      rename(chrom = chr, chr_position = pos, gvc = !!enquo(col_gvc)) %>%
      left_join(top_trait_pop %>% select(chrom,chr_position,ref,alt,starts_with("ld_score_")),
                by=c("chrom","chr_position","ref","alt")) %>%
      select(-starts_with(c("af_","pval_","beta_","se_", "low_confidence_","neglog10_pval_")))
    # loops through both unadjusted and adjusted LD score calculations
    for (adj_status in c("","_adj")) {
      col_ldscore <- paste0("ld_score",adj_status,"_",pop)
      # renames columns and removes rows with NA LD scores
      sf_top_pop2 <- sf_top_pop %>%
        select(gvc,ld_score = !!as.name(col_ldscore)) %>% drop_na()
      n_LD_data <- nrow(sf_top_pop2)
      # computes gvc-weighted mean LD score for trait's top pop SNPs 
      traitLD <- weighted.mean(sf_top_pop2$ld_score, sf_top_pop2$gvc)
      # log10 transforms unadjusted traitLD score
      if (adj_status == "") {
        col_traitLD <- paste0("traitLD_unadj_",pop)
        #traitLD <- log10(traitLD) # no longer log10'ing unadjusted traitLD values
      } else {
        col_traitLD <- paste0("traitLD_adj_",pop)
      }
      # adds traitLD info to traits table
      traits_table[i,col_traitLD] <- traitLD
      print(paste(i, code, n_LD_data, pop, col_traitLD, round(traitLD,2)))
    }
  }
}
### computes meta-metrics related to traitLD across population
LD_table <- traits_table %>%
  select(prive_code, starts_with("traitLD_unadj_")) %>%
  pivot_longer(
    cols = paste0("traitLD_unadj_",pops),
    names_to = "pop",
    names_prefix = "traitLD_unadj_",
    values_to="traitLD_unadj") %>%
  group_by(prive_code) %>%
  summarize(traitLD_unadj_mean = mean(traitLD_unadj, na.rm=TRUE),
            traitLD_unadj_sd = sd(traitLD_unadj, na.rm=TRUE),
            traitLD_unadj_max = max(traitLD_unadj, na.rm=TRUE),
            traitLD_unadj_min = min(traitLD_unadj, na.rm=TRUE),
            traitLD_unadj_range = traitLD_unadj_max - traitLD_unadj_min,
            traitLD_unadj_ratio = traitLD_unadj_max / traitLD_unadj_min,
            traitLD_unadj_CoV = traitLD_unadj_sd / traitLD_unadj_mean)
# joins LD meta-metrics to traits_table
traits_table <- traits_table %>% left_join(LD_table, by="prive_code")
# looks at how trait groups (in quantitative traits) differ in LD meta-metrics
traits_table %>% filter(GWAS_trait_type=="quantitative") %>%
  group_by(group_consolidated) %>%
  summarize(traitLD_unadj_mean = mean(traitLD_unadj_mean, na.rm=TRUE),
            traitLD_unadj_range = mean(traitLD_unadj_range, na.rm=TRUE),
            traitLD_unadj_CoV = mean(traitLD_unadj_CoV, na.rm=TRUE))
# writes traits_table to system
fwrite(traits_table, loc_traits_table, sep="\t")
