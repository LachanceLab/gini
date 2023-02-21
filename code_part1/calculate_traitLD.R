# calculate_traitLD - calculates a weighted average LD score for SNPs associated with a trait


## libraries and directories

library(tidyverse)
library(data.table)


dir_generated_data <- "../generated_data/"
dir_sf <- paste0(dir_generated_data,"panUKB_sf/")
dir_input_data <- "../input_data/"
dir_ldscores <- paste0(dir_input_data, "ldscores/")

pops <- c("EUR","AFR","AMR","CSA","EAS","MID")
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
loc_traits_table <- paste0(dir_generated_data,"traits_table.txt")
traits_table <- as_tibble(fread(loc_traits_table))

loc_top_indep_SNPs <- paste0(dir_generated_data,"top_indep_SNPs.txt")
top_indep_SNPs <- as_tibble(fread(loc_top_indep_SNPs)) %>%
  filter(chrom != "X") %>%
  mutate(chrom = as.numeric(chrom))


if (!top_SNPs_already_appended) {
  for (pop in pops) {
    
    loc_ldscores <- paste0(dir_input_data, "ldscores/UKBB.",pop,".ldscores_FULL_adj.txt")
    
    print(paste0(pop, ": Reading ldscores"))
    ldscores <- as_tibble(fread(loc_ldscores))
    
    if (!ldscores_already_ajdusted) {
      print(paste0(pop, ": AF-adjusting ldscores"))
      ldscores <- ldscores %>% adjust_LDscores()
      
      print(paste0(pop, ": Writing adjusted ldscores to system"))
      fwrite(ldscores, loc_ldscores, sep="\t")
    }
    
    print(paste0(pop, ": Joining ldscores to top_indep_SNPs"))
    top_indep_SNPs <- top_indep_SNPs %>%
      left_join(ldscores %>% select(CHR,BP, A0, A1, ld_score, ld_score_adj), 
                by=c("chrom"="CHR","chr_position"="BP", "A1"="A0","A2"="A1"))
    
    colnames(top_indep_SNPs) <- c(colnames(top_indep_SNPs)[1:(ncol(top_indep_SNPs)-2)],
                                  paste0("ld_score_",pop),paste0("ld_score_adj_",pop))
    
  }
  
  #fwrite(top_indep_SNPs, loc_top_indep_SNPs, sep="\t")
}

for (i in 1:nrow(traits_table)) {
  code <- traits_table$prive_code[i]
  
  filename <- paste0(code,"_sf_indep.txt")
  loc_summary_file <- paste0(dir_sf,filename)
  sf <- as_tibble(fread(loc_summary_file)) %>%
    group_by(SNP) %>% filter(pval == min(pval)) %>% ungroup()
  
  top_indep_SNPs_trait <- top_indep_SNPs %>% filter(prive_code==code)
  
  sf_top <- sf %>%
    filter(SNP %in% top_indep_SNPs_trait$SNP) %>%
    rename(chrom = chr, chr_position = pos, A1 = ref, A2 = alt) %>%
    distinct() %>%
    group_by(SNP) %>% filter(gvc == max(gvc)) %>% ungroup() %>%
    filter(chrom != "X") %>%
    mutate(chrom = as.numeric(chrom)) %>%
    left_join(top_indep_SNPs_trait %>%
                select(chrom, chr_position, A1, A2, starts_with("ld_score_")),
              by=c("chrom","chr_position","A1","A2")) %>%
    select(-starts_with(c("af_","pval_","beta_","se_", "low_confidence_")))
  
  for (pop in pops) {
    for (adj_status in c("","_adj")) {
      col_ldscore <- paste0("ld_score",adj_status,"_",pop)
      sf_top_pop <- sf_top %>%
        select(gvc,ld_score = !!as.name(col_ldscore)) %>% drop_na()
      traitLD <- weighted.mean(sf_top_pop$ld_score, sf_top_pop$gvc)
      
      if (adj_status == "") {
        col_traitLD <- paste0("traitLD_unadj_",pop)
        traitLD <- log10(traitLD)
      } else {
        col_traitLD <- paste0("traitLD_adj_",pop)
      }
      
      traits_table[i,col_traitLD] <- traitLD
      print(paste(i, code, pop, col_traitLD, round(traitLD,2)))
    }
  }
  
}
fwrite(traits_table, loc_traits_table, sep="\t")


#traits_table <- traits_table %>% select(-starts_with("traitLD_"))