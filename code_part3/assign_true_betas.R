# 2 assign_true_betas.R

# Libraries and paths ####
library(tidyverse)
library(data.table)
source("../code_part1/helper_functions/helper_functions.R")

dir_sims <- "../generated_data/simulations/"

# Functions ####
sample_causal_SNPs <- function(LD_all, Mc, seed=as.numeric(NA)) {
  if (!is.na(seed)) {set.seed(seed)}
  
  causal_i <- sample(1:nrow(LD_all), Mc, replace=FALSE)
  causal_SNPs <- LD_all[causal_i,] %>% arrange(CHR, BP, A0, A1)
  
  return(causal_SNPs)
}
get_true_betas <- function(causal_SNPs, Mc, h2, shape, seed=as.numeric(NA)) {
  if (!is.na(seed)) {set.seed(seed)}
  causal_SNPs <- causal_SNPs %>%
    mutate(beta_true = as.numeric(NA),
           gvc_EUR = as.numeric(NA))
  
  scale <- h2 / (Mc * shape)
  
  for (i in 1:Mc) {
    # adapted from: https://www.biorxiv.org/content/10.1101/2022.12.29.522270v2.full
    AF <- causal_SNPs$AF_EUR[i]
    causal_SNPs$beta_true[i] <- sqrt( rgamma(1, shape = shape, scale = scale) /
                                        (2 * AF * (1-AF)) )
  }
  causal_SNPs <- causal_SNPs %>%
    mutate(beta_true = beta_true * sample(c(-1,1),Mc, replace=TRUE),
           gvc_EUR = 2 * beta_true^2 * AF_EUR * (1-AF_EUR))
  #print(sum(causal_SNPs$gvc_EUR))
  
  return(causal_SNPs)
}

# Code ####
LD_all <- as_tibble(fread(paste0(dir_sims,"LD_all_SNPs.txt")))

# simulation parameters
set.seed(1)
threshold <- 500
Mcs <- c(50,100,500,1000,5000)
h2s <- c(0.1, 0.3, 0.5)
shape <- 0.5
trials <- 5

all_causal_SNPs <- LD_all[0,"varid"]

simphenos_tbl <- tibble(ID = 0, Mc = 0, h2 = 0, trial = 0,
                        traitname = "", gini_true = 0)[0,]
# generates true betas
for (Mc in Mcs) {
  for (h2 in h2s) {
    for (i in 1:trials) {
      print(paste(Mc, h2, i))
      
      causal_SNPs <- LD_all %>%
        sample_causal_SNPs(Mc) %>%
        get_true_betas(Mc, h2, shape) %>%
        select(varid, beta_true)

      col_beta <- paste0("tb_h",h2*10,"_Mc",Mc,"_t",i)
      colnames(causal_SNPs)[ncol(causal_SNPs)] <- col_beta


      all_causal_SNPs <- all_causal_SNPs %>%
        full_join(causal_SNPs, by="varid")
      
      simphenos_tbl <- simphenos_tbl %>% add_row(
        Mc = Mc, h2 = h2, trial = i,
        traitname = substring(col_beta,4)
      )
    }
  }
}

# sets non-causal SNP effect sizes to zero
all_causal_SNPs[is.na(all_causal_SNPs)] <- 0
# joins LD data
all_causal_SNPs <- LD_all %>% inner_join(all_causal_SNPs, by="varid")

# saves true betas
loc_out <- paste0(dir_sims, "all_true_betas.txt")
fwrite(all_causal_SNPs, loc_out, sep="\t")

# computes true gini
AFs <- all_causal_SNPs$AF_EUR
for (i in 1:nrow(simphenos_tbl)) {
  # extracts betas and computes gvcs
  betas <- all_causal_SNPs[[11+i]]
  gvcs <- 2 * betas^2 * AFs * (1 - AFs)
  
  gvcs <- gvcs[gvcs != 0]
  # pads with zeros if necessary
  gvc_list <- pad_zeros(gvcs, threshold)
  # computes gini
  gini <- get_gini(gvc_list)
  
  simphenos_tbl$gini_true[i] <- gini
}

# saves simulation phenotypes table
simphenos_tbl$ID <- 1:nrow(simphenos_tbl)
loc_out <- paste0(dir_sims, "simphenos_tbl.txt")
fwrite(simphenos_tbl, loc_out, sep="\t")