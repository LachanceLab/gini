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
get_true_betas <- function(causal_SNPs, Mc, h2, seed=as.numeric(NA)) {
  if (!is.na(seed)) {set.seed(seed)}
  causal_SNPs <- causal_SNPs %>%
    mutate(beta_true = as.numeric(NA),
           gvc_EUR = as.numeric(NA))
  
  for (i in 1:Mc) {
    # adapted from: https://www.biorxiv.org/content/10.1101/2022.12.29.522270v2.full
    AF <- causal_SNPs$AF_EUR[i]
    prior_beta_var <- h2 / (2 * AF * (1-AF) * Mc)
    
    
    # draws from normal distribution
    beta_true <- rnorm(1, mean=0, sd=sqrt(prior_beta_var))
    causal_SNPs$beta_true[i] <- beta_true
    causal_SNPs$gvc_EUR[i] <- 2 * beta_true^2 * AF * (1-AF)
  }
  #print(sum(causal_SNPs$gvc_EUR))
  
  return(causal_SNPs)
}

# Code ####
LD_all <- as_tibble(fread(paste0(dir_sims,"LD_all_SNPs.txt")))

# simulation parameters
set.seed(1)
threshold <- 500
Mcs <- c(10,100,500,1000,5000)
h2s <- c(0.1, 0.3, 0.5)
trials <- 5

all_causal_SNPs <- LD_all[0,"varid"]

simphenos_tbl <- tibble(ID = as.numeric(),
                        Mc = as.numeric(),
                        h2 = as.numeric(),
                        trial = as.numeric(),
                        traitname = as.character())
# generates true betas
for (Mc in Mcs) {
  for (h2 in h2s) {
    for (i in 1:trials) {
      print(paste(Mc, h2, i))
      
      causal_SNPs <- LD_all %>%
        sample_causal_SNPs(Mc) %>%
        get_true_betas(Mc, h2) %>%
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

# saves simulation phenotypes table
simphenos_tbl$ID <- 1:nrow(simphenos_tbl)
loc_out <- paste0(dir_sims, "simphenos_tbl.txt")
fwrite(simphenos_tbl, loc_out, sep="\t")