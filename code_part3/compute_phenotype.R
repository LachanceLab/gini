# 4 compute_phenotype.R

# Libraries and paths ####
library(tidyverse)
library(data.table)

dir_sims <- "../generated_data/simulations/"
loc_PCs <- "/storage/home/hcoda1/1/ncarvalho6/scratch/03-06_UKB/ukb_IID_16PCs.txt"

# Compute genetic_liabilities ####
dir_genliabs <- paste0(dir_sims,"genetic_liabilities/")

# loads simmed phenotype table
simphenos_tbl <- as_tibble(fread(paste0(dir_sims,"simphenos_tbl.txt")))

# loops through each chromosome
for (i in 1:22) {
  print(i)
  loc_gl <- paste0(dir_genliabs, "genetic_liabilities_chr",i,".sscore")
  chr_gl <- as_tibble(fread(loc_gl))
  
  if (i==1) {
    gl_matrix <- matrix(0, nrow=nrow(chr_gl), ncol=ncol(chr_gl)-2)
  }
  gl_matrix <- gl_matrix + as.matrix(chr_gl[1:nrow(chr_gl), 3:ncol(chr_gl)])
}

# makes genetic liabilities table
GLs <- as_tibble(chr_gl %>% select(FID=`#FID`,IID) %>% cbind(gl_matrix))
colnames(GLs)[-c(1:2)] <- paste0("gl_",simphenos_tbl$traitname)

# saves genetic liabilities
loc_out <- paste0(dir_sims, "genetic_liabilities.txt")
fwrite(GLs, loc_out, sep="\t")


# Compute phenotype ####
pheno <- GLs
for (i in 1:nrow(simphenos_tbl)) {
  print(i)
  # extracts trait genetic liabilities
  gl <- GLs[[2+i]]
  # generates environmental noise of variance ~= 1 - h2
  env_noise <- rnorm(nrow(pheno), mean=0, sd=sqrt(1 - var(gl)))
  
  # centers phenotype and adds environmental noise
  pheno[[2+i]] <- gl - mean(gl) + env_noise
  # variance of phenotype should be ~1
}
# renames columns
colnames(pheno)[-c(1:2)] <- paste0("ph_",simphenos_tbl$traitname)

# saves phenotype data
loc_out <- paste0(dir_sims, "pheno_values.txt")
fwrite(pheno, loc_out, sep="\t")


# makes PCs/covariate file
PCs <- as_tibble(fread(loc_PCs))
colnames(PCs) <- c("IID",paste0("PC",1:40))

PCs <- PCs %>% mutate(FID = IID) %>% select(FID, everything())
loc_out <- paste0(dir_sims, "covar_PCs.txt")
fwrite(PCs, loc_out, sep="\t")

# gets training cohort for GWAS
pop_all <- as_tibble(fread("../generated_data/pop_ALLrc_IIDs.txt", fill=TRUE)) %>% select(FID=V1,IID=V1,pop=V3)
pop_PGS <- as_tibble(fread("../generated_data/pop_sampled_IIDs.txt", fill=TRUE)) %>% select(FID=V1,IID=V1,pop=V3)
# filters to non-PGS UK individuals
pop_GWAS <- pop_all %>% filter(!(IID %in% pop_PGS$IID), pop == "United")
N <- 20000
pop_GWAS2 <- pop_GWAS[sample(1:nrow(pop_GWAS),N,replace = FALSE),]

loc_out <- paste0(dir_sims,"pop_GWAS2.txt")
fwrite(pop_GWAS2, loc_out, sep=" ", col.names = FALSE)