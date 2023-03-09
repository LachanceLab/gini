# 6 - calculate_divergence.R

# Combines the chromosome-specific PGS data from calculate_divergence.sh and
# computes an ANOVA by ancestry on PGS for each trait

### Libraries and directories ####
library(tidyverse)
library(data.table)

# sets working directory
setwd("./")

# sets location of traits_table
loc_table <- "../generated_data/traits_table.txt"
# sets directory for generated_data 
dir_generated_data <- "../generated_data/"

loc_IIDs <- paste0(dir_generated_data,"pop_ALLrc_IIDs.txt")
IIDs <- as_tibble(fread(loc_IIDs, fill=TRUE)) %>% select(FID=V1,IID=V1,pop=V3)

for (i in 1:22) {
  print(i)
  loc_PGS <- paste0(dir_generated_data, "polygenic_scores/ALL_traits-PGS_chr",i,".sscore")
  chr_PGS <- as_tibble(fread(loc_PGS))
  
  if (i==1) {
    PGS_matrix <- matrix(0, nrow=nrow(chr_PGS), ncol=ncol(chr_PGS)-2)
  }
  PGS_matrix <- PGS_matrix + as.matrix(chr_PGS[1:nrow(chr_PGS), 3:ncol(chr_PGS)])
}
PGSs <- as_tibble(chr_PGS %>% select(FID=`#FID`,IID) %>%
                    left_join(IIDs, by=c("FID","IID")) %>%
                    cbind(PGS_matrix)) %>% filter(pop != "Ashkenazi")
colnames(PGSs)[-c(1:3)] <- sapply(colnames(PGSs)[-c(1:3)],
                                  function(x) substring(x, 1, nchar(x)-4)) %>% unname()

loc_out <- paste0(dir_generated_data, "pop_sampled_PGSs.txt")
fwrite(PGSs, loc_out, sep="\t")

traits_table <- as_tibble(fread(loc_table))
pop_vec <- PGSs[["pop"]]
for (i in 1:nrow(traits_table)) {
  code <- traits_table$prive_code[i]
  print(code)
  PGS_vec <- PGSs[[code]]
  aov1 <- aov(PGS_vec ~ pop_vec, data=NULL)
  f_stat <- summary(aov1)[[1]][1,4]
  p_value <- summary(aov1)[[1]][1,5]
  
  traits_table$log_F[i] <- log10(f_stat)
}

fwrite(traits_table, loc_table, sep="\t")