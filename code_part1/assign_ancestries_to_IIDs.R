# 1 - assign_ancestries_to_IIDs.R

# This script runs the ancestral grouping code found on the homepage of Prive's
# github page for the UKBB-PGS project: https://github.com/privefl/UKBB-PGS
# It produces a tsv file containing the IIDs of everyone in UKBB, their ancestry
# (according to Prive's method), and whether they self identify as white-British)

### Packages and directories ----

library(tidyverse)
library(data.table)

# sets location of phenotype file containing UKB IIDs and their first 16
# principal components (or more)
loc_PC <- "/directory/ukb_IID_16PCs.txt"
loc_output <- "/storage/coda1/p-jlachance6/0/shared/gini/UKB/pop_ALL_iids.txt"
# sets location of file with list of IIDs to remove (due to withdrawing consent
# from study) and of directory where ancestry-IID files will be created
loc_remove <- "/directory/IIDS_to_remove.csv"
dir_output <- "/directory/output/"


### Code ----

# reads PC data and selects for just first 16 PCs (if more are present)
raw_data <- as_tibble(loc_PC)
PC_UKBB <- raw_data %>% select(paste0("V", 26:41))

# Prive's code from github
all_centers <- read.csv(
  "https://raw.githubusercontent.com/privefl/UKBB-PGS/main/pop_centers.csv",
  stringsAsFactors = FALSE)
all_sq_dist <- apply(all_centers[-1], 1, function(one_center) {
  rowSums(sweep(PC_UKBB, 2, one_center, '-')^2)
})
thr_sq_dist <- max(dist(all_centers[-1])^2) * 0.002 / 0.16
group <- apply(all_sq_dist, 1, function(x) {
  grp <- NA
  ind <- which.min(x)
  if (isTRUE(x[ind] < thr_sq_dist)) {
    grp <- all_centers$Ancestry[ind]
    # We used a more stringent cutoff for the Ashkenazi group
    if (grp == "Ashkenazi" && x[ind] > 12.5^2) grp <- NA
  }
  grp
})
table(group, exclude = NULL)

# creates table containing IIDs and ancestry (keep order of original data)
pop <- tibble(IDD=raw_data$IID,
              ancestry=group)

# reads the list of IIDs to remove (who withdrew consent) and filters them out
IIDs_to_remove <- fread(loc_remove)$V1
pop <- pop %>% filter(!(IID %in% IIDs_to_remove))

# Produces file containing FID (which is just IID), IID, and their ancestry
pop <- pop %>% mutate(FID=IID) %>% select(FID,IID,ancestry)
# r stands for 'removed those who withdrew consent
# c stands for 'clean', meaning it contains FID, IID, and ancestry
loc_output <- paste0(dir_output,"pop_ALLrc_iids.txt")
write.table(pop,loc_output,row.names=FALSE,col.names=FALSE,quote=FALSE)