pheno <- commandArgs(trailingOnly = TRUE)[1]
#sink(paste0("filter_GWAS_independent.Rout"))

library(data.table)

pval_cutoff <- 0.0001
dir_sf <- "../generated_data/panUKB_sf/"

loc_sf <- paste0(dir_sf,pheno,"_sf.tsv")
loc_clumped <- paste0(dir_sf,pheno,"_sf.clumped")

full_sf <- fread(loc_sf)
full_clumped <- fread(loc_clumped)
independent_SNPs <- full_clumped$SNP

sf <- full_sf[full_sf$pval < pval_cutoff,]
indep_sf <- sf[sf$SNP %in% independent_SNPs]

loc_out <- paste0(dir_sf,pheno,"_sf_indep.txt")
write.table(indep_sf, loc_out, quote=FALSE, row.names=FALSE)
