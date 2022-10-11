# gini_PLR_vs_LDpred2

# compares Gini calculated from using PLR effect sizes and LDpred2 effect sizes

### Libraries and directories ####
library(tidyverse)
library(data.table)
library(ggrepel)
library(ggpubr)
source("../code_part1/helper_functions.R")

# sets working directory
setwd("./")

# sets location to traits_table.txt
loc_table <- "../generated_data/traits_table.txt"

# location to a file we generated that vastly speeds up the process of binning
# can be obtained from our github under ~/generated_data/
loc_chr_max_bps <- "../code_part1/chr_max_bps.txt"

# sets location to allele frequencies generated in previous PLINK script
dir_AFs <- "../generated_data/allele_frequencies/"

# contains list of LDpred2 effect sizes obtained from Prive's figshare
# must be unzipped and left as .csv: 'PGS-effects-ldpred2.csv'
loc_LDP_betas <- "../prive_data/PGS-effects-ldpred2.csv"

# sets directory of generated_data
dir_generated_data <- "../generated_data/"

# sets directory of generated figures to output to
dir_out <- "../generated_figures/"


### Code ####
traits_table <- as_tibble(fread(loc_table))

# copied code section from "../code_part1/append_betas_AFs.R"
# reads the 'PGS-effects-PLR.csv' file and extract the traits (codes)
betas_wide <- as_tibble(fread(loc_LDP_betas))
codes <- colnames(betas_wide)[6:length(colnames(betas_wide))]

# creates tibble with correctly formatted column names
betas <- tibble(
  rsid=as.character(),
  chr=as.numeric(),
  pos=as.numeric(),
  a0=as.character(),
  a1=as.character(),
  A1=as.character(),
  A2=as.character()
)
for (code in codes) {
  betas[,as.character(code)] <- as.numeric()
}
temp <- as_tibble(fread(paste0(dir_AFs,"pop_ALL_AFs_chr22.frq.strat")))
col_ancestries <- paste0("VarFreq_",levels(as.factor(temp$CLST)))
for (col_ancestry in col_ancestries) {
  betas[,col_ancestry] <- as.numeric()
}

# defines function for cases where the reference alleles are mismatched between
# PLINK's allele frequencies and Prive et al.'s betas
flip_AF <- function(AF) {return(1-AF)}

# appends ancestry-specific allele frequencies to list of all SNPs
for (chr_i in 1:22) {
  loc_AFs_chr <- paste0(dir_AFs,"pop_ALL_AFs_chr",chr_i,".frq.strat")
  
  AFs_chr <- as_tibble(fread(loc_AFs_chr, drop=c("CHR"))) %>%
    select(SNP, A1, A2, CLST, MAF) %>%
    pivot_wider(names_from=CLST,values_from=MAF,names_prefix="VarFreq_")
  
  betas_wide_chr <- betas_wide %>% filter(chr == chr_i) %>%
    left_join(AFs_chr, by=c("rsid"="SNP"))
  
  betas <- betas %>%
    add_row(betas_wide_chr %>% filter(a0==A1)) %>%
    add_row(betas_wide_chr %>%
              filter(a0!=A1) %>%
              mutate_at(col_ancestries,flip_AF))
  
  print(paste("Loaded allele frequencies for chromosome",chr_i))
}

# I recommend saving this file since it takes a while to generate
# (WARNING: LARGE FILE, approx 3.25 GB)
loc_out <- paste0(dir_generated_data,"LDpred2-betasAFs.txt")
write.table(betas,loc_out,sep="\t",quote=FALSE,row.names=FALSE)

# calculates ginis using LDpred2 effect sizes
LDP_ginis <- tibble(
  prive_code = as.character(),
  gini_United_LDP = as.numeric(),
)
bin_size <- 100000
col_AF <- "VarFreq_United"
bin_summary_method <- "sum"
threshold <- 100

for (code in codes) {
  if (!(code %in% traits_table$prive_code)) {next}
  sf <- betas[which(betas[code] != 0),] %>%
    select(chrom=chr,rsID=rsid,chr_position=pos,A1=a0,A2=a1,
           effect_weight=!!as.name(code),"VarFreq_United") %>%
    bin_snps(bin_size) %>% # bin_size
    select(effect_weight, !!as.name(col_AF), bin_ID) %>%
    drop_na() %>%
    get_h2("effect_weight", col_AF) %>%
    get_data_binned(bin_summary_method)
  n_significant_bins <- nrow(sf)
  sf_filtered <- sf %>% filter(rank <= threshold)
  
  if (threshold <= n_significant_bins) {h2_list <- sf_filtered$h2
  } else {h2_list <- c(rep(0, (threshold - n_significant_bins)), sf$h2)}
  gini <- get_gini(h2_list)
  
  LDP_ginis <- LDP_ginis %>% add_row(
    prive_code = code,
    gini_United_LDP = gini
  )
  print(paste("Calculated LDP gini for", code))
}

### comparing with gini PLR ####
traits_table <- traits_table %>% left_join(LDP_ginis, by="prive_code")
write.table(traits_table,loc_table,sep="\t",quote=FALSE,row.names=FALSE)

### Printing function ####
print_plot <- function(gg, loc_out, print_mode, plot_width, plot_height, sf) {
  if (print_mode == "png") {
    png(loc_out, width = plot_width*sf, height = plot_height*sf)
  } else if (print_mode == "pdf") {
    pdf(loc_out, width = plot_width*sf / 75, height = plot_height*sf / 75)
  }
  print(gg)
  dev.off()
}
sf <- 2
print_mode <- "pdf"

gini_p_theme <- theme(
  plot.subtitle = element_text(size=20*sf),
  axis.title = element_text(size=30*sf),
  axis.text = element_text(size=22*sf),
  legend.title = element_text(size=30*sf),
  legend.text = element_text(size=25*sf),
  plot.margin = unit(1*sf*c(1,1,1,1), "cm")
)
low_prevalence <- (traits_table %>% filter(prevalence < 0.01))$prive_code
traits_table2 <- traits_table %>% filter(!(prive_code %in% low_prevalence))

cor1 <- cor.test(traits_table2$gini_United, traits_table2$gini_United_LDP)
r <- cor1$estimate[[1]]
p <- cor1$p.value
gini_PLR_LDP_p <- ggplot(traits_table2, aes(x=gini_United,y=gini_United_LDP)) +
  geom_abline(slope=1,intercept=0, size = 1*sf) +
  geom_point(alpha = 0.75, size = 4*sf, aes(color=trait_type)) +
  scale_x_continuous(expand=c(0.01,0.01), limits=c(0,1)) +
  scale_y_continuous(expand=c(0.01,0.01), limits=c(0,1)) +
  xlab("Gini (using PLR effect sizes)") +
  ylab("Gini (using LDpred2 effect sizes)") +
  labs(color="Trait Type") +
  scale_color_manual(labels=c("Binary","Quantitative"),
                     breaks=c("binary", "quantitative"),
                     values = c("binary"="#F8766D", "quantitative"="#00BFC4")) +
  theme_light() +
  gini_p_theme +
  annotate("text", x=0.05, y=0.95, vjust=1, hjust=0, size = 10*sf,
           #label = paste0("r = ", round(r,4), "\np = ", formatC(p,format="E",digits=2)))
           label = paste0("r==", round(r,4)), parse=TRUE)
loc_out <- paste0(dir_out,"gini_PLR_vs_LDP.", print_mode)
print_plot(gini_PLR_LDP_p, loc_out, print_mode, 1200, 1000, sf)
print(loc_out)
