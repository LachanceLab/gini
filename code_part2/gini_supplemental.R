# gini_supplemental

# contains code pertaining to supplemental analyses of the Gini coefficient

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

# sets directory of betas with appended allele frequencies
dir_betas <- "../generated_data/betas_and_AFs/"

# sets directory of generated_data
dir_generated_data <- "../generated_data/"

# sets directory of generated figures to output to
dir_out <- "../generated_figures/"


### Code ####
traits_table <- as_tibble(fread(loc_table))

# loads a file that contains the max base pair position for each chromosome
chr_max_bps <- as_tibble(fread(loc_chr_max_bps))

# will loop through each of these
codes <- traits_table$prive_code
ancestries <- sapply(str_split(colnames(traits_table)[substring(colnames(traits_table),1,5)=="gini_"],"gini_"),"[",2)
bin_sizes <- c(1, 1E5)
bin_summary_methods <- c("sum","max")
thresholds <- c(10, 25, 50, 100, 200, 400, 800, 1600)

# used to estimate the progress bar in the console output
max_i <- length(codes) * length(ancestries) * length(bin_sizes) * length(bin_summary_methods) * (length(thresholds)+1)

# Helper function
record_gini <- function(ginis, i, prive_code, n_snps, ancestry, bin_size, n_significant_bins, bin_summary_method, threshold, threshold_zero_padding, gini) {
  ginis <- ginis %>% add_row(
    prive_code = prive_code,
    n_snps = n_snps,
    ancestry = ancestry,
    bin_size = bin_size,
    n_significant_bins = n_significant_bins,
    bin_summary_method = bin_summary_method,
    threshold = threshold,
    threshold_zero_padding = threshold_zero_padding,
    gini = gini
  )
  pct_done <- (i / max_i) * 100
  print(paste(prive_code, bin_size, ancestry, bin_summary_method, threshold, threshold_zero_padding,
              ":: About",round(pct_done,3),"% done."))
  
  return(ginis)
}

ginis <- tibble(
  prive_code = as.character(),
  n_snps = as.numeric(),
  ancestry = as.character(),
  bin_size = as.numeric(),
  n_significant_bins = as.numeric(),
  bin_summary_method = as.character(),
  threshold = as.numeric(),
  threshold_zero_padding = as.logical(),
  gini = as.numeric()
)

i <- 0
# The Hypernested Loop
for (code in codes) {
  
  loc_betas <- paste0(dir_betas,code,"-betasAFs.txt")
  if (!file.exists(loc_betas)) {next}
  
  print(paste("Reading summary file for",code))
  sf <- as_tibble(fread(loc_betas))
  n_snps <- nrow(sf)
  
  for (bin_size in bin_sizes) {
    
    if (bin_size > 1) {sf_binID <- bin_snps(sf, bin_size)}
    else {sf_binID <- sf %>% mutate(bin_ID = row_number())}
    
    for (ancestry in ancestries) {
      
      col_AF <- paste0("VarFreq_", ancestry)
      sf_h2 <- sf_binID %>%
        select(effect_weight, !!as.name(col_AF), bin_ID) %>%
        drop_na() %>%
        get_h2("effect_weight", col_AF)
      
      for (bin_summary_method in bin_summary_methods) {
        
        temp <- sf_h2
        if (bin_size > 1) {sf_binned <- get_data_binned(sf_h2, bin_summary_method)}
        else {sf_binned <- sf_h2}
        n_significant_bins <- nrow(sf_binned)
        
        threshold_zero_padding <- FALSE
        
        for (threshold in thresholds) {
          
          if (threshold < n_significant_bins) {
            sf_filtered <- sf_binned %>% filter(rank <= threshold)
            h2_list <- sf_filtered$h2
            gini <- get_gini(h2_list)
            i <- i + 1
            ginis <- record_gini(ginis, i, code, n_snps, ancestry, bin_size, n_significant_bins, bin_summary_method, threshold, threshold_zero_padding, gini)
          }
          
          if (threshold > n_significant_bins & !threshold_zero_padding) {
            threshold_zero_padding <- TRUE
            h2_list <- c(rep(0, (threshold - n_significant_bins)), sf_binned$h2)
            gini2 <- get_gini(h2_list)
            
            h2_list <- sf_binned$h2
            gini1 <- get_gini(h2_list)
            i <- i + 1
            ginis <- record_gini(ginis, i, code, n_snps, ancestry, bin_size, n_significant_bins, bin_summary_method, n_significant_bins, FALSE, gini1)
            i <- i + 1
            ginis <- record_gini(ginis, i, code, n_snps, ancestry, bin_size, n_significant_bins, bin_summary_method, threshold, TRUE, gini2)
          } else if (threshold > n_significant_bins & threshold_zero_padding) {
            h2_list <- c(rep(0, (threshold - n_significant_bins)), sf_binned$h2)
            gini <- get_gini(h2_list)
            i <- i + 1
            ginis <- record_gini(ginis, i, code, n_snps, ancestry, bin_size, n_significant_bins, bin_summary_method, threshold, threshold_zero_padding, gini)
          }
        }
      }
    }
  }
}

# Uncomment if you want to save results for future analysis
# loc_out <- paste0(dir_generated_data,"ginis_extra.txt")
# write.table(ginis,loc_out,sep="\t",row.names=FALSE, quote=FALSE)

#### VISUALIZING DATA ####

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
print_mode <- "png"

# common theme for supplemental plots
gini_p_theme <- theme(
  plot.subtitle = element_text(size=20*sf),
  axis.title = element_text(size=30*sf),
  axis.text = element_text(size=22*sf),
  legend.title = element_text(size=30*sf),
  legend.text = element_text(size=25*sf),
  plot.margin = unit(1*sf*c(1,1,1,1), "cm")
)

# reads supplemental ginis generated (if saved to system)
ginis <- as_tibble(fread("../generated_data/ginis_extra.txt"))
# reads traits table and gets list of low prevalence binary traits
traits_table <- as_tibble(fread(loc_table))
low_prevalence <- (traits_table %>% filter(prevalence < 0.01))$prive_code

# converts variables to factors and makes unique character string for parameter combination
ginis <- ginis %>%
  filter(!(prive_code %in% low_prevalence)) %>%
  mutate(
    bin_size = as.factor(bin_size),
    bin_summary_method = as.factor(bin_summary_method),
    settings = paste(bin_size, bin_summary_method, threshold, threshold_zero_padding),
)

# filters supplemental ginis to allow bin size comparison
ginis_binsize <- ginis %>% filter(
  bin_summary_method == "sum",
  threshold_zero_padding = TRUE,
  threshold == 100,
  ancestry == "United"
) %>%
  select(prive_code, ancestry, bin_size, gini) %>%
  pivot_wider(
    names_from = bin_size,
    values_from = gini,
    names_prefix = "gini_"
  ) %>%
  left_join(traits_table %>% select(prive_code,short_label), by="prive_code")

# gets correlation between bin size 1 gini and bin size 100k gini
cor1 <- cor.test(ginis_binsize$gini_1, ginis_binsize$`gini_1e+05`)
r <- cor1$estimate[[1]]
# plots bin size comparison
gini_bz_p <- ggplot(ginis_binsize, aes(x=gini_1,y=`gini_1e+05`)) +
  geom_abline(slope=1,intercept=0, size = 1*sf) +
  geom_point(alpha = 0.75, size = 4*sf) +
  scale_x_continuous(expand=c(0.01,0.01), limits=c(0,1)) +
  scale_y_continuous(expand=c(0.01,0.01), limits=c(0,1)) +
  xlab("Gini (bin size = 1; i.e. no binning)") +
  ylab("Gini (bin size = 100kb)") +
  theme_light() +
  gini_p_theme +
  annotate("text", x=0.05, y=0.95, vjust=1, hjust=0, size = 10*sf,
           label = paste0("r==", round(r,4)), parse=TRUE)
# prints bin size comparison plot
loc_out <- paste0(dir_out,"gini_bz_scatterplot.", print_mode)
print_plot(gini_bz_p, loc_out, print_mode, 1000, 1000, sf)
print(loc_out)

# comparing rank order by threshold
ginis_t <- ginis %>%
  filter(ancestry=="United", bin_size == 100000, bin_summary_method == "sum", (threshold_zero_padding == TRUE)|(n_significant_bins > threshold) ) %>%
  filter(threshold %in% c(50,100,200,400,800)) %>%
  select(prive_code, threshold, gini, n_significant_bins) %>%
  group_by(threshold) %>%
  mutate(gini_t_rank = order(order(gini, decreasing=FALSE))) %>%
  left_join(traits_table %>% select(prive_code,group_consolidated), by="prive_code")

# sets color legend
gc_scale <- list(
  labels = c("Biological Measures","Diseases","Lifestyle/Psychological","Physical Measures"),
  breaks = c("biological measures","diseases","lifestyle/psychological","physical measures"),
  values = c("biological measures"="#F8766D", "diseases"="#A3A500","lifestyle/psychological"="#00BF7D","physical measures"="#00B0F6","psychological"="#E76BF3")
)

# plots gini rank comparison by threshold
gini_t_p1 <- ggplot(ginis_t,
       aes(x=as.factor(threshold), y=gini_t_rank, group=prive_code)) +
  geom_line(aes(color = group_consolidated), size=0.75*sf, key_glyph = "rect") +
  scale_x_discrete(expand=c(0.01,0.01)) +
  scale_y_continuous(expand=c(0.01,0.01)) +
  scale_color_manual(labels=gc_scale[["labels"]],
                     breaks=gc_scale[["breaks"]],
                     values=gc_scale[["values"]]) +
  xlab("Threshold") +
  ylab("Trait Rank by Gini (low rank = low gini)") +
  labs(color = "Trait Group") +
  theme_light() +
  gini_p_theme

# plots gini comparison by threshold
gini_t_p2 <- ggplot(ginis_t,
       aes(x=as.factor(threshold), y=gini, group=prive_code)) +
  geom_line(aes(color = group_consolidated), size=0.75*sf, key_glyph = "rect") +
  scale_x_discrete(expand=c(0.01,0.01)) +
  scale_y_continuous(expand=c(0.01,0.01), limits=c(0,1)) +
  scale_color_manual(labels=gc_scale[["labels"]],
                     breaks=gc_scale[["breaks"]],
                     values=gc_scale[["values"]]) +
  xlab("Threshold") +
  ylab("Gini") +
  labs(color = "Trait Group") +
  theme_light() +
  gini_p_theme

gini_t_plot <- ggarrange(plotlist = list(gini_t_p1, gini_t_p2),
          common.legend = TRUE, legend="bottom", labels="AUTO",
          font.label = list(size=32*sf))
loc_out <- paste0(dir_out,"gini_threshold_rank.", print_mode)
print_plot(gini_t_plot, loc_out, print_mode, 2*1000, 1000, sf)
print(loc_out)

## comparing ginis by ancestry to UK
ginis_pop <- ginis %>% filter(
  threshold == 100,
  (threshold_zero_padding == TRUE | threshold < n_significant_bins),
  bin_size == 100000,
  bin_summary_method == "sum"
) %>% select(prive_code, ancestry, gini) %>%
  pivot_wider(
    names_from = ancestry,
    values_from = gini
  )
pop_centers <- read.csv("https://raw.githubusercontent.com/privefl/UKBB-PGS/main/pop_centers.csv",
                        stringsAsFactors = FALSE)
distances <- tibble(
  ancestry = pop_centers$Ancestry,
  prive_dist_to_UK = as.matrix(dist(pop_centers[2:17]))[,1]
) %>% arrange(prive_dist_to_UK)
distances$ancestry <- str_replace(distances$ancestry,"United Kingdom","UK")
ancestries_sorted <- distances$ancestry[c(2:4,6:9)]
gini_UK_plots <- list()
for (i in 1:length(ancestries_sorted)) {
  pop2 <- ancestries_sorted[i]
  cor <- cor.test(ginis_pop[,"United"][[1]], ginis_pop[,pop2][[1]])
  p<-ggplot(data=ginis_pop, aes(x = !!as.name(pop2), y = United)) +
    geom_abline(slope=1,intercept=0, size=1*sf) +
    geom_point(alpha=0.75, size=4*sf) +
    scale_x_continuous(expand=c(0.01,0.01),limits=c(0,1)) +
    scale_y_continuous(expand=c(0.01,0.01),limits=c(0,1)) +
    xlab(bquote(Gini[100][','][.(pop2)])) +
    ylab(bquote(Gini[100][','][UK])) +
    theme_light() +
    theme(aspect.ratio = 1) +
    annotate("text", x=0.05, y=0.95,hjust=0,vjust=1, parse=TRUE, size=10*sf,
             label = paste0("r==",round(cor$estimate[[1]],4))) +
    gini_p_theme
  gini_UK_plots[[i]] <- p
}
ukginis <- ggarrange(plotlist = gini_UK_plots, ncol = 4, nrow = 2)

loc_out <- paste0(dir_out,"gini_ancestry.", print_mode)
print_plot(ukginis, loc_out, print_mode, 4*500, 2*500, sf)
print(loc_out)
