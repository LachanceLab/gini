# plot_lorenz_curve.R

# This script takes the summary files of a trait, calculates SNP gvc, bins SNPs,
# and plots a Lorenz curve for the top 100 bins. Loops for all traits and makes
# Figure 2 consisting of two Lorenz curves.

### Libraries and directories ####
library(tidyverse)
library(ggpubr)
library(data.table)
source("../code_part1/helper_functions.R")

# sets working directory
setwd("./")

# Sets the directory of the summary files with appended allele frequencies
dir_sfs <- "../generated_data/betas_and_AFs/"
# sets the location of the traits table
loc_table <- "../generated_data/traits_table.txt"
# location to a file we generated that vastly speeds up the process of binning
# can be obtained from our github under ~/generated_data/
loc_chr_max_bps <- "../code_part1/chr_max_bps.txt"
# sets directory of outputted figures
dir_out <- "../generated_figures/"

### Code ####

# reads traits table
traits_table <- as_tibble(fread(loc_table))

# loads a file that contains the max base pair position for each chromosome
chr_max_bps <- as_tibble(fread(loc_chr_max_bps))

# function that reads a trait's summary file and extracts the needed columns
cleanup_data_lorenz <- function(code, ancestry="United", threshold=100, threshold_padding=TRUE, bin_size=100000, bin_summary_method="sum") {
  col_AF <- paste0("VarFreq_",ancestry)
  
  loc_sf <- paste0(dir_sfs,code,"-betasAFs.txt")
  sf <- as_tibble(fread(loc_sf)) %>%
    bin_snps(bin_size) %>%
    get_h2("effect_weight", col_AF) %>%
    get_data_binned(bin_summary_method) %>%
    arrange(-h2) %>%
    filter(row_number() <= threshold) %>%
    arrange(h2) %>%
    mutate(h2_csum = cumsum(h2)) %>%
    select(h2,h2_csum)
  sf$h2_cshare <- sf$h2_csum / sf$h2_csum[nrow(sf)]
  if ((nrow(sf) < threshold) & (threshold_padding)) {
    sf <- sf %>%
      add_row(
        h2 = rep(0, threshold - nrow(sf)),
        h2_csum = rep(0, threshold - nrow(sf)),
        h2_cshare = rep(0, threshold - nrow(sf)),
        ) %>% arrange(h2)
  }
  sf <- sf %>% mutate(percentile = row_number() / threshold)
  sf
}
# function that actually plots Lorenz curve
plot_lorenz <- function(code, sf, ancestry="United") {
  
  slice <- traits_table %>% filter(prive_code == code)
  description <- slice$description
  gini <- slice[1,paste0("gini_",ancestry)]
  
  if (ancestry=="United") {ancestry <- "UK"}
  
  title <- paste0("Lorenz Curve: ", description)
  subtitle <- paste0("Gini (",ancestry,") = ", round(gini,3))
  
  gg<-ggplot(sf, aes(x=100*percentile, y=h2_cshare)) +
    geom_col(position = position_nudge(-0.5), fill="gray20", width=0.8) +
    geom_abline(slope=1/100,color="dodgerblue1") +
    labs(title=title, subtitle=subtitle) +
    xlab(paste("Percentile of bin genetic variance contribution within top", threshold,"bins")) +
    ylab(paste("Proportion of top", threshold,"bins' genetic variance contribution")) +
    theme_light() +
    theme(aspect.ratio = 1) +
    xlim(0,100) +
    ylim(0,1) +
    scale_x_continuous(expand = c(0.0,0.0)) +
    scale_y_continuous(expand = c(0,0))
  gg
}
# settings used for plotting. Deviation from these (other than ancestry) will
# lead to the displayed Gini score on the plot being incorrect
ancestry <- "United"
threshold <- 100
threshold_padding <- TRUE
bin_size <- 100000
bin_summary_method <- "sum"

# plot save settings for each plot (in pixels)
width <- 2200
height <- 2200

# Plots all traits' divergence plot and saves to [dir_out]/divergence_plots/
dir.create(paste0(dir_out,"lorenz_curves/"), recursive = TRUE)
for (code in traits_table$prive_code) {
  sf <- cleanup_data_lorenz(code, ancestry, threshold, threshold_padding, bin_size, summary_method)
  gg <- plot_lorenz(code, sf, ancestry)
  
  loc_out <- paste0(dir_out,"lorenz_curve/lorenz_curve_",code,".png")
  ggsave(loc_out,plot=gg,width=width,height=height,units="px")
  print(paste("Saved lorenz curve for",code))
}

# Makes Figure 2: Lorenz Curve
figure_codes <- c("geek_time","277.4")
plots <- list()
for (i in 1:length(figure_codes)) {
  code <- figure_codes[i]
  gg <- plot_divergence(code)
  plots[[i]] <- gg
}

gg_fig <- ggarrange(plots[[1]],plots[[2]])

loc_out <- paste0(dir_out,"figure_2_lorenz_curves.png")
ggsave(loc_out,plot=gg_fig,width=2*width,height=height,units="px")
