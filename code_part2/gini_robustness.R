# gini_robustness.R
# Calculates gini using different parameters to check robustness

## Libraries and directories ####
library(tidyverse)
library(data.table)
library(GGally)
library(ggpubr)
library(extrafont)
loadfonts(device = "win", quiet=TRUE)
source("../code_part1/helper_functions/helper_functions.R")

# Check whether robustness test should be redone or not
run_robustness_loop <- FALSE
# sets working directory
setwd("./")

# directory where panUKB GWAS summary statistics were saved to
dir_sf <- "../generated_data/panUKB_sf/"
# sets location to traits_table.txt
loc_table <- "../generated_data/traits_table.txt"

# reads traits_table and extract quantitative trait codes
traits_table <- as_tibble(fread(loc_table))
qcodes <- (traits_table %>% filter(GWAS_trait_type=="quantitative",
                                   PGS_trait_type=="quantitative"))$prive_code

if (run_robustness_loop) {
  ## Computes Gini at different parameters ####
  # parameter options
  thresholds <- c(50, 100, 250, 500, 750, 1000)
  # makes empty table for collecting data
  robustness_tbl <- tibble(
    prive_code = as.character(),
    pop = as.character(),
    threshold = as.numeric(),
    gini = as.numeric(),
    n_sig_SNPs = as.numeric(),
    n_sig_SNPs_pop = as.numeric(),
    sum_gvc_top = as.numeric(),
    sum_gvc_all = as.numeric()
  )
  
  ii <- 0
  # mega-loop
  for (code in qcodes) {
    # reads summary file
    #slice <- traits_table %>% filter(prive_code==code)
    filename <- paste0(code,"_sf_indep.txt")
    loc_summary_file <- paste0(dir_sf,filename)
    sf <- as_tibble(fread(loc_summary_file))
    # gets number of SNPs
    n_sig_SNPs <- nrow(sf)
    
    # loops through available gvc columns
    cols_gvc <- colnames(sf %>% select(starts_with("gvc_")))
    pops2use <- sapply(cols_gvc, function(x) substring(x, 5, nchar(x))) %>% unname()
    for (pop in pops2use) {
      col_gvc <- paste0("gvc_",pop)
      sf_pop <- sf %>%
        select(gvc_pop = !!enquo(col_gvc)) %>%
        arrange(-gvc_pop) %>%
        drop_na()
      # gets number of SNPs and sum_gvc (for this population)
      n_sig_SNPs_pop <- nrow(sf_pop)
      sum_gvc_all <- sum(sf_pop$gvc_pop)
      
      for (threshold in thresholds) {
        # gets top n SNPs, padding if necessary
        gvc_list <- (sf_pop %>% filter(row_number() <= threshold))$gvc_pop %>%
          pad_zeros(threshold)
        
        # calculates gini
        gini <- get_gini(gvc_list)
        # gets the sum of gvc among top SNPs
        sum_gvc_top <- sum(gvc_list)
        
        # records data
        robustness_tbl <- robustness_tbl %>% add_row(
          prive_code = code,
          pop = pop,
          threshold = threshold,
          gini = gini,
          n_sig_SNPs = n_sig_SNPs,
          n_sig_SNPs_pop = n_sig_SNPs_pop,
          sum_gvc_top = sum_gvc_top,
          sum_gvc_all = sum_gvc_all
        )
        
        print(paste(ii,code, pop, threshold, round(gini,3), sep="    "))
        ii <- ii + 1
      }
    }
  }
  robustness_tbl <- robustness_tbl %>%
    mutate(prop_gvc_top = sum_gvc_top / sum_gvc_all)
  # saves results
  loc_out <- paste0("../generated_data/gini_robustness_data.txt")
  fwrite(robustness_tbl, loc_out, sep="\t")
}

## VISUALIZATION ####
loc_out <- paste0("../generated_data/gini_robustness_data.txt")
robustness_tbl <- as_tibble(fread(loc_out))

### Printing function ####
print_plot <- function(gg, loc_out, print_mode, plot_width, plot_height, sf) {
  if (print_mode == "png") {
    png(loc_out, width = plot_width*sf, height = plot_height*sf)
  } else if (print_mode == "pdf") {
    #pdf(loc_out, width = plot_width*sf / 75, height = plot_height*sf / 75)
    cairo_pdf(file = loc_out,
              width = plot_width*sf / 75,
              height = plot_height*sf / 75)
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

# color legend for plots
gc_scale <- scale_color_manual(
  labels = c("Biological Measures","Lifestyle/Psychological","Physical Measures"),
  breaks = c("biological measures","lifestyle/psychological","physical measures"),
  values = c("biological measures"="#F8766D","lifestyle/psychological"="#00BF7D","physical measures"="#00B0F6","psychological"="#E76BF3")
)

### Actual plots ####

#### compares gini by population ####
# pop_threshold_tbl <- robustness_tbl %>% group_by(pop, threshold) %>%
#   dplyr::summarize(gini = mean(gini),
#             prop_gvc_top = mean(prop_gvc_top)) %>%
#   filter(!pop %in% c("meta","meta_hq"))
# # plots gini by population and threshold
# ggplot(pop_threshold_tbl, aes(x= as.factor(threshold), y=gini)) +
#   geom_line(aes(color = pop, group = pop), size=0.5*sf) +
#   xlab("Top # of SNPs used in Gini calculation") +
#   ylab("Gini") +
#   labs(title="Mean Gini for different top # of SNPs and populations",
#        color = "Population") +
#   theme_light()
# # plots prop_gvc_top by population and threshold
# ggplot(pop_threshold_tbl, aes(x=as.factor(threshold), y=prop_gvc_top)) +
#   geom_line(aes(color = pop, group = pop), size=0.5*sf) +
#   xlab("Top # of SNPs used in Gini calculation") +
#   ylab("Proportion of total gvc captured within top SNPs") +
#   labs(title="Mean Proportion of gvc in top SNPs for different top # of SNPs and populations",
#        color = "Population") +
#   theme_light()

#### compares prop_gvc_top against gini ####
prop_vs_gini_tbl <- robustness_tbl %>%
  filter(pop == "meta2use", threshold == 500) %>%
  left_join(traits_table[c("prive_code","group_consolidated")], by="prive_code")
gg <- ggplot(prop_vs_gini_tbl, aes(x=prop_gvc_top, y=gini)) +
  geom_point(aes(color=group_consolidated), alpha=0.75, size=6*sf) +
  xlim(0,1) + ylim(0,1) +
  coord_fixed() +
  labs(#title = expression(paste("Gini vs Proportion of ",italic(gvc)," among top SNPs for quantitative traits")),
    x = expression(paste("Proportion of total ", italic(gvc), " captured within top SNPs")),
    y = expression(paste(italic("G"[list(500,Meta)]))),
    color = "Trait Group") +
  theme_light() +
  gini_p_theme +
  theme(legend.position = "top",
        axis.title.y = element_text(family="Georgia")) +
  gc_scale
# prints out plot
loc_out <- paste0("../generated_figures/gini_vs_prop-gvc-top")
print(loc_out)
#print_plot(gg, loc_out, print_mode, 1200, 1200, sf)
print_plot(gg, paste0(loc_out,".png"), "png", 1200, 1200, sf)
print_plot(gg, paste0(loc_out,".pdf"), "pdf", 1200, 1200, sf)

#### gini at different thresholds ####
gini_rank_tbl <- robustness_tbl %>%
  filter(pop == "meta2use") %>%
  group_by(threshold) %>%
  mutate(gini_rank = rank(gini)) %>%
  left_join(traits_table[c("prive_code","group_consolidated")], by="prive_code")
gg <- ggplot(gini_rank_tbl, aes(x = as.factor(threshold), y = gini)) +
  geom_line(aes(group=prive_code, color=group_consolidated), size=0.75*sf, key_glyph = "rect") +
  labs(#title = 'Gini vs top # of SNPs for 96 quantitative traits,
       x = 'Number of SNPs used in Genomic Inequality (Gini) calculations',
       y = expression(paste(italic("G"[list(500,Meta)]))),
       color = "Trait Group") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_continuous(expand=c(0.01,0.01), limits=c(0,1)) +
  #scale_y_reverse() +
  theme_light() + gini_p_theme +
  theme(legend.position = "top",
        axis.title.y = element_text(family="Georgia")) +
  gc_scale
# prints out plot
loc_out <- paste0("../generated_figures/gini_vs_threshold")
print(loc_out)
#print_plot(gg, loc_out, print_mode, 1200, 1200, sf)
print_plot(gg, paste0(loc_out,".png"), "png", 1200, 1200, sf)
print_plot(gg, paste0(loc_out,".pdf"), "pdf", 1200, 1200, sf)

#### gini_pop vs gini_meta
# makes table of ginis
pop_ginis2 <- traits_table %>%
  select(prive_code, GWAS_trait_type, PGS_trait_type, group_consolidated, #portability_index,
         starts_with("gini_")) %>%
  filter(GWAS_trait_type == "quantitative", PGS_trait_type=="quantitative")

pops <- c("EUR","AMR","CSA","AFR","EAS") # "MID" ommitted

gini_UK_plots <- list()
for (i in 1:length(pops)) {
  # computes correlation
  pop2 <- pops[i]
  col_pop2 <- paste0("gini_",pop2)
  cor <- cor.test(pop_ginis2[["gini_panUKB"]], pop_ginis2[[col_pop2]])
  # makes actual plot
  p<-ggplot(data=pop_ginis2, aes(x = !!as.name(col_pop2), y = gini_panUKB, color=group_consolidated)) +
    geom_abline(slope=1,intercept=0, size=1*sf) +
    geom_point(alpha=0.75, size=4*sf) +
    scale_x_continuous(expand=c(0.01,0.01),limits=c(0,1)) +
    scale_y_continuous(expand=c(0.01,0.01),limits=c(0,1)) +
    xlab(bquote(italic(G[500][','][.(pop2)]))) +
    ylab(bquote(italic(G[500][','][Meta]))) +
    labs(color = "Trait Group") +
    coord_fixed() +
    theme_light() +
    gini_p_theme +
    gc_scale +
    theme(legend.position = "top",
          axis.title = element_text(family="Georgia")) +
    annotate("text", x=0.025, y=0.975,hjust=0,vjust=1, parse=TRUE, size=10*sf,
             label = paste0("r==",round(cor$estimate,3)))
  gini_UK_plots[[i]] <- p
}
ukginis <- ggarrange(plotlist = gini_UK_plots, ncol = 3, nrow = 2,
                     common.legend = TRUE,
                     legend.grob = get_legend(gini_UK_plots, position="top"))

loc_out <- paste0("../generated_figures/gini_by_pop")
print(loc_out)
#print_plot(ukginis, loc_out, print_mode, 3*600, 2*600, sf)
print_plot(ukginis, paste0(loc_out,".png"), "png", 3*600, 2*600, sf)
print_plot(ukginis, paste0(loc_out,".pdf"), "pdf", 3*600, 2*600, sf)