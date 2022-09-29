# generate_extra_ginis.R

# An early version of this script was used to generate extra ginis using
# different settings in the gini calculation in order to compare their results.
# This data was used to select an ideal combination of settings

### Libraries and directories ####
library(tidyverse)
library(data.table)
library(ggrepel)
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


### Code ###
traits_table <- as_tibble(fread(loc_table))

# loads a file that contains the max base pair position for each chromosome
chr_max_bps <- as_tibble(fread(loc_chr_max_bps))

# will loop through each of these
codes <- traits_table$prive_code
ancestries <- sapply(str_split(colnames(traits_table)[substring(colnames(traits_table),1,5)=="gini_"],"gini_"),"[",2)
bin_sizes <- c(1, 1E4, 1E5, 1E6)
bin_summary_methods <- c("sum","max")
thresholds <- c(10, 25, 50, 100, 150, 250, 500, 1000, 5000, 20000, 50000, 200000)

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
ginis <- as_tibble(fread("../generated_data/ginis_extra.txt"))

ginis <- ginis %>% mutate(
  bin_size = as.factor(bin_size),
  bin_summary_method = as.factor(bin_summary_method),
  settings = paste(bin_size, bin_summary_method, threshold, threshold_zero_padding),
)

# Overall
ginis_grouped <- ginis %>% filter(threshold %in% thresholds) %>%
  group_by(bin_size, bin_summary_method, threshold, threshold_zero_padding) %>%
  summarize(mean_gini = mean(gini),
            sd_gini = sd(gini),
            n = n()) %>%
  mutate(
    settings = paste(bin_size, bin_summary_method, threshold, threshold_zero_padding)
    )
summary(ginis_grouped$sd_gini)
nrow(ginis_grouped %>% filter(sd_gini >= 0.2))
summary((ginis_grouped %>% filter(sd_gini >= 0.2))$n)
ginis2 <- ginis %>% filter(settings %in% (ginis_grouped %>% filter(sd_gini >= 0.2))$settings)

ggplot(ginis2, aes(x = gini)) +
  geom_histogram(fill = "dodgerblue1") +
  scale_x_continuous(limits=c(0,1), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Gini") +
  ylab("Count") +
  labs(title = "Distribution of Ginis",
       subtitle = "Filtered to settings combinations with standard deviation >= 0.2") +
  theme_light()

## Bin Size ##
ginis2 %>% group_by(bin_size) %>%
  summarize(mean_gini = mean(gini),
            sd_gini = sd(gini),
            n = n())
summary(aov(data=ginis2, gini ~ bin_size))
ggplot(ginis2, aes(x=gini, y=bin_size, fill=bin_size)) +
  geom_violin() +
  stat_summary(fun=mean, geom="point") +
  stat_summary(fun.data=mean_sdl) +
  xlab("Gini") +
  ylab("Bin Size") +
  labs(title = "Distribution of Ginis by Bin Size",
       subtitle = "Filtered to settings combinations with standard deviation >= 0.2") +
  theme_light() +
  theme(legend.position = "none")

## Bin Summary Method ##
ginis2 %>% group_by(bin_summary_method) %>%
  summarize(mean_gini = mean(gini),
            sd_gini = sd(gini),
            n = n())
summary(aov(data=ginis2, gini ~ bin_summary_method))
ggplot(ginis2, aes(x=gini, y=bin_summary_method, fill=bin_summary_method)) +
  geom_violin() +
  stat_summary(fun=mean, geom="point") +
  stat_summary(fun.data=mean_sdl) +
  xlab("Gini") +
  ylab("Bin Summary Method") +
  labs(title = "Distribution of Ginis by Bin Summary Method",
       subtitle = "Filtered to settings combinations with standard deviation >= 0.2") +
  theme_light() +
  theme(legend.position = "none")

## Threshold ##
ginis2 %>%
  group_by(as.factor(threshold)) %>%
  summarize(mean_gini = mean(gini),
            sd_gini = sd(gini),
            n = n())
summary(lm(ginis2$gini ~ log2(ginis2$threshold)))
ggplot(ginis2, aes(x=gini, y=as.factor(threshold), fill=as.factor(threshold))) +
  geom_violin() +
  stat_summary(fun=mean, geom="point") +
  stat_summary(fun.data=mean_sdl) +
  xlab("Gini") +
  ylab("Threshold") +
  labs(title = "Distribution of Ginis by Threshold",
       subtitle = "Filtered to settings combinations with standard deviation >= 0.2") +
  theme_light() +
  theme(legend.position = "none")

## Threshold Zero Padding ##
n_snps_tbl <- ginis %>% group_by(prive_code) %>% summarise(n_snps = mean(n_snps))
cor.test(traits_table$pcor_United, n_snps_tbl$n_snps)


ginis_tzp <- ginis %>%
  filter(bin_size==100000,
         bin_summary_method=="sum",
         (threshold==100) | ((threshold==n_significant_bins) & threshold < 100))
ginis_tzp %>%
  group_by(threshold_zero_padding) %>%
  summarize(mean_gini = mean(gini),
            sd_gini = sd(gini),
            n = n())
ginis_tzp_pw <- ginis_tzp %>%
  select(-settings, -threshold) %>%
  filter(n_significant_bins < 100) %>%
  pivot_wider(
    values_from = gini,
    names_from = threshold_zero_padding,
    names_prefix = "gini_TZP_"
  ) %>%
  mutate(gini_TZP_diff = gini_TZP_TRUE - gini_TZP_FALSE)

t.test(x=ginis_tzp_pw$gini_TZP_TRUE, y=ginis_tzp_pw$gini_TZP_FALSE, paired=TRUE)

(cor1 <- cor.test(ginis_tzp_pw$gini_TZP_diff, ginis_tzp_pw$n_significant_bins))
r <- cor1$estimate[[1]]
p <- cor1$p.value

ginis_tzp_pw_labels <- ginis_tzp_pw %>%
  group_by(prive_code) %>%
  summarize(gini_TZP_diff = median(gini_TZP_diff),
            n_significant_bins = mean(n_significant_bins)) %>%
  left_join(traits_table %>% select(prive_code,short_label), by="prive_code")
ggplot(ginis_tzp_pw, aes(x=gini_TZP_diff, y = n_significant_bins)) +
  geom_point(aes(color=ancestry)) +
  geom_smooth(method="lm") +
  geom_text_repel(data=ginis_tzp_pw_labels, aes(label=short_label)) +
  ylim(0,100) +
  xlab("Difference between Gini with TZP=TRUE and TZP=FALSE") +
  ylab("Number of Significant Bins") +
  theme_light() +
  labs(title="The relationship between Threshold Zero Padding and the Number of Significant Bins",
       subtitle=paste0("r = ", round(r,4), " :: p = ", formatC(p,format="E",digits=2),
                       "\nFiltered to bin size = 100kb, bin summary method = sum, threshold = 100"))


## Ancestry ##

# quick helper_function

ginis_ancestry <- traits_table %>%
  select(prive_code, short_label, starts_with("gini_"), f_stat) %>%
  pivot_longer(
    cols = starts_with("gini_"),
    names_to = "ancestry",
    names_prefix = "gini_",
    values_to = "gini"
  ) %>%
  left_join(traits_table %>% select(prive_code, gini_United), by="prive_code") %>%
  mutate(gini_sqdiff_UK = (gini - gini_United)**2)
ginis_ancestry[ginis_ancestry$ancestry=="United","ancestry"] <- "UK"

ginis_ancestry_grouped <- ginis_ancestry %>%
  group_by(prive_code, short_label, f_stat) %>%
  summarize(mean_gini = mean(gini),
            median_gini = median(gini),
            sd_gini = sd(gini)) %>%
  arrange(mean_gini)
ginis_ancestry <- ginis_ancestry %>%
  select(-short_label) %>%
  left_join(ginis_ancestry_grouped, by="prive_code") %>%
  arrange(mean_gini) %>%
  mutate(short_label = factor(short_label, levels=ginis_ancestry_grouped$short_label),
         gini_diff_sq = (gini - mean_gini)**2)

summary(aov(data=ginis_ancestry, gini ~ anestry))
ggplot(ginis_ancestry, aes(x=short_label, y=gini, color=ancestry)) +
  geom_point(alpha=0.5) +
  geom_point(aes(y=mean_gini), color="black", size=0.75) +
  geom_label_repel(data=ginis_ancestry %>% filter(gini_diff_sq > 0.010),
            aes(label = short_label)) +
  coord_flip() +
  xlab("Trait") +
  ylab("Gini") +
  labs(title = "Differences in Gini by Ancestry",
       subtitle = "Black dots represent mean gini across nine ancestries. Highlighted points deviate more than 0.1 units from mean gini.",
       color = "Ancestry") +
  theme_light() + 
  theme(axis.text.y = element_blank(), #element_text(size = 5),
        axis.ticks.y = element_blank()
        )

pop_centers <- read.csv("https://raw.githubusercontent.com/privefl/UKBB-PGS/main/pop_centers.csv",
                        stringsAsFactors = FALSE)
distances <- tibble(
  ancestry = pop_centers$Ancestry,
  prive_dist_to_UK = as.matrix(dist(pop_centers[2:17]))[,1]
) %>% arrange(ancestry)
distances$ancestry <- str_replace(distances$ancestry,"United Kingdom","UK")
distances <- distances %>%
  left_join(ginis_ancestry %>%
              group_by(ancestry) %>%
              summarize(mean_gini=mean(gini), mean_gini_sqdiff_UK = mean(gini_sqdiff_UK)) %>%
              arrange(-mean_gini_sqdiff_UK))

cor1 <- cor.test(distances$mean_gini_sqdiff_UK,distances$prive_dist_to_UK)
r <- cor1$estimate[[1]]
p <- cor1$p.value
ggplot(distances, aes(x=prive_dist_to_UK, y=mean_gini_sqdiff_UK)) +
  geom_smooth(method="lm") +
  geom_point(aes(color=ancestry)) +
  geom_label_repel(aes(label=ancestry, color=ancestry)) +
  xlab("PC distance to UK") +
  ylab("Mean squared difference from Gini(UK)") +
  labs(title = "Relationship between Gini differences from UK and genetic distance from UK",
       subtitle=paste0("r = ", round(r,4), " :: p = ", formatC(p,format="E",digits=2)),
       color = "Ancestry") +
  theme_light() +
  theme(legend.position = "none")



#### Verifying robustness of parameters ####

### bin size ####
ginis_binsize <- ginis %>% filter(
  bin_summary_method == "sum",
  threshold_zero_padding = TRUE,
  threshold == 100
  ) %>%
  select(prive_code, ancestry, bin_size, gini) %>%
  group_by(prive_code, bin_size) %>%
  summarize(gini = mean(gini)) %>%
  pivot_wider(
    names_from = bin_size,
    values_from = gini,
    names_prefix = "gini_"
  ) %>%
  left_join(traits_table %>% select(prive_code,short_label), by="prive_code")
cor1 <- cor.test(ginis_binsize$gini_1, ginis_binsize$`gini_1e+05`)
r <- cor1$estimate[[1]]
p <- cor1$p.value
ggplot(ginis_binsize, aes(x=gini_1,y=`gini_1e+05`)) +
  geom_abline(slope=1,intercept=0) +
  #geom_text(aes(label=short_label), size=2) +
  geom_point() +
  xlim(0,1) +
  ylim(0,1) +
  xlab("Gini (bin size = 1)") +
  ylab("Gini (bin size = 100kb)") +
  labs(title = "Gini is robust by bin size",
       subtitle = paste0("r = ", round(r,4), " :: p = ", formatC(p,format="E",digits=2),
                         "\nFiltered to bin summary method = sum, threshold = 100, TZP = true, Gini averaged across ancestries")) +
  theme_light()

### bin summary method ####
ginis_bsm <- ginis %>% filter(
  bin_size == 100000,
  threshold_zero_padding = TRUE,
  threshold == 100
) %>%
  select(prive_code, ancestry, bin_summary_method, gini) %>%
  group_by(prive_code, bin_summary_method) %>%
  summarize(gini = mean(gini)) %>%
  pivot_wider(
    names_from = bin_summary_method,
    values_from = gini,
    names_prefix = "gini_"
  ) %>%
  left_join(traits_table %>% select(prive_code,short_label), by="prive_code")
cor1 <- cor.test(ginis_bsm$gini_sum, ginis_bsm$gini_max)
r <- cor1$estimate[[1]]
p <- cor1$p.value
ggplot(ginis_bsm, aes(x=gini_sum,y=gini_max)) +
  geom_abline(slope=1,intercept=0) +
  #geom_text(aes(label=short_label), size=2) +
  geom_point() +
  xlim(0,1) +
  ylim(0,1) +
  xlab("Gini (bin summary method = sum)") +
  ylab("Gini (bin summary method = max)") +
  labs(title = "Gini is robust by bin summary method",
       subtitle = paste0("r = ", round(r,4), " :: p = ", formatC(p,format="E",digits=2),
                         "\nFiltered to bin size = 100kb, threshold = 100, TZP = true, Gini averaged across ancestries")) +
  theme_light()

# comparing rank order by threshold
ginis_t <- ginis %>%
  filter(ancestry=="United", bin_size == 100000, bin_summary_method == "sum", (threshold_zero_padding == TRUE)|(n_significant_bins > threshold) ) %>%
  filter(threshold %in% c(50,100,150, 250)) %>%
  select(prive_code, threshold, gini, n_significant_bins) %>%
  group_by(threshold) %>%
  mutate(gini_t_rank = order(order(gini, decreasing=FALSE))) %>%
  ungroup() %>% group_by(prive_code) %>%
  mutate(rank_sd = sd(gini_t_rank)) %>%
  left_join(traits_table %>% select(prive_code,group_consolidated, trait_type) %>% filter(row_number() %% 1 == 0), by="prive_code") %>%
  ungroup() %>% group_by(threshold) %>%
  mutate(t_mean_gini = mean(gini))

ggplot(ginis_t ,
       aes(x=as.factor(threshold), y=gini_t_rank, group=prive_code)) +
  geom_line(aes(color=group_consolidated)) +
  facet_wrap(~ rank_sd > summary(unique(ginis_t$rank_sd))[[5]])
  
ggplot(ginis_t ,
       aes(x=as.factor(threshold), y=gini - t_mean_gini, group=prive_code)) +
  geom_line(aes(color=group_consolidated)) +
  facet_wrap(~ rank_sd > summary(unique(ginis_t$rank_sd))[[5]])

##### sandbox
traits_table <- traits_table %>%
  mutate(N_total = pmax(N, N_case, na.rm=TRUE))
ggplot(traits_table, aes(x=log10(N_case), y=gini_United)) +
  geom_point(aes(color=group_consolidated)) +
  labs(title="Binary traits only")

ggplot(traits_table, aes(x=log10(N_case/(N_case+N_control)), y=gini_United)) +
  geom_point(aes(color=group_consolidated)) +
  labs(title="Binary traits only")

ggplot(traits_table, aes(x=log10(N_case), y=ldpred2_h2)) +
  geom_point(aes(color=group_consolidated))
cor.test(log10(traits_table$N_case), traits_table$gini_United)
ggplot(traits_table, aes(x=log10(N_case), y=pcor_United)) +
  geom_point(aes(color=group_consolidated))

ggplot(traits_table, aes(x=N, y=gini_United)) +
  geom_point(aes(color=group_consolidated)) +
  geom_smooth(method="lm")
####
traits_table <- traits_table %>%
  mutate(K = N_case / (N_case + N_control)) %>%
  rowwise() %>%
  mutate(
    ldpred2_h2_raw = ldpred2_h2 / bigsnpr::coef_to_liab(K,K),
    ldsc_h2_raw = ldsc_h2 / bigsnpr::coef_to_liab(K,K)
  )
ggplot(traits_table, aes(x=pmax(ldpred2_h2,ldpred2_h2_raw,na.rm=TRUE), y=gini_United)) +
  geom_point(aes(color=trait_type))

ggplot(traits_table %>% filter(!(prive_code %in% low_prevalence)), aes(x=ldpred2_h2, pcor_United)) +
  geom_text(aes(color=trait_type, label=paste(short_label,N_case)), size=2)

ggplot(traits_tab)

cor.test((traits_table%>%filter(trait_type=="binary"))$ldpred2_h2, (traits_table%>%filter(trait_type=="binary"))$pcor_United)
cor.test((traits_table%>%filter(trait_type=="quantitative"))$ldpred2_h2, (traits_table%>%filter(trait_type=="quantitative"))$pcor_United)
