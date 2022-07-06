#!/usr/bin/env Rscript
#Example input: Rscript calc_fstat.R -f <path_to_f_stat.csv> -m <path_to_master_file> -p <path_to_pop_iids.csv> -t1 <> -t2 <>
library(rio)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(optparse)

option_list = list(
  make_option(c("-f", "--f_stat"), type="character", default=NULL, 
              help="csv for f_stat"),
  make_option(c("-m", "--master"), type="character", default=NULL, 
              help="master table"),
  make_option(c("-p", "--pop_iids"), type="character", default=NULL, 
              help="csv with pop_idds, generated from downsample.py"),
  make_option(c("-t1", "--trait1"), type="character", default=NULL, 
              help="csv of scores for individuals for low f-stat trait"),
  make_option(c("-t2", "--trait2"), type="character", default=NULL, 
              help="csv of scores for individuals for high f-stat trait")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

f_df <- import(opt$f_stat)
#f_df <- import('~/Desktop/Gini-PGS/f_stat.csv')

master_df <- import(opt$master)
#master_df <- import('~/Desktop/Gini-PGS/table_PLR_100k_sum_100_zp.txt')
master_df <- master_df %>% select(c(prive_code, trait, description))
f_df <- merge(f_df, master_df, by ='prive_code')

divergence_df <- import(opt$trait1)
#divergence_df <- import('~/Desktop/final/250.1-1KG_PLR.csv')
divergence_df <- divergence_df[,2:ncol(divergence_df)] #drop index column
divergence_df$sum <- rowSums(divergence_df[,3:ncol(divergence_df)])
divergence_df <- divergence_df %>% select(c(FID, IID, sum))

divergence_df <- import(opt$trait2)
#divergence_df2 <- import('~/Desktop/final/darker_skin0-1KG_PLR.csv')
divergence_df2 <- divergence_df2[,2:ncol(divergence_df2)] #drop index column
divergence_df2$sum <- rowSums(divergence_df2[,3:ncol(divergence_df2)])
divergence_df2 <- divergence_df2 %>% select(c(FID, IID, sum))

pop_iids_df <- import(opt$pop_iids)
#pop_iids_df <- import('/Users/adrianharris/Desktop/Gini-PGS/UKB/randomized_pop_iids.csv')
pop_iids_df <- pop_iids_df[,1:2]

divergence_df <- merge(divergence_df, pop_iids_df, by = 'IID')
divergence_df2 <- merge(divergence_df2, pop_iids_df, by = 'IID')

#anova to calculate f-stat
one.way <- aov(as.numeric(sum) ~ ancestry, data = divergence_df)
f_val <- round(summary(one.way)[[1]]["F value"][[1]][[1]], digits = 3)
f_val <- round(log(f_val), digits = 3)
#p_val <- summary(one.way)[[1]]["Pr(>F)"][[1]][[1]]

one.way2 <- aov(as.numeric(sum) ~ ancestry, data = divergence_df2)
f_val2 <- round(summary(one.way2)[[1]]["F value"][[1]][[1]], digits = 3)
f_val2 <- round(log(f_val2), digits = 3)
#p_val <- summary(one.way)[[1]]["Pr(>F)"][[1]][[1]]

histogram1 <- ggplot(divergence_df, aes(x=as.numeric(sum), fill=ancestry)) + geom_density(color='#e9ecef', alpha=0.6, position='identity') + theme_bw() + labs(x='Polygenic Score per UKBB Individual', y='Density') + ggtitle('Type 1 Diabetes') + annotate("text",x=2,y=1,label=paste0("log(F) =", f_val, sep=" "), size=6)
histogram1 <- histogram1 + theme(legend.title=element_blank())
print(histogram1)

histogram2 <- ggplot(divergence_df2, aes(x=as.numeric(sum), fill=ancestry)) + geom_density(color='#e9ecef', alpha=0.6, position='identity') + theme_bw() + labs(x='Polygenic Score per UKBB Individual', y='Density') + ggtitle('Skin Colour') + annotate("text",x=0.75,y=2.2,label=paste0("log(F) =", f_val2, sep=" "), size = 6)
histogram2 <- histogram2 + theme(legend.title=element_blank()) 

print(histogram2)

ggarrange(histogram1, histogram2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")



