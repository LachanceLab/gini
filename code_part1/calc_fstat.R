#!/usr/bin/env Rscript
#Example input: Rscript calc_fstat.R -p <path_to_pop_iids.csv> -d <directory_of_calc_effect_dose.py_results> -o <output_file_path>

library(rio)
library(ggplot2)
library(tidyverse)
library(optparse)

option_list = list(
  make_option(c("-p", "--pop_iids"), type="character", default=NULL, 
              help="csv with pop_idds, generated from downsample.py"),
  make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="directory with resulting files from calc_effect_dose.py"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
             help="output file path (with csv extension)")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

f_df <- data.frame(prive_code = character(), f_stat = double(), p_value = double())
iteration = 1

pop_iids_df <- import(opt$pop_iids)
pop_iids_df <- pop_iids_df[,1:2]

calc_stat <- function(filename) {
  print(iteration)
  df <- import(filename)
  df <- df[,2:ncol(df)]
  final_df <- merge(pop_iids_df, df, by = 'IID')
  final_df['sum'] <- rowSums(final_df[,4:ncol(final_df)])
  
  final_df <- final_df %>% select(c(FID, IID, ancestry, sum))
  
  #Acquire stats
  one.way <- aov(as.numeric(sum) ~ ancestry, data = final_df)
  f_val <- summary(one.way)[[1]]["F value"][[1]][[1]]
  p_val <- summary(one.way)[[1]]["Pr(>F)"][[1]][[1]]
  
  rm(final_df)
  rm(df)
  filename <- unlist(strsplit(filename, "/"))
  filename <- rev(filename)[1]
  filename <- unlist(strsplit(filename, "-"))
  filename <- filename[1]
  return(c(filename, f_val, p_val))
}

files <- list.files(path=opt$dir, pattern="*.csv", full.names=TRUE, recursive=FALSE)

for (file in files) {
  f_df[nrow(f_df) + 1,] <- calc_stat(file)
  iteration = iteration + 1
}

#Pull first, mutate remaining rows, and append column with names
row.names(f_df) <- f_df[,1]
f_df <- mutate_all(f_df[,2:ncol(f_df)], function(x) as.numeric(as.character(x)))
f_df['prive_code'] <- row.names(f_df)
f_df <- f_df %>% select(c(prive_code, f_stat, p_value))

write.csv(f_df, opt$out, row.names = FALSE)

