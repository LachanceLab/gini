# print_tables.R

# This script produces the tables found in the publication and supplemental

# Libraries and directories ####
library(tidyverse)
library(webshot)
library(gt)
library(data.table)

# sets working directory
setwd("./")

# sets the location of the traits table
loc_table <- "../generated_data/traits_table.txt"
# sets directory of outputted figures
dir_out <- "../generated_figures/"

# Shared Code ####

gt_theme <- function(data) {
  tab_options(
    data = data,
    table.font.name = "Helvetica"
  ) %>%
    cols_align(
      align = "center",
      columns = everything()
    )
}

# reads traits table
traits_table <- as_tibble(fread(loc_table))
rounding_decimals <- 3

# Big Traits Table ####

## tab-delimited .txt S1 ####
TS1 <- traits_table %>%
  select(prive_code,
         shortname = short_label,
         description,
         panUKB_code = filename,
         panUKB_trait_type = GWAS_trait_type,
         prive_trait_type = PGS_trait_type,
         group = group_consolidated,
         n_sig_ind_SNPs = n_sig_SNPs,
         heritability = ldpred2_h2,
         LD_variability = traitLD_unadj_CoV,
         gini = gini_panUKB,
         PGS_accuracy = pcor_United,
         portability = portability_index,
         PGS_divergence = log_F
         )
TS1$panUKB_code <- sapply(TS1$panUKB_code, function(x) substring(x, 1, nchar(x)-8))
loc_out <- "../generated_figures/traits_table_S1.txt"
fwrite(TS1,loc_out, sep="\t")


## GT Table S1 (deprecated) ####

#vars <- c("ldpred2_h2","traitLD_unadj_CoV","gini_panUKB","pcor_United","portability_index", "log_F")
# big_table <- traits_table %>%
#   filter(!(description %in% low_prevalence)) %>%
#   select(description, trait_type, group_consolidated, all_of(vars)) %>%
#   mutate(trait_type = ifelse(trait_type=="binary","Binary","Quantitative"),
#          description = ifelse(description %in% low_prevalence, paste0(description,"*"),description),
#          f_stat = log10(f_stat),
#          portability_index = ifelse(portability_index>0,0,portability_index)) %>%
#   group_by(trait_type) %>%
#   arrange(group_consolidated, description)
# 
# # converts Table S1 into a GT table
# big_table_gt <- gt(big_table) %>%
#   gt_theme() %>%
#   summary_rows(
#     columns = vars,
#     groups = TRUE,
#     fns = list("Mean" = ~mean(.)),
#     formatter = fmt_number,
#     decimals = rounding_decimals
#   ) %>%
#   cols_label(
#     description = "Trait",
#     group_consolidated = "Trait Group",
#     ldpred2_h2 = md("h<sup>2</sup><sub>SNP</sup>"),
#     cMperMb = "R",
#     gini_United = md("G<sub>100,UK</sub>"),
#     pcor_United = md("&rho;<sub>UK</sub>"),
#     portability_index = "m",
#     f_stat = md("D")
#   ) %>%
#   fmt_number(
#     columns = vars[vars !="portability_index"],
#     suffixing = F,
#     decimals = rounding_decimals
#   ) %>%
#   fmt_number(
#     columns = portability_index,
#     suffixing = F,
#     decimals = 5
#   ) %>%
#   cols_align(
#     align = "left",
#     columns = c("description")
#   ) %>%
#   tab_header(
#     title = paste("The Six Summary Statistics for", nrow(big_table),"Traits"),
#     subtitle = "Full table provided in Github repository"
#   )
# 
# # saves Table S1 to system
# big_table_gt
# gtsave(big_table_gt, "table_big_table.png", dir_out, vwidth = 1000, zoom = 1, cliprect=c(0,0,1000,6000))
# # for some reason, text gets cutoff when printed as a PDF
# big_table_gt <- big_table_gt %>% tab_options(table.font.size = "90%") 
# gtsave(big_table_gt, "table_big_table.pdf", dir_out, vwidth = 1000, zoom = 1, cliprect=c(0,0,1000,6000))


# H/L Gini Table ####
n_show <- 5 # shows top n_show highest and top n_show lowest

# temporary table with just the gini, trait, and rank 
temp <- traits_table %>%
  arrange(gini_panUKB) %>%
  mutate(gini_panUKB = as.character(formatC(gini_panUKB,rounding_decimals,format="f")),
         rank = as.character(row_number())) %>%
  select(rank, short_label,GWAS_trait_type, gini_panUKB)

# total number of traits (should be 163)
N_total <- nrow(temp)

# makes Table 1
table1 <- temp %>%
  filter(as.numeric(rank) <= n_show) %>%
  add_row(
    rank = "...", short_label = "...", GWAS_trait_type = "...",
    gini_panUKB = "..."
  ) %>%
  add_row(
    temp %>% filter(row_number() > (N_total - n_show))
  )

# converts Table 1 into a GT table
table1gt <- gt(table1) %>%
  gt_theme() %>%
  cols_label(
    rank = "Rank",
    short_label = "Trait",
    GWAS_trait_type = "Type",
    gini_panUKB = md("G<sub>500,Meta</sub>")
  ) %>%
  cols_align(
    align = "left",
    columns = c("short_label")
  ) %>%
  tab_header(
    title = "Traits with the lowest and highest Gini coefficients"
  )

# saves Table 1 to system
table1gt
gtsave(table1gt, paste0("table1_gini.png"), dir_out, vwidth = 2160, vheight=1620)
gtsave(table1gt, paste0("table1_gini.rtf"), dir_out)
print(paste("Made gini table"))


# H/L Divergence Table ####
n_show <- 5 # shows top n_show highest and top n_show lowest

# temporary table with just the divergence, trait, and rank 
temp <- traits_table %>%
  arrange(log_F) %>%
  mutate(log_F = as.character(formatC(log_F,rounding_decimals,format="f")),
         rank = as.character(row_number())) %>%
  select(rank, short_label,PGS_trait_type, log_F)

# total number of traits (should be 163)
N_total <- nrow(temp)

# makes Table 2
table2 <- temp %>%
  filter(as.numeric(rank) <= n_show) %>%
  add_row(
    rank = "...", short_label = "...", PGS_trait_type = "...",
    log_F = "..."
  ) %>%
  add_row(
    temp %>% filter(row_number() > (N_total - n_show))
  )

# converts Table 2 into a GT table
table2gt <- gt(table2) %>%
  gt_theme() %>%
  cols_label(
    rank = "Rank",
    short_label = "Trait",
    PGS_trait_type = "Type",
    log_F = md("D")
  ) %>%
  cols_align(
    align = "left",
    columns = c("short_label")
  ) %>%
  tab_header(
    title = "Traits with the lowest and highest PGS divergence"
)

# saves Table 2 to system
table2gt
gtsave(table2gt, paste0("table2_divergence.png"), dir_out, vwidth = 2160, vheight=1620)
gtsave(table2gt, paste0("table2_divergence.rtf"), dir_out)
print(paste("Made divergence table"))
