# print_tables.R

# This script produces the tables found in the publication

### Libraries and directories ####
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

### Code ####

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

#### Makes big traits table ####

# reads traits table
traits_table <- as_tibble(fread(loc_table))
vars <- c("ldpred2_h2","cMperMb","gini_United","pcor_United","portability_index", "f_stat")
rounding_decimals <- 3

big_table <- traits_table %>%
  select(description, trait_type, group_consolidated, all_of(vars)) %>%
  mutate(trait_type = ifelse(trait_type=="binary","Binary","Quantitative"),
         f_stat = log10(f_stat),
         portability_index = ifelse(portability_index>0,0,1000*portability_index)) %>%
  group_by(trait_type) %>%
  arrange(group_consolidated, description)

big_table_gt <- gt(big_table) %>%
  gt_theme() %>%
  summary_rows(
    columns = vars,
    groups = TRUE,
    fns = list("Mean" = ~mean(.)),
    formatter = fmt_number,
    decimals = rounding_decimals
  ) %>%
  cols_label(
    description = "Trait",
    group_consolidated = "Trait Group",
    ldpred2_h2 = md("Heritability (h<sup>2</sup>)"),
    cMperMb = "Recombination Rate (R)",
    gini_United = md("Gini<sub>100,UK</sub>"),
    pcor_United = md("PGS Efficacy (&rho;<sub>UK</sub>)"),
    portability_index = "Portability (1000*m)",
    f_stat = md("Divergence (D)")
  ) %>%
  fmt_number(
    columns = vars,
    suffixing = F,
    decimals = rounding_decimals
  ) %>%
  cols_align(
    align = "left",
    columns = c("description")
  ) %>%
  tab_header(
    title = "The Six Summary Statistics for 211 Traits",
    subtitle = "Full table provided in Github repository"
  )

big_table_gt
gtsave(big_table_gt, "table_big_table.png", dir_out, vwidth = 2160, vheight=1620)


#### Makes High and Low Gini Table ####
n_show <- 5 # shows top n_show highest and top n_show lowest
trait_types <- c("binary", "quantitative")

for (the_trait_type in trait_types[1]) {
  #N_total <- nrow(traits_table %>% filter(trait_type == the_trait_type))
  N_total <- nrow(traits_table)
  
  temp <- traits_table %>%
    arrange(gini_United) %>%
    #filter(trait_type == the_trait_type) %>%
    mutate(gini_United = as.character(round(gini_United,rounding_decimals)),
           rank = row_number()) %>%
    select(rank, description, 
           #trait_type, group_consolidated, 
           gini_United)
  
  table1 <- temp %>%
    filter(rank <= n_show) %>%
    add_row(
      rank = 0, description = "...",
      #group_consolidated = "", trait_type = "",
      gini_United = ""
    ) %>%
    add_row(
      temp %>% filter(row_number() > (N_total - n_show))
    ) %>%
    mutate(rank = ifelse(rank==0,"",as.character(rank)))
  
  table1gt <- gt(table1) %>%
    gt_theme() %>%
    cols_label(
      rank = "#",
      description = "Trait",
      #trait_type = "Trait Type", group_consolidated = "Trait Group",
      gini_United = md("Gini<sub>100,UK</sub>")
    ) %>%
    cols_align(
      align = "left",
      columns = c("description")
    ) %>%
    tab_header(
      #subtitle = paste("Only",the_trait_type,"traits"),
      title = paste0(n_show," Lowest and ",n_show," Highest Gini Traits")
    )
  
  table1gt
  #gtsave(table1gt, paste0("table1",the_trait_type,".png"), dir_out)
  gtsave(table1gt, paste0("table1_ALL",n_show,".png"), dir_out)
  print(paste("Made gini table for",the_trait_type,"traits"))
}

#### Makes High and Low Divergence Table ####
n_show <- 10 # shows top n_show highest and top n_show lowest

N_total <- nrow(traits_table)

temp <- traits_table %>%
  mutate(f_stat = log10(f_stat)) %>%
  arrange(f_stat) %>%
  mutate(f_stat = as.character(round(f_stat,rounding_decimals)),
         rank = row_number()) %>%
  select(rank, description,
         #trait_type, group_consolidated,
         f_stat)

table2 <- temp %>%
  filter(rank <= n_show) %>%
  add_row(
    rank = 0, description = "...",
    #group_consolidated = "", trait_type = "",
    f_stat = ""
  ) %>%
  add_row(
    temp %>% filter(row_number() > (N_total - n_show))
  ) %>%
  mutate(rank = ifelse(rank==0,"",as.character(rank)))

table2gt <- gt(table2) %>%
  gt_theme() %>%
  cols_label(
    rank = "#",
    description = "Trait",
    #trait_type = "Trait Type", group_consolidated = "Trait Group",
    f_stat = md("Divergence (D)")
  ) %>%
  cols_align(
    align = "left",
    columns = c("description")
  ) %>%
  tab_header(
    #subtitle = paste("Only",the_trait_type,"traits"),
    title = paste0(n_show," Lowest and ",n_show," Highest Divergence Traits")
  )

table2gt
gtsave(table2gt, paste0("table2_ALL",n_show,".png"), dir_out, vwidth = 1160, vheight=1620)
print(paste("Made divergence table"))
