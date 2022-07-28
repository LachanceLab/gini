# plot_statistics_PCA.R

# Creates PCA plots of the traits, using the six summary statistics as the initial dimensions

### Libraries and directories ####
library(tidyverse)
library(data.table)
source("custom_ggbiplot.R") # loads custom 'ggbiplot' function from 'ggbiplot' package

# sets working directory
setwd("./")

# sets location of trait_table generated in part1
loc_table <- "../generated_data/traits_table.txt"
# sets directory for generated figures
dir_out <- "../generated_figures/"

# sets scaling factors for image output. Default = 2
sf <- 3
print_mode <- "png" # set to either "png" or "pdf"

### Code ###
vars <- c("Heritability"="ldpred2_h2",
          "Recombination Rate"="cMperMb",
          "Polygenicity"="gini_United",
          "Prediction"="pcor_United",
          "Portability"="portability_index",
          "Divergence"="f_stat")

traits_table <- as_tibble(fread(loc_table))

# plots PCA for all traits, comparing binary vs quantitative traits
matrix <- traits_table %>% select(prive_code, short_label, group_consolidated, trait_type, all_of(unname(vars)))
colnames(matrix)[which(colnames(matrix) %in% unname(vars))] <- names(vars)
matrix$Divergence <- log10(matrix$Divergence)
matrix.pca <- prcomp(matrix[names(vars)], center=TRUE, scale. = TRUE)

gg <- custom_ggbiplot(matrix.pca, groups = matrix$trait_type, ellipse=TRUE, labels=matrix$short_label,
                varname.adjust = 1.75, varname.size = 6*sf, var.color="gray20", ell.size = 0.6*sf,
               labels.size = 3*sf, var.scale = 1, obs.scale = 1, arrow.size = 0.75*sf) +
  geom_text(aes(label="", color=matrix$trait_type), key_glyph = "rect") + # empty geom
  theme_light() +
  theme(legend.position = "bottom",
        text = element_text(size = 20*sf)) +
  labs(title="PCA of summary statistics for all traits",
       subtitle = NULL,
       color = "Trait Type") +
  scale_color_manual(labels=c("Binary","Quantitative"),
                     breaks=c("binary", "quantitative"),
                     values = c("binary"="#F8766D", "quantitative"="#00BFC4"))
# Saves image onto system
plot_width <- 1000
plot_height <- 1000
loc_out <- paste0(dir_out,"summary_stats_PCA_ALL.", print_mode)
if (print_mode == "png") {
  png(loc_out, width = plot_width*sf, height = plot_height*sf)
} else if (print_mode == "pdf") {
  pdf(loc_out, width = plot_width*sf / 75, height = plot_height*sf / 75)
}
print(gg)
dev.off()
print(paste0("Saved PCA plot of all traits"))


# plots PCA for each of binary and quantitative trait, comparing groups
gg_pca_scale <- list(
  labels = c("Biological Measures","Diseases","Lifestyle/Psychological","Physical Measures"),
  breaks = c("biological measures","diseases","lifestyle/psychological","physical measures"),
  values = c("biological measures"="#F8766D", "diseases"="#A3A500","lifestyle/psychological"="#00BF7D","physical measures"="#00B0F6","psychological"="#E76BF3")
)
trait_types <- c("binary","quantitative")
for (the_trait_type in trait_types) {
  overlap_fix = FALSE
  if (the_trait_type == "binary") {gg_pca_scale_subset <- c(2,3,4)}
  else if (the_trait_type == "quantitative") {
    gg_pca_scale_subset <- c(1,3,4)
    overlap_fix = TRUE
  }
  
  matrix_filtered <- matrix %>%
    filter(trait_type == the_trait_type)
  
  matrix.pca <- prcomp(matrix_filtered[names(vars)], center=TRUE, scale. = TRUE)
  gg <- custom_ggbiplot(matrix.pca, groups = matrix_filtered$group_consolidated, ellipse=TRUE, labels=matrix_filtered$short_label,
                        varname.adjust = 1.25, varname.size = 6*sf, var.color="gray20", ell.size = 0.6*sf,
                        labels.size = 3*sf, var.scale = 1, obs.scale = 1, arrow.size = 0.75*sf, overlap_fix = overlap_fix) +
    geom_text(aes(label="", color=matrix_filtered$group_consolidated), key_glyph = "rect") + # empty geom
    theme_light() +
    theme(legend.position = "bottom",
          text = element_text(size = 20*sf)) +
    labs(title=paste("PCA of summary statistics for",the_trait_type,"traits"),
         subtitle = NULL,
         color = "Trait Group")  +
    scale_color_manual(labels=gg_pca_scale[["labels"]][gg_pca_scale_subset],
                       breaks=gg_pca_scale[["breaks"]][gg_pca_scale_subset],
                       values=gg_pca_scale[["values"]][gg_pca_scale_subset])
  
  loc_out <- paste0(dir_out,"summary_stats_PCA_",the_trait_type,".", print_mode)
  if (print_mode == "png") {
    png(loc_out, width = plot_width*sf, height = plot_height*sf)
  } else if (print_mode == "pdf") {
    pdf(loc_out, width = plot_width*sf / 75, height = plot_height*sf / 75)
  }
  print(gg)
  dev.off()
  print(paste("Saved PCA plot of",the_trait_type,"traits"))
}
