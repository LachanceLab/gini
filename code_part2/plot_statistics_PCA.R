# plot_statistics_PCA.R

# Creates PCA plots of the traits, using the six summary statistics as the initial dimensions

# Libraries and directories ####
library(tidyverse)
library(data.table)
source("custom_ggbiplot.R") # loads custom 'ggbiplot' function from 'ggbiplot' package

# sets working directory
setwd("./")

# sets location of trait_table generated in part1
loc_table <- "../generated_data/traits_table.txt"
# sets directory for generated figures
dir_out <- "../generated_figures/"

# sets scaling factors for image output. Default = 3
sf <- 3
print_mode <- "png" # set to either "png" or "pdf"

# Printing function ####
print_plot <- function(gg, loc_out, print_mode, plot_width, plot_height, sf) {
  if (print_mode == "png") {
    png(loc_out, width = plot_width*sf, height = plot_height*sf)
  } else if (print_mode == "pdf") {
    pdf(loc_out, width = plot_width*sf / 75, height = plot_height*sf / 75)
  }
  print(gg)
  dev.off()
}


# Shared Code ####
vars <- c("SNP Heritability"="ldpred2_h2",
          "LD Variability"="traitLD_unadj_CoV",
          "Genomic Inequality"="gini_panUKB",
          "PGS Accuracy"="PGS_R2_United",
          "Portability"="portability_index",
          "Divergence"="log_F")

# reads traits table, filters out low prevalence binary traits,
traits_table <- as_tibble(fread(loc_table)) %>%
  filter(PGS_trait_type == GWAS_trait_type)

# color legend for plots
gg_pca_scale <- list(
  labels = c("Biological Measures","Diseases","Lifestyle/Psychological","Physical Measures"),
  breaks = c("biological measures","diseases","lifestyle/psychological","physical measures"),
  values = c("biological measures"="#F8766D", "diseases"="#A3A500","lifestyle/psychological"="#00BF7D","physical measures"="#00B0F6","psychological"="#E76BF3")
)
trait_types <- c("binary","quantitative")

PCA_theme <- theme_light() +
  theme(legend.position = "top",
        text = element_text(size = 20*sf),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(size=0.75*sf),
        axis.ticks.length = unit(4*sf, "points"))
plot_width <- 1000
plot_height <- 1000

# makes PCA for all traits, comparing binary vs quantitative traits
matrix <- traits_table %>% select(prive_code, short_label, group_consolidated, GWAS_trait_type, all_of(unname(vars)))
colnames(matrix)[which(colnames(matrix) %in% unname(vars))] <- names(vars)
matrix.pca <- prcomp(matrix[names(vars)], center=TRUE, scale. = TRUE)

# Binary vs Quant ####
# plots PCA for all traits, comparing binary vs quantitative traits
gg <- custom_ggbiplot(matrix.pca, groups = matrix$GWAS_trait_type, ellipse=TRUE, #labels=matrix$short_label,
                varname.adjust = 1.25, varname.size = 6*sf, var.color="gray20", ell.size = 0.6*sf,
               labels.size = 5*sf, var.scale = 1, obs.scale = 1, arrow.size = 0.75*sf,
               point.size=2*sf, #plot.outside = FALSE
               ) +
  geom_text(aes(label="", color=matrix$GWAS_trait_type), key_glyph = "rect") + # empty geom
  PCA_theme +
  labs(title = NULL, subtitle = NULL, color = "Trait Type") +
  scale_color_manual(labels=c("Binary","Quantitative"),
                     breaks=c("binary", "quantitative"),
                     values = c("binary"="#F8766D", "quantitative"="#00BFC4")) +
  scale_x_continuous(expand = c(0.1,0.075))
ratio <- diff(layer_scales(gg)$y$range$range) / diff(layer_scales(gg)$x$range$range)

# Saves image onto system
loc_out <- paste0(dir_out,"summary_stats_PCA_ALL")
#print_plot(gg, loc_out, print_mode, plot_width, plot_height*ratio, sf)
print_plot(gg, paste0(loc_out,".png"), "png", plot_width, plot_height*ratio, sf)
print_plot(gg, paste0(loc_out,".pdf"), "pdf", plot_width, plot_height*ratio, sf)
print(paste0("Saved PCA plot of all traits"))


# Trait Groups ####
# makes and plots PCA for each of binary and quantitative trait, comparing groups
# makes PCA by group for each trait type
for (the_trait_type in trait_types) {
  overlap_fix = FALSE
  # changes plotting settings depending on trait type
  if (the_trait_type == "binary") {gg_pca_scale_subset <- c(2,3,4)
  } else if (the_trait_type == "quantitative") {
    gg_pca_scale_subset <- c(1,3,4)
    #overlap_fix = TRUE
  }
  
  # selects just traits of trait type
  matrix_filtered <- matrix %>% filter(GWAS_trait_type == the_trait_type)
  
  # makes PCA
  matrix.pca <- prcomp(matrix_filtered[names(vars)], center=TRUE, scale. = TRUE)
  
  # plots PCA
  gg <- custom_ggbiplot(matrix.pca, groups = matrix_filtered$group_consolidated, ellipse=TRUE, #labels=matrix_filtered$short_label,
                        varname.adjust = 1.25, varname.size = 6*sf, var.color="gray20", ell.size = 0.6*sf,
                        labels.size = 5*sf, var.scale = 1, obs.scale = 1, arrow.size = 0.75*sf, overlap_fix = overlap_fix,
                        point.size=2*sf, #plot.outside = FALSE
                        ) +
    geom_text(aes(label="", color=matrix_filtered$group_consolidated), key_glyph = "rect") + # empty geom
    PCA_theme +
    labs(title = NULL, subtitle = NULL, color = "Trait Group") +
    scale_color_manual(labels=gg_pca_scale[["labels"]][gg_pca_scale_subset],
                       breaks=gg_pca_scale[["breaks"]][gg_pca_scale_subset],
                       values=gg_pca_scale[["values"]][gg_pca_scale_subset]) +
    scale_x_continuous(expand = c(0.075,0.075))
  if(the_trait_type=="binary") {gg <- gg + scale_y_continuous(expand = c(0.085,0.085))}
  
  ratio <- diff(layer_scales(gg)$y$range$range) / diff(layer_scales(gg)$x$range$range)
  
  # saves image to system
  loc_out <- paste0(dir_out,"summary_stats_PCA_",the_trait_type)
  print_plot(gg, paste0(loc_out,".png"), "png", plot_width, plot_height*ratio, sf)
  print_plot(gg, paste0(loc_out,".pdf"), "pdf", plot_width, plot_height*ratio, sf)
  print(paste("Saved PCA plot of",the_trait_type,"traits"))
}
