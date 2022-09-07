# plot_divergence.R

# Plots the PRS distribution of different ancestries, highlighting their divergence


### Libraries and directories ####
library(tidyverse)
library(ggpubr)
library(data.table)

# sets working directory
setwd("./")

# Sets the location of the sampled individuals' PRSs
loc_PRSs <- "../generated_data/pop_sampled_PRSs.txt"
# sets the location of the traits table
loc_table <- "../generated_data/traits_table.txt"
# sets directory of outputted figures
dir_out <- "../generated_figures/"

### Code ####

# reads traits table
traits_table <- as_tibble(fread(loc_table))
# reads the PRSs and changes "United" to "UK"
PRSs <- as_tibble(fread(loc_PRSs)) %>% filter(ancestry != "Ashkenazi")
PRSs[PRSs$ancestry=="United","ancestry"] <- "UK"

# function that generates divergence plot
plot_divergence <- function(code) {
  
  # extracts info about trait
  PRS_trait <- PRSs %>% select(Ancestry = ancestry, PRS = all_of(code))
  slice <- traits_table %>% filter(prive_code == code)
  description <- slice$description[1]
  f_stat <- slice$f_stat[1]
  p_value <- slice$p_value_f[1]
  logfstat <- round(log10(f_stat),2)
  
  # converts p-value to more legible text
  if (p_value < 1E-320) {
    p_text <- "< 1E-320"
    p_text <- bquote(ANOVA:~log[10](F)==.(logfstat)~~~~~p-value<10^{-320})
  } else {
    p_text <- formatC(p_value,format="E", digits=2)
    p_text_stem <- as.numeric(substr(p_text,1,4))
    p_text_exp <- as.numeric(substr(p_text,6,10))
    p_text <- bquote(ANOVA:~log[10](F)==.(logfstat)~~~~~p-value==.(p_text_stem)%*%10^{.(p_text_exp)})
  }
  
  #ff <- bquote(10^{.(threshold)})
  #subtitle <- paste0("ANOVA: log10(F-stat) = ", round(log10(f_stat),2),". p-value ",p_text )
  
  # plots divergence
  gg<-ggplot(PRS_trait, aes(x=PRS, fill=Ancestry)) +
    geom_density(color='#e9ecef', alpha=0.6, position='identity') +
    theme_light() +
    theme(legend.position = "bottom",
          legend.title=element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab("Polygenic Score per UKBB Individual") +
    ylab("Density") +
    scale_y_continuous(expand=expansion(mult = c(0, .1))) +
    labs(title=paste0("PRS Divergence: ",description),
         subtitle = p_text)
  gg
}
# plot save settings for each plot (in pixels)
width <- 2200
height <- 2200

# Plots all traits' divergence plot and saves to [dir_out]/divergence_plots/
dir.create(paste0(dir_out,"divergence_plots/"), recursive = TRUE)
for (code in traits_table$prive_code) {
  gg <- plot_divergence(code)
  loc_out <- paste0(dir_out,"divergence_plots/divergence_plot_",code,".png")
  ggsave(loc_out,plot=gg,width=width,height=height,units="px")
  print(paste("Saved divergence plot for",code))
}

# Makes Figure 2: Divergence plot
figure_codes <- c("250.1","darker_skin0")
plots <- list()
for (i in 1:length(figure_codes)) {
  code <- figure_codes[i]
  gg <- plot_divergence(code)
  plots[[i]] <- gg
}

gg_fig <- ggarrange(plots[[1]],plots[[2]],
                    common.legend = TRUE,
                    legend = "bottom")

loc_out <- paste0(dir_out,"figure_2_divergence_plots.png")
ggsave(loc_out,plot=gg_fig,width=2*width,height=height,units="px")
