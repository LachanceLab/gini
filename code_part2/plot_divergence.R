library(tidyverse)
library(ggpubr)
library(data.table)

loc_PRSs <- "D:/Genee_local/divergence_plots/pop_sampled_PRSs.txt"
loc_table <- "D:/Genee_local/divergence_plots/traits_table.txt"
dir_out <- "D:/Genee_local/divergence_plots/"

traits_table <- as_tibble(fread(loc_table))
PRSs <- as_tibble(fread(loc_PRSs))
PRSs[PRSs$ancestry=="United","ancestry"] <- "UK"

plot_divergence <- function(code) {
  PRS_trait <- PRSs %>% select(Ancestry = ancestry, PRS = all_of(code))
  slice <- traits_table %>% filter(prive_code == code)
  description <- slice$description[1]
  f_stat <- slice$f_stat[1]
  p_value <- slice$p_value_f[1]
  if (p_value < 1E-320) {
    p_text <- "< 1E-320"
  } else {
    p_text <- paste0("= ", formatC(p_value, format="E", digits=2))
  }
  
  
  gg<-ggplot(PRS_trait, aes(x=PRS, fill=Ancestry)) +
    geom_density(color='#e9ecef', alpha=0.6, position='identity') +
    theme_light() +
    theme(legend.position = "bottom",
          legend.title=element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    xlab("Polygenic Score per UKBB Individual") +
    ylab("Density") +
    labs(title=description,
         subtitle = paste0("ANOVA: log10(F-stat) = ", round(log10(f_stat),2),
                           ". p-value ",p_text ))
  gg
}

# Plots all traits' divergence plot
for (code in traits_table$prive_code) {
  gg <- plot_divergence(code)
  loc_out <- paste0(dir_out,"divergence_plot_",code,".png")
  ggsave(loc_out,plot=gg,width=2200,height=2200,units="px")
  print(paste("Saved divergence plot for",code))
}

# Makes figure 2
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
ggsave(loc_out,plot=gg_fig,width=4400,height=2200,units="px")
