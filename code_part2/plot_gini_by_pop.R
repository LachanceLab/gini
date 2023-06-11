# plot_gini_by_pop.R

# plots comparison of Gini values using different population allele frequencies

### Libraries and directories ####
library(tidyverse)
library(data.table)
library(ggpubr)

# sets working directory
setwd("./")

# sets location to traits_table.txt
loc_table <- "../generated_data/traits_table.txt"
# sets directory of generated figures to output to
dir_out <- "../generated_figures/"


### Printing function ####
print_plot <- function(gg, loc_out, print_mode, plot_width, plot_height, sf) {
  if (print_mode == "png") {
    png(loc_out, width = plot_width*sf, height = plot_height*sf)
  } else if (print_mode == "pdf") {
    pdf(loc_out, width = plot_width*sf / 75, height = plot_height*sf / 75)
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

### Code ####
traits_table <- as_tibble(fread(loc_table))

pop_ginis2 <- traits_table %>%
  select(prive_code, GWAS_trait_type, PGS_trait_type, group_consolidated, #portability_index,
         starts_with("gini_"), pops, pops_pass_qc) %>%
  filter(GWAS_trait_type == "quantitative", PGS_trait_type=="quantitative")
#pops <- substring(colnames(traits_table %>% select(starts_with("gini_"))),6)[-1]
pops <- c("EUR","AMR","CSA","AFR","EAS") # "MID" ommitted

gini_UK_plots <- list()
for (i in 1:length(pops)) {
  pop2 <- pops[i]
  col_pop2 <- paste0("gini_",pop2)
  cor <- cor.test(pop_ginis2[,"gini_panUKB"][[1]], pop_ginis2[,col_pop2][[1]])
  p<-ggplot(data=pop_ginis2, aes(x = !!as.name(col_pop2), y = gini_panUKB, color=group_consolidated)) +
    geom_abline(slope=1,intercept=0, size=1*sf) +
    geom_point(alpha=0.75, size=4*sf) +
    scale_x_continuous(expand=c(0.01,0.01),limits=c(0,1)) +
    scale_y_continuous(expand=c(0.01,0.01),limits=c(0,1)) +
    xlab(bquote(Gini[100][','][.(pop2)])) +
    ylab(bquote(Gini[100][','][meta])) +
    labs(color = "Trait Type") +
    theme_light() +
    theme(aspect.ratio = 1,
          legend.position = c(0.99,0.01),
          legend.title.align = 0.5,
          legend.background = element_rect(fill = "transparent", color="gray90"),
          legend.justification = c(1,0)) +
    # scale_color_manual(labels=c("Binary","Quantitative"),
    #                    breaks=c("binary", "quantitative"),
    #                    values = c("binary"="#F8766D", "quantitative"="#00BFC4")) +
    annotate("text", x=0.05, y=0.95,hjust=0,vjust=1, parse=TRUE, size=10*sf,
             label = paste0("r==",round(cor$estimate[[1]],4))) +
    gini_p_theme
  if (i != length(pops)+1) {p <- p + theme(legend.position = "none")}
  gini_UK_plots[[i]] <- p
}
ukginis <- ggarrange(plotlist = gini_UK_plots, ncol = 3, nrow = 2)

loc_out <- paste0(dir_out,"gini_by_pop.", print_mode)
print_plot(ukginis, loc_out, print_mode, 4*500, 2*500, sf)
print(loc_out)

# ###
# 
# cor_ginis <- tibble(pop1 = character(), pop2 = character(), r = numeric(), p = numeric())
# cols_ginis <- colnames(pop_ginis2 %>% select(starts_with("gini_")))
# for (col1 in cols_ginis) {
#   pop1 <- substr(col1, 6, nchar(col1))
#   for (col2 in cols_ginis) {
#     pop2 <- substr(col2, 6, nchar(col2))
#     
#     cor1 <- cor.test(pop_ginis2[[col1]], pop_ginis2[[col2]])
#     cor_ginis <- cor_ginis %>% add_row(
#       pop1=pop1, pop2=pop2, r=cor1$estimate, p=cor1$p.value
#     )
#   }
# }
# cor_ginis <- cor_ginis %>%
#   filter(pop1 != "panUKB", pop2 != "panUKB", pop1 != pop2)
# mean((cor_ginis %>% filter(pop1=="meta", pop2 != "EUR"))$r)

####
library(GGally)
popxpop <- ggpairs(pop_ginis2, mapping=aes(color=group_consolidated),
                   columns = paste0("gini_",pops)) +
  labs(title="not finalized plot")

loc_out <- paste0(dir_out,"gini_by_pop2.", print_mode)
print_plot(popxpop, loc_out, print_mode, 1.5*500, 1.5*500, sf)
ggsave(loc_out, popxpop,device="png", width=3500,height=3500, units="px")
print(loc_out)
