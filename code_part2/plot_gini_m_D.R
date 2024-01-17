# plot_gini_m_D

# Makes figure 2, which consists of showing gini (Lorenz Curve), portability,
# and PGS Divergence

### Libraries and directories ####
library(tidyverse)
library(data.table)
library(scales)
library(ggpubr)
library(ggrepel)
library(extrafont)
loadfonts(device = "win", quiet=TRUE)
loadfonts(device = "pdf", quiet=TRUE)
source("../code_part1/helper_functions/helper_functions.R")

# sets working directory
setwd("./")

# Sets the directory of the summary files with appended allele frequencies
dir_sfs <- "../generated_data/panUKB_sf/"
# Sets the location of the sampled individuals' PGSs
loc_PGSs <- "../generated_data/pop_sampled_PGSs.txt"
# sets the location of the traits table
loc_table <- "../generated_data/traits_table.txt"
# sets directory of outputted figures
dir_out <- "../generated_figures/"

# sets scaling factors for image output. Default = 2
sf <- 2
print_mode <- "pdf" # set to either "png" or "pdf"

### Functions ####

# function that sets a common theme for plots:
common_theme <- theme_light() +
  theme(
  aspect.ratio = 1,
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(colour = "black", fill=NA, size=1),
  plot.title = element_text(hjust = 0.5, size=15*sf),
  axis.title = element_text(size=14*sf),
  axis.text = element_text(size=13*sf),
  legend.title = element_text(size=13*sf),
  legend.text = element_text(size=11*sf),
  plot.margin = unit(0.25*sf*c(1.2,1.4,1.2,1.2), "cm")
)

# function that reads a trait's summary file and extracts the needed columns for plotting
cleanup_data_lorenz <- function(code, threshold=500) {
  loc_sf <- paste0(dir_sfs,code,"_sf_indep.txt")
  
  sf <- as_tibble(fread(loc_sf)) %>%
    select(gvc = gvc_meta2use) %>% drop_na() %>%
    arrange(-gvc) %>%
    filter(row_number() <= threshold) %>% arrange(gvc) %>%
    mutate(gvc_csum = cumsum(gvc))
  sf$gvc_cshare <- sf$gvc_csum / sf$gvc_csum[nrow(sf)]
  if (nrow(sf) < threshold) {
    sf <- sf %>%
      add_row(
        gvc = rep(0, threshold - nrow(sf)),
        gvc_csum = rep(0, threshold - nrow(sf)),
        gvc_cshare = rep(0, threshold - nrow(sf)),
      ) %>% arrange(h2)
  }
  sf <- sf %>% mutate(percentile = dplyr::row_number() / nrow(sf))
}
# function that actually plots Lorenz curve
plot_lorenz <- function(code, sfile) {
  
  slice <- traits_table %>% filter(prive_code == code)
  description <- slice$description
  gini <- slice$gini_panUKB
  
  title <- paste0(description)
  # sets proper math formatting for plot annotation
  gini_text <- formatC(gini[[1]],digits=3, format="f")
  text <- paste0("italic(G[list(500,Meta)])==",gini_text)
  
  # makes Lorenz curve plot
  gg <- ggplot(sfile, aes(x=100*percentile, y=gvc_cshare)) +
    geom_col(position = position_nudge(-0.5), fill="gray20", width=0.8) +
    geom_abline(slope=1/100,color="dodgerblue1", size=0.5*sf) +
    labs(title=title) +
    xlab(expression(paste("Percentile of SNP ",italic(gvc)))) +
    ylab(expression(paste("Cumulative sum of SNP ",italic(gvc)))) +
    common_theme +
    scale_x_continuous(limits=c(0,100), expand = c(0.0,0.0)) +
    scale_y_continuous(limits=c(0,1), expand = c(0,0)) +
    annotate("text", x = 0.5 * 100, y = 0.975*(1), label = text,
             parse=TRUE, vjust=1, hjust=0.5, size=5*sf, family = "Georgia")
  gg
}
# function that plots portability
plot_portability <- function(code) {
  
  # properly formats tibble with PGS accuracy to compute relative PGS accuracy
  pcor_data <- traits_table %>% filter(prive_code == code) %>%
    select(starts_with("PGS_R2_")) %>%
    pivot_longer(
      cols = starts_with("PGS_R2_"),
      names_prefix = "PGS_R2_",
      names_to = "ancestry",
      values_to = "PGS_R2"
    ) %>% mutate(
      ancestry = factor(str_replace(ancestry, 'United','UK'), levels = distances$ancestry), 
      relative_PGS_R2 = PGS_R2 / (traits_table %>% filter(prive_code == code))$PGS_R2_United[1],
    ) %>% left_join(distances, by="ancestry")
  
  # extracts trait information (including portability slope)
  slice <- traits_table %>% filter(prive_code == code)
  m <- slice$portability_index
  description <- slice$description
  if (code == "haemoglobin") {description <- "Hemoglobin concentration"}
  # properly formats annotation
  if ((m < 0.001) & (m != 0)) {
    exponent <- floor(log10(abs(m)))
    base <- signif(m, digits = 3) / 10^exponent
    text <- paste0("italic(m)==",base,"%*%~10^",exponent)
  } else {
    text <- paste0("italic(m)==",round(m,3))
  }
  
  gg<-ggplot(pcor_data, aes(x=prive_dist_to_UK)) +
    geom_segment(aes(x=0,y=1, xend=max(prive_dist_to_UK), yend=1 + m * max(prive_dist_to_UK)),
                 size=0.5*sf, color="dodgerblue1") +
    geom_point(aes(y = relative_PGS_R2, color=ancestry), size = 2.5*sf) +
    scale_x_continuous(expand=expansion(mult = c(0.02, .02))) +
    common_theme +
    theme(legend.position = "none") +
    xlab("Genetic PC Distance to UK") +
    ylab("PGS Accuracy Relative to UK") +
    labs(title = description)
  
  # annotation position and scales depend on existing plot scales
  xrange <- layer_scales(gg)$x$range$range
  yrange <- layer_scales(gg)$y$range$range
  if (yrange[1] > 0) {yrange[1] <- 0}
  yrange[2] <- 1.15 #manually pads y-max of plot
  gg <- gg +
    scale_y_continuous(limits = c(yrange[1],yrange[2]),
                       expand=expansion(mult = c(0, .01)),
                       breaks = c(0,0.25,0.5,0.75,1),
                       label = format(c(0,0.25,0.5,0.75,1))) +
    annotate("text",
             x=0.5*(xrange[2]-xrange[1]),
             y = 0.985*(yrange[2]-0), label = text, family = "Georgia",
             parse=TRUE, vjust=1, hjust=0.5, size=5*sf)
  gg
}
# function that generates divergence plot
plot_divergence <- function(code) {
  
  # extracts info about trait
  PGS_trait <- PGSs %>% select(Ancestry = "pop", PGS = all_of(code))
  slice <- traits_table %>% filter(prive_code == code)
  description <- slice$description[1]
  if (code == "darker_skin0") {description <- "Skin color"}
  log_F <- slice$log_F[1]
  logfstat <- formatC(log_F,digits=2, format="f")
  text <- paste0("italic(D)==",logfstat)
  
  # plots divergence
  gg<-ggplot(PGS_trait, aes(x=PGS, fill=Ancestry)) +
    geom_density(color='#e9ecef', alpha=0.6, position='identity') +
    common_theme +
    theme(legend.position = "bottom",
          legend.title.align = 0.5,
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    #xlab("Polygenic Score per UKBB Individual") +
    xlab("Polygenic Score") +
    ylab("Density") +
    scale_x_continuous(expand=expansion(mult = c(0, 0))) +
    labs(title=paste0(description))
  
  # annotation position and scales depend on existing plot scales
  xrange <- layer_scales(gg)$x$range$range
  yrange <- layer_scales(gg)$y$range$range
  yrange[2] <- yrange[2] * 1.1
  gg <- gg +
    scale_y_continuous(limits = c(0,yrange[2]), expand=expansion(mult = c(0, 0.01))) +
    annotate("text",
             x=(xrange[2]+xrange[1])/2, y = 0.985*(yrange[2]-0), label = text,
             parse=TRUE, vjust=1, hjust=0.5, size=5*sf, family = "Georgia")
  gg
}
# function that generates traitLD bar plot
plot_traitLD <- function(code) {
  # extracts info about trait
  slice <- traits_table %>% filter(prive_code==code)
  description <- slice$description
  traitLD_tbl <- pop_LDs %>% filter(prive_code == code)
  
  ylims <- c(slice$traitLD_unadj_mean - 2 * slice$traitLD_unadj_mean * max_CoV,
             slice$traitLD_unadj_mean + 2 * slice$traitLD_unadj_mean * max_CoV)
  
  # adds text for traitLD_unadj_CV
  text <- paste0("LDCV==",round(slice$traitLD_unadj_CoV,3))
  gg <- ggplot(traitLD_tbl, aes(x = factor(pop, levels=pops), y = traitLD_unadj)) +
    common_theme +
    geom_hline(yintercept = slice$traitLD_unadj_mean, color="gray", size=2) +
    geom_col(aes(fill=pop)) +
    scale_y_continuous(limits=ylims,oob = rescale_none, expand = c(0,0)) +
    xlab("Continental Population") +
    ylab("LD Score (trait mean)") +
    labs(title = description, fill = "Ancestry") +
    scale_fill_manual(#values = hue_pal()(5),
                      values = c('#F8766D','#00BFC4',"#00A9FF","#C77CFF",'#FF61CC'),
                      breaks = pops) +
    theme(axis.ticks.x = element_blank()) +
    annotate("text",
             x=pops[3], y = 0.975*diff(ylims) + ylims[1], label = text,
             parse=TRUE, vjust=1, hjust=0.5, size=5*sf, family = "Georgia")
  
  gg
}
### Code ####

# reads traits table
traits_table <- as_tibble(fread(loc_table))

#pops <- c("AFR","CSA","EAS","AMR","EUR")
pops <- c('EUR','CSA','EAS','AMR','AFR')
pop_LDs <- traits_table %>% filter(GWAS_trait_type=="quantitative",PGS_trait_type=="quantitative")
max_CoV <- max(pop_LDs$traitLD_unadj_CoV)
pop_LDs <- pop_LDs %>%
  select(prive_code,any_of(paste0("traitLD_unadj_",pops))) %>% pivot_longer(
    cols = starts_with("traitLD_unadj_"),
    names_to = "pop",
    names_prefix = "traitLD_unadj_",
    values_to = "traitLD_unadj"
  )

# reads pop_centers, which contains PC distance information
pop_centers <- read.csv(
  "https://raw.githubusercontent.com/privefl/UKBB-PGS/main/pop_centers.csv",
  stringsAsFactors = FALSE)
ancestries <- sort(pop_centers$Ancestry)
ancestries[9] <- "UK" # in order to match other data
# calculates average PCA distance between ancestries
prive_PC <- pop_centers %>% select(PC1:PC16)
prive_dist_to_UK <- as.matrix(dist(prive_PC))[,1]
distances <- tibble(
  ancestry = pop_centers$Ancestry,
  prive_dist_to_UK = prive_dist_to_UK
) %>% arrange(prive_dist_to_UK)
distances$ancestry <- str_replace(distances$ancestry,"United Kingdom","UK")
distances$ancestry <- factor(distances$ancestry, levels = distances$ancestry)

# reads the PGSs and changes "United" to "UK"
PGSs <- as_tibble(fread(loc_PGSs)) %>% filter(pop != "Ashkenazi")
PGSs[PGSs$pop=="United","pop"] <- "UK"
PGSs$pop <- factor(PGSs$pop, levels = distances$ancestry)

## makes Lorenz plots
low_gini_code <- "height"
low_gini_sf <- cleanup_data_lorenz(low_gini_code)
low_gini_plot <- plot_lorenz(low_gini_code, low_gini_sf)

high_gini_code <- "log_potassium_urine"
high_gini_sf <- cleanup_data_lorenz(high_gini_code)
high_gini_plot <- plot_lorenz(high_gini_code, high_gini_sf)

## makes portability plots
high_m_code <- "cholesterol"
high_m_plot <- plot_portability(high_m_code)

low_m_code <- "haemoglobin"
low_m_plot <- plot_portability(low_m_code)

## makes divergence plots
low_D_code <- "log_platelet_crit"
low_D_plot <- plot_divergence(low_D_code) +
  theme(legend.justification = c(0,1),
        legend.position = c(0.01, 0.99),
        legend.background = element_rect(fill="transparent"))

high_D_code <- "darker_skin0"
high_D_plot <- plot_divergence(high_D_code) +
  theme(legend.position = "none")

# makes traitLD plots
low_CV_code <- "log_glucose"
low_CV_plot <- plot_traitLD(low_CV_code) +
  theme(legend.position = "none")
  # theme(legend.justification = c(0,1),
  #       legend.position = c(0.01, 0.99),
  #       legend.background = element_rect(fill="transparent"))

high_CV_code <- "neuroticism"
high_CV_plot <- plot_traitLD(high_CV_code) +
  theme(legend.position = "none")

plots <- list(low_CV_plot, NULL, low_gini_plot , NULL, low_m_plot, NULL, low_D_plot, 
              high_CV_plot, NULL, high_gini_plot, NULL, high_m_plot, NULL, high_D_plot)

ncol = 7 # NULL plots used for extra spacing
nrow = 2
gg <- ggarrange(plotlist = plots, ncol = ncol, nrow = nrow,
                widths = c(1, 0.075, 1, 0.075, 1, 0.075, 1))

# plot save settings for each plot (in pixels)
plot_width <- 1200
plot_height <- 1200

# saves as png and pdf
loc_out <- paste0(dir_out,"gini_m_D")
ggsave(paste0(loc_out,".png"),
       width=4*plot_width*sf,
       height=2*plot_height*sf,units="px")
cairo_pdf(file = paste0(loc_out,".pdf"),
          width = plot_width*sf / 75,
          height = 0.5*plot_height*sf / 75)
gg
dev.off()
