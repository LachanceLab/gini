# plot_scatterplot_matrix_density.R

# Plots the scatterplot matrix, the dual density plots, and the scatterplots
# between the six summary statistics and prevalence

### Libraries and directories ####
library(tidyverse)
library(data.table)
library(GGally)
library(rlang)

# sets working directory
setwd("./")

# sets location of trait_table generated in part1
loc_table <- "../generated_data/traits_table.txt"
# sets directory for generated figures
dir_out <- "../generated_figures/"
# sets scaling factors for image output. Default = 2
sf <- 2
p_adjust_method <- "fdr" # used in p.adjust()
print_mode <- "png" # set to either "png" or "pdf"
# columns to plot
vars <- c("ldpred2_h2","traitLD_unadj_EUR","gini_panUKB","pcor_United","portability_index", "log_F")

### Printing function ####
print_plot <- function(gg, loc_out, print_mode, plot_width, plot_height, sf) {
  # can print as pdf or as png
  if (print_mode == "png") {
    png(loc_out, width = plot_width*sf, height = plot_height*sf)
  } else if (print_mode == "pdf") {
    pdf(loc_out, width = plot_width*sf / 75, height = plot_height*sf / 75)
  }
  print(gg)
  dev.off()
  # I don't really know what this does but it fixes a bug I was running into
  dev.set(dev.next())
  dev.set(dev.next())
}

### Code ####

# reads traits table, filters out low prevalence traits, and defines lifestyle traits
traits_table <- as_tibble(fread(loc_table)) #%>% mutate(cMperMb = log10(cMperMb))
traits_table2 <- traits_table %>%
  #filter(PGS_trait_type != GWAS_trait_type) %>%
  filter(PGS_trait_type == "quantitative",
         GWAS_trait_type == "quantitative") %>%
  select(prive_code, description, PGS_trait_type, GWAS_trait_type, group,
         group_consolidated, prevalence, all_of(vars)) %>%
  mutate(lifestyle = group_consolidated == "lifestyle/psychological")
  

# Caps maximum portability to 0
#traits_table2[which(traits_table2$portability_index > 0),"portability_index"] <- 0
# Log10 transforms F-statistic (D statistic)
#traits_table2[,"f_stat"] <- log10(traits_table2[,"f_stat"])

# Helper function for writing p-values onto plots
color_p_significant <- "gray5"
color_p_nonsignificant <- "gray60"
p_value_to_text <- function(p_value) {
  p_text <- round(p_value,3)
  sci <- FALSE
  # uses scientific notation when < 0.001
  if (p_value < 0.001) {
    p_text <- formatC(p_value,format="E", digits=2)
    p_text_stem <- as.numeric(substr(p_text,1,4))
    p_text_exp <- as.numeric(substr(p_text,6,10))
    p_text <- paste0("adjusted~p==",p_text_stem,"%*%10^",p_text_exp)
  } else {
    p_text <- paste0("adjusted~p==",p_text)
  }
  
  # changes text to darker color if significant
  if (p_value < 0.05) {text_color = color_p_significant
  } else { text_color = color_p_nonsignificant}
  
  # returns p-value text and significance color 
  list(p_text,text_color)
}

### Scatterplot Matrix ####

# calculates adjusted p-values for correlation measurement between variables
p_values_cor <- tibble(
  var1 = as.character(),
  var2 = as.character(),
  cor = as.numeric(),
  unadj_p_value = as.numeric(),
  adj_p_value = as.numeric()
)
# loops through unique pairs of summary statistics
for (i in 1:(length(vars) - 1)) {
  x <- vars[[i]]
  for (j in (i+1):(length(vars))) {
    y <- vars[[j]]
    cor1 <- cor.test(as.data.frame(traits_table2)[,x],
                     as.data.frame(traits_table2)[,y])
    cor_value <- cor1$estimate[[1]]
    p_value <- cor1$p.value
    p_values_cor <- p_values_cor %>% add_row(
      var1 = x,
      var2 = y,
      cor = cor_value,
      unadj_p_value = p_value,
      adj_p_value = NA
    )
  }
}
# adjust p-values for multiple-testing
adj_p_values_cor <- p.adjust(p_values_cor$unadj_p_value,p_adjust_method)
p_values_cor$adj_p_value <- adj_p_values_cor

## Matrix subplot functions

# Top-right plots: correlation and p-value between measurements
upper_corr_p <- function(data,mapping) {
  # extracts x and y from ggpairs data argument
  x <- as.character(quo_get_expr(mapping[[1]]))
  y <- as.character(quo_get_expr(mapping[[2]]))
  
  # extracts r- and p-value from previous computations
  slice <- p_values_cor %>% filter( (var1==x & var2==y) | (var1==y & var2==x) )
  cor_value <- slice$cor
  p_value <- slice$adj_p_value
  
  p_text_list <- p_value_to_text(p_value)
  p_text <- p_text_list[[1]]
  
  # determines full text to display
  cor_text <- formatC(cor_value,digits=3, format="f")
  text_rline <- paste0("r==", cor_text)
  
  
  
  # makes GGally textplot using custom text
  p <- ggally_text(label="") +
    geom_text(aes(x=0.5,y=0.55),hjust=0.5,vjust=0,size=6*sf,color=p_text_list[[2]],
              label = text_rline, parse=TRUE ) + # correlation value
    geom_text(aes(x=0.5,y=0.45),hjust=0.5,vjust=1,size=6*sf,color=p_text_list[[2]],
              label = p_text, parse=TRUE ) + # p-value
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_rect(linetype = "solid", 
                                      color = theme_get()$panel.background$fill,
                                      fill = "transparent"))
  p
}

# Diagonal plots: display name and symbol of variable
var_labels <- list(
  "ldpred2_h2" = c("Heritability","({h^{2}}[SNP])"),
  #"cMperMb" = c("Recombination","Rate~(R)"),
  "traitLD_unadj_EUR" = c("Trait~LD","Scores~(L[EUR])"),
  "gini_panUKB" = c("Gini","(G[list(100,UK)])"),
  "pcor_United" = c("PGS~Accuracy","(symbol(r)[UK])"),
  "portability_index" = c("Portability","(m)"),
  "log_F" = c("Divergence","(D)"))
diag_label <- function(data, mapping) {
  variable <- as.character(quo_get_expr(mapping[[1]]))
  
  # makes GGally textplot using custom text
  p <- ggally_text(label="") +
    geom_text(aes(x=0.5,y=0.55), size=9*sf, hjust=0.5, vjust=0, color="black",
              label = var_labels[[variable]][1], parse = TRUE) + # top line
    geom_text(aes(x=0.5,y=0.45), size=9*sf, hjust=0.5, vjust=1, color="black",
              label = var_labels[[variable]][2], parse = TRUE) + # bottom line
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_rect(linetype = "solid", 
                                      color = "black",
                                      fill = "transparent"))
  p
}

# Bottom-left plots: scatterplot

# sets limits for each scatterplot
get_axis_lims <- function(vector,hard_min=NA,hard_max=NA) {
  vector <- vector[!is.na(vector)]
  lim_range <- diff(range(vector)) * 1.05
  lims <- mean(range(vector)) + c(-1,1) * lim_range/2
  if (!is.na(hard_min)) {lims[1] <- hard_min}
  if (!is.na(hard_max)) {lims[2] <- hard_max}
  lims
}
axis_lims <- list(
  "log_F" = get_axis_lims((traits_table %>% filter(PGS_trait_type=="quantitative"))$log_F),
  #"cMperMb"= get_axis_lims((traits_table %>% filter(GWAS_trait_type=="quantitative"))$cMperMb),
  "traitLD_unadj_EUR"= get_axis_lims((traits_table %>% filter(GWAS_trait_type=="quantitative"))$traitLD_unadj_EUR),
  "ldpred2_h2" = c(0,1),
  "pcor_United"= get_axis_lims((traits_table %>% filter(PGS_trait_type=="quantitative"))$pcor_United,0),
  "portability_index"= get_axis_lims((traits_table %>% filter(PGS_trait_type=="quantitative"))$portability_index,NA,0),
  "gini_panUKB"=c(0,1))
lm_scatterplot <- function(data, mapping) {
  # determines two variables being plotted
  x <- as.character(quo_get_expr(mapping[[1]]))
  y <- as.character(quo_get_expr(mapping[[2]]))
  
  # gets axis limits
  xlims <- axis_lims[[x]]
  ylims <- axis_lims[[y]]
  
  # extracts adjusted p-value and sets significance color
  p_value <- (p_values_cor %>% filter( (var1==x & var2==y) | (var1==y & var2==x) ))$adj_p_value
  if (p_value < 0.05) {linealpha <- 0.9}
  else {linealpha <- 0.4}
  
  p <- ggplot(data=data, mapping=mapping) +
    geom_line(stat="smooth", method="lm", color="dodgerblue1", formula=y~x, size=1*sf, alpha=linealpha) +
    geom_smooth(method="lm", linetype=0, formula=y~x, size=1*sf, alpha=linealpha/2) +
    geom_point(#aes(color=GWAS_trait_type),
               alpha=0.75,shape=19, size=1.75*sf, ) +
    xlim(xlims) +
    ylim(ylims) +
    theme_light()
  
  # log10 scales for h^2
  if (x == "ldpred2_h2") {p <- p + scale_x_log10(limits = c(0.0096,1),
                                                 breaks = c(0.01,0.1,1),
                                                 labels = c("0.01","0.10","1.00"))}
  if (y == "ldpred2_h2") {p <- p + scale_y_log10(limits = c(0.0096,1),
                                                 breaks = c(0.01,0.1,1),
                                                 labels = c("0.01","0.10","1.00"))}
  
  p
}
scplot_textsize <- 20

# makes actual scatterplot matrix
p_sc <- ggpairs(data = traits_table2,
                columns=vars,
                lower = list(continuous = lm_scatterplot),
                diag = list(continuous = diag_label),
                upper = list(continuous = upper_corr_p),
                axisLabels = "show"
) +
  theme(strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        text = element_text(size = scplot_textsize*sf),
        axis.text.x = element_text(size=(scplot_textsize*0.5)*sf),
        axis.text.y = element_text(size=(scplot_textsize*0.5)*sf))

# Saves image onto system
smplot_width <- 1200
smplot_height <- 1150
loc_out <- paste0(dir_out,"scatterplot_matrix.", print_mode)
print_plot(p_sc, loc_out, print_mode, smplot_width, smplot_height, sf)
print(paste0("Saved ",length(vars),"x",length(vars)," scatterplot matrix"))

### Dual Density Plots ####

column_labels <- c("Heritability","Trait LD Score","Gini","PGS Accuracy","Portability","Divergence")

# uses Wilcoxon-ranked test to compare means differences between types and groups
# for each of the 6 measurements
p_values_WRT <- tibble(
  var_measurement = as.character(),
  var_comparison = as.character(),
  unadj_p_value = as.numeric(),
  adj_p_value = as.numeric()
)
for (i in 1:2) {
  if (i==1) {
    var_comparison <- "group"
    # restricts to just quantitative to demonstrate relationships accurately
    subtable1 <- traits_table2 %>% filter(lifestyle) %>% select(vars)
    subtable2 <- traits_table2 %>% filter(!lifestyle) %>% select(vars)
  } else if (i==2) {
    var_comparison <- "type"
    subtable1 <- traits_table %>% filter(PGS_trait_type=="binary", GWAS_trait_type=="binary") %>% select(vars)
    subtable2 <- traits_table2 %>% select(vars)
  }
  
  for (j in 1:length(vars)) {
    # uses Wilcoxon signed-rank test to determine significant differences
    var_measurement <- colnames(subtable1)[j]
    WRT1 <- wilcox.test(subtable1[[j]], subtable2[[j]], paired=FALSE)
    p_value <- WRT1$p.value
    p_values_WRT <- p_values_WRT %>% add_row(
      var_measurement = var_measurement,
      var_comparison = var_comparison,
      unadj_p_value = p_value,
      adj_p_value = NA
    )
  }
  # Adjusts p-values for multiple-testing. Adjusts within 6x1 plot
  adj_p_values_WRT <- p.adjust(p_values_WRT$unadj_p_value[(length(vars)*i-(length(vars)-1)):(length(vars)*i)],p_adjust_method)
  p_values_WRT$adj_p_value[(length(vars)*i-(length(vars)-1)):(length(vars)*i)] <- adj_p_values_WRT
}

## Dual Density Plot function ####
dual_density <- function(data, mapping, the_var_comparison, the_var_measurement) {
  # extracts name of variable and xlim from scatterplot matrix
  x <- as.character(quo_get_expr(mapping[[1]]))
  xlims <- axis_lims[[x]]
  
  # gets adjusted p-values from Wilcoxon Ranked Test already done
  if (the_var_comparison=="group") {
    col_var <- "lifestyle"
    #data <- data %>% filter(GWAS_trait_type=="quantitative")
  } else if (the_var_comparison=="type") {col_var <- "GWAS_trait_type"}
  
  adj_p_value <- (p_values_WRT %>%
                    filter(var_measurement==the_var_measurement,
                           var_comparison==the_var_comparison))$adj_p_value
  
  # converts p-value to text version
  p_text_list <- p_value_to_text(adj_p_value)
  text <- p_text_list[[1]]
  
  # Adds density plots
  p <- ggplot(data=data, mapping=aes(x=!!as.name(x) )) +
    geom_density(mapping=aes(fill=!!as.name(col_var)),alpha=0.5,size=0.5*sf) +
    theme_light()
  
  # Adjusts legend to match the data
  if (the_var_comparison == "group") {
    p <- p +
      labs(fill="Trait Group") +
      scale_fill_manual(labels=c("Lifestyle/Psychological","Non-lifestyle/psychological"),
                        breaks=c(TRUE, FALSE),
                        values = c("TRUE"="#00BF7D", "FALSE"="gray20"))
  } else if (the_var_comparison == "type") {
    p <- p +
      labs(fill="Trait Type") +
      scale_fill_manual(labels=c("Binary","Quantitative"), 
                        breaks=c("binary", "quantitative"),
                        values = c("binary"="#F8766D", "quantitative"="#00BFC4"))
  }
  # log10 scales for h^2
  if (x == "ldpred2_h2") {p <- p + scale_x_log10(limits = c(0.0096,1))}
  else {p <- p + xlim(xlims)}
  # Adds p-value to plot
  # pads top to allow space for p-value
  yrange <- layer_scales(p)$y$range$range
  padding <- 1.20
  yrange[2] <- (yrange[2]-yrange[1]) * padding + yrange[1]
  p <- p +
    scale_y_continuous(limits = yrange, expand=expansion(mult = c(0, .05))) +
    annotate("text",
             x = ifelse(x=="ldpred2_h2",0.1,mean(xlims)),
             y=diff(yrange)*((0.95+padding)/(2*padding)),
             label = text, parse=TRUE, vjust=0, hjust=0.5, size=4*sf,
             color=p_text_list[[2]])
  
  p
}
# Actually generates and saves dual density plots to system
ddplot_width <- 1200
ddplot_height <- 200
ddplot_textsize <- 20
for (var_comparison in c("group","type")) {
  density_plots <- list()
  for (i in 1:length(vars)) {
    var_measurement <- vars[i]
    #print(paste(var_comparison, var_measurement))
    
    if (var_comparison=="group") {
      density_plot <- dual_density(data=traits_table2,
                                   mapping = aes(x=!!as.name(var_measurement)),
                                   var_comparison, var_measurement)
    } else {
      density_plot <- dual_density(data=traits_table %>% filter(PGS_trait_type==GWAS_trait_type),
                                   mapping = aes(x=!!as.name(var_measurement)),
                                   var_comparison, var_measurement)
    }
    
    density_plots[[i]] <- density_plot
  }
  
  ddp <- ggmatrix(plots=density_plots,
                  nrow=1,
                  ncol=length(vars),
                  legend=1,
                  xAxisLabels = column_labels) +
    theme(legend.position="top",
          legend.key.size = unit(1,"cm"),
          legend.text = element_text(size = ddplot_textsize*sf),
          axis.text.x = element_text(size=(ddplot_textsize*0.5)*sf),
          axis.text.y = element_blank(),
          text = element_text(size = ddplot_textsize*sf),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  loc_out <- paste0(dir_out,"dual_density_plot_",var_comparison,".", print_mode)
  print_plot(ddp, loc_out, print_mode, ddplot_width, ddplot_height, sf)
  
  print(paste("Saved dual density plots for",var_comparison))
}


### Prevalence Plots ###

# calculates adjusted p-values for correlation measurement between prevalence 
# and summary statistics
p_values_cor2 <- tibble(
  var1 = as.character(),
  var2 = as.character(),
  cor = as.numeric(),
  unadj_p_value = as.numeric(),
  adj_p_value = as.numeric()
)
for (j in 1:(length(vars))) {
  y <- vars[[j]]
  cor1 <- cor.test(log10(as.data.frame(traits_table)[,"prevalence"]),
                   as.data.frame(traits_table)[,y])
  cor_value <- cor1$estimate[[1]]
  p_value <- cor1$p.value
  p_values_cor2 <- p_values_cor2 %>% add_row(
    var1 = "log10_prevalence",
    var2 = y,
    cor = cor_value,
    unadj_p_value = p_value,
    adj_p_value = NA
  )
}
# uses False Discovery Rate to adjust p-values
adj_p_values_cor2 <- p.adjust(p_values_cor2$unadj_p_value,p_adjust_method)
p_values_cor2$adj_p_value <- adj_p_values_cor2

# Actually generates and saves dual density plots to system
prevalence_plots <- list()
for (i in 1:length(vars)) {
  var_measurement <- vars[i]
  xlims <- axis_lims[[var_measurement]]
  
  adj_p_value <- (p_values_cor2 %>% filter(var2==var_measurement))$adj_p_value
  if (adj_p_value < 0.05) {linealpha <- 0.9
  } else {linealpha <- 0.4}
  
  # converts p-value to text version
  p_text_list <- p_value_to_text(adj_p_value)
  text <- p_text_list[[1]]
  
  p <- ggplot(traits_table %>% filter(PGS_trait_type == "binary", GWAS_trait_type=="binary"),
              aes(x = !!as.name(var_measurement), y = log10(prevalence))) +
    geom_line(stat="smooth", method="lm", color="#F8766D", formula=y~x, size=1*sf, alpha=linealpha) +
    geom_smooth(method="lm", linetype=0, formula=y~x, size=1*sf, alpha=linealpha/2) +
    geom_point(alpha=0.75,shape=19, size=1.75*sf) +
    xlab(var_labels[[var_measurement]]) +
    ylab(bquote(Log[10](Prevalence))) +
    theme_light()
  if ( var_measurement == "ldpred2_h2") {
    p <- p + scale_x_log10(limits = c(0.0096,1),
                           breaks = c(0.01,0.1,1),
                           labels = c("0.01","0.10","1.00"))
  } else {p <- p + xlim(xlims)}
  # Adds p-value to plot
  yrange <- layer_scales(p)$y$range$range
  padding <- 1.15
  yrange[2] <- (yrange[2]-yrange[1]) * padding + yrange[1] # pads top to allow space for p-value
  p <- p +
    scale_y_continuous(limits = yrange) +
    annotate("text",
             x = ifelse(var_measurement=="ldpred2_h2",0.1,mean(xlims)),
             y=(diff(yrange))*((1+padding)/(2*padding)) + yrange[1],
             label = text, parse=TRUE, vjust=0, hjust=0.5, size=4*sf,
             color=p_text_list[[2]])
  p
  prevalence_plots[[i]] <- p
}

# combines plots into one figure
pscp <- ggmatrix(plots=prevalence_plots,
                nrow=1,
                ncol=length(vars),
                xAxisLabels = NULL,
                ylab = bquote(Log[10](Prevalence))) +
  theme(legend.position="top",
        legend.key.size = unit(1,"cm"),
        legend.text = element_text(size = ddplot_textsize*sf),
        axis.text = element_text(size=(ddplot_textsize*0.5)*sf),
        axis.title.y = element_text(size=(ddplot_textsize*0.75)*sf),
        text = element_text(size = ddplot_textsize*sf),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# prints figure
loc_out <- paste0(dir_out,"prevalence_scatterplot.", print_mode)
print_plot(pscp, loc_out, print_mode, ddplot_width+100, ddplot_height, sf)
print("Saved prevalence scatterplot")

