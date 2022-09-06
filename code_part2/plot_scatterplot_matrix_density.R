# plot_scatterplot_matrix_density.R

# Plots the scatterplot matrix and the dual density plots


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
vars <- c("ldpred2_h2","cMperMb","gini_United","pcor_United","portability_index", "f_stat")

### Code ####
traits_table <- as_tibble(fread(loc_table)) %>%
  select(prive_code, description, trait_type, group, group_consolidated,
         all_of(vars)) %>%
  mutate(lifestyle = group_consolidated == "lifestyle/psychological") %>%
  drop_na()


# Caps maximum portability to 0
traits_table[which(traits_table$portability_index > 0),"portability_index"] <- 0
# Log10 transforms F-statistic
traits_table[,"f_stat"] <- log10(traits_table[,"f_stat"])

# Helper function for writing p-values onto plots
digits <- 3 # how many digits to round p-values to (non-scientific notation)
color_p_significant <- "gray20"
color_p_nonsignificant <- "gray50"
p_value_to_text <- function(p_value) {
  p_text <- round(p_value,digits)
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
  if (p_value < 0.05) {
    text_color = color_p_significant
  } else {
    text_color = color_p_nonsignificant
  }
  
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
for (i in 1:(length(vars) - 1)) {
  x <- vars[[i]]
  for (j in (i+1):(length(vars))) {
    y <- vars[[j]]
    cor1 <- cor.test(as.data.frame(traits_table)[,x],
                     as.data.frame(traits_table)[,y])
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
# uses False Discovery Rate to adjust p-values
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
  text_rline <- paste0("r = ", round(cor_value,digits),"\n")
  
  
  # makes GGally textplot using custom text
  p <- ggally_text(label=text_rline, color=p_text_list[[2]], size=10*sf) +
    geom_text(aes(x=0.5,y=0.42),hjust=0.5,vjust=1,size=6*sf,color=p_text_list[[2]],
              label = p_text, parse=TRUE ) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_rect(linetype = "solid", 
                                      color = theme_get()$panel.background$fill, fill = "transparent"))
  p
}

# Diagonal plots: display name and symbol of variable
var_labels <- list(
  "ldpred2_h2" = c("Heritability","({h^{2}}[SNP])"),
  "cMperMb" = c("Recombination","Rate~(R)"),
  "gini_United" = c("Polygenicity","(Gini[list(100,UK)])"),
  "pcor_United" = c("Prediction","(symbol(r)[UK])"),
  "portability_index" = c("Portability","(m)"),
  "f_stat" = c("Divergence","(log[10](F))"))
diag_label <- function(data, mapping) {
  variable <- as.character(quo_get_expr(mapping[[1]]))
  
  # makes GGally textplot using custom text
  p <- ggally_text(label="") +
    geom_text(aes(x=0.5,y=0.55), size=9*sf, hjust=0.5, vjust=0, color="black",
              label = var_labels[[variable]][1], parse = TRUE) +
    geom_text(aes(x=0.5,y=0.45), size=9*sf, hjust=0.5, vjust=1, color="black",
              label = var_labels[[variable]][2], parse = TRUE) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_rect(linetype = "solid", 
                                      color = "black",
                                      fill = "transparent"))
  p
}

# Bottom-left plots: scatterplot
get_axis_lims <- function(vector,hard_min=NA,hard_max=NA) {
  vector <- vector[!is.na(vector)]
  lim_range <- diff(range(vector)) * 1.05
  lims <- mean(range(vector)) + c(-1,1) * lim_range/2
  if (!is.na(hard_min)) {lims[1] <- hard_min}
  if (!is.na(hard_max)) {lims[2] <- hard_max}
  lims
}
axis_lims <- list(
  "f_stat" = get_axis_lims(traits_table$f_stat),
  "cMperMb"= get_axis_lims(traits_table$cMperMb),
  "ldpred2_h2" = c(0,1),
  "pcor_United"= get_axis_lims(traits_table$pcor_United,0),
  "portability_index"= get_axis_lims(traits_table$portability_index,NA,0),
  "gini_United"=c(0,1))
# format_axis_sci <- function(x) {
#   formatC(x,format="E", digits=1, drop0trailing=TRUE)
# }
lm_scatterplot <- function(data, mapping) {
  x <- as.character(quo_get_expr(mapping[[1]]))
  y <- as.character(quo_get_expr(mapping[[2]]))
  
  xlims <- axis_lims[[x]]
  ylims <- axis_lims[[y]]
  
  p <- ggplot(data=data, mapping=mapping) +
    geom_smooth(method="lm", color="dodgerblue1", formula=y~x, size=1*sf) +
    geom_point(alpha=0.5,shape=19, size=1.75*sf) +
    xlim(xlims) +
    ylim(ylims) +
    theme_light()
  # log10 scales for h^2 and portability
  if (x == "ldpred2_h2") {p <- p + scale_x_log10()}
  if (y == "ldpred2_h2") {p <- p + scale_y_log10()}
  #if (x == "portability_index") {p <- p + scale_x_continuous(labels = format_axis_sci)}
  #if (y == "portability_index") {p <- p + scale_y_continuous(labels = format_axis_sci)}
  
  p
}
scplot_textsize <- 20
# Plots actual scatterplot matrix
p_sc <- ggpairs(data = traits_table,
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
if (print_mode == "png") {
  png(loc_out, width = smplot_width*sf, height = smplot_height*sf)
} else if (print_mode == "pdf") {
  pdf(loc_out, width = smplot_width*sf / 75, height = smplot_height*sf / 75)
}
print(p_sc)
dev.off()
print(paste0("Saved ",length(vars),"x",length(vars)," scatterplot matrix"))

### Dual Density Plots ####

column_labels <- c("Heritability","Recombination Rate","Polygenicity","Prediction","Portability","Divergence")

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
    subtable1 <- traits_table %>% filter(trait_type=="quantitative") %>% filter(lifestyle) %>% select(vars)
    subtable2 <- traits_table %>% filter(trait_type=="quantitative") %>% filter(!lifestyle) %>% select(vars)
  } else if (i==2) {
    var_comparison <- "type"
    subtable1 <- traits_table %>% filter(trait_type=="binary") %>% select(vars)
    subtable2 <- traits_table %>% filter(trait_type=="quantitative") %>% select(vars)
  }
  
  for (j in 1:length(vars)) {
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
  # Uses False Discovery Rate to adjust p-values. Adjusts within 5x1 plot
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
    data <- data %>% filter(trait_type=="quantitative")}
  else if (the_var_comparison=="type") {col_var <- "trait_type"}
  
  adj_p_value <- (p_values_WRT %>%
                    filter(var_measurement==the_var_measurement,
                           var_comparison==the_var_comparison))$adj_p_value
  
  # converts p-value to text version
  p_text_list <- p_value_to_text(adj_p_value)
  text <- p_text_list[[1]]
  
  # Adds density plots
  p <- ggplot(data=data, mapping=aes(x=!!as.name(x) )) +
    geom_density(mapping=aes(fill=!!as.name(col_var)),alpha=0.5,size=0.5*sf) +
    xlim(xlims) +
    theme_light()
  
  # Adjusts legend to match the data
  if (the_var_comparison == "group") {
    p <- p +
      labs(fill="Trait Group") +
      scale_fill_manual(labels=c("Lifestyle/Psychological","Non-lifestyle/psychological"),
                        breaks=c(TRUE, FALSE),
                        values = c("TRUE"="dodgerblue1", "FALSE"="gray20"))
  } else if (the_var_comparison == "type") {
    p <- p +
      labs(fill="Trait Type") +
      scale_fill_manual(labels=c("Binary","Quantitative"), 
                        breaks=c("binary", "quantitative"),
                        values = c("binary"="gray70", "quantitative"="gray10"))
  }
  # log10 scales for h^2 and portability
  if (x == "ldpred2_h2") {p <- p + scale_x_log10()}
  #if (x == "portability_index") {p <- p + scale_x_continuous(labels = format_axis_sci)}
  # Adds p-value to plot
  is_h2 <- x=="ldpred2_h2"
  xlims <- axis_lims[[x]]
  yrange <- layer_scales(p)$y$range$range
  padding <- 1.20
  yrange[2] <- yrange[2] * padding # pads top to allow space for p-value
  p <- p +
    ylim(yrange) +
    geom_text(data=NULL, label=text, parse=TRUE,
              aes(x=ifelse(x=="ldpred2_h2",0.01,min(xlims)), y=max(yrange)*((0.95+padding)/(2*padding))),
              vjust=0, hjust=0, color=p_text_list[[2]], size = 4*sf)
  
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
    
    density_plot <- dual_density(data=traits_table,
                                 mapping = aes(x=!!as.name(var_measurement)),
                                 var_comparison, var_measurement)
    
    density_plots[[i]] <- density_plot
  }
  
  ddp <- ggmatrix(plots=density_plots,
                  nrow=1,
                  ncol=length(vars),
                  legend=1,
                  xAxisLabels = column_labels) +
    theme(legend.position="bottom",
          axis.text.y = element_blank(),
          legend.key.size = unit(1,"cm"),
          legend.text = element_text(size = ddplot_textsize*sf),
          text = element_text(size = ddplot_textsize*sf),
          axis.text.x = element_text(size=(ddplot_textsize*0.5)*sf))
  
  loc_out <- paste0(dir_out,"dual_density_plot_",var_comparison,".", print_mode)
  if (print_mode == "png") {
    png(loc_out, width = ddplot_width*sf, height = ddplot_height*sf)
  } else if (print_mode =="pdf") {
    pdf(loc_out, width = ddplot_width*sf / 75, height = ddplot_height*sf / 75)
  }
  print(ddp)
  dev.off()
  
  print(paste("Saved dual density plots for",var_comparison))
}