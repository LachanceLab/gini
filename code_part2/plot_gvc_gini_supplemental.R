#Summed gvc for top100 vs. total 
rm(list = ls())
dev.off()

#Load libraries 
library(ggplot2)
library(rio)
library(tidyverse)
library(data.table)

# sets working directory
setwd("./")

#Location of the traits table 
loc_table <- "./Desktop/Gini-PGS/traits_table.txt" #Remove in final version
loc_table <- "../generated_data/traits_table.txt" 

# reads traits table
traits_table <- as.data.frame(fread(loc_table))

#Filter out problematic binary traits from traits_table
remove <- c("Malignant neoplasm of testis", 
            "Cancer of bladder", 
            "Cancer of brain", 
            "Thyroid cancer", 
            "Polycythemia vera", 
            "Nontoxic multinodular goiter", 
            "Thyrotoxicosis with or without goiter", 
            "Type 1 diabetes", 
            "Diabetic retinopathy", 
            "Hypoglycemia", 
            "Gout", 
            "Disorders of iron metabolism", 
            "Disorders of bilirubin excretion", 
            "Congenital deficiency of other clotting factors (including factor VII)", 
            "Dementias", 
            "Alzheimer's disease", 
            "Multiple sclerosis", 
            "Retinal detachments and defects", 
            "Macular degeneration (senile) of retina NOS", 
            "Corneal dystrophy", 
            "Peripheral vascular disease, unspecified", 
            "Nasal polyps", 
            "Appendiceal conditions", 
            "Celiac disease", 
            "Other chronic nonalcoholic liver disease", 
            "Rhesus isoimmunization in pregnancy", 
            "Lupus (localized and systemic)", 
            "Psoriasis", 
            "Sarcoidosis", 
            "Rheumatoid arthritis", 
            "Ankylosing spondylitis", 
            "Polymyalgia Rheumatica", 
            "Ganglion and cyst of synovium, tendon, and bursa", 
            "Contracture of palmar fascia [Dupuytren's disease]")

traits_table <- traits_table[!(traits_table$description %in% remove),]

#Revise the trait groups 
groups_consolidated <- list(
  "Diseases" = c("circulatory system","dermatologic","digestive",
                 "endocrine/metabolic","genitourinary","hematopoietic",
                 "musculoskeletal","neoplasms","neurological",
                 "psychiatric disorders","respiratory","sense organs",         
                 "symptoms"),
  "Biological measures" = c("biological measures"),
  "Lifestyle/psychological" = c("lifestyle and environment", "psychiatric disorders"),
  "Physical measures" = c("injuries & poisonings","physical measures","sex-specific factors")
)
for (i in 1:nrow(traits_table)) {
  for (group_consolidated in names(groups_consolidated)) {
    if (traits_table$group[i] %in% groups_consolidated[[group_consolidated]]) {
      traits_table$group_consolidated[i] <- group_consolidated
      break
    }
  }
}
rm(i, group_consolidated)

#Remove duplicate column 
traits_table$group <- traits_table$group_consolidated

#Slight edit to traits_table 
traits_table <- traits_table %>% mutate(trait_type = ifelse(as.character(trait_type) == "binary", "Binary", as.character(trait_type))) %>% mutate(trait_type = ifelse(as.character(trait_type) == "quantitative", "Quantitative", as.character(trait_type)))

#Load gvc file 
traits_table$summed_gvc100_raw <- (traits_table$summed_gvc100_raw / traits_table$total_summed_gvc)

x1 <- traits_table$gini_United
y1 <- traits_table$summed_gvc100_raw

binary <- traits_table[(traits_table$trait_type == 'Binary'),]
quant <- traits_table[(traits_table$trait_type == 'Quantitative'),]

bin_x <- binary$gini_United
bin_y <- binary$summed_gvc100_raw

quant_x <- quant$gini_United
quant_y <- quant$summed_gvc100_raw

bin_lm <- lm(bin_y ~ bin_x)
bin_r2 <- summary(bin_lm)$adj.r.squared

quant_lm <- lm(quant_y ~ quant_x)
quant_r2 <- summary(quant_lm)$adj.r.squared


#plot percentage of gvc vs. gini_UK
plot1 <- ggplot(data=traits_table, mapping = aes(x = x1, y = y1, color = trait_type)) + geom_smooth(method = "lm") + labs(color = 'Trait Types') + geom_point(size = 2) + theme_light() + labs(x=expression('Gini'[UK]), y="Proportion of Total gvc") + scale_x_continuous(name=expression('Gini'[UK]), breaks=c(0, 0.25, 0.5, 0.75, 1.0, 1.25)) + scale_y_continuous(name="Proportion of Total gvc", breaks=c(0, 0.25, 0.5, 0.75, 1.0, 1.25)) 
plot1 <-  plot1 + xlim(0, 1.0) + ylim(0,1.0) + annotate("text",  x=0.25, y = 1, label = paste("Binary Adj. R-squared: ", round(bin_r2, 3), sep = "")) + annotate("text",  x=0.30, y = 0.95, label = paste("Quantitative Adj. R-squared: ", round(quant_r2, 3), sep = ""))
print(plot1)

ggsave(file = "FigureS3.pdf", units = c("in"), width=6, height=5, dpi=300, plot1)


