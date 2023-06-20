#Generate overlap histogram 
rm(list = ls())
dev.off()

#Unload all packages 
library(pacman)
p_unload(all)

#Loading libraries
library(tidyverse)
library(rio)
library(data.table)
library(ggplot2)
library(igraph)

# sets working directory
setwd("./")

#Location of the traits table 
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

#Filter tibble to include prive_code, description, trait_type, group
traits_table <- traits_table %>% select(c(prive_code, description, trait_type, group))

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
traits_table <- traits_table[, 1:(ncol(traits_table)-1)]

#Load the table with top bins for each trait 
top_bins_loc <- "./Desktop/Gini-PGS/bin_overlap.csv" #Remove this in final version
top_bins_loc <- "../generated_data/bin_overlap.csv" 
top_bins_table <- import(file=top_bins_loc, header = TRUE)

#Load trait names into a vector
traits <- top_bins_table$prive_code

#Filter down the traits vector using a column of the traits_table
filtered_traits <- traits_table$prive_code
traits <- traits[traits %in% filtered_traits]

#Filter down top_bins_table using traits vector
top_bins_table <- top_bins_table[top_bins_table$prive_code %in% traits,]

#Match order of traits_table with top_bins_traits
traits_table <- traits_table[match(traits, traits_table$prive_code),]

#Prepare top_bins_table for construction of similarity matrix
top_bins_table <- top_bins_table[, 3:ncol(top_bins_table)]

#Remove unnecessary values 
rm(filtered_traits, loc_table, top_bins_loc, remove, groups_consolidated)

#Creation of similarity matrix 
similarity_matrix <- data.frame(matrix(ncol=177, nrow=177))

#Fill similarity matrix
for (row in 1:nrow(top_bins_table)) {
  vector <- c(unlist(top_bins_table[row,]))
  vector <- vector[!is.na(vector)]
  for (row2 in 1:nrow(top_bins_table)) {
    vector2 <- c(unlist(top_bins_table[row2,]))
    vector2 <- vector2[!is.na(vector2)]
    sim <- length(intersect(vector2, vector))
    similarity_matrix[row, row2] <- sim
  }
}
rm(row, row2, sim, vector, vector2)

#Set the row and column names to the traits vector
row.names(similarity_matrix) <- traits
colnames(similarity_matrix) <- traits

#Keep variable with original values for matrix prior to applying threshold
sim <- similarity_matrix

plot_network <- function(threshold) {
  similarity_matrix[similarity_matrix<threshold] <- 0
  sim_threshold <- similarity_matrix
  similarity_matrix <- sim
  
  #Generate network
  network_threshold <- graph_from_adjacency_matrix(as.matrix(sim_threshold), weighted=T, mode="undirected", diag=F)
  
  #Color palette for plot 
  coul <- c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6")
  
  # Map the color to cylinders
  my_color <- coul[as.numeric(as.factor(traits_table$group))]
  
  set.seed(1)
  title <- paste("Overlap: >= ", threshold, sep = "")
  plot(network_threshold, vertex.color = my_color, vertex.size=3, edge.color="black", vertex.label = NA, edge.curved=0, edge.width = 2, main=paste(title, " bins", sep = ""))
  
}


#Determine final output location - Generate Sup. Figure 2
pdf(file = "FigureS2.pdf", width = 8, height = 8)
par(mfrow = c(2, 2))

#Loop through thresholds - Can be changed to desired thresholds
threshold_vector <- c(5, 10, 15, 20)
for (num in threshold_vector) {
  plot_network(num)
}

dev.off()


