#plot network graphs from similarity matrix for 211 traits 

#Makes figure 1b
rm(list = ls())

#Unload all packages prior to starting 
library(pacman)
p_unload(all)

#Loading libraries
library(igraph)
library(tidyverse)
library(rio)
library(data.table)
library(stringr)

# sets working directory
setwd("./")

#Location of the traits table 
loc_table <- "../Desktop/Gini-PGS/traits_table.csv" #Remove in final version
loc_table <- "../generated_data/traits_table.txt" 
# reads traits table
traits_table <- as.data.frame(fread(loc_table))

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
top_bins_loc <- "../Desktop/traits_overlap.csv" #Remove this in final version
top_bins_loc <- "../generated_data/traits_overlap.csv" 
top_bins_table <- import(file=top_bins_loc, header = TRUE)

#Load trait names into a vector
traits <- top_bins_table$V1

#Remove  ".0" from  numbers in the traits vector - format correction
traits_vector <- c()
for (item in traits) {
  end_char <- substr(item, nchar(item)-1, nchar(item))
  if (end_char == '.0') {
    new_item <- substr(item, 1, nchar(item)-2)
    traits_vector <- c(traits_vector, new_item)
  } else {
    traits_vector <- c(traits_vector, item)
  }
}
traits <- traits_vector
rm(traits_vector, item, new_item, end_char)

#Set row names using traits vector
row.names(top_bins_table) <- traits

#Prepare top_bins_table for construction of similarity matrix
top_bins_table <- top_bins_table[, 2:ncol(top_bins_table)]

#Creation of similarity matrix 
similarity_matrix <- data.frame(matrix(ncol=211, nrow=211))

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

#Change row names of traits_table, match order of traits vector and convert to description
row.names(traits_table) <- traits_table$prive_code
traits_table <- traits_table[match(traits,rownames(traits_table)),]

#Set the row and column names to the traits vector
row.names(similarity_matrix) <- traits
colnames(similarity_matrix) <- traits

#Optional - Write to CSV
#write.csv(similarity_matrix, file = "~/Desktop/similarity_matrix.csv", row.names = TRUE)

#Set threshold 
similarity_matrix[similarity_matrix<10] <- 0

#Generate undirected network graph
network <- graph_from_adjacency_matrix(as.matrix(similarity_matrix), weighted=T, mode="undirected", diag=F)

#Color palette for plot 
coul <- c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6")

# Map the color to cylinders
my_color <- coul[as.numeric(as.factor(traits_table$group))]

#Determine final output location 
pdf(file = "~/Desktop/Figure1.pdf", width = 8, height = 8)
set.seed(5)
#Plot with labels 
#plot(network, vertex.color = my_color, vertex.size=3, vertex.label.color="black", edge.color="black", vertex.label = traits_table$description, vertex.label.cex=0.5, edge.curved=0, edge.width = 2)

#Plot without labels
plot(network, vertex.color = my_color, vertex.size=3, edge.color="black", vertex.label = NA, edge.curved=0, edge.width = 2)

legend(x=-1.35, y=1.0, legend=levels(as.factor(traits_table$group)), fill = coul, border = "black")

dev.off()

