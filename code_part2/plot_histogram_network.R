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
loc_table <- "./Desktop/Gini-PGS/traits_table.csv" #Remove in final version
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
top_bins_loc <- "./Desktop/traits_overlap.csv" #Remove this in final version
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

#Section specific to the histogram 
#Empty vectors for overlap in trait_groups - Maybe just append to generate_network_graph script
all_vector <- c()
lifestyle_vector <- c()
biological_vector <- c()
disease_vector <- c()
physical_vector <- c()

#Iterate though similarity matrix and append to corresponding vectors
for (row in 1:nrow(similarity_matrix)) {
  trait <- rownames(similarity_matrix)[row]
  trait_row <- traits_table[(traits_table$prive_code == trait),]
  trait_group <- trait_row$group
  # print(trait)
  # print(trait_group)
  for (col in 1:ncol(similarity_matrix)) {
    if (row > col ) {
      trait2 <- colnames(similarity_matrix)[col]
      trait_col <- traits_table[(traits_table$prive_code == trait2),]
      trait2_group <- trait_col$group
      overlap <- similarity_matrix[row, col]
      all_vector <- c(all_vector, overlap)
      if (trait_group == "Biological measures" & trait2_group == "Biological measures") {
        biological_vector <- c(biological_vector, overlap)
      } else if (trait_group == "Diseases" & trait2_group == "Diseases") {
        disease_vector <- c(disease_vector, overlap)
      } else if (trait_group == "Physical measures" & trait2_group == "Physical measures") {
        physical_vector <- c(physical_vector, overlap)
      } else if (trait_group == "Lifestyle/psychological" & trait2_group == "Lifestyle/psychological") {
        lifestyle_vector <- c(lifestyle_vector, overlap)
      }
    }
  }
}
rm(col, row, trait, trait_group, trait2, trait2_group, overlap)

#Create variables for the mean and median 
mean_all <- mean(all_vector)
median_all <- median(all_vector)
mean_lifestyle <- mean(lifestyle_vector)
median_lifestyle <- median(lifestyle_vector)
mean_biological <- mean(biological_vector)
median_biological <- median(biological_vector)
mean_disease <- mean(disease_vector)
median_disease <- median(disease_vector)
mean_physical <- mean(physical_vector)
median_physical <- median(physical_vector)

#Generate dataframe from vectors
df <- data.frame(matrix(ncol=2, nrow=0))
makedf <- function(group, vectors) {
  vector_df <- data.frame(matrix(ncol=2, nrow=length(vectors)))
  vector_df[,1] <- group
  vector_df[,2] <- as.numeric(vectors)
  return(vector_df)
}

#Appending to dataframe
df <- rbind(df, makedf('all', all_vector))
df <-rbind(df, makedf('Lifestyle/psychological', lifestyle_vector))
df <-rbind(df, makedf('Biological measures', biological_vector))
df <-rbind(df, makedf('Diseases', disease_vector))
df <-rbind(df, makedf('Physical measures', physical_vector))
colnames(df) <- c("trait_group", "overlap")

#Subset the dataframe to overlap of all traits
df_all <- subset(df, trait_group == 'all')
df_all$trait_group <- ""

#Plot histogram in ggplot
plot1 <- ggplot(data=df_all, aes(x=overlap)) +
  geom_histogram(fill = "#000000", alpha = 1.0, binwidth=1) + labs(x="Number of Overlapping Bins", y="Count") + geom_label(label="Mean = 1.63", x=75, y = 6000, colour="black", label.size= NA, size=7) + geom_label(label="Max = 94", x=75, y = 4700, colour="black", label.size= NA, size=7) + theme(axis.text=element_text(size=20), axis.title=element_text(size=20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=0.25, fill=NA))
plot1 <- plot1 + scale_x_continuous(name="Number of Overlapping Bins", breaks=c(0, 20, 40, 60, 80, 100))

ggsave(file = "Figure1A.pdf", units = c("in"), width=10, height=3.5, dpi=300, plot1)

#Section specific to the network graph 
#Set threshold 
similarity_matrix[similarity_matrix<10] <- 0

#Generate undirected network graph
network <- graph_from_adjacency_matrix(as.matrix(similarity_matrix), weighted=T, mode="undirected", diag=F)

#Color palette for plot 
coul <- c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6")

# Map the color to trait_groups
my_color <- coul[as.numeric(as.factor(traits_table$group))]

#Determine final output location 
pdf(file = "Figure1B.pdf", width = 8, height = 8)
set.seed(5)
#Plot with labels 
#plot(network, vertex.color = my_color, vertex.size=3, vertex.label.color="black", edge.color="black", vertex.label = traits_table$description, vertex.label.cex=0.5, edge.curved=0, edge.width = 2)

#Plot without labels
plot(network, vertex.color = my_color, vertex.size=3, edge.color="black", vertex.label = NA, edge.curved=0, edge.width = 2)

legend(x=-1.35, y=1.0, legend=levels(as.factor(traits_table$group)), fill = coul, border = "black")

dev.off()


