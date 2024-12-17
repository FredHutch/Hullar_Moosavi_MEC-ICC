################################################################################
#ICC Analysis_Functional_Genes
#Author: Darya Moosavi
#Date: 7/10/2023
#Repository: https://github.com/FredHutch/Hullar_Moosavi_MEC-ICC.git
################################################################################

# Loading required packages
library(vegan)
library(lme4)
library(sjstats)
library(compositions) 
library(lme4)
library(dplyr)
library(performance)
library(ggplot2)
library(lattice)

setwd("/Users/daryamoosavi/Desktop")

#Load Genus data with rows being observations and columns being taxa. 
data <- read.csv("MCT_metagenomics_filtered_gene_in_pathway_20240617.csv", header = TRUE)
metadata <- read.csv("MCT_metadata.csv", header = TRUE) 

######################PREPROCESSING################################
rownames(data) <- data$ID
data <- as.matrix(data[, -1])

# Filter columns based on mean values greater than 20% 
Data_filtered <-data[,colMeans(data>0)>=.2] 

#Add a Pseudocount to handle zeros
smallest_nonzero <- min(Data_filtered[Data_filtered > 0])
data_adjusted <- (Data_filtered + smallest_nonzero)

#Apply the CLR transformation:
data_transformed <- t(apply(data_adjusted, 1, clr)) #Number 2 indicates that the function should be applied column-wise. 

#shift data to eliminate negative values
min_value <- min(data_transformed)
if (min_value < 0) {
  data_transformed <- data_transformed - min_value + 1
}


# Write the final data to a CSV file:
write.csv(data_transformed, "Transformed_genes_in_pathways.csv", row.names = TRUE)
######################################################################################


##########################DISTANCE MATRIX AND PCOA ANALYSIS######################################
######1. Calculate distance matrix 
#load the data with IDs that count for the repeated measures
datadis <- read.csv("Transformed_genes_in_pathways.csv", header = TRUE)
# Exclude the first column (sample IDs) for the dissimilarity calculation
data_for_vegdist <- datadis[,-1]

dist_matrix_bray <- vegdist(data_for_vegdist, method = "bray")
dist_matrix_bray <- as.matrix(dist_matrix_bray)
print(dist_matrix_bray)

#Principal Coordinate Analysis (PCoA)
cmdscale_result <- cmdscale(dist_matrix_bray)


############ICC Analysis##############
dataICC <- read.csv("Transformed_genes_in_pathways.csv", header = TRUE)
icc_results <- lapply(dataICC[, 2:ncol(dataICC)], function(x) 
  performance::icc(lmer(x ~ 1 + (1 | Obesity_id), data = dataICC)))


# Extract ICC values
extracted_iccs <- lapply(icc_results, function(result) {
  if (is.list(result) && !is.null(result$ICC_adjusted) && !is.null(result$ICC_conditional)) {
    return(c(ICC_adjusted = result$ICC_adjusted, ICC_conditional = result$ICC_conditional))
  } else {
    return(c(ICC_adjusted = NA, ICC_conditional = NA))
  }
})

# Convert the list to a data frame
icc_dataframe <- do.call(rbind.data.frame, extracted_iccs)

# Set the rownames
rownames(icc_dataframe) <- names(icc_results)

# Change the column names
colnames(icc_dataframe) <- c("ICC_adjusted", "ICC_conditional")

# Write the ICC results to a CSV file
write.csv(icc_dataframe, "ICC_Results_Genes_in_pathways.csv", row.names = TRUE)

# Mean abundance and ICC visualization
mean_abundance <- as.data.frame(colMeans(dataICC[, 2:ncol(dataICC)])) #computes the mean of each column from the 2nd column
icc_visualization <- cbind(icc_dataframe, mean_abundance) #bind
# Change the column names
colnames(icc_visualization) <- c("ICC_adjusted", "ICC_conditional", "mean_abundance")

# Write the ICC results to a CSV file
write.csv(icc_visualization, "ICC_Results_Genes_in_pathways.csv", row.names = TRUE)



#################Visualization##########################################
# Scatterplot 
ggplot(icc_visualization, aes(x = `mean_abundance`, y = ICC_adjusted)) +
  geom_point() +
  labs(x = "Mean CLR-transformed abundance", y = "ICC") +
  geom_hline(yintercept = 0.4, linetype = 2) +
  ylim(0, 1) +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size = 1, fill = NA), axis.line = element_line(colour = "black"))

#############################################
ggplot(icc_visualization, aes_string(x = "mean_abundance", y = "ICC_adjusted")) + 
  geom_point(aes(color = ICC_adjusted)) +
  geom_hline(yintercept = 0.4, linetype = 2) +
  ylim(0, 1) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "Mean CLR-transformed abundance", y = "ICC", color = "ICC Score") + 
  theme_minimal() +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size = 1, fill = NA),
        axis.line = element_line(colour = "black"))

#################################################################################



