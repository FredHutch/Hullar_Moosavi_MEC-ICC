################################################################################
#ICC Analysis
#Author: Darya Moosavi
#Date: 12-5-2024
#Repository: https://github.com/FredHutch/Hullar_Moosavi_MEC-ICC.git
################################################################################

# Load required libraries
library(vegan)
library(lme4)
library(compositions) # For CLR transformation
library(performance)
library(dplyr)
library(ggplot2)
library(readr)

setwd("/Users/daryamoosavi/Desktop")

#Load Genus data with rows being observations and columns being taxa. 
data <- read.csv("MCT_metaphlan_metagenomics_L7_NEW.csv", header = TRUE)
metadata <- read.csv("MCT_metadata.csv", header = TRUE) 

################################################################################
# Data Preprocessing
################################################################################
rownames(data) <- data$taxa
data <- as.matrix(data[, -1])

# Filter columns based on mean values greater than 20% 
species_filtered <-data[, colMeans(data > 0) >= .2] 

#Add a Pseudocount to handle zeros
smallest_nonzero <- min(species_filtered[species_filtered > 0])
species_adjusted <- (species_filtered + smallest_nonzero)

#Apply CLR transformation
data_transformed <- t(apply(species_adjusted, 1, clr)) 

#shift data to eliminate negative values
min_value <- min(data_transformed)
if (min_value < 0) {
  data_transformed <- data_transformed - min_value + 1
}

# Write the final data to a CSV file:
write.csv(data_transformed, "Species_Transformed_data_NEW.csv", row.names = TRUE)




################################################################################
# Distance Matrix and Principal Coordinate Analysis (PCoA)
################################################################################
datadis <- data_transformed
# Exclude the first column (sample IDs) for the dissimilarity calculation
data_for_vegdist <- datadis[,-1]
dist_matrix <- vegdist(data_for_vegdist, method = "bray")
dist_matrix <- as.matrix(dist_matrix)

# Save the distance matrix
write.csv(dist_matrix, "dist_matrix.csv", row.names = FALSE)


#Perform PCoA
pcoa_result <- cmdscale(dist_matrix)
write.csv(pcoa_result, "pcoa_result.csv", row.names = FALSE)



################################################################################
# Intraclass Correlation Coefficient (ICC)
################################################################################
data_transformed2 <- read.csv("Species_Transformed.csv", header = TRUE)

icc_results <- lapply(data_transformed2[, 2:ncol(data_transformed2)], function(x) 
  performance::icc(lmer(x ~ 1 + (1 | Obesity_id), data = metadata)))


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
write.csv(icc_dataframe, "ICC_Results_Species_NEW.csv", row.names = TRUE)
#########################################################################



##############################################################
#Mean Abundance and Visualization
##############################################################
mean_abundance <- as.data.frame(colMeans(data_transformed2[, 2:ncol(data_transformed2)])) 
icc_visualization <- cbind(icc_dataframe, mean_abundance) 
colnames(icc_visualization) <- c("ICC_adjusted", "ICC_conditional", "mean_abundance")

# Write the ICC results to a CSV file
write.csv(icc_visualization, "ICC_Results_Species_NEW.csv", row.names = TRUE)


#Plot
ggplot(icc_visualization, aes(x = mean_abundance, y = ICC_adjusted)) + 
  # Scatter points with colors
  geom_point(aes(color = ICC_adjusted)) +
  
  # Horizontal line
  geom_hline(yintercept = 0.4, linetype = 2) +
  
  # Limits
  ylim(0, 1) +
  
  # Color scale
  scale_color_gradient(low = "blue", high = "red") +
  
  # Labels
  labs(
    x = "Mean CLR-transformed abundance",
    y = "ICC",
    color = "ICC Scores"
  ) + 
  
  # Theme
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size = 1, fill = NA),
    axis.line = element_line(colour = "black")
  )
