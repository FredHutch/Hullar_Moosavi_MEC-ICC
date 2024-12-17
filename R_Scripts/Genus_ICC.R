################################################################################
# ICC Analysis - Genus Level
# Author: Darya Moosavi
# Date: 3-28-2024
################################################################################

# Load required packages
library(vegan)
library(lme4)
library(compositions) 
library(dplyr)
library(performance)
library(ggplot2)
library(ggrepel)
library(plotly)
library(MASS)


#Load data 
data <- read.csv("MCT_metaphlan_metagenomics_L6_NEW.csv", header = TRUE)
metadata <- read.csv("MCT_metadata.csv", header = TRUE) 

################################################################################
# Data Preprocessing
################################################################################
rownames(data) <- data$taxa
genus_data_matrix <- as.matrix(data[, -1])

# Filter taxa columns based on mean values > 20%
filtered_data <-genus_data_matrix[,colMeans(genus_data_matrix > 0) >= .2] 
dim(genus_data_matrix)

# #Add a Pseudocount to the whole data set to get rid of all zeros
smallest_nonzero <- min(filtered_data[filtered_data > 0])
adjusted_data <- (filtered_data + smallest_nonzero)

#Apply the CLR transformation:
data_transformed <- t(apply(adjusted_data, 1, clr)) # Number 1 indicates that the function should be applied row-wise.

# Shift data to eliminate negative values
min_val <- min(data_transformed)
if (min_val < 0) {
  data_transformed <- data_transformed - min_val + 1
}


# Write the final data to a CSV file:
write.csv(data_transformed, "Genus_Transformed_data.csv", row.names = TRUE)
######################################################################################

 
################################################################################
# Distance Matrix and Principal Coordinate Analysis (PCoA)
################################################################################
# Calculate Bray-Curtis distance matrix
datadis <- read.csv("Transformed_data_Genus.csv", header = TRUE)

# Exclude the first column (Obesity_id) for the dissimilarity calculation
data_for_vegdist <- datadis[, -1]

# Calculate Bray-Curtis distance matrix
dist_matrix_bray <- vegdist(data_for_vegdist, method = "bray")
dist_matrix_bray <- as.matrix(dist_matrix_bray)


#Principal Coordinate Analysis (PCoA)
cmdscale_result <- cmdscale(dist_matrix_bray)
pcoa_df <- as.data.frame(cmdscale_result)
colnames(pcoa_df) <- c("PC1", "PC2")  



################################################################################
# Intraclass Correlation Coefficient (ICC) Calculation
################################################################################
#Load the transformed data and metadata
dataICC <- read.csv("Transformed_data_Genus.csv", header = TRUE)

#ICC
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
write.csv(icc_dataframe, "ICC_Results_data_NEW.csv", row.names = TRUE)
#########################################################################




################################################################################
# Visualization of ICC and mean Abundance
################################################################################
mean_abundance <- as.data.frame(colMeans(dataICC[, 2:ncol(dataICC)])) 
icc_visualization <- cbind(icc_dataframe, mean_abundance) 
colnames(icc_visualization) <- c("ICC_adjusted", "ICC_conditional", "mean_abundance")

# Write the ICC results to a CSV file
write.csv(icc_visualization, "Genus-ICC_Results-20%_CLR_data_NEW.csv", row.names = TRUE)


#Plot ICC against mean abundance
ggplot(icc_visualization, aes(x = mean_abundance, y = ICC_adjusted)) + 
  geom_point(aes(color = ICC_adjusted)) +
  geom_hline(yintercept = 0.4, linetype = 2) +
  ylim(0, 1) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    x = "Mean CLR-transformed abundance",
    y = "ICC",
    color = "ICC Scores"
  ) + 
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size = 1, fill = NA),
    axis.line = element_line(colour = "black")
  )
############################################################################

