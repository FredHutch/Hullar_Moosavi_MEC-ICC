################################################################################
#ICC Analysis
#Author: Darya Moosavi
#Date: 3/28/2024
#Repository: https://github.com/FredHutch/Hullar_Moosavi_MEC-ICC.git
################################################################################

#Load required libraries
library(vegan)
library(lme4)
library(compositions)
library(performance)
library(dplyr)
library(ggplot2)
library(readr)


#Load data
data <- read.csv("MCT_metaphlan_metagenomics_L2_NEW.csv", header = TRUE)
metadata <- read.csv("MCT_metadata.csv", header = TRUE) 

################################################################################
# Data Preprocessing
################################################################################
rownames(data) <- data$taxa
data <- as.matrix(data[, -1])

# Filter columns based on mean values greater than 20% 
Data_filtered <-data[,colMeans(data>0)>=.2] 

#Add a Pseudocount to handle zeros
smallest_nonzero <- min(Data_filtered[Data_filtered > 0])
phyla_adjusted <- (Data_filtered + smallest_nonzero)

#Apply CLR transformation
data_transformed <- t(apply(phyla_adjusted, 1, clr)) 

#shift data to eliminate negative values
min_value <- min(data_transformed)
if (min_value < 0) {
  data_transformed <- data_transformed - min_value + 1
}

# Write the final data to a CSV file:
write.csv(data_transformed, "Phylum_transformed_data.csv", row.names = TRUE)




################################################################################
# Distance Matrix and Principal Coordinate Analysis (PCoA)
################################################################################
#ICC
dataICC <- read.csv("Phylum_transformed_data_NEW.csv", header = TRUE)
data_for_vegdist <- dataICC[,-1]

dist_matrix_bray <- vegdist(data_for_vegdist, method = "bray")
dist_matrix_bray <- as.matrix(dist_matrix_bray)
print(dist_matrix_bray)

#PCoA
cmdscale_result <- cmdscale(dist_matrix_bray)
dist_matrix_bray



################################################################################
# Intraclass Correlation Coefficient (ICC)
################################################################################
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
rownames(icc_dataframe) <- names(icc_results)
colnames(icc_dataframe) <- c("ICC_adjusted", "ICC_conditional")

#save ICC results to a CSV file
write.csv(icc_dataframe, "ICC_Results_Phylum.csv", row.names = TRUE)




##############################################################
#Mean Abundance and Visualization
##############################################################
mean_abundance <- as.data.frame(colMeans(dataICC[, 2:ncol(dataICC)])) 
icc_visualization <- cbind(icc_dataframe, mean_abundance) 
colnames(icc_visualization) <- c("ICC_adjusted", "ICC_conditional", "mean_abundance")

# Write the ICC results to a CSV file
write.csv(icc_visualization, "ICC_Results_Phylum.csv", row.names = TRUE)



ggplot(icc_visualization, aes_string(x = "mean_abundance", y = "ICC_adjusted")) + 
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





