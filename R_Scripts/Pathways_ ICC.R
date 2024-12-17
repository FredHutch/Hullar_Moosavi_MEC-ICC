################################################################################
#ICC Analysis_Pathways Level
#Author: Darya Moosavi
#Date: 3/30/2024
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


#Load data
data <- read.csv("MCT_metagenomics_filtered_pathways_20240325_NEW.csv", header = TRUE)
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
pathways_adjusted <- (Data_filtered + smallest_nonzero)

#Apply the CLR transformation:
data_transformed <- t(apply(pathways_adjusted, 1, clr)) #Number 2 indicates that the function should be applied column-wise. 

##shift data to eliminate negative values
min_value <- min(data_transformed)
if (min_value < 0) {
  data_transformed <- data_transformed - min_value + 1
}


# Write the final data to a CSV file:
write.csv(data_transformed, "Metagenomics_Filtered_Pathways_Transformed.csv", row.names = TRUE)
######################################################################################



##########################DISTANCE MATRIX AND PCOA ANALYSIS######################################
datadis <- read.csv("Metagenomics_Filtered_Pathways_Transformed_NEW.csv", header = TRUE)
data_for_vegdist <- datadis[,-1]
 
 dist_matrix_bray <- vegdist(data_for_vegdist, method = "bray")
 dist_matrix_bray <- as.matrix(dist_matrix_bray)
 print(dist_matrix_bray)
 
#Save
write.csv(dist_matrix_bray, "dist_matrix_Pathway_NEW.csv", row.names = FALSE)
 

#PCA
cmdscale_result <- cmdscale(dist_matrix_bray)
 
 
############ICC Analsis##############
#ICC
icc_results <- lapply(datadis[, 2:ncol(datadis)], function(x) 
performance::icc(lmer(x ~ 1 + (1 | ID), data = datadis)))
 
#Extract ICC values
extracted_iccs <- lapply(icc_results, function(result) {
   if (is.list(result) && !is.null(result$ICC_adjusted) && !is.null(result$ICC_conditional)) {
     return(c(ICC_adjusted = result$ICC_adjusted, ICC_conditional = result$ICC_conditional))
   } else {
     return(c(ICC_adjusted = NA, ICC_conditional = NA))
   }
 })
 
# Convert the list to a data frame
icc_dataframe <- do.call(rbind.data.frame, extracted_iccs)
 
#Set the rownames
rownames(icc_dataframe) <- names(icc_results)
 
#Change the column names
colnames(icc_dataframe) <- c("ICC_adjusted", "ICC_conditional")
 
#Write the ICC results to a CSV file
write.csv(icc_dataframe, "ICC_Results_Pathways.csv", row.names = TRUE)
 
#Mean abundance and ICC visualization
mean_abundance <- as.data.frame(colMeans(datadis[, 2:ncol(datadis)])) 
icc_visualization <- cbind(icc_dataframe, mean_abundance) 
colnames(icc_visualization) <- c("ICC_adjusted", "ICC_conditional", "mean_abundance")
 
#Write the ICC results to a CSV file
write.csv(icc_visualization, "ICC_Results_Pathways_NEW.csv", row.names = TRUE)
 
#################Visualization##########################################
ggplot(icc_visualization, aes_string(x = "mean_abundance", y = "ICC_adjusted")) + 
   geom_point(aes(color = ICC_adjusted)) +
   geom_hline(yintercept = 0.4, linetype = 2) +
   
   ylim(0, 1) +
   
   scale_color_gradient(low = "blue", high = "red") +
   
   labs(
     x = "Mean CLR-transformed abundance",
     y = "ICC",
     color = "ICC Score"
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
#################################################################################
