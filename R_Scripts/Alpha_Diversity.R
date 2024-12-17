################################################################################
# Alpha Diversity Analysis at Species Level
# Author: Darya Moosavi
# Date: 12-5-2024
# Repository: https://github.com/FredHutch/Hullar_Moosavi_MEC-ICC.git
################################################################################

# Load required libraries
library(vegan)
library(survival)
library(coin)
library(ggplot2)
library(performance)
library(lme4) 
library(dplyr)
library(RColorBrewer)


#load datasets
data <- read.csv("data/MCT_metaphlan_metagenomics_L7_NEW.csv", header = TRUE, row.names = 1)
metadata <- read.csv("data/MCT_Metadata.csv", header = TRUE)

# Calculate diversity indices
shannon_diversity <- diversity(data, index = "shannon") #Shannon diversity
simpson_index <- diversity(data, index = "simpson") #Simpson diversity
species_richness <- specnumber(data) #species richness

# Combine results into a data frame
diversity_results <- data.frame(
  taxa = rownames(data),
  Shannon_Diversity = shannon_diversity,
  Simpson_Diversity = simpson_index,
  Species_Richness = species_richness
)

write.csv(diversity_results, file = "output/Alpha_diversity.csv", row.names = FALSE)


# Merge diversity results with metadata
merged_data <- merge(diversity_results, metadata, by = "taxa")

# List of diversity indices to check for normality
normality_columns <- c("Shannon_Diversity", "Simpson_Diversity", "Species_Richness")

# Function to check normality and visualize distribution
check_normality <- function(data, var_name) {
  # Generate histogram with density plot
  plot <- ggplot(data, aes_string(var_name)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
    geom_density(color = "red") +
    labs(title = paste("Histogram of", var_name)) +
    theme_minimal()
  
  # Print the plot
  print(plot)
  
  # Perform Shapiro-Wilk Test
  shapiro <- shapiro.test(data[[var_name]])
  cat(sprintf("Shapiro-Wilk Test for %s: W = %.3f, p = %.3f\n", var_name, shapiro$statistic, shapiro$p.value))
}

# Apply the normality check to each variable
lapply(normality_columns, function(var) check_normality(merged_data, var))



# Function to calculate ICC and return results
calculate_icc <- function(data, response, random_effect) {
  model <- lmer(as.formula(paste(response, "~ 1 + (1|", random_effect, ")")), data = data)
  icc_values <- performance::icc(model)
  return(data.frame(
    Variable = response,
    ICC_adjusted = icc_values$ICC_adjusted,
    ICC_unadjusted = icc_values$ICC_unadjusted
  ))
}

# Calculate ICC for multiple variables
icc_results <- do.call(rbind, lapply(
  c("Shannon_Diversity", "Simpson_Diversity", "Species_Richness"),
  function(var) calculate_icc(merged_data, var, "Obesity_id")
))

# Save ICC results
write.csv(icc_results, "output/ICC_results.csv", row.names = FALSE)


# Perform Repeated Measures ANOVA for All Diversity Indices
# List of diversity indices to test
diversity_indices <- c("Shannon_Diversity", "Simpson_Diversity", "Species_Richness")

# Function to perform ANOVA for a given diversity index
run_anova <- function(data, index, subject, time) {
  formula <- as.formula(paste(index, "~", time, "+ Error(", subject, "/", time, ")"))
  anova_result <- aov(formula, data = data)
  cat("\nRepeated Measures ANOVA -", index, ":\n")
  print(summary(anova_result))
}

# List to store p-values
p_values <- list()

# Run ANOVA for each diversity index
for (index in diversity_indices) {
  p_values[[index]] <- run_anova(merged_data, index, "Obesity_id", "TP")
}

p_values



# Visualizations
plot_data <- merged_data %>%
  select(taxa, Shannon_Diversity, TP) %>%
  mutate(Timepoint = as.factor(TP))

p_val <- 0.34

# Boxplot of Shannon Diversity for Each Timepoint
boxplot_data <- ggplot(plot_data, aes(x = Timepoint, y = Shannon_Diversity, group = Timepoint)) +
  geom_boxplot() +
  scale_fill_grey(start = 0.3, end = 0.34) +  # Use grayscale for the fill
  labs(#title = "Shannon diversity at each timepoints",
    x = "Timepoints (months)", 
    y = "Shannon Diversity") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16,  hjust = 0.5),  # Increase plot title size
    axis.title = element_text(size = 14),  # Increase axis titles size
    axis.text = element_text(size = 14),   # Increase axis text size
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 14),   # Increase legend text size
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add a black border around the plot
    plot.background = element_rect(color = "black", fill = NA, size = 1)  # Add a black border around the entire plot area
  ) +
  annotate("text", x = max(as.numeric(plot_data$Timepoint)), y = max(plot_data$Shannon_Diversity), 
           label = paste("p =", p_val), vjust = 1, size = 6, colour = "black")  # Increase annotation text size


# Display the box plot
print(boxplot_data)

