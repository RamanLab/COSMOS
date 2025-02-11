rm(list = ls())
library(readxl)
library(dplyr)
library(tidyr)
library(dunn.test)
library(mgcv)
library(stats)
library(writexl)

file_path <- 'C:/Users/lavan/OneDrive/Desktop/CommvsMono/TotalData.xls'
data <- read_excel(file_path)

# Ensure 'Environment' and 'Product' are categorical
data$Environment <- as.factor(data$Environment)
data$Product <- as.factor(data$Product)
data$InteractionA <- as.factor(data$InteractionA)
data$InteractionB <- as.factor(data$InteractionB)

# Ensure 'Productivity' is numeric
data$Productivity_Ratio <- as.numeric(data$Productivity_Ratio)

#Perform Shapiro-Wilk test
split_data <- split(data$Productivity_Ratio, data$Environment)
# Perform Shapiro-Wilk test for each environment
shapiro_results <- lapply(split_data, function(x) {shapiro.test(x)})
for (int in names(shapiro_results)) {cat("\nEnvironment:", int, "\n")
  print(shapiro_results[[int]])
}

# Fit a generalized linear model to remove the effect of Product
residual_model <- gam(Productivity_Ratio ~ Product+Community+InteractionA+InteractionB, data = data)

# Extract residuals
data <- data %>%
  mutate(Residuals = residuals(residual_model))

# Perform a Kruskal-Wallis test to check for differences in residuals by Environment
kruskal_result <- kruskal.test(Residuals ~ Environment, data = data)

# Print Kruskal-Wallis test result
print(kruskal_result)

# If significant, perform a Dunn test for pairwise comparisons
dunn_result <- dunn.test(data$Residuals, data$Environment, method = "BH")
  
# Print Dunn test results
print(dunn_result)
  
# Save results to a data frame for further analysis
dunn_df <- as.data.frame(dunn_result)
  
# Specify your desired FDR threshold
FDR_threshold <- 0.05

# Add a Significance column based on the condition
dunn_df$Significance <- dunn_df$P.adjusted < FDR_threshold

# View the updated data frame
head(dunn_df)

# Write the results
write_xlsx(as.data.frame(dunn_df), "Env_PrdtRatio_StatTest.xlsx")
