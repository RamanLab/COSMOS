rm(list = ls())
library(readxl)
library(dplyr)
library(tidyr)
library(dunn.test)
library(stats)
library(writexl)
library(stringr)
library(mgcv)

file_path <- 'C:/Users/lavan/OneDrive/Desktop/CommvsMono/TotalData.xls'
data <- read_excel(file_path)
head(data)

# Create a new column with characters before the underscore if it is present
data$InteractionComm <- sub("_.*", "", data$InteractionA)

# Ensure 'Interaction' and 'Product' are categorical
data$InteractionComm <- as.factor(data$InteractionComm)
data$Product <- as.factor(data$Product)
data$Environment <- as.factor(data$Environment)

# Ensure 'Productivity' is numeric
data$Productivity_Ratio <- as.numeric(data$Productivity_Ratio)

# Fit a linear model to remove the effect of Community and Product
residual_model <- gam(Productivity_Ratio ~ Product+Community+Environment, data = data)

# Extract residuals and add them to the dataset
data <- data %>%
  mutate(Residuals = residuals(residual_model))

# Perform a Kruskal-Wallis test to check for differences in residuals by Environment
kruskal_result <- kruskal.test(Residuals ~ InteractionComm, data = data)

# Print Kruskal-Wallis test result
print(kruskal_result)

# Perform the Dunn test
dunn_result <- dunn.test(data$Residuals, data$InteractionComm, method = "BH")
dunn_df <- as.data.frame(dunn_result)

# Specify your desired FDR threshold
FDR_threshold <- 0.05

# Add a Significance column based on the condition
dunn_df$Significance <- dunn_df$P.adjusted < FDR_threshold

# View the updated data frame
head(dunn_df)

# Write the results
write_xlsx(as.data.frame(dunn_df), "Int_PrdtRatio_StatTest.xlsx")