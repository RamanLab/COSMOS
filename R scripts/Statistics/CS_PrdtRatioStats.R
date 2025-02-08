rm(list = ls())
library(readxl)
library(tidyr)
library(dplyr)
library(dunn.test)
library(writexl)
library(mgcv)
library(stats)

# Load your data
file_path <- 'C:/Users/lavan/OneDrive/Desktop/CommvsMono/CSData.xls'
data <- read_excel(file_path)

#Perform Shapiro-Wilk test
split_data <- split(data$Productivity_Ratio, data$CarbonSource)
# Perform Shapiro-Wilk test for each environment
shapiro_results <- lapply(split_data, function(x) {shapiro.test(x)})
for (int in names(shapiro_results)) {cat("\nCarbonSource:", int, "\n")
  print(shapiro_results[[int]])
}

# Fit a linear model to remove the effect of Environment and Product
model <- gam(Productivity_Ratio ~ Product+Environment+Community, data = data)

# Extract the residuals
data$residuals <- residuals(model)

# Perform Kruskal-Wallis H test on the residuals
kruskal_test <- kruskal.test(residuals ~ CarbonSource, data = data)

# Display the results
print(kruskal_test)

# If significant, perform a Dunn test for pairwise comparisons
dunn_result <- dunn.test(data$residuals, data$CarbonSource, method = "BH")

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
write_xlsx(as.data.frame(dunn_df), "CS_PrdtRatio_Stat.xlsx")
