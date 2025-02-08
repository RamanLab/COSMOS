rm(list = ls())
library(readxl)
library(tidyr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(extrafont)
library(showtext)

# Load your data
file_path <- 'C:/Users/lavan/OneDrive/Desktop/CommvsMono/TotalData.xls'
data <- read_excel(file_path)

# Extract relevant columns
data <- data %>%
  select(Environment = 1, Community = 2, Product = 5, Productivity_Ratio = 6)

# Ensure 'Environment' and 'Product' are categorical
data$Environment <- as.factor(data$Environment)
data$Product <- as.factor(data$Product)

# Ensure 'Productivity' is numeric
data$Productivity_Ratio <- as.numeric(data$Productivity_Ratio)

# Group by Environment and Product, then compute the mean Productivity_Ratio
averaged_data <- data %>%
  group_by(Environment, Product) %>%
  summarize(mean_Productivity_Ratio = mean(Productivity_Ratio, na.rm = TRUE)) %>%
  ungroup()

# Rename columns
colnames(averaged_data) <- c('Environment', 'Product', 'Productivity_Ratio')

# Function to calculate the threshold for top 2% values
calculate_threshold <- function(data_column) {
  top5_threshold <- quantile(data_column, 0.98, na.rm = TRUE)
  return(top5_threshold)
}

# Calculate threshold for Productivity_Ratio
threshold <- calculate_threshold(averaged_data$Productivity_Ratio)

# Pivot the data to create a matrix with Communities as rows and Products as columns
heatmap_data <- pivot_wider(averaged_data, names_from = Product, values_from = Productivity_Ratio)

# Convert to a data frame and set row names
heatmap_data <- as.data.frame(heatmap_data)
rownames(heatmap_data) <- heatmap_data$Environment
heatmap_matrix <- as.matrix(heatmap_data[,-1])  # Exclude the Environment column

# Handle missing values by replacing them with 0 (or another suitable value)
heatmap_matrix[is.na(heatmap_matrix)] <- 0

# Mark values above the threshold as 'high' and replace them with a specific value (e.g., 100)
heatmap_matrix[heatmap_matrix > threshold] <- threshold+1

# Transpose the matrix
heatmap_matrix <- t(heatmap_matrix)

# Reorder columns
desired_order <- c("Aer_Rich", "Aer_Min", "Anaer_Rich", "Anaer_Min")
heatmap_matrix <- heatmap_matrix[, desired_order, drop = FALSE]

# Annotation data
annotation <- data.frame(Product_Class = c("Diol", "Organic acid", "Carboxylic acid", "Alcohol", "SCFA", "Phenol", 
                                           "Organic acid", "Alcohol", "Organic acid", "Organic acid", "Organic acid", 
                                           "Tripeptide", "Triol", "Gas", "Organic acid", "Alcohol", "Diol", "Diol", 
                                           "Carboxylic acid", "Polyamine", "Organic acid", "Sugar alcohol", "Polyamine", 
                                           "Organic acid", "Sugar alcohol"))
annotation$Product_Class <- factor(annotation$Product_Class, levels = unique(annotation$Product_Class))
rownames(annotation) <- rownames(heatmap_matrix)  # Adjust to match the rows of the transposed matrix

# Assigning colors to each Product_Class
annotation_colors <- list(Product_Class = c(
  "Organic acid" = "#E69F00",  # Yellow-Orange
  "Carboxylic acid" = "#F0E442",  # Bright Yellow
  "Sugar alcohol" = "#9B59B6",  # Purple
  "Alcohol" = "#0072B2",  # Strong Blue
  "Diol" = "#D55E00",  # Dark Red-Orange
  "Triol" = "#CC79A7",  # Pink
  "Gas" = "#56B4E9",  # Light Blue
  "SCFA" = "#E41A1C",  # Bright Red
  "Polyamine" = "#999999",  # Gray
  "Phenol" = "#FF69B4",  # Hot Pink
  "Tripeptide" = "#8B4513"  # Saddle Brown
))

# Create heatmap annotation
ha <- rowAnnotation(Product_Class = annotation$Product_Class, 
                    col = annotation_colors,
                    show_annotation_name = FALSE,
                    annotation_legend_param = list(title = "Product Class"))

# Define custom color function
col_fun <- function(x) {
  ifelse(x > threshold, "#000080", colorRamp2(c(min(x), 0, threshold), c("red", "white", "blue"))(x))
}

# Rotate column labels for better readability
column_names_rotated <- gpar(fontsize = 10, rot = 45)

# Adjust row label settings
row_names_adjusted <- gpar(fontsize = 10)

# Enhanced Heatmap
Heatmap(
  heatmap_matrix, 
  name = "Productivity_Ratio",
  left_annotation = ha,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = row_names_adjusted,      # Adjust row label formatting
  column_names_gp = column_names_rotated, # Rotate column labels
  column_names_rot = 45, 
  col = col_fun,                          # Use custom color palette
  heatmap_legend_param = list(
    title = "Mean Productivity Ratio",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 10),
    at = c(min(heatmap_matrix), 0, threshold, threshold+1),  # Add custom color for high values
    labels = c(round(min(heatmap_matrix),2), "0", round(threshold,2), "High")
  ),
  column_title = "Environments",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"), # Adjust column title formatting
  row_title = "Products",
  row_title_gp = gpar(fontsize = 14, fontface = "bold"),    # Adjust row title formatting
  height = unit(4 * nrow(heatmap_matrix), "mm"),
  width = unit(6 * ncol(heatmap_matrix), "mm"),
  row_dend_width = unit(2, "cm"),
  column_dend_height = unit(0.5, "cm")
)
