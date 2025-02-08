rm(list = ls())
library(readxl)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)

file_path <- 'C:/Users/lavan/OneDrive/Desktop/CommvsMono/TotalData.xls'
data <- read_excel(file_path)
head(data)

# Function to calculate the threshold for top 2% values
calculate_threshold <- function(data_column) {
  top2_threshold <- quantile(data_column, 0.95, na.rm = TRUE)
  return(top2_threshold)
}

# Create a new column with characters before the underscore if it is present
data$InteractionComm <- sub("_.*", "", data$InteractionA)

data <- data%>%
  select(InteractionComm,Product,Productivity_Ratio)

# Count the frequency of each interaction
interaction_freq <- data %>%
  group_by(InteractionComm) %>%
  summarize(Freq = n()) %>%
  ungroup()

# Join frequency data with the main data
data_with_freq <- data %>%
  left_join(interaction_freq, by = "InteractionComm")

# Calculate weighted productivity ratio
data_with_freq <- data_with_freq %>%
  mutate(Weighted_Productivity_Ratio = Productivity_Ratio * Freq)

# Aggregate data to handle duplicates by taking the mean of weighted productivity
aggregated_data <- data_with_freq %>%
  group_by(Product, InteractionComm) %>%
  summarize(AvgProductivity_Ratio = mean(Weighted_Productivity_Ratio, na.rm = TRUE)) %>%
  ungroup()

# Proceed with the rest of your analysis as before
# Calculate threshold for Productivity
threshold <- calculate_threshold(aggregated_data$AvgProductivity_Ratio)

# Ensure the data has the necessary columns: Community, Product, Productivity
# Pivot the data to create a matrix with Communities as rows and Products as columns
heatmap_data <- pivot_wider(aggregated_data, names_from = Product, values_from = AvgProductivity_Ratio)

# Convert to a data frame and set row names
heatmap_data <- as.data.frame(heatmap_data)
rownames(heatmap_data) <- heatmap_data$InteractionComm
heatmap_matrix <- as.matrix(heatmap_data[,-1])  

# Handle missing values by replacing them with 0 (or another suitable value)
heatmap_matrix[is.na(heatmap_matrix)] <- 0

# Mark values above the threshold as 'high' and replace them with a specific value (threshold + 1)
heatmap_matrix[heatmap_matrix > threshold] <- threshold + 1
heatmap_matrix <- t(heatmap_matrix)

# Annotation data
annotation <- data.frame(Product_Class = c("Diol", "Organic acid", "Carboxylic acid", "Alcohol", "SCFA", "Phenol",
                                           "Organic acid", "Alcohol", "Organic acid", "Organic acid", "Organic acid",
                                           "Tripeptide", "Triol", "Gas", "Organic acid", "Alcohol", "Diol", "Diol",
                                           "Carboxylic acid", "Polyamine", "Organic acid", "Sugar alcohol", "Polyamine",
                                           "Organic acid", "Sugar alcohol"))
annotation$Product_Class <- factor(annotation$Product_Class, levels = unique(annotation$Product_Class))
rownames(annotation) <- rownames(heatmap_matrix)

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

# Create row annotation
ha <- rowAnnotation(Product_Class = annotation$Product_Class,
                    col = annotation_colors,
                    show_annotation_name = FALSE,
                    annotation_legend_param = list(title = "Product Class"))

# Define custom color function
col_fun <- colorRamp2(c(min(heatmap_matrix), 0, threshold, threshold + 1), c("red", "white", "blue", "#000080"))

# Generate the heatmap with ComplexHeatmap
Heatmap(
  heatmap_matrix,
  name = "Productivity Ratio",
  left_annotation = ha,  # Use rowAnnotation() instead of HeatmapAnnotation()
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  col = col_fun,
  heatmap_legend_param = list(
    title = "Weighted Productivity Ratio",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 10),
    at = c(min(heatmap_matrix), 0, threshold, threshold + 1),  # Add custom color for high values
    labels = c("Low", "0", "High", "Very High")
  ),
  column_title = "Interactions",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"), # Adjust column title formatting
  row_title = "Products",
  row_title_gp = gpar(fontsize = 14, fontface = "bold"),    # Adjust row title formatting
  height = unit(4 * nrow(heatmap_matrix), "mm"),
  width = unit(4 * ncol(heatmap_matrix), "mm"),
  row_dend_width = unit(2, "cm"),
  column_dend_height = unit(0.5, "cm")
)