rm(list = ls())
library(readxl)
library(tidyr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)


# Load your data
file_path <- 'C:/Users/lavan/OneDrive/Desktop/CommvsMono/CSData.xls'
data <- read_excel(file_path)

# Group by the first three columns and compute the mean of duplicates
data <- data %>%
  group_by(Environment, CarbonSource, Product) %>%
  summarise(AvgProductivity_Ratio = mean(Productivity_Ratio), .groups = "drop")

# Reorder rows according to the desired order
data <- data %>%
  mutate(Environment = factor(Environment, levels = c("Aer_Rich", "Aer_Min", "Anaer_Rich", "Anaer_Min"))) %>%
  arrange(Environment)

df <- as.data.frame(data)

# Split the data by environment
env_list <- split(df, df$Environment)

# Function to calculate the threshold for top 2% values
calculate_threshold <- function(data_column) {
  top2_threshold <- quantile(data_column, 0.95, na.rm = TRUE)
  return(top2_threshold)
}


# Function to create heatmap
create_heatmap <- function(data, env) {
  # Pivot the data to create a matrix for the heatmap
  heatmap_data <- acast(data, CarbonSource ~ Product, value.var = "AvgProductivity_Ratio")
  
  # Calculate the threshold for top 2% values
  threshold <- calculate_threshold(as.vector(heatmap_data))
  
  # Set values greater than the threshold to the threshold + 1 (to mark them as 'high')
  heatmap_data[heatmap_data > threshold] <- threshold + 1
  
  # Apply log transformation to the data (adding 1 to avoid log(0))
  heatmap_data <- t((heatmap_data))
  
  # Annotation data
  annotation <- data.frame(Product_Class = c("Diol", "Organic acid", "Carboxylic acid","Alcohol", "SCFA", "Phenol","Organic acid", "Alcohol",
                                             "Organic acid","Organic acid","Organic acid","Tripeptide","Triol","Gas","Organic acid","Alcohol",
                                             "Diol","Diol","Carboxylic acid","Polyamine","Organic acid","Sugar alcohol","Polyamine","Organic acid",
                                             "Sugar alcohol"))
  annotation$Product_Class <- factor(annotation$Product_Class, levels = unique(annotation$Product_Class))
  rownames(annotation) <- rownames(heatmap_data)
  
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
  
  # Create the heatmap with red, white, and blue colors
  heatmap <- Heatmap(heatmap_data,  
                     col = colorRamp2(c(min(heatmap_data, na.rm = TRUE), 0, threshold, threshold + 1), c("red", "white", "blue", "#000080")),
                     left_annotation = ha,
                     column_title = paste(env), # Set environment name as the title of each heatmap
                     column_title_gp = gpar(fontsize = 12), # Adjust column title formatting
                     column_names_gp = gpar(fontsize = 10), # Increase font size for column labels
                     row_names_gp = gpar(fontsize = 10),    # Increase font size for row labels
                     heatmap_legend_param = list(title = paste("MPR:",env), # Legend with environment info
                                                 at = c(min(heatmap_data), 0, threshold, threshold + 1),  # Add custom color for high values
                                                 labels = c(round(min(heatmap_data), 2), "0", round(threshold, 2), "High"),
                                                 labels_gp = gpar(fontsize = 10, rot = 90),
                                                 legend_direction = "horizontal",  # Horizontal legend placement
                                                 legend_width = unit(3, "cm"),     # Adjust legend width
                                                 title_position = "topcenter",     # Centered title
                                                 legend_gp = gpar(fontsize = 10),
                                                 legend_gap = unit(5, "cm")),      # Increase space between legend bars
                     row_title = "Products",
                     row_title_gp = gpar(fontsize = 14, fontface = "bold"),    # Adjust row title formatting
                     height = unit(4 * nrow(heatmap_data), "mm"),
                     width = unit(4 * ncol(heatmap_data), "mm"),
                     row_dend_width = unit(2, "cm"),  
                     column_dend_height = unit(0.5, "cm"))  
  
  return(heatmap)
}

# Create a list of heatmaps with correct environment names
heatmap_list <- lapply(names(env_list), function(env) create_heatmap(env_list[[env]], env))

# Combine heatmaps into a single figure using ComplexHeatmap draw()
ht_list <- Reduce(`+`, heatmap_list)  # Combine heatmaps using the `+` operator

# Draw the heatmap with proper annotation
draw(ht_list, 
     heatmap_legend_side = "right",  # Place heatmap legends on the right
     annotation_legend_side = "right",  # Place Product Class legend on the right
     merge_legend = FALSE,  # Do not merge legends
     padding = unit(c(0, 0, 0, 0), "cm"))  # Increase space between legends

# Add the overall title
grid.text("Carbon sources", x = unit(0.45, "npc"), y = unit(0.85, "npc"), gp = gpar(fontsize = 14, fontface = "bold"))