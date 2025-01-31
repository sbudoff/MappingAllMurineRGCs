library(tidyverse)
library(uwot)
library(scales)

root <- '/media/sam/Data2/baysor_rbpms_consolidated'
out_dir <- '/home/sam/FinalRGC_xenium/'
# Load RGC data
cat("Loading RGC data...")
rgc_df <- read_csv(paste0(root,'/OldModel/merged_rgc_prediction_expmat.csv'))

# Store metadata columns to remove
metadata_cols <- c("slide", "slice", "dapi_max", "dapi_min", "dapi_mean", 
                   "dapi_sd", "sample", "retina", "Class", "dapi_max_norm",
                   "dapi_min_norm", "dapi_mean_norm", "dapi_sd_norm", 
                   "volume_hull", "nn_dist", "nn_id", "cell", "x", "y", "z",
                   "volume", "x_range", "y_range", "z_range", "rect_vol",
                   "elongation", "flatness", "sphericity", "compactness")

# Select only gene columns and Prediction
gene_data <- rgc_df %>%
  select(-all_of(metadata_cols))

# Separate Prediction column from gene data
prediction_col <- gene_data$Prediction

# Create UMAP embedding
set.seed(18)  # for reproducibility
# 
# find_discriminative_genes <- function(gene_data, min_var_percentile = 0.75, max_pval = 0.01) {
#   # Separate prediction and gene data
#   prediction <- gene_data$Prediction
#   gene_matrix <- gene_data %>% select(-Prediction)
#   
#   # Calculate variance for each gene
#   gene_vars <- apply(gene_matrix, 2, var)
#   var_cutoff <- quantile(gene_vars, min_var_percentile)
#   
#   # First filter: Keep high variance genes
#   high_var_genes <- names(gene_vars)[gene_vars >= var_cutoff]
#   
#   # Second filter: ANOVA for each gene
#   anova_pvals <- sapply(high_var_genes, function(gene) {
#     # Create temporary dataframe for ANOVA
#     temp_df <- data.frame(
#       value = gene_matrix[[gene]],
#       group = prediction
#     )
#     # Perform one-way ANOVA
#     fit <- aov(value ~ group, data = temp_df)
#     # Extract p-value
#     summary(fit)[[1]]$`Pr(>F)`[1]
#   })
#   
#   # Get significant genes
#   significant_genes <- high_var_genes[anova_pvals <= max_pval]
#   
#   # Return filtered gene matrix
#   return(list(
#     genes = significant_genes,
#     filtered_data = gene_data %>% 
#       select(Prediction, all_of(significant_genes))
#   ))
# }
# 
# # Usage:
# filtered_results <- find_discriminative_genes(gene_data, min_var_percentile = 0.5)
# print(paste("Original number of genes:", ncol(gene_data) - 1))
# print(paste("Number of discriminative genes:", length(filtered_results$genes)))

# Use the filtered data for UMAP
gene_matrix <- gene_data %>%#filtered_results$filtered_data %>%
  select(-Prediction) %>%
  as.matrix()

normalize_gene_expression <- function(count_matrix) {
  # Calculate the total counts per cell (row sums)
  total_counts_per_cell <- rowSums(count_matrix)
  
  # Calculate the median total counts across all cells
  median_counts <- median(total_counts_per_cell)
  
  # Create size factors for normalization
  size_factors <- median_counts / total_counts_per_cell
  
  # Create normalized counts matrix by multiplying each row by its size factor
  normalized_counts_matrix <- t(t(count_matrix) * size_factors)
  
  # Log transform: add 1 and take natural log
  log_normalized_matrix <- log1p(normalized_counts_matrix)
  
  return(list(
    normalized_counts = normalized_counts_matrix,
    log_normalized = log_normalized_matrix,
    median_depth = median_counts,
    size_factors = size_factors
  ))
}


# Normalize the gene expression data
normalized_results <- normalize_gene_expression(gene_matrix)
log_normalized_data <- normalized_results$log_normalized

# Scale the normalized data
scaled_data <- log_normalized_data
scaled_data <- scale(log_normalized_data)

######################################################################################
# Take random 20% of rows
subset_idx <- sample(1:nrow(scaled_data), size = round(0.2 * nrow(scaled_data)))
subset_scaled_data <- scaled_data[subset_idx, ]
subset_prediction_col <- prediction_col[subset_idx]
###################################################################################


umap_coords <- umap(subset_scaled_data, #scaled_data, #
                    n_neighbors = 50,     # Reduced from 15
                    min_dist = 0.01,       # Increased from 0.1
                    spread = 5#,         # Added spread parameter
                    #repulsion_strength = 1 # Added repulsion
)

# Create dataframe for plotting
umap_df <- data.frame(
  UMAP1 = umap_coords[,1],
  UMAP2 = umap_coords[,2],
  Prediction = subset_prediction_col #prediction_col# 
)

# umap_df <- read_csv(paste0(out_dir, "xenium_umap_coordinates.csv"))

# Create a custom color palette for 45 classes
n_colors <- length(unique(prediction_col))
colors <- colorRampPalette(
  c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
    "#FFFF33", "#A65628", "#F781BF", "#999999")
)(n_colors)

# Create the plot with modified aesthetics for many classes
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Prediction)) +
  geom_point(size = 1, alpha = 0.3) +
  theme_minimal() +
  scale_color_manual(values = colors) +
  labs(title = "UMAP visualization of gene expression data",
       x = "UMAP1",
       y = "UMAP2") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    # Adjust legend to handle many categories
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.8, "lines")
  ) +
  guides(color = guide_legend(ncol = 2,  # Display legend in 2 columns
                              override.aes = list(size = 3)))  # Larger legend points

# Base layer with all points in black
ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(color = "black", size = 1, alpha = 0.1) +
  # Filtered overlay
  geom_point(data = umap_df[umap_df$Prediction %in% c("21_Tbr1_S2","16_ooDS_DV"), ], # c("21_Tbr1_S2","16_ooDS_DV")  # c("01_W3D1.1","02_W3D1.2", "03_FminiON")
             aes(color = Prediction), size = 1, alpha = 0.5) +
  theme_minimal() +
  # scale_color_manual(values = colors) +
  labs(title = "UMAP visualization of gene expression data",
       x = "UMAP1", y = "UMAP2") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.8, "lines")
  ) +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))


# To save with appropriate dimensions for many classes
ggsave(paste0(out_dir, "xenium_umap_plot.pdf"), width = 12, height = 8, dpi = 300)
write.csv(paste0(out_dir,umap_coords), "xenium_umap_coordinates.csv")

umap_df_subset <- umap_df %>%
 # filter(UMAP1 >= -17.5 & UMAP1 <= -10 & UMAP2 >= -2.5 & UMAP2 <= 8) # for 20% subset
  filter(UMAP1 >= 12.5 & UMAP1 <= 18 & UMAP2 >= -7.5 & UMAP2 <= -1.25)

umap_df_subset %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = Prediction)) +
    geom_point(size = 1, alpha = 0.3) +
    theme_minimal() +
    scale_color_manual(values = colors) +
    labs(title = "UMAP visualization of gene expression data",
         x = "UMAP1",
         y = "UMAP2") +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right",
      # Adjust legend to handle many categories
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 8),
      legend.key.size = unit(0.8, "lines")
    ) +
    guides(color = guide_legend(ncol = 2,  # Display legend in 2 columns
                                override.aes = list(size = 3))) 
summary_stats <- umap_df %>%
  # Get full dataset counts and proportions
  group_by(Prediction) %>%
  mutate(total_type = n()) %>%
  ungroup() %>%
  # Get subset data and calculate all metrics
  left_join(
    umap_df_subset %>%
      group_by(Prediction) %>%
      summarise(
        N = n(),
        .groups = 'drop'
      )
  ) %>%
  # Replace NA with 0 for types not in subset
  mutate(N = replace_na(N, 0)) %>%
  # Calculate proportions
  mutate(
    p_blob = N / nrow(umap_df_subset),
    p_full = N / total_type
  ) %>%
  # Select unique rows with desired columns
  select(Prediction, N, p_blob, p_full) %>%
  distinct() %>%
  # Sort by N descending
  arrange(desc(N))

# View results
ggplot(summary_stats %>% arrange(desc(p_full)), 
       aes(y = reorder(Prediction, p_full), x = p_full)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Proportion in subset relative to full dataset",
       y = "Prediction",
       title = "Distribution of cell types in UMAP subset") +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  )

ggplot(summary_stats %>% arrange(desc(Prediction)), 
       aes(y = reorder(Prediction, desc(Prediction)), x = p_full)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Proportion in subset relative to full dataset",
       y = "Prediction",
       title = "Distribution of cell types in UMAP subset") +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  )


write.csv(umap_df, file = paste0(out_dir, "xenium_umap_coordinates.csv"))
write.csv(umap_df_subset, file = paste0(out_dir, "xenium_umap_coordinates_subset.csv"))
write.csv(summary_stats, file = paste0(out_dir, "xenium_umap_subset_stats.csv"))
