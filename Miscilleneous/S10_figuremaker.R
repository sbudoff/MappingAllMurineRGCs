# Required libraries
library(tidyverse)
library(patchwork)
library(ggrepel)
library(deldir)

# Configuration
output_dir <- "/home/sam/MappingAllMurineRGCs/Data/Figure4_Outputs"
LABEL_SIZE <- 30  # Hyperparameter for label font size

# Custom color scheme matching the cartoon
cluster_colors <- c(
  "1" = "#A5B6D4",  # Pink for binocular region
  "2" = "#90F7B5",  # Light yellow for peripheral region
  "3" = "#F8691B",  # Light blue for sky region
  "4" = "#1967FB"   # Light green for ground region
)

# Function to create enhanced UMAP visualization with mask regions
# Function to create enhanced UMAP visualization with optional mask regions
create_enhanced_umap_visualization <- function(umap_coords, cluster_assignments,
                                               mask_projections, 
                                               bkgrnd_tessellation = FALSE) {
  # Combine cluster assignments with UMAP coordinates
  plot_data <- merge(umap_coords, cluster_assignments, by = "Label")
  
  # Initialize the plot
  p <- ggplot()
  
  # Add tessellated background if enabled
  if (bkgrnd_tessellation) {
    # Create a fine grid for the background tessellation
    grid_points <- 100
    x_range <- range(plot_data$UMAP1)
    y_range <- range(plot_data$UMAP2)
    margin <- 0.1 * c(diff(x_range), diff(y_range))
    x_seq <- seq(x_range[1] - margin[1], x_range[2] + margin[1], length.out = grid_points)
    y_seq <- seq(y_range[1] - margin[2], y_range[2] + margin[2], length.out = grid_points)
    grid_df <- expand.grid(UMAP1 = x_seq, UMAP2 = y_seq)
    
    # Function to find nearest point for grid cells
    find_nearest_point <- function(x, y, points) {
      dists <- sqrt((points$UMAP1 - x)^2 + (points$UMAP2 - y)^2)
      return(points$Cluster[which.min(dists)])
    }
    
    # Assign cluster to each grid point
    grid_df$Cluster <- vapply(1:nrow(grid_df), function(i) {
      find_nearest_point(grid_df$UMAP1[i], grid_df$UMAP2[i], plot_data)
    }, FUN.VALUE = numeric(1))
    
    # Add tessellation layer
    p <- p + geom_tile(data = grid_df, 
                       aes(x = UMAP1, y = UMAP2, fill = factor(Cluster)),
                       alpha = 0.4)
  }
  
  # Add remaining layers
  p <- p +
    # Add points
    geom_point(data = plot_data,
               aes(x = UMAP1, y = UMAP2, color = factor(Cluster)),
               size = 3) +
    # Add point labels
    geom_text_repel(data = plot_data,
                    aes(x = UMAP1, y = UMAP2, label = Label),
                    size = 3, 
                    max.overlaps = Inf,
                    force = 10) +
    # Add mask points
    geom_point(data = mask_projections,
               aes(x = UMAP1, y = UMAP2),
               size = 5, shape = 17, color = "red", alpha=0.3) +
    # Add mask labels with stronger repulsion
    geom_text_repel(data = mask_projections,
                    aes(x = UMAP1, y = UMAP2, label = Label),
                    size = 4, 
                    color = "red", 
                    fontface = "bold", 
                    alpha = 0.6,
                    force = 20,
                    max.overlaps = Inf,
                    box.padding = 1,
                    point.padding = 0.5,
                    direction = "both") +
    # Customize appearance
    scale_fill_manual(name = "Cluster", values = cluster_colors) +
    scale_color_manual(name = "Cluster", values = cluster_colors) +
    theme_minimal() +
    theme(legend.position = "right",
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0, vjust = -1, size = LABEL_SIZE)) +
    labs(title = "b")
  
  return(p)
}

# Read the required files
silhouette_scores <- read.csv(file.path(output_dir, "silhouette_scores_all_methods.csv"))
umap_coords <- read.csv(file.path(output_dir, "umap_coordinates.csv"))
cluster_assignments <- read.csv(file.path(output_dir, "cluster_assignments.csv"))
mask_projections <- read.csv(file.path(output_dir, "mask_projections.csv"))

# Find consensus k (modal optimal k)
consensus_k <- as.numeric(names(sort(table(silhouette_scores$k[
  silhouette_scores$silhouette_score == ave(silhouette_scores$silhouette_score,
                                            silhouette_scores$method,
                                            FUN = max)
]), decreasing = TRUE)[1]))

# Create silhouette plot
sil_plot <- ggplot(silhouette_scores, 
                   aes(x = k, y = silhouette_score, 
                       color = method, group = method)) +
  geom_line() +
  geom_point(size = 2) +
  geom_vline(xintercept = consensus_k, linetype = "dashed", color = "red", alpha = 0.5) +
  theme_minimal() +
  scale_color_viridis_d() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0, vjust = -1, size = LABEL_SIZE)) +
  labs(title = "a",
       x = "Number of clusters (k)",
       y = "Silhouette Score")

# Create UMAP visualization
consensus_plot <- create_enhanced_umap_visualization(
  umap_coords = umap_coords,
  cluster_assignments = cluster_assignments,
  mask_projections = mask_projections
)

# Create combined plot
(combined_plot <- sil_plot + consensus_plot +
  plot_layout(ncol=2, widths=c(1, 1)))

# Save the plot
ggsave(file.path(output_dir, "S10_clustering_summary.png"),
       combined_plot, width=14.5, height=6, dpi=300)
