# Load all required libraries
library(tidyverse)
library(proxy)
library(clusterSim)
library(ggdendro)
library(patchwork)

root <- '/home/sam/MappingAllMurineRGCs/Data/Figure6_Outputs/'

# First, let's read the raw data
data <- read.delim(paste0(root, 'mouse_subtype_normalized_count_densities.txt'), header = T, sep = "\t")

# Add group identifier
data$group <- 1:nrow(data)

# Convert to long format and extract coordinates
data_long <- data %>%
  pivot_longer(
    cols = -group,
    names_to = "coordinate",
    values_to = "density"
  ) %>%
  # Extract x and y from the coordinate names using the actual format
  mutate(
    x = as.numeric(gsub("x\\.(\\d+)\\.y\\.\\d+", "\\1", coordinate)),
    y = as.numeric(gsub("x\\.\\d+\\.y\\.(\\d+)", "\\1", coordinate))
  )

# Create heatmaps for groups 1-16
# Create heatmaps for groups 1-16 with group-wise normalization
data_long %>%
  filter(group <= 16) %>%
  group_by(group) %>%
  rename(
    normalized_density = density#(density - min(density)) / (max(density) - min(density))
  ) %>%
  ggplot(aes(x = x, y = y, fill = normalized_density)) +
  geom_tile() +
  facet_wrap(~group, nrow = 4) +
  scale_fill_viridis_c(limits = c(0, 1)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 8),
    strip.text = element_text(size = 8),
    panel.spacing = unit(0.5, "lines")
  ) +
  coord_fixed() +
  scale_y_reverse()


################################################################################
library(proxy)

# Reshape the data so each row is a grid position and each column is a group
clustering_matrix_positions <- data_long %>%
  pivot_wider(
    id_cols = c(x, y),
    names_from = group,
    values_from = density
  ) %>%
  mutate(position = paste0("x", x, "_y", y)) %>%
  dplyr::select(-x, -y) %>%
  column_to_rownames("position") %>%
  as.matrix()

cosine_sim <- function(A) {
  # Normalize each row
  A_norm <- t(apply(A, 1, function(x) x/norm_vec(x)))
  # Calculate similarity matrix
  sim <- A_norm %*% t(A_norm)
  return(sim)
}
  
max_k <- 6
results <- list()
for (dist.method in c("euclidean", "manhattan", "canberra")){

  if (dist.method != "cosine") {
    # This is depreciated and not implemented but produces the same result. Removed to avoid hand coded-math where possible
    dist_pos <- dist(clustering_matrix_positions, method = dist.method) 
  } else {
    # Calculate cosine similarity manually
    norm_vec <- function(x) sqrt(sum(x^2))
    # Calculate similarity and convert to distance
    sim_matrix <- cosine_sim(clustering_matrix_positions)
    # Convert similarity to distance (1 - similarity)
    dist_matrix <- 1 - sim_matrix
    # Convert to dist object for hclust
    dist_pos <- as.dist(dist_matrix)
  }
  # Perform hierarchical clustering
  
  for (link.method in c("ward.D", "ward.D2", "complete")) {
    hc_pos <- hclust(dist_pos, method = link.method)
    
    # Calculate silhouette scores and DB index
    silhouette_scores <- numeric()
    db_scores <- numeric()
    k_values <- 2:max_k
    
    for(k in k_values) {
      clusters <- cutree(hc_pos, k = k)
      sil <- silhouette(clusters, dist_pos)
      silhouette_scores[k-1] <- mean(sil[,3])
      db_score <- index.DB(x = clustering_matrix_positions, cl = clusters)
      db_scores[k-1] <- db_score
    }
    
    # Get optimal clusters for both methods
    optimal_k_sil <- k_values[which.max(silhouette_scores)]
    optimal_k_db <- k_values[which.min(db_scores)]  # DB index should be minimized
    
    # Get cluster assignments for both methods
    optimal_clusters_sil <- cutree(hc_pos, k = optimal_k_sil)
    optimal_clusters_db <- cutree(hc_pos, k = optimal_k_db)
    
    # Create visualization dataframes for both clustering results
    clustering_results_sil <- data.frame(
      position = rownames(clustering_matrix_positions),
      cluster = factor(optimal_clusters_sil)
    ) %>%
      separate(position, into = c("x", "y"), sep = "_") %>%
      mutate(
        x = as.numeric(gsub("x", "", x)),
        y = as.numeric(gsub("y", "", y))
      )
    
    clustering_results_db <- data.frame(
      position = rownames(clustering_matrix_positions),
      cluster = factor(optimal_clusters_db)
    ) %>%
      separate(position, into = c("x", "y"), sep = "_") %>%
      mutate(
        x = as.numeric(gsub("x", "", x)),
        y = as.numeric(gsub("y", "", y))
      )
    
    # Create plots
    # 1. Dendrogram
    dendr <- ggdendrogram(hc_pos, rotate = FALSE) +
      theme_minimal() +
      labs(title = "Hierarchical Clustering of Grid Positions") +
      theme(axis.text.x = element_blank())
    
    # 2. Silhouette plot
    silhouette_df <- data.frame(
      k = 2:max_k,  
      silhouette_score = silhouette_scores
    )
    
    sil_plot <- ggplot(silhouette_df, aes(x = k, y = silhouette_score)) +
      geom_line() +
      geom_point() +
      theme_minimal() +
      labs(
        title = "Silhouette Analysis",
        x = "Number of clusters (k)",
        y = "Average Silhouette Score"
      )
    
    # Create DB dataframe correctly
    db_df <- data.frame(
      k = k_values,
      db_score = as.numeric(db_scores)  # Ensure it's a single numeric vector
    )
    
    # Update db_plot
    db_plot <- ggplot(db_df, aes(x = k, y = db_score)) +
      geom_line() +
      geom_point() +
      theme_minimal() +
      labs(
        title = "Davies-Bouldin Index",
        x = "Number of clusters (k)",
        y = "DB Index"
      )
    
    db_plot <- ggplot(db_df, aes(x = k, y = db_score)) +
      geom_line() +
      geom_point() +
      theme_minimal() +
      labs(
        title = "Davies-Bouldin Index",
        x = "Number of clusters (k)",
        y = "DB Index"
      )
    
    # 4. Spatial cluster plots
    spatial_plot_sil <- ggplot(clustering_results_sil, aes(x = x, y = y, fill = cluster)) +
      geom_tile() +
      scale_fill_viridis_d() +
      theme_minimal() +
      coord_fixed() +
      scale_y_reverse() +
      labs(title = "Clusters (Silhouette)",
           fill = "Cluster")
    
    spatial_plot_db <- ggplot(clustering_results_db, aes(x = x, y = y, fill = cluster)) +
      geom_tile() +
      scale_fill_viridis_d() +
      theme_minimal() +
      coord_fixed() +
      scale_y_reverse() +
      labs(title = "Clusters (DB Index)",
           fill = "Cluster")
    
    # Combine plots using patchwork
    combined_plot <- (dendr + sil_plot + db_plot) / 
      (plot_spacer() + spatial_plot_sil + spatial_plot_db) +
      plot_layout(heights = c(1, 1.5),
                  widths = c(1, 1, 1))
    
    # Print the combined plot
    print(combined_plot)
    
    results[[paste0(dist.method,'_',link.method)]] <- list(
      'dist' = dist.method,
      'link' = link.method,
      'plot' = combined_plot,
      'silhouette_df' = silhouette_df,
      'clusters' = clustering_results_sil,
      'db_clusters' = clustering_results_db,
      'db_df' = db_df
    )
  }
}




# Initialize with the first clustering result
first_key <- names(results)[1]
cluster_df <- results[[first_key]]$clusters %>%
  dplyr::select(x, y) %>%
  mutate(!!first_key := results[[first_key]]$clusters$cluster)

# Add all other clustering results
for (method_name in names(results)[-1]) {
  cluster_df <- cluster_df %>%
    left_join(
      results[[method_name]]$clusters %>%
        dplyr::select(x, y, cluster) %>%
        rename(!!method_name := cluster),
      by = c("x", "y")
    )
}

# Function to calculate mode
get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Add consensus column
cluster_df <- cluster_df %>%
  rowwise() %>%
  mutate(
    consensus = get_mode(c_across(matches("_")))
  ) %>%
  ungroup()

# Create combined silhouette plot
silhouette_combined <- map_dfr(
  names(results),
  ~bind_rows(
    results[[.x]]$silhouette_df %>%
      mutate(
        method = .x
      )
  )
)

# Create the visualization
p1 <- ggplot(silhouette_combined, 
         aes(x = k, y = silhouette_score, 
             color = method, group = method)) +
  geom_line() +
  geom_point(size = 2) +
  geom_vline(xintercept = max(as.numeric(cluster_df$consensus)), linetype = "dashed", color = "red", alpha = 0.5) +
  theme_minimal() +
  scale_color_viridis_d() +
  labs(
    title = "a",
    x = "Number of clusters (k)",
    y = "Silhouette Score",
    color = "method"
  ) +
  theme(legend.position = "bottom")

# Create consensus cluster plot
p2 <- ggplot(cluster_df, aes(x = x, y = y, fill = factor(consensus))) +
  geom_tile() +
  scale_fill_viridis_d(name = "Consensus\nCluster") +
  theme_minimal() +
  coord_fixed() +
  theme_void() +
  labs(title = "b") +
  theme(legend.position = "right")

# Combine plots using patchwork
combined_plot <- p1 + p2 +
  plot_layout(widths = c(1, 1))

# Display the plot
print(combined_plot)

write.csv(cluster_df, file=paste0(root,'hc_consensus.csv'))
write.csv(silhouette_combined, file=paste0(root,'hc_silhouettes.csv'))


# Save the plot
ggsave(file.path(root, "S13_clustering_summary.png"),
       combined_plot, width=14.5, height=6, dpi=300)


