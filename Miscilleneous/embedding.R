library(progress)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(umap)

library(cluster)      # For k-medoids (PAM)
library(dbscan)       # For DBSCAN and HDBSCAN
library(e1071)        # For SVM
library(fastcluster)  # For efficient hierarchical clustering

library(spatstat)  # For Ripley's K and other spatial statistics


root <- '/media/sam/Data2/baysor_rbpms_consolidated'
rgc_path <- paste0(root,'/OldModel/merged_rgc_prediction_expmat.csv')
out_path <- '/home/sam/FinalRGC_xenium/'
alignment_path <- paste0(out_path, "Alon_registration.csv")

rgc_df <- read_csv(rgc_path)
alignment <- read_csv(alignment_path, col_names = c('Prediction', 'X', 'Y'))%>% #, 'row_0')) %>%
  mutate('row_0' = row_number()) %>%
  select(-X4) %>%
  drop_na()


################################################################################
############################## FUNCTIONS #######################################
################################################################################
comprehensive_embedding <- function(data, group_var, x_col = "X", y_col = "Y", 
                                    n_bins = 100, 
                                    k_min_dist=0, k_max_dist=0.5,
                                    method = "composite") {
  
  groups <- unique(data[[group_var]])
  n_groups <- length(groups)
  
  # Initialize lists for different statistics
  kde_features <- matrix(NA, nrow = n_groups, ncol = n_bins * n_bins)
  ripley_features <- matrix(NA, nrow = n_groups, ncol = 50)
  quadrant_features <- matrix(NA, nrow = n_groups, ncol = 4)
  peak_features <- matrix(NA, nrow = n_groups, ncol = 4)
  
  rownames(kde_features) <- groups
  rownames(ripley_features) <- groups
  rownames(quadrant_features) <- groups
  rownames(peak_features) <- groups
  
  for (i in seq_along(groups)) {
    group <- groups[i]
    subset_data <- data %>%
      dplyr::filter(!!sym(group_var) == group)
    
    # 1. KDE features
    kde_result <- MASS::kde2d(
      subset_data[[x_col]], 
      subset_data[[y_col]], 
      n = sqrt(n_bins)
    )
    kde_vec <- as.vector(kde_result$z)
    kde_features[i,] <- kde_vec / sum(kde_vec)  # Normalize to sum to 1
    
    # 2. Ripley's K
    pp <- ppp(subset_data[[x_col]], subset_data[[y_col]], 
              window = owin(c(-1.5, 1.5), c(-1.5, 1.5)))
    k_func <- Kest(pp, r = seq(k_min_dist, k_max_dist, length.out = 50), correction = "isotropic")
    k_values <- sqrt(k_func$iso/pi) - k_func$r
    ripley_features[i,] <- replace(k_values, is.infinite(k_values), NA)
    
    # 3. Quadrant analysis
    quadrant_counts <- c(
      sum(subset_data[[x_col]] >= 0 & subset_data[[y_col]] >= 0),
      sum(subset_data[[x_col]] < 0 & subset_data[[y_col]] >= 0),
      sum(subset_data[[x_col]] < 0 & subset_data[[y_col]] < 0),
      sum(subset_data[[x_col]] >= 0 & subset_data[[y_col]] < 0)
    )
    quadrant_features[i,] <- quadrant_counts / sum(quadrant_counts)
    
    # 4. Peak analysis from KDE
    z_values <- kde_result$z
    peak_features[i,] <- c(
      max(z_values) / sum(z_values),  # Normalized maximum intensity
      mean(z_values) / sum(z_values),  # Normalized average intensity
      sd(as.vector(z_values)) / sum(z_values),  # Normalized spread
      sum(z_values > quantile(z_values, 0.75)) / length(z_values)  # Proportion high intensity (already normalized)
    )
  }
  
  # Safe scaling function
  safe_scale <- function(m) {
    # Replace any remaining Inf or NaN with NA
    m[!is.finite(m)] <- NA
    
    # For any columns that are all NA, replace with zeros
    na_cols <- colSums(is.na(m)) == nrow(m)
    m[, na_cols] <- 0
    
    # For remaining columns with some NAs, replace NAs with column mean
    m <- apply(m, 2, function(x) {
      if(any(is.na(x))) x[is.na(x)] <- mean(x, na.rm = TRUE)
      return(x)
    })
    
    # Scale
    return(scale(m, center = TRUE, scale = TRUE))
  }
  
  # Combine features with safe scaling
  z_matrix <- cbind(
    safe_scale(kde_features),
    safe_scale(ripley_features),
    safe_scale(quadrant_features),
    safe_scale(peak_features)
  )
  
  rownames(z_matrix) <- groups
  
  # Final check for any remaining NaN/Inf
  z_matrix[!is.finite(z_matrix)] <- 0
  
  return(list(
    density_list = list(),
    z_matrix = z_matrix,
    x_grid = NULL,
    y_grid = NULL,
    method = "composite"
  ))
}
create_density_embedding <- function(data, group_var, x_col = "X", y_col = "Y",
                                     n_bins = 100, method = "kde2d",
                                     k_max_dist = NULL,  # For Ripley's K
                                     n_neighbors = 5) {   # For NND
  # Input validation
  if (!group_var %in% colnames(data)) {
    stop(sprintf("Group variable '%s' not found in data", group_var))
  }
  if (!x_col %in% colnames(data) || !y_col %in% colnames(data)) {
    stop("X or Y columns not found in data")
  }
  if (!method %in% c("kde2d", "density", "ripley", "quadrat", "nnd")) {
    stop("Method must be one of: 'kde2d', 'density', 'ripley', 'quadrat', 'nnd'")
  }

  # Initialize lists
  density_list <- list()

  # Get study area bounds
  x_range <- range(data[[x_col]])
  y_range <- range(data[[y_col]])
  window_area <- diff(x_range) * diff(y_range)

  if (method %in% c("kde2d", "density")) {
    # Original code for kde2d and density methods remains the same
    if (method == "kde2d") {
      # KDE2D approach
      for (group in unique(data[[group_var]])) {
        subset_data <- data %>%
          dplyr::filter(!!sym(group_var) == group)

        density_list[[as.character(group)]] <- MASS::kde2d(
          subset_data[[x_col]],
          subset_data[[y_col]],
          n = n_bins
        )
      }

      # Extract common grid points
      x_grid <- density_list[[1]]$x
      y_grid <- density_list[[1]]$y

      # Create matrix to store z-values
      z_matrix <- matrix(
        NA,
        nrow = length(density_list),
        ncol = length(x_grid) * length(y_grid),
        dimnames = list(names(density_list), NULL)
      )

      # Populate the matrix
      for (i in seq_along(density_list)) {
        z_matrix[i, ] <- as.vector(density_list[[i]]$z)
      }

    } else {
      # Regular density-based approach
      # Calculate range for consistent binning across groups
      x_range <- range(data[[x_col]])
      y_range <- range(data[[y_col]])

      # Create common grid
      x_grid <- seq(x_range[1], x_range[2], length.out = n_bins)
      y_grid <- seq(y_range[1], y_range[2], length.out = n_bins)

      # Initialize z_matrix
      z_matrix <- matrix(
        NA,
        nrow = length(unique(data[[group_var]])),
        ncol = n_bins * n_bins,
        dimnames = list(unique(data[[group_var]]), NULL)
      )

      # Calculate density for each group
      for (group in unique(data[[group_var]])) {
        subset_data <- data %>%
          dplyr::filter(!!sym(group_var) == group)

        # Create 2D histogram
        hist_2d <- hist2d(
          x = subset_data[[x_col]],
          y = subset_data[[y_col]],
          nbins = n_bins,
          return.breaks = TRUE
        )

        # Store density information
        density_list[[as.character(group)]] <- list(
          x = hist_2d$x,
          y = hist_2d$y,
          z = hist_2d$counts / sum(hist_2d$counts)  # Normalize to get density
        )

        # Populate z_matrix
        z_matrix[as.character(group), ] <- as.vector(density_list[[as.character(group)]]$z)
      }
    }

  } else if (method == "ripley") {
    # Set default k_max_dist
    if (is.null(k_max_dist)) {
      k_max_dist <- 0.5  # Conservative default
    }

    # Create r sequence
    r_values <- seq(0, k_max_dist, length.out = n_bins)

    groups <- unique(data[[group_var]])
    z_matrix <- matrix(
      NA,
      nrow = length(groups),
      ncol = length(r_values),
      dimnames = list(groups, NULL)
    )

    for (group in groups) {
      subset_data <- data %>%
        dplyr::filter(!!sym(group_var) == group)

      pp <- ppp(subset_data[[x_col]], subset_data[[y_col]],
                window = owin(c(-1, 1), c(-1, 1)))

      k_func <- Kest(pp, r = r_values, correction = "isotropic")
      l_func <- sqrt(k_func$iso/pi) - r_values

      z_matrix[group, ] <- l_func
    }

    x_grid <- r_values
    y_grid <- NULL

  } else if (method == "quadrat") {
    # Calculate optimal number of quadrats based on number of points
    n_quadrats <- ceiling(sqrt(nrow(data)/20))  # Rule of thumb

    # Create quadrat counts matrix
    z_matrix <- matrix(
      NA,
      nrow = length(unique(data[[group_var]])),
      ncol = n_quadrats * n_quadrats,
      dimnames = list(unique(data[[group_var]]), NULL)
    )

    for (group in unique(data[[group_var]])) {
      subset_data <- data %>%
        dplyr::filter(!!sym(group_var) == group)

      # Create ppp object
      pp <- ppp(subset_data[[x_col]], subset_data[[y_col]],
                window = owin(x_range, y_range))

      # Calculate quadrat counts
      q_counts <- quadratcount(pp, nx = n_quadrats, ny = n_quadrats)

      # Calculate standardized residuals
      expected <- sum(q_counts) / length(q_counts)
      std_resids <- (as.vector(q_counts) - expected) / sqrt(expected)

      # Store in matrix
      z_matrix[as.character(group), ] <- std_resids
    }

    x_grid <- seq(x_range[1], x_range[2], length.out = n_quadrats)
    y_grid <- seq(y_range[1], y_range[2], length.out = n_quadrats)

  } else if (method == "nnd") {
    # Initialize matrix for nearest neighbor distances
    z_matrix <- matrix(
      NA,
      nrow = length(unique(data[[group_var]])),
      ncol = n_bins,
      dimnames = list(unique(data[[group_var]]), NULL)
    )

    for (group in unique(data[[group_var]])) {
      subset_data <- data %>%
        dplyr::filter(!!sym(group_var) == group)

      # Create ppp object
      pp <- ppp(subset_data[[x_col]], subset_data[[y_col]],
                window = owin(x_range, y_range))

      # Calculate nearest neighbor distances
      nnd <- nndist(pp, k = 1:n_neighbors)

      # Create distribution of NND
      nnd_dist <- density(as.vector(nnd), n = n_bins,
                          from = 0, to = max(nnd))

      # Store in matrix
      z_matrix[as.character(group), ] <- nnd_dist$y
    }

    x_grid <- seq(0, max(nnd), length.out = n_bins)
    y_grid <- NULL
  }

  # Return results
  list(
    density_list = density_list,
    z_matrix = z_matrix,
    x_grid = x_grid,
    y_grid = y_grid,
    method = method
  )
}



Zmat2umap <- function(z_matrix, visualize = TRUE) {
  # Perform UMAP embedding
  umap_result <- umap(z_matrix)
  
  # Create a data frame for plotting
  umap_df <- data.frame(
    UMAP1 = umap_result$layout[,1],
    UMAP2 = umap_result$layout[,2],
    Subtype = rownames(z_matrix)
  )
  
  if (visualize){
    p<-ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Subtype, label = Subtype)) +
      geom_point(size = 3) +
      geom_text(vjust = 1.5, hjust = 1) +
      theme_minimal() +
      labs(
        title = "UMAP Embedding of Density Estimates by Subtype",
        x = "UMAP 1",
        y = "UMAP 2"
      ) +
      theme(legend.position = "none")
    print(p)
  }
  return(umap_df)
}


plot_cluster_densities <- function(alignment_data, 
                                   normalize = FALSE, n_bins = 100) {
  # Calculate cluster proportions
  proportions <- alignment_data %>%
    group_by(Cluster) %>%
    summarize(p = n()/nrow(alignment_data))
  
  # Get density embeddings
  embeddings <- create_density_embedding(
    data = alignment_data,
    group_var = "Cluster",
    n_bins = n_bins
  )
  
  # Normalize and scale densities by cluster proportion
  density_list <- list()
  for (clust in names(embeddings$density_list)) {
    # Get original density values
    kde_result <- embeddings$density_list[[clust]]
    
    if (normalize){
      # Normalize z-values to max of 1
      normalized_z <- kde_result$z / max(kde_result$z)
      
      # Scale by cluster proportion
      cluster_prop <- proportions$p[proportions$Cluster == clust]
      scaled_z <- normalized_z * cluster_prop
      
      # Store in density list
      density_list[[clust]] <- list(
        x = kde_result$x,
        y = kde_result$y,
        z = scaled_z
      )
    } else {
      # Store in density list
      density_list[[clust]] <- list(
        x = kde_result$x,
        y = kde_result$y,
        z = kde_result$z
      )
    }
  }
  
  # Prepare data for plotting
  plot_data <- map_dfr(names(density_list), function(clust) {
    data.frame(
      x = density_list[[clust]]$x,
      y = rep(density_list[[clust]]$y, each = length(density_list[[clust]]$x)),
      z = as.vector(density_list[[clust]]$z),
      Cluster = as.factor(clust)
    )
  })
  
  # Create plot
  global_max <- max(plot_data$z)
  
  ggplot(plot_data, aes(x = x, y = y, fill = z)) +
    geom_raster() +
    facet_wrap(~Cluster, ncol = 2) +
    scale_fill_viridis_c(limits = c(0, global_max)) +
    theme_minimal() +
    labs(
      title = "Density Plots by Cluster",
      x = "X",
      y = "Y"
    ) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "lightgray"),
      strip.text = element_text(color = "black")
    )
}



find_optimal_clusters <- function(data, method, max_clusters, seed = 18) {
  wss <- numeric(max_clusters)
  
  if (method %in% c("kmeans", "kmedoids", "hierarchical")) {
    set.seed(seed)
    for (i in 2:max_clusters) {
      if (method == "kmeans") {
        temp_clust <- kmeans(data, centers = i)
        wss[i] <- temp_clust$tot.withinss
      } else if (method == "kmedoids") {
        temp_clust <- pam(data, k = i)
        wss[i] <- temp_clust$objective
      } else if (method == "hierarchical") {
        # Calculate distance matrix once
        dist_mat <- dist(data)
        hc <- hclust(dist_mat, method = "ward.D2")
        clusters <- cutree(hc, k = i)
        # Calculate within-cluster sum of squares
        wss[i] <- sum(sapply(1:i, function(cl) {
          cluster_points <- data[clusters == cl, , drop = FALSE]
          if (nrow(cluster_points) > 1) {
            sum((cluster_points - colMeans(cluster_points))^2)
          } else {
            0
          }
        }))
      }
    }
    
    elbow_df <- data.frame(
      Clusters = 2:max_clusters,
      WSS = wss[2:max_clusters]
    )
    optimal_k <- find_elbow(elbow_df$Clusters, elbow_df$WSS)
    
  } else {
    # For density-based methods, use different optimization criteria
    if (method == "dbscan") {
      # Use kNN distance plot to find optimal eps
      kNNdist <- dbscan::kNNdist(data, k = 4)
      optimal_k <- which.max(diff(sort(kNNdist))) # Find elbow in kNN distances
      optimal_k <- max(optimal_k, 2) # Ensure at least 2 clusters
    } else if (method == "hdbscan") {
      # HDBSCAN automatically determines the number of clusters
      optimal_k <- NULL
    } else if (method == "svm") {
      # For SVM, use silhouette score to determine optimal number of clusters
      silhouette_scores <- numeric(max_clusters - 1)
      for (i in 2:max_clusters) {
        svm_model <- svm(data, y = NULL, type = "one-classification", kernel = "radial", nu = 0.1)
        pred <- predict(svm_model, data)
        silhouette_scores[i-1] <- silhouette(as.numeric(pred), dist(data))
      }
      optimal_k <- which.max(silhouette_scores) + 1
    }
  }
  
  return(optimal_k)
}

find_optimal_clusters_multi <- function(data, method, max_clusters, seed = 18) {
  # Initialize metrics
  metrics <- data.frame(
    k = 2:max_clusters,
    wss = NA,
    silhouette = NA,
    calinski = NA
  )
  
  # Calculate distance matrix once for reuse
  dist_mat <- dist(data)
  
  for (k in 2:max_clusters) {
    # Perform clustering based on method
    if (method == "kmeans") {
      set.seed(seed)
      clust <- kmeans(data, centers = k)
      clusters <- clust$cluster
      centers <- clust$centers
      metrics$wss[k-1] <- clust$tot.withinss
      
    } else if (method == "kmedoids") {
      set.seed(seed)
      clust <- pam(data, k = k)
      clusters <- clust$clustering
      centers <- clust$medoids
      metrics$wss[k-1] <- clust$objective
      
    } else if (method == "hierarchical") {
      hc <- hclust(dist_mat, method = "ward.D2")
      clusters <- cutree(hc, k = k)
      # Calculate centers
      centers <- t(sapply(1:k, function(i) {
        colMeans(data[clusters == i, , drop = FALSE])
      }))
      metrics$wss[k-1] <- sum(sapply(1:k, function(i) {
        sum((data[clusters == i, ] - centers[i, ])^2)
      }))
    }
    
    # Calculate Silhouette score
    sil <- silhouette(clusters, dist_mat)
    metrics$silhouette[k-1] <- mean(sil[, 3])
    
    # Calculate Calinski-Harabasz Index
    between_ss <- sum(sapply(1:k, function(i) {
      ni <- sum(clusters == i)
      if (ni > 0) {
        ni * sum((centers[i,] - colMeans(data))^2)
      } else {
        0
      }
    }))
    within_ss <- metrics$wss[k-1]
    n <- nrow(data)
    metrics$calinski[k-1] <- (between_ss/(k-1)) / (within_ss/(n-k))
  }
  
  # Normalize metrics to 0-1 scale
  metrics_norm <- metrics
  metrics_norm$wss <- 1 - scale(metrics$wss, center = min(metrics$wss), 
                                scale = diff(range(metrics$wss)))
  metrics_norm$silhouette <- scale(metrics$silhouette, center = min(metrics$silhouette), 
                                   scale = diff(range(metrics$silhouette)))
  metrics_norm$calinski <- scale(metrics$calinski, center = min(metrics$calinski), 
                                 scale = diff(range(metrics$calinski)))
  
  # Calculate composite score
  metrics_norm$composite <- rowMeans(metrics_norm[, c("wss", "silhouette", "calinski")])
  
  # Find optimal k
  optimal_k <- metrics_norm$k[which.max(metrics_norm$composite)]
  
  # Create visualization
  metrics_long <- tidyr::pivot_longer(metrics_norm, 
                                      cols = c(wss, silhouette, calinski, composite),
                                      names_to = "metric", 
                                      values_to = "value")
  
  p <- ggplot(metrics_long, aes(x = k, y = value, color = metric)) +
    geom_line() +
    geom_point() +
    geom_vline(xintercept = optimal_k, linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(title = paste("Cluster Quality Metrics -", method),
         subtitle = paste("Optimal number of clusters:", optimal_k),
         x = "Number of clusters (k)",
         y = "Normalized Score",
         color = "Metric")
  
  print(p)
  
  return(list(
    optimal_k = optimal_k,
    metrics = metrics,
    normalized_metrics = metrics_norm
  ))
}

# Update the cluster_data function to use the new method
cluster_data <- function(umap_df, max_clusters = 10, method = "kmeans", 
                         visualize = TRUE, seed = 18, SubtypeLabel = TRUE) {
  supported_methods <- c("kmeans", "kmedoids", "hierarchical", "dbscan", "hdbscan", "svm")
  if (!method %in% supported_methods) {
    stop(sprintf("Method must be one of: %s", paste(supported_methods, collapse = ", ")))
  }
  
  clustering_data <- as.matrix(umap_df[, c("UMAP1", "UMAP2")])
  
  # Use new optimal cluster finding for applicable methods
  if (method %in% c("kmeans", "kmedoids", "hierarchical")) {
    optimal_result <- find_optimal_clusters_multi(clustering_data, method, max_clusters, seed)
    optimal_clusters <- optimal_result$optimal_k
  } else {
    # Use existing logic for density-based methods
    optimal_clusters <- find_optimal_clusters(clustering_data, method, max_clusters, seed)
  }
  
  # Perform clustering with optimal parameters
  set.seed(seed)
  if (method == "kmeans") {
    cluster_result <- kmeans(clustering_data, centers = optimal_clusters)
    cluster_assignments <- cluster_result$cluster
  } else if (method == "kmedoids") {
    cluster_result <- pam(clustering_data, k = optimal_clusters)
    cluster_assignments <- cluster_result$clustering
  } else if (method == "hierarchical") {
    dist_mat <- dist(clustering_data)
    hc <- hclust(dist_mat, method = "ward.D2")
    cluster_assignments <- cutree(hc, k = optimal_clusters)
  } else if (method == "dbscan") {
    eps <- sort(dbscan::kNNdist(clustering_data, k = 4))[optimal_clusters]
    cluster_result <- dbscan::dbscan(clustering_data, eps = eps, minPts = 4)
    cluster_assignments <- cluster_result$cluster
  } else if (method == "hdbscan") {
    cluster_result <- hdbscan(clustering_data, minPts = 4)
    cluster_assignments <- cluster_result$cluster
    optimal_clusters <- length(unique(cluster_assignments))
  } else if (method == "svm") {
    svm_model <- svm(clustering_data, y = NULL, type = "one-classification", 
                     kernel = "radial", nu = 0.1)
    cluster_assignments <- predict(svm_model, clustering_data)
  }
  
  # Add cluster information to the dataframe
  umap_df$Cluster <- as.factor(cluster_assignments)
  
  if (visualize) {
    # Create elbow plot if applicable
    if (method %in% c("kmeans", "kmedoids", "hierarchical")) {
      wss <- numeric(max_clusters)
      for (i in 2:max_clusters) {
        if (method == "kmeans") {
          temp <- kmeans(clustering_data, centers = i)
          wss[i] <- temp$tot.withinss
        } else if (method == "kmedoids") {
          temp <- pam(clustering_data, k = i)
          wss[i] <- temp$objective
        } else if (method == "hierarchical") {
          temp <- cutree(hc, k = i)
          wss[i] <- sum(sapply(1:i, function(cl) {
            points <- clustering_data[temp == cl, , drop = FALSE]
            if (nrow(points) > 1) sum((points - colMeans(points))^2) else 0
          }))
        }
      }
      
      elbow_df <- data.frame(
        Clusters = 2:max_clusters,
        WSS = wss[2:max_clusters]
      )
      
      p1 <- ggplot(elbow_df, aes(x = Clusters, y = WSS)) +
        geom_line() +
        geom_point() +
        theme_minimal() +
        labs(
          title = paste("Elbow Plot for", toupper(method)),
          x = "Number of Clusters",
          y = "Total Within-Cluster Sum of Squares"
        )
      print(p1)
    }
    
    # Plot clustering results
    if (SubtypeLabel){
      p2 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, 
                                color = Cluster,  
                                label = Subtype)) +
        geom_point(size = 3) +
        geom_text(vjust = 1.5, hjust = 1, show.legend = FALSE) +
        theme_minimal() +
        labs(
          title = paste("UMAP Embedding with", method, 
                        ifelse(method != "hdbscan", 
                               paste("(", optimal_clusters, "clusters)"), "")),
          x = "UMAP 1",
          y = "UMAP 2"
        ) +
        theme(legend.position = "right")
    } else {
      p2 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, 
                                color = Cluster)) +
        geom_point(size = 3) +
        theme_minimal() +
        labs(
          title = paste("UMAP Embedding with", method, 
                        ifelse(method != "hdbscan", 
                               paste("(", optimal_clusters, "clusters)"), "")),
          x = "UMAP 1",
          y = "UMAP 2"
        ) +
        theme(legend.position = "right")
    }
    print(p2)
  }
  
  return(list(
    "umap_df" = umap_df, 
    "optimal_clusters" = if(method != "hdbscan") optimal_clusters else NA,
    "algorithm" = method,
    "cluster_object" = if(exists("cluster_result")) cluster_result else NULL
  ))
}




################################################################################
############################## Spatial Embedding  Clusters #####################
################################################################################

# For your original use case
embeddings <- comprehensive_embedding(
  data = alignment,
  group_var = "Prediction",
  n_bins = 100,
  k_min_dist=0, k_max_dist=1,
)

# extract embeddings matrix components
z_matrix <- embeddings$z_matrix

# Compute umap
umap_df <- Zmat2umap(z_matrix)


# Using k-means
# clust_list <- cluster_data(umap_df, method = "kmeans", max_clusters=10) # basically gives a random answer every time
# Using k-medoids
clust_list <- cluster_data(umap_df, method = "kmedoids", max_clusters=10) # 3 clusters found with ripley
# # Using hierarchical clustering
# clust_list <- cluster_data(umap_df, method = "hierarchical", max_clusters=10) # 2 clusters found
# # Using DBSCAN
# clust_list <- cluster_data(umap_df, method = "dbscan", max_clusters=10) # 5 clusters found
# # Using HDBSCAN
# clust_list <- cluster_data(umap_df, method = "hdbscan", max_clusters=10) # 4 promising clusters with quadrant
# # Using SVM
# clust_list <- cluster_data(umap_df, method = "svm", max_clusters=10)# 2 clusters found


umap_df <- clust_list[['umap_df']]
optimal_clusters <- clust_list[['optimal_clusters']]

# Print cluster composition
cluster_composition <- table(umap_df$Cluster, umap_df$Subtype)
# print(cluster_composition)

# Add cluster IDs to the alignment dataframe
alignment <- alignment %>%
  group_by(Prediction) %>%
  mutate(Cluster = filter(umap_df, Subtype == first(Prediction))$Cluster) %>%
  ungroup()

# Visualize these clusters
(cluster_plot <- plot_cluster_densities(alignment))

# Write the output
write.csv(alignment, file=paste0(out_path, "/embedding_clusters_k",optimal_clusters,".csv"))
write.csv(cluster_composition,  file=paste0(out_path, "/cluster_composition_k",optimal_clusters,".csv"))
write.csv(umap_df,  file=paste0(out_path, "/clustered_umap_k",optimal_clusters,".csv"), row.names = F)



#####################################################################################

# Find the numerical row indices to match alignment
rgc_df$row_0 <- seq_len(nrow(rgc_df))

# Perform the join using row_0 and Prediction as conditions
# Using left_join to keep all rgc_df rows and add NA where no matches exist
rgc_df <- rgc_df %>%
  left_join(
    alignment,
    by = c("row_0", "Prediction")
  ) %>%
  drop_na()

meta_data <- c("cell", "row_0", # unique identifiers
               "Cluster", "Class", "Prediction", # Grouping Factors
               "X", "Y", # Registered XY position
               "x", "y", "z", # raw xyz position
               "nn_dist", "nn_id", # nearest neighbor
               "slide","slice", "sample", "retina", # experimental metadata
               "dapi_max", "dapi_min", "dapi_mean", "dapi_sd", "dapi_max_norm", 
                    "dapi_min_norm", "dapi_mean_norm", "dapi_sd_norm", # dapi values     
               "volume_hull", "volume", "x_range", "y_range", "z_range", "rect_vol",      
                    "elongation", "flatness", "sphericity", "compactness" # geometric values
               )

rgc_exp_mat_z <- rgc_df %>%
  distinct(across(all_of(c("cell", "x", "y"))), .keep_all = TRUE) %>%
  column_to_rownames("cell") %>%
  dplyr::select(-any_of(meta_data[-c(3,5, 6, 7)])) %>%  # Remove Cluster Prediction and XY from meta_data
  group_by(Cluster) %>%
  mutate(
    across(
      where(is.numeric) & !c(Prediction, X, Y),
      ~scale(.) %>% as.vector()
    )
  ) %>%
  ungroup()

#################################################################################
######################## Gene Embedding Functions ###############################
#################################################################################

library(progress)

create_gene_spatial_embedding <- function(data, n_bins = 50, 
                                          pos_threshold = 0.5, 
                                          neg_threshold = -0.5) {
  # Get unique clusters and genes
  clusters <- unique(data$Cluster)
  genes <- colnames(data)[!colnames(data) %in% c("Cluster", "Prediction", "X", "Y")]
  
  # Calculate total size of embedding vector
  grid_size <- n_bins * n_bins
  stats_size <- 6  # peak stats size
  radial_size <- 4  # radial stats size
  total_features <- 2 * grid_size +  # positive and negative grids
    2 * stats_size +   # positive and negative peak stats
    2 * radial_size    # positive and negative radial stats
  
  # Create progress bar
  total_steps <- length(clusters) * length(genes)
  pb <- progress_bar$new(
    format = "Processing cluster :cluster_id [:bar] :percent eta: :eta",
    total = total_steps,
    clear = FALSE,
    width = 80
  )
  
  # Initialize list to store results for each cluster
  cluster_embeddings <- list()
  
  # Rest of the helper functions remain the same...
  compute_grid_stats <- function(grid_values) {
    if (all(is.na(grid_values)) || length(grid_values) == 0) {
      return(rep(0, 6))
    }
    peak_val <- max(grid_values, na.rm = TRUE)
    peak_idx <- which.max(grid_values)
    peak_x <- (peak_idx - 1) %% n_bins / n_bins
    peak_y <- (peak_idx - 1) %/% n_bins / n_bins
    mean_val <- mean(grid_values, na.rm = TRUE)
    sd_val <- sd(grid_values, na.rm = TRUE)
    if (is.na(sd_val)) sd_val <- 0
    
    c(peak_val, peak_x, peak_y, mean_val, sd_val, 
      sum(!is.na(grid_values)) / length(grid_values))
  }
  
  compute_radial_symmetry <- function(grid_values) {
    if (all(is.na(grid_values)) || length(grid_values) == 0) {
      return(rep(0, 4))
    }
    grid_matrix <- matrix(grid_values, nrow = n_bins, ncol = n_bins)
    x_coords <- rep(1:n_bins, n_bins)
    y_coords <- rep(1:n_bins, each = n_bins)
    total_mass <- sum(grid_values, na.rm = TRUE)
    
    if (total_mass == 0) return(rep(0, 4))
    
    com_x <- sum(x_coords * grid_values, na.rm = TRUE) / total_mass
    com_y <- sum(y_coords * grid_values, na.rm = TRUE) / total_mass
    distances <- sqrt((x_coords - com_x)^2 + (y_coords - com_y)^2)
    weighted_distances <- distances * grid_values
    mean_dist <- mean(weighted_distances, na.rm = TRUE)
    sd_dist <- sd(weighted_distances, na.rm = TRUE)
    if (is.na(sd_dist)) sd_dist <- 0
    
    if (sd_dist == 0) {
      skew_dist <- 0
      kurt_dist <- 0
    } else {
      skew_dist <- mean((weighted_distances - mean_dist)^3, na.rm = TRUE) / (sd_dist^3)
      kurt_dist <- mean((weighted_distances - mean_dist)^4, na.rm = TRUE) / (sd_dist^4) - 3
    }
    
    stats <- c(mean_dist, sd_dist, skew_dist, kurt_dist)
    stats[!is.finite(stats)] <- 0
    stats
  }
  
  safe_kde2d <- function(x, y, n) {
    valid <- is.finite(x) & is.finite(y)
    x <- x[valid]
    y <- y[valid]
    
    if (length(x) < 2 || length(y) < 2) {
      return(list(
        x = seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = n),
        y = seq(min(y, na.rm = TRUE), max(y, na.rm = TRUE), length.out = n),
        z = matrix(0, n, n)
      ))
    }
    
    tryCatch({
      MASS::kde2d(x, y, n = n)
    }, error = function(e) {
      return(list(
        x = seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = n),
        y = seq(min(y, na.rm = TRUE), max(y, na.rm = TRUE), length.out = n),
        z = matrix(0, n, n)
      ))
    })
  }
  
  for (cluster_id in clusters) {
    pb$tick(0, tokens = list(cluster_id = cluster_id))
    
    cluster_data <- data[data$Cluster == cluster_id, ]
    
    # Initialize matrix with correct dimensions
    gene_embeddings <- matrix(0, 
                              nrow = length(genes), 
                              ncol = total_features)
    rownames(gene_embeddings) <- genes
    
    for (gene in genes) {
      gene_data <- cluster_data[, c("X", "Y", gene)]
      gene_data <- gene_data[is.finite(gene_data[[gene]]), ]
      
      # Process positive values
      pos_data <- gene_data[gene_data[[gene]] > pos_threshold, ]
      if (nrow(pos_data) > 0) {
        pos_grid <- safe_kde2d(pos_data$X, pos_data$Y, n = n_bins)
        pos_grid_vec <- as.vector(pos_grid$z)
      } else {
        pos_grid_vec <- rep(0, n_bins * n_bins)
      }
      
      # Process negative values
      neg_data <- gene_data[gene_data[[gene]] < neg_threshold, ]
      if (nrow(neg_data) > 0) {
        neg_grid <- safe_kde2d(neg_data$X, neg_data$Y, n = n_bins)
        neg_grid_vec <- as.vector(neg_grid$z)
      } else {
        neg_grid_vec <- rep(0, n_bins * n_bins)
      }
      
      # Compute statistics
      pos_peak_stats <- compute_grid_stats(pos_grid_vec)
      pos_radial_stats <- compute_radial_symmetry(pos_grid_vec)
      neg_peak_stats <- compute_grid_stats(neg_grid_vec)
      neg_radial_stats <- compute_radial_symmetry(neg_grid_vec)
      
      # Normalize grid vectors
      pos_grid_vec <- pos_grid_vec / sum(pos_grid_vec + .Machine$double.eps)
      neg_grid_vec <- neg_grid_vec / sum(neg_grid_vec + .Machine$double.eps)
      
      # Combine all features
      start_idx <- 1
      gene_embeddings[gene, start_idx:(start_idx + grid_size - 1)] <- pos_grid_vec
      start_idx <- start_idx + grid_size
      gene_embeddings[gene, start_idx:(start_idx + grid_size - 1)] <- neg_grid_vec
      start_idx <- start_idx + grid_size
      gene_embeddings[gene, start_idx:(start_idx + stats_size - 1)] <- pos_peak_stats
      start_idx <- start_idx + stats_size
      gene_embeddings[gene, start_idx:(start_idx + radial_size - 1)] <- pos_radial_stats
      start_idx <- start_idx + radial_size
      gene_embeddings[gene, start_idx:(start_idx + stats_size - 1)] <- neg_peak_stats
      start_idx <- start_idx + stats_size
      gene_embeddings[gene, start_idx:(start_idx + radial_size - 1)] <- neg_radial_stats
      
      pb$tick(1)
    }
    
    cluster_embeddings[[as.character(cluster_id)]] <- gene_embeddings
  }
  
  return(cluster_embeddings)
}

# Function to perform clustering on gene embeddings
cluster_gene_embeddings <- function(gene_embeddings, method = "hierarchical", max_clusters = 10) {
  # Initialize lists to store results
  cluster_assignments <- list()
  umap_dfs <- list()
  optimal_clusters <- list()
  clustering_results <- list()
  
  # Process each original cluster separately
  for (cluster_id in names(gene_embeddings)) {
    # Get features for this cluster
    cluster_mat <- gene_embeddings[[cluster_id]]
    
    # Create UMAP embedding
    set.seed(42)
    umap_result <- umap(cluster_mat)
    
    # Create UMAP dataframe
    umap_df <- data.frame(
      UMAP1 = umap_result$layout[,1],
      UMAP2 = umap_result$layout[,2],
      Gene = rownames(cluster_mat)
    )
    
    # Perform clustering
    clust_list <- cluster_data(umap_df, method = method, 
                               max_clusters = max_clusters, SubtypeLabel = FALSE)
    
    # Store results
    cluster_assignments[[cluster_id]] <- setNames(
      as.character(clust_list$umap_df$Cluster),
      clust_list$umap_df$Gene
    )
    umap_dfs[[cluster_id]] <- clust_list$umap_df
    optimal_clusters[[cluster_id]] <- clust_list$optimal_clusters
    clustering_results[[cluster_id]] <- clust_list
  }
  
  return(list(
    cluster_assignments = cluster_assignments,
    umap_dfs = umap_dfs,
    optimal_clusters = optimal_clusters,
    clustering_results = clustering_results
  ))
}

# Updated visualization function to handle per-cluster assignments
create_cluster_umaps <- function(gene_embeddings, color_labels = NULL) {
  # Initialize list to store plots and UMAP ranges
  umap_plots <- list()
  all_umap1 <- numeric()
  all_umap2 <- numeric()
  
  # First pass: create all UMAP embeddings and collect ranges
  umap_data <- list()
  for (cluster_id in names(gene_embeddings)) {
    # Get embedding matrix for this cluster
    embedding_matrix <- gene_embeddings[[cluster_id]]
    
    # Create UMAP embedding
    set.seed(42)  # For reproducibility
    umap_result <- umap(embedding_matrix)
    
    # Store UMAP data
    umap_data[[cluster_id]] <- list(
      coords = umap_result,
      genes = rownames(embedding_matrix)
    )
    
    # Collect all UMAP coordinates
    all_umap1 <- c(all_umap1, umap_result$layout[,1])
    all_umap2 <- c(all_umap2, umap_result$layout[,2])
  }
  
  # Calculate global ranges with some padding
  padding <- 0.05  # 5% padding
  x_range <- range(all_umap1)
  y_range <- range(all_umap2)
  x_padding <- diff(x_range) * padding
  y_padding <- diff(y_range) * padding
  x_limits <- c(x_range[1] - x_padding, x_range[2] + x_padding)
  y_limits <- c(y_range[1] - y_padding, y_range[2] + y_padding)
  
  # Second pass: create plots with consistent limits
  for (cluster_id in names(gene_embeddings)) {
    # Create data frame for plotting
    umap_df <- data.frame(
      UMAP1 = umap_data[[cluster_id]]$coords$layout[,1],
      UMAP2 = umap_data[[cluster_id]]$coords$layout[,2],
      Gene = umap_data[[cluster_id]]$genes
    )
    
    # Add color labels if provided
    if (!is.null(color_labels)) {
      umap_df$Cluster <- color_labels[[cluster_id]][umap_df$Gene]
    }
    
    # Create base plot
    p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2))
    
    # Add points with or without color
    if (!is.null(color_labels)) {
      p <- p + geom_point(aes(color = factor(Cluster)), size = 1, alpha = 0.6) +
        scale_color_brewer(palette = "Set1")
    } else {
      p <- p + geom_point(size = 1, alpha = 0.6)
    }
    
    # Add common theme elements
    p <- p + theme_minimal() +
      theme(
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 11, face = "bold")
      ) +
      labs(title = paste("Cluster", cluster_id, 
                         "\n(", length(unique(umap_df$Cluster)), " gene groups)")) +
      coord_fixed(
        xlim = x_limits,
        ylim = y_limits
      )
    
    umap_plots[[cluster_id]] <- p
  }
  
  # Combine plots using patchwork
  combined_plot <- wrap_plots(umap_plots, ncol = 2) +
    plot_annotation(
      title = "Gene Expression Pattern UMAP Embeddings by Cluster",
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
      )
    )
  
  # Return both the plot and the UMAP data for further use
  return(list(
    plot = combined_plot,
    umap_data = umap_data
  ))
}
# Function to create gene cluster assignment table
create_gene_cluster_table <- function(cluster_assignments) {
  # Get all unique genes across all clusters
  all_genes <- unique(unlist(lapply(cluster_assignments, names)))
  
  # Create an empty data frame with genes as row names
  assignment_df <- data.frame(
    row.names = all_genes
  )
  
  # Fill in cluster assignments for each original cluster
  for (cluster_id in names(cluster_assignments)) {
    # Get assignments for this cluster
    cluster_vec <- cluster_assignments[[cluster_id]]
    
    # Create column name
    col_name <- paste0("Cluster_", cluster_id)
    
    # Add to data frame
    assignment_df[[col_name]] <- cluster_vec[rownames(assignment_df)]
  }
  
  # Write to CSV
  write.csv(assignment_df, 
            file = paste0(out_path, "/gene_cluster_assignments.csv"))
  
  return(assignment_df)
}
#################################################################################
#################################################################################
#################################################################################



# Create embeddings
gene_embeddings <- create_gene_spatial_embedding(rgc_exp_mat_z,
                                                 pos_threshold = 0.01, 
                                                 neg_threshold = -0.01,
                                                 n_bins = 100)

# Perform clustering on each UMAP separately
gene_clusters <- cluster_gene_embeddings(gene_embeddings, 
                                         method = "dbscan", 
                                         max_clusters = 10)

# Create visualization with cluster colors
cluster_viz <- create_cluster_umaps(gene_embeddings, 
                                    color_labels = gene_clusters$cluster_assignments)

# Display the plot
print(cluster_viz$plot)

# Examine cluster composition for each original cluster
for(cluster_id in names(gene_clusters$cluster_assignments)) {
  cat("\nCluster", cluster_id, "gene group sizes:\n")
  print(table(gene_clusters$cluster_assignments[[cluster_id]]))
}




# Create and save the table
gene_cluster_table <- create_gene_cluster_table(gene_clusters$cluster_assignments)

# Print first few rows to verify
head(gene_cluster_table)



library(tidyverse)
library(patchwork)

visualize_gene_patterns <- function(data, gene_cluster_table, 
                                    ncol = 5, # number of genes per row
                                    point_size = 0.5,
                                    max_z = 3) { # cap z-scores for visualization
  
  # Get original cluster numbers from the table columns
  original_clusters <- sub("Cluster_", "", colnames(gene_cluster_table))
  
  # Initialize list for plots
  plots_list <- list()
  
  # Process each original cluster
  for(orig_clust in original_clusters) {
    # Get gene cluster assignments for this original cluster
    gene_clusters <- gene_cluster_table[[paste0("Cluster_", orig_clust)]]
    names(gene_clusters) <- rownames(gene_cluster_table)
    
    # Get data for this original cluster
    cluster_data <- data %>%
      dplyr::filter(Cluster == orig_clust) %>%
      dplyr::select(X, Y, all_of(names(gene_clusters)))
    
    # Get gene columns
    gene_cols <- names(gene_clusters)
    
    # Create a long format dataframe for plotting
    plot_data <- cluster_data %>%
      tidyr::pivot_longer(
        cols = all_of(gene_cols),
        names_to = "Gene",
        values_to = "Expression"
      ) %>%
      dplyr::mutate(
        # Cap z-scores for visualization
        Expression = pmin(pmax(Expression, -max_z), max_z),
        # Add gene cluster information
        GeneCluster = factor(gene_clusters[Gene])
      )
    
    # Create plot
    p <- ggplot(plot_data, aes(x = X, y = Y, color = Expression)) +
      geom_point(size = point_size) +
      scale_color_gradient2(
        low = "blue", mid = "white", high = "red",
        midpoint = 0, limits = c(-max_z, max_z)
      ) +
      facet_wrap(
        GeneCluster ~ Gene,
        ncol = ncol,
        scales = "free",
        labeller = labeller(Gene = label_wrap_gen(width = 10))
      ) +
      theme_minimal() +
      theme(
        strip.text = element_text(size = 8),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        panel.spacing = unit(0.1, "lines")
      ) +
      labs(
        title = paste("Original Cluster", orig_clust),
        subtitle = paste("Genes grouped by cluster assignment"),
        color = "Z-score"
      )
    
    plots_list[[orig_clust]] <- p
  }
  
  # Combine all plots using patchwork
  combined_plot <- wrap_plots(plots_list, ncol = 1) +
    plot_annotation(
      title = "Gene Expression Patterns by Cluster",
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
      )
    )
  
  return(combined_plot)
}

# Create visualization
gene_patterns <- visualize_gene_patterns(
  data = rgc_exp_mat_z,
  gene_cluster_table = gene_cluster_table,
  ncol = 5,  # adjust based on your preference
  point_size = 0.5,
  max_z = 3
)

# Display the plot
print(gene_patterns)

# Optionally save the plot
ggsave(
  paste0(out_path, "/gene_patterns.pdf"),
  gene_patterns,
  width = 20,
  height = 5 * length(unique(rgc_exp_mat_z$Cluster)),
  limitsize = FALSE
)



#################################################################################
# generate_grid_embeddings <- function(data, 
#                                      group_var, 
#                                      x_col = "X", 
#                                      y_col = "Y", 
#                                      genes = NULL,
#                                      n_bins = 100) {
#   # Input validation
#   if (!group_var %in% colnames(data)) {
#     stop(sprintf("Group variable '%s' not found in data", group_var))
#   }
#   if (!x_col %in% colnames(data) || !y_col %in% colnames(data)) {
#     stop("X or Y columns not found in data")
#   }
#   
#   # If genes not specified, use all numeric columns except X, Y, and grouping variables
#   if (is.null(genes)) {
#     genes <- data %>%
#       dplyr::select(where(is.numeric)) %>%
#       dplyr::select(-any_of(c(x_col, y_col))) %>%
#       colnames()
#   }
#   
#   # Create grid breaks
#   x_breaks <- seq(min(data[[x_col]]), max(data[[x_col]]), length.out = n_bins + 1)
#   y_breaks <- seq(min(data[[y_col]]), max(data[[y_col]]), length.out = n_bins + 1)
#   
#   # Get grid midpoints for output
#   x_grid <- x_breaks[-1] - diff(x_breaks)/2
#   y_grid <- y_breaks[-1] - diff(y_breaks)/2
#   
#   # Get unique groups
#   groups <- unique(data[[group_var]])
#   
#   # Initialize list to store all grid averages
#   all_grid_averages <- list()
#   
#   # Set up progress bar
#   total_iterations <- length(groups) * length(genes)
#   pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)
#   counter <- 0
#   
#   # Generate grid averages for each group-gene combination
#   for (group in groups) {
#     group_data <- data %>%
#       dplyr::filter(!!sym(group_var) == group)
#     
#     for (gene in genes) {
#       counter <- counter + 1
#       setTxtProgressBar(pb, counter)
#       message("\rProcessing cluster: ", group, "                    ", appendLF = FALSE)
#       
#       # Cut data into bins and compute mean gene expression
#       binned_data <- group_data %>%
#         dplyr::mutate(
#           x_bin = cut(!!sym(x_col), breaks = x_breaks, labels = FALSE),
#           y_bin = cut(!!sym(y_col), breaks = y_breaks, labels = FALSE)
#         ) %>%
#         dplyr::group_by(x_bin, y_bin) %>%
#         dplyr::summarise(
#           avg_expression = mean(!!sym(gene), na.rm = TRUE),
#           .groups = 'drop'
#         )
#       
#       # Create complete grid including empty bins
#       full_grid <- expand.grid(
#         x_bin = 1:n_bins,
#         y_bin = 1:n_bins
#       ) %>%
#         dplyr::left_join(binned_data, by = c("x_bin", "y_bin")) %>%
#         dplyr::mutate(avg_expression = replace(avg_expression, is.na(avg_expression), 0))
#       
#       # Convert to matrix
#       grid_matrix <- matrix(
#         full_grid$avg_expression,
#         nrow = n_bins,
#         ncol = n_bins
#       )
#       
#       # Store with informative name
#       map_name <- paste(group, gene, sep = "_")
#       all_grid_averages[[map_name]] <- grid_matrix
#     }
#   }
#   
#   # Close progress bar and add newline
#   close(pb)
#   message("")  # Add final newline
#   
#   # Create matrix to store grid values
#   z_matrix <- matrix(
#     NA,
#     nrow = length(all_grid_averages),
#     ncol = n_bins * n_bins,
#     dimnames = list(names(all_grid_averages), NULL)
#   )
#   
#   # Populate the matrix - identical to original kde2d version
#   for (i in seq_along(all_grid_averages)) {
#     z_matrix[i, ] <- as.vector(all_grid_averages[[i]])
#   }
#   
#   # Return both the raw grid averages and the matrix
#   return(list(
#     grid_averages = all_grid_averages,
#     embedding_matrix = z_matrix,
#     x_grid = x_grid,
#     y_grid = y_grid
#   ))
# }
# 
# library(uwot)
# library(ggplot2)
# 
# # Function to perform UMAP embedding
# create_umap_embedding <- function(z_matrix, n_neighbors = 15, min_dist = 0.1) {
#   set.seed(18)
#   umap_result <- umap(z_matrix, n_neighbors = n_neighbors, min_dist = min_dist)
#   
#   data.frame(
#     UMAP1 = umap_result$layout[,1],
#     UMAP2 = umap_result$layout[,2],
#     Subtype = rownames(z_matrix)
#   )
# }
# 
# # Function to calculate elbow curve
# calculate_elbow_curve <- function(umap_df, max_clusters = 10) {
#   wss <- numeric(max_clusters)
#   
#   for (i in 2:max_clusters) {
#     kmeans_temp <- kmeans(umap_df[, c("UMAP1", "UMAP2")], centers = i)
#     wss[i] <- kmeans_temp$tot.withinss
#   }
#   
#   data.frame(
#     Clusters = 2:max_clusters,
#     WSS = wss[2:max_clusters]
#   )
# }
# 
# # Function to find elbow point
# find_elbow_point <- function(elbow_df) {
#   x <- elbow_df$Clusters
#   y <- elbow_df$WSS
#   
#   # Calculate the angles between consecutive points
#   angles <- atan((y[-1] - y[-length(y)]) / (x[-1] - x[-length(x)]))
#   
#   # Find the point where the angle changes most significantly
#   elbow_index <- which.max(diff(angles)) + 1
#   
#   x[elbow_index]
# }
# 
# # Function to perform clustering and return results
# perform_clustering <- function(umap_df, n_clusters) {
#   set.seed(18)
#   kmeans_result <- kmeans(umap_df[, c("UMAP1", "UMAP2")], centers = n_clusters)
#   
#   # Add cluster information
#   umap_df$Cluster <- as.factor(kmeans_result$cluster)
#   
#   # Calculate cluster composition
#   cluster_composition <- table(umap_df$Cluster, umap_df$Subtype)
#   
#   list(
#     umap_df = umap_df,
#     cluster_composition = cluster_composition
#   )
# }
# 
# # Main function to analyze gene embeddings for a single cluster
# analyze_cluster_embeddings <- function(z_matrix, cluster_name, max_clusters = 10, 
#                                        make_plots = TRUE, out_path = NULL) {
#   # Create UMAP embedding
#   umap_df <- create_umap_embedding(z_matrix)
#   
#   # Calculate elbow curve
#   elbow_df <- calculate_elbow_curve(umap_df, max_clusters)
#   
#   # Find optimal number of clusters
#   optimal_clusters <- find_elbow_point(elbow_df)
#   
#   # Perform clustering
#   clustering_results <- perform_clustering(umap_df, optimal_clusters)
#   
#   if (make_plots) {
#     # Elbow plot
#     p1 <- ggplot(elbow_df, aes(x = Clusters, y = WSS)) +
#       geom_line() +
#       geom_point() +
#       theme_minimal() +
#       labs(
#         title = paste("Elbow Plot for Cluster", cluster_name),
#         x = "Number of Clusters",
#         y = "Total Within-Cluster Sum of Squares"
#       )
#     
#     # UMAP plot
#     p2 <- ggplot(clustering_results$umap_df, 
#                  aes(x = UMAP1, y = UMAP2, 
#                      color = Cluster, 
#                      label = Subtype)) +
#       geom_point(size = 3) +
#       geom_text(vjust = 1.5, hjust = 1, show.legend = FALSE) +
#       theme_minimal() +
#       labs(
#         title = paste("UMAP Embedding for Cluster", cluster_name, 
#                       "\nOptimal Clusters:", optimal_clusters),
#         x = "UMAP 1",
#         y = "UMAP 2"
#       ) +
#       theme(legend.position = "none")
#     
#     if (!is.null(out_path)) {
#       ggsave(paste0(out_path, "/elbow_plot_cluster_", cluster_name, ".pdf"), p1)
#       ggsave(paste0(out_path, "/umap_plot_cluster_", cluster_name, ".pdf"), p2)
#     }
#   }
#   
#   list(
#     optimal_clusters = optimal_clusters,
#     umap_df = clustering_results$umap_df,
#     cluster_composition = clustering_results$cluster_composition,
#     elbow_df = elbow_df
#   )
# }
# 
# analyze_all_clusters <- function(embeddings, out_path = NULL) {
#   # Get unique clusters
#   clusters <- unique(gsub("_.*", "", names(embeddings$grid_averages)))
#   
#   # Initialize results list
#   results <- list()
#   
#   # Process each cluster
#   for (cluster in clusters) {
#     message("Processing Cluster: ", cluster)
#     
#     # Filter z_matrix for current cluster
#     cluster_rows <- grepl(paste0("^", cluster, "_"), rownames(embeddings$embedding_matrix))
#     cluster_z_matrix <- embeddings$embedding_matrix[cluster_rows, ]
#     
#     # Analyze embeddings for this cluster
#     results[[cluster]] <- analyze_cluster_embeddings(
#       z_matrix = cluster_z_matrix,
#       cluster_name = cluster,
#       out_path = out_path
#     )
#     
#     message("Cluster ", cluster, " optimal number of subclusters: ", 
#             results[[cluster]]$optimal_clusters)
#   }
#   
#   # Create summary
#   summary_df <- data.frame(
#     Cluster = names(results),
#     Optimal_Subclusters = sapply(results, function(x) x$optimal_clusters)
#   )
#   
#   if (!is.null(out_path)) {
#     write.csv(summary_df, file = paste0(out_path, "/cluster_summary.csv"))
#   }
#   
#   list(
#     results = results,
#     summary = summary_df
#   )
# }
# 
# 
# 
# # Basic usage with default parameters
# embeddings <- generate_grid_embeddings(
#   data = rgc_exp_mat_z,
#   group_var = "Cluster",
#   n_bins = 100  # Adjust resolution as needed
# )
# 
# 
# test <- embeddings[['embedding_matrix']][['3_Isl2']]
# 
# # Assuming you have your embeddings from generate_grid_embeddings()
# results <- analyze_all_clusters(
#   embeddings = embeddings,
#   out_path = out_path
# )
# 
# # View summary of optimal clusters for each main cluster
# print(results$summary)
# 
# # Access detailed results for a specific cluster
# cluster1_results <- results$results[["1"]]
# 






################################################################################
# # Perform PCA on the z-matrix
# pca_result <- prcomp(z_matrix, scale. = TRUE)
# 
# # Create a data frame for plotting
# pca_df <- data.frame(
#   PC1 = pca_result$x[,1],
#   PC2 = pca_result$x[,2],
#   Subtype = rownames(z_matrix)
# )
# 
# 
# ggplot(pca_df, aes(x = PC1, y = PC2, color = Subtype, label = Subtype)) +
#   geom_point(size = 3) +
#   geom_text(vjust = 1.5, hjust = 1) +
#   theme_minimal() +
#   labs(
#     title = "PCA Embedding of Density Estimates by Subtype",
#     x = "First Principal Component",
#     y = "Second Principal Component"
#   ) +
#   theme(legend.position = "none")