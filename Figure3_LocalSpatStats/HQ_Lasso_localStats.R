library(spatstat)
library(sf)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

# Analysis hyperparameters
set.seed(18)
MIN_CELLS_PER_REGION <- 3  # Minimum number of cells required for analysis in a region

# Add debugging function
print_debug_info <- function(msg, data) {
  cat("\nDEBUG:", msg, "\n")
  if (is.data.frame(data)) {
    cat("Dimensions:", dim(data), "\n")
    if (nrow(data) > 0) {
      cat("First few rows:\n")
      print(head(data))
    }
  } else if (is.list(data)) {
    cat("List length:", length(data), "\n")
    if (length(data) > 0) {
      cat("First element type:", class(data[[1]]), "\n")
    }
  }
}



####################################################
# Function to calculate all pairwise distances for a subtype
calculate_neighbor_distances <- function(subtype_data, minimum_distance = 10) {
  if (nrow(subtype_data) < 2) {
    return(data.frame(
      distance = numeric(0)
    ))
  }
  
  # Create distance matrix
  dist_mat <- as.matrix(dist(subtype_data[, c("x", "y")]))
  
  # Get lower triangle (excluding diagonal)
  distances <- dist_mat[lower.tri(dist_mat)]
  
  # Filter out distances below minimum_distance if minimum_distance > 0
  if (minimum_distance > 0) {
    distances <- distances[distances >= minimum_distance]
    
    # If no distances remain after filtering, return empty dataframe
    if (length(distances) == 0) {
      return(data.frame(
        distance = numeric(0)
      ))
    }
  }
  
  return(data.frame(
    distance = distances
  ))
}
####################################################


# Function to pre-filter points within window
filter_points_in_window <- function(points_data, window_obj) {
  # Convert points to ppp first
  temp_ppp <- ppp(
    x = points_data$x,
    y = points_data$y,
    window = window_obj,
    check = FALSE
  )
  
  # Find which points are inside the window
  inside <- inside.owin(temp_ppp$x, temp_ppp$y, window_obj)
  
  # Return filtered dataframe
  points_data[inside, ]
}

# Function to load lasso regions with filtering for HQ regions
load_lasso_regions <- function(filepath) {
  cat(sprintf("\nLoading lasso regions from %s\n", basename(filepath)))
  
  # Read and print first few lines of raw file for debugging
  raw_lines <- readLines(filepath, n = 5)
  cat("\nFirst 5 lines of file:\n")
  print(raw_lines)
  
  # Skip the first two metadata rows, read from row 3
  coords_df <- read.csv(filepath, skip = 2)
  print_debug_info("Loaded coordinates dataframe", coords_df)
  
  # Filter for HQ_Lasso regions
  coords_df <- coords_df[grep("^HQ_Lasso_\\d+$", coords_df$Selection), ]
  print_debug_info("Filtered for HQ_Lasso regions", coords_df)
  
  # Split into list of polygons based on Selection ID
  lasso_list <- split(coords_df, coords_df$Selection) %>%
    lapply(function(region) {
      # Ensure the polygon is closed (first and last points match)
      if (!identical(region[1, c("X","Y")], region[nrow(region), c("X","Y")])) {
        region <- rbind(region, region[1,])
      }
      coords_matrix <- as.matrix(region[, c("X","Y")])
      st_polygon(list(coords_matrix))
    })
  
  print_debug_info("Created lasso region list", lasso_list)
  return(lasso_list)
}

# Utility function to convert sf polygon to owin
sf_to_owin <- function(sf_poly) {
  # Get coordinates from sf polygon
  coords <- st_coordinates(sf_poly)
  
  # Create owin polygon
  owin(poly = list(x = coords[, "X"], y = coords[, "Y"]))
}

# Test function to plot a lasso region with data points
plot_lasso_region <- function(lasso_region, points_data, region_name) {
  # Convert lasso region to sf object and owin
  lasso_sf <- st_sf(geometry = st_sfc(lasso_region))
  lasso_window <- sf_to_owin(lasso_region)
  
  # Filter points to only those within the lasso region
  filtered_points <- filter_points_in_window(points_data, lasso_window)
  
  # Create plot
  p <- ggplot() +
    geom_sf(data = lasso_sf, fill = NA, color = "black") +
    geom_point(data = filtered_points, aes(x = x, y = y, color = Prediction), size = 1) +
    theme_minimal()  +
    theme(
      legend.position = 'None',
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(color = "black")
    ) +
    labs(title = paste("Lasso Region:", region_name),
         subtitle = sprintf("Points in region: %d", nrow(filtered_points))) +
    coord_sf() +
    theme(legend.position = 'None')
  
  print(p)
  return(p)
}

# Function to cluster nearby points of the same subtype
cluster_nearby_points <- function(points_data, minimum_distance = 10) {
  if (nrow(points_data) < 2) {
    return(points_data)
  }
  
  # Initialize cluster assignments
  cluster_id <- 1
  clusters <- rep(NA, nrow(points_data))
  processed <- rep(FALSE, nrow(points_data))
  
  # Create distance matrix once
  dist_mat <- as.matrix(dist(points_data[, c("x", "y")]))
  
  # Find clusters
  for (i in 1:nrow(points_data)) {
    if (!processed[i]) {
      # Find all points within minimum_distance of current point
      nearby <- which(dist_mat[i,] < minimum_distance)
      
      # Only cluster points of the same subtype
      same_subtype <- points_data$Prediction[nearby] == points_data$Prediction[i]
      nearby <- nearby[same_subtype]
      
      if (length(nearby) > 1) {  # If we found points to cluster
        clusters[nearby] <- cluster_id
        processed[nearby] <- TRUE
        cluster_id <- cluster_id + 1
      } else {
        clusters[i] <- cluster_id
        processed[i] <- TRUE
        cluster_id <- cluster_id + 1
      }
    }
  }
  
  # Calculate mean positions for each cluster
  clustered_points <- data.frame()
  for (i in unique(clusters)) {
    cluster_data <- points_data[clusters == i,]
    
    # If cluster has multiple points, use mean position
    if (nrow(cluster_data) > 1) {
      mean_point <- data.frame(
        x = mean(cluster_data$x),
        y = mean(cluster_data$y),
        Prediction = cluster_data$Prediction[1]  # All points in cluster have same prediction
      )
      clustered_points <- rbind(clustered_points, mean_point)
    } else {
      # If point is alone in its cluster, keep original position
      clustered_points <- rbind(clustered_points, cluster_data[, c("x", "y", "Prediction")])
    }
  }
  
  return(clustered_points)
}

create_ppp <- function(points_data, window_obj, 
                       cluster_points =TRUE,
                       minimum_distance = 10 ) {
  # Filter points to only those within the window
  temp_ppp <- ppp(
    x = points_data$x,
    y = points_data$y,
    window = window_obj,
    check = FALSE
  )
  
  inside <- inside.owin(temp_ppp$x, temp_ppp$y, window_obj)
  
  if (cluster_points) {
    points_in_window <- points_data[inside,]
    
    if(nrow(points_in_window) < 1) {
      return(NULL)
    }
    
    # Apply clustering to points within the window
    clustered_points <- cluster_nearby_points(points_in_window, minimum_distance)
    
    # Create final ppp object with clustered points
    ppp_obj <- ppp(
      x = clustered_points$x,
      y = clustered_points$y,
      window = window_obj,
      check = TRUE
    )
    
    return(ppp_obj)
  } else {
    if(sum(inside) < 1) {
      return(NULL)
    }
    
    ppp_obj <- ppp(
      x = temp_ppp$x[inside],
      y = temp_ppp$y[inside],
      window = window_obj,
      check = TRUE
    )
    
    return(ppp_obj)
  }
}

# Function to calculate spatial statistics for a point pattern
calculate_spatial_stats <- function(ppp_obj, minimum_distance = 10) {
  if (npoints(ppp_obj) < 3) {
    return(list(nnri = NA, vdri = NA))
  }
  
  # Calculate edge-corrected nearest neighbor distances using border method
  all_nn <- nndist(ppp_obj) # border correction is now default
  
  # Filter out nearest neighbors that are too close
  nn <- all_nn[all_nn >= minimum_distance]
  
  # If no valid distances remain after filtering, return NA
  if(length(nn) < 3) {
    return(list(nnri = NA, vdri = NA))
  }
  
  # NNRI with edge correction
  nnri <- mean(nn) / sd(nn)
  
  # Calculate VDRI using tessellation with edge handling
  # Create Voronoi tessellation
  voro <- dirichlet(ppp_obj)
  
  # Create buffer to identify edge cells
  win <- Window(ppp_obj)
  # Use area.owin() instead of area()
  buffer_distance <- sqrt(area.owin(win)) * 0.05  # 10% of square root of area as buffer
  
  # Identify points that are far enough from the boundary
  inner_points <- bdist.points(ppp_obj) > buffer_distance
  
  if(sum(inner_points) < 3) {
    return(list(nnri = nnri, vdri = NA))
  }
  
  # Get areas only for cells not on the edge
  voro_areas <- tile.areas(voro)[inner_points]
  
  # Calculate VDRI using only interior cells
  vdri <- mean(voro_areas) / sd(voro_areas)
  
  return(list(
    nnri = nnri,
    vdri = vdri,
    n_edge = sum(!inner_points),  # number of edge cells
    n_interior = sum(inner_points),  # number of interior cells
    nn = nn # nearest neighbor distances                 
  ))
}

analyze_subtype_in_region <- function(subtype_data, all_region_data, region_window, 
                                      n_bootstrap, minimum_distance = 10) {
  # Create observed point pattern
  obs_ppp <- ppp(
    x = subtype_data$x,
    y = subtype_data$y,
    window = region_window,
    check = TRUE
  )
  
  # Calculate observed statistics
  obs_stats <- calculate_spatial_stats(obs_ppp)
  
  # Get the filtered nearest neighbor distances
  nn_distances <- data.frame(
    distance = obs_stats$nn,
    Subtype = unique(subtype_data$Prediction)
  )
  
  # Define the columns we want in our results
  result_columns <- c(
    "NNRI", "VDRI", "N", "N_region", "Observed", "Subtype",
    "NNRI_min", "NNRI_max", "NNRI_CI_95_lower", "NNRI_CI_95_upper",
    "NNRI_CI_99_lower", "NNRI_CI_99_upper",
    "VDRI_min", "VDRI_max", "VDRI_CI_95_lower", "VDRI_CI_95_upper",
    "VDRI_CI_99_lower", "VDRI_CI_99_upper"
  )
  
  # Initialize results dataframe with observed values
  results <- data.frame(
    NNRI = obs_stats$nnri,
    VDRI = obs_stats$vdri,
    N = npoints(obs_ppp),
    N_region = nrow(all_region_data),
    Observed = TRUE,
    Subtype = unique(subtype_data$Prediction),
    NNRI_min = NA_real_,
    NNRI_max = NA_real_,
    NNRI_CI_95_lower = NA_real_,
    NNRI_CI_95_upper = NA_real_,
    NNRI_CI_99_lower = NA_real_,
    NNRI_CI_99_upper = NA_real_,
    VDRI_min = NA_real_,
    VDRI_max = NA_real_,
    VDRI_CI_95_lower = NA_real_,
    VDRI_CI_95_upper = NA_real_,
    VDRI_CI_99_lower = NA_real_,
    VDRI_CI_99_upper = NA_real_
  )[, result_columns]
  
  # Initialize storage for bootstrap results
  bootstrap_results <- list()
  attempts <- 0
  max_attempts <- n_bootstrap * 100
  
  # Create progress bar
  pb <- txtProgressBar(min = 0, max = n_bootstrap, style = 3)
  
  # While loop for bootstrap
  while(length(bootstrap_results) < n_bootstrap && attempts < max_attempts) {
    attempts <- attempts + 1
    
    # Shuffle the subtype labels
    shuffled_data <- all_region_data
    shuffled_data$Prediction <- sample(shuffled_data$Prediction)
    
    # Get points for this subtype from shuffled data
    sampled_points <- shuffled_data %>% 
      filter(Prediction == unique(subtype_data$Prediction))
    
    # Create point pattern
    boot_ppp <- ppp(
      x = sampled_points$x,
      y = sampled_points$y,
      window = region_window,
      check = TRUE
    )
    
    # Calculate statistics
    boot_stats <- calculate_spatial_stats(boot_ppp)
    
    # Only add valid results
    if(!is.na(boot_stats$nnri) || !is.na(boot_stats$vdri)) {
      boot_result <- data.frame(
        NNRI = boot_stats$nnri,
        VDRI = boot_stats$vdri,
        N = npoints(boot_ppp),
        N_region = nrow(all_region_data),
        Observed = FALSE,
        Subtype = unique(subtype_data$Prediction),
        NNRI_min = NA_real_,
        NNRI_max = NA_real_,
        NNRI_CI_95_lower = NA_real_,
        NNRI_CI_95_upper = NA_real_,
        NNRI_CI_99_lower = NA_real_,
        NNRI_CI_99_upper = NA_real_,
        VDRI_min = NA_real_,
        VDRI_max = NA_real_,
        VDRI_CI_95_lower = NA_real_,
        VDRI_CI_95_upper = NA_real_,
        VDRI_CI_99_lower = NA_real_,
        VDRI_CI_99_upper = NA_real_
      )[, result_columns]
      
      bootstrap_results[[length(bootstrap_results) + 1]] <- boot_result
      
      # Update progress bar
      setTxtProgressBar(pb, length(bootstrap_results))
    }
  }
  
  # Close progress bar
  close(pb)
  
  # Warning if we didn't get enough samples
  if(length(bootstrap_results) < n_bootstrap) {
    warning(sprintf(
      "Could only generate %d valid bootstrap samples out of %d requested after %d attempts for subtype %s",
      length(bootstrap_results), n_bootstrap, attempts, unique(subtype_data$Prediction)
    ))
  }
  
  # Combine results and compute CIs for bootstrap data
  if(length(bootstrap_results) > 0) {
    boot_df <- do.call(rbind, bootstrap_results)
    
    # Compute confidence intervals
    boot_nnri <- boot_df$NNRI[!is.na(boot_df$NNRI)]  # Remove NAs
    boot_vdri <- boot_df$VDRI[!is.na(boot_df$VDRI)]  # Remove NAs
    
    # Update the first row of boot_df (bootstrap summary)
    boot_df[1, c("NNRI_min", "NNRI_max")] <- c(min(boot_nnri), max(boot_nnri))
    boot_df[1, c("NNRI_CI_95_lower", "NNRI_CI_95_upper")] <- quantile(boot_nnri, c(0.025, 0.975))
    boot_df[1, c("NNRI_CI_99_lower", "NNRI_CI_99_upper")] <- quantile(boot_nnri, c(0.005, 0.995))
    
    boot_df[1, c("VDRI_min", "VDRI_max")] <- c(min(boot_vdri), max(boot_vdri))
    boot_df[1, c("VDRI_CI_95_lower", "VDRI_CI_95_upper")] <- quantile(boot_vdri, c(0.025, 0.975))
    boot_df[1, c("VDRI_CI_99_lower", "VDRI_CI_99_upper")] <- quantile(boot_vdri, c(0.005, 0.995))
    
    # Remove CI columns from all but first row of bootstrap results
    boot_df[-1, grep("_CI_|_min|_max", names(boot_df))] <- NA
    
    results <- rbind(results, boot_df)
  }
  
  # Add diagnostic information
  attr(results, "attempts") <- attempts
  attr(results, "valid_samples") <- length(bootstrap_results)
  
  return(list(
    results = results,
    nn_distances = nn_distances
  ))
}


analyze_region <- function(region_data, region_window, n_bootstrap = 1000) {
  subtypes <- unique(region_data$Prediction)
  all_results <- list()
  all_distances <- list()
  all_nn_distances <- list()
  
  for(subtype in subtypes) {
    cat(sprintf("Analyzing subtype %s...\n", subtype))
    
    # Get data for this subtype
    subtype_data <- region_data %>% filter(Prediction == subtype)
    
    # Only analyze if we have enough cells
    if(nrow(subtype_data) >= 3) {
      # Original spatial statistics analysis
      analysis <- analyze_subtype_in_region(
        subtype_data = subtype_data,
        all_region_data = region_data,
        region_window = region_window,
        n_bootstrap = n_bootstrap
      )
      
      # Add metadata to nn_distances
      if (!is.null(analysis$nn_distances) && nrow(analysis$nn_distances) > 0) {
        # Metadata should be available in region_data
        analysis$nn_distances$slide <- unique(region_data$slide)
        analysis$nn_distances$region <- unique(region_data$region)
      }
      
      # Add metadata
      analysis$Subtype <- subtype
      all_results[[length(all_results) + 1]] <- analysis$results
      all_nn_distances[[length(all_nn_distances) + 1]] <- analysis$nn_distances
      
      
      # Calculate neighbor distances
      distances <- calculate_neighbor_distances(subtype_data)
      if (nrow(distances) > 0) {
        distances$Subtype <- subtype
        all_distances[[length(all_distances) + 1]] <- distances
      }
    } else {
      cat(sprintf("Skipping subtype %s: insufficient cells (%d)\n", 
                  subtype, nrow(subtype_data)))
    }
  }
  
  # Return both results and distances
  return(list(
    results = do.call(rbind, all_results),
    distances = do.call(rbind, all_distances),
    nn_distances = do.call(rbind, all_nn_distances)
    
  ))
}

create_summary_statistics <- function(results_df) {
  # First group by Subtype and Observed
  summary_stats <- results_df %>%
    group_by(Subtype, Observed) %>%
    summarize(
      # Cell counts
      total_cells = sum(N, na.rm = TRUE),
      total_region_cells = sum(N_region, na.rm = TRUE),
      
      # NNRI statistics
      NNRI_mean = mean(NNRI, na.rm = TRUE),
      NNRI_median = median(NNRI, na.rm = TRUE),
      NNRI_min = min(c(NNRI, NNRI_min), na.rm = TRUE),  # Include computed min
      NNRI_max = max(c(NNRI, NNRI_max), na.rm = TRUE),  # Include computed max
      NNRI_CI_95_lower = quantile(NNRI, 0.025, na.rm = TRUE),
      NNRI_CI_95_upper = quantile(NNRI, 0.975, na.rm = TRUE),
      NNRI_CI_99_lower = quantile(NNRI, 0.005, na.rm = TRUE),
      NNRI_CI_99_upper = quantile(NNRI, 0.995, na.rm = TRUE),
      
      # VDRI statistics
      VDRI_mean = mean(VDRI, na.rm = TRUE),
      VDRI_median = median(VDRI, na.rm = TRUE),
      VDRI_min = min(c(VDRI, VDRI_min), na.rm = TRUE),  # Include computed min
      VDRI_max = max(c(VDRI, VDRI_max), na.rm = TRUE),  # Include computed max
      VDRI_CI_95_lower = quantile(VDRI, 0.025, na.rm = TRUE),
      VDRI_CI_95_upper = quantile(VDRI, 0.975, na.rm = TRUE),
      VDRI_CI_99_lower = quantile(VDRI, 0.005, na.rm = TRUE),
      VDRI_CI_99_upper = quantile(VDRI, 0.995, na.rm = TRUE),
      
      # Add number of observations for reference
      n_observations = n(),
      .groups = 'drop'
    )
  
  return(summary_stats)
}
# Modified analyze_all_regions function
analyze_all_regions <- function(rgc_df, lasso_regions, slide_id) {
  all_results <- list()
  all_distances <- list()
  all_nn_distances <- list() 
  
  for(i in seq_along(lasso_regions)) {
    region_name <- names(lasso_regions)[i]
    cat(sprintf("\nAnalyzing region %s...\n", region_name))
    
    # Convert region to owin
    region_window <- sf_to_owin(lasso_regions[[i]])
    
    # Filter points for this region
    region_data <- filter_points_in_window(rgc_df, region_window)
    
    if(nrow(region_data) >= 3) {
      # Analyze region
      analysis <- analyze_region(region_data, region_window)
      
      # Add metadata to results
      analysis$results$slide <- slide_id
      analysis$results$region <- region_name
      
      # Add metadata to distances
      if (!is.null(analysis$distances) && nrow(analysis$distances) > 0) {
        analysis$distances$slide <- slide_id
        analysis$distances$region <- region_name
      }
      
      if (!is.null(analysis$nn_distances) && nrow(analysis$nn_distances) > 0) {
        analysis$nn_distances$slide <- slide_id
        analysis$nn_distances$region <- region_name
      }
      
      all_results[[length(all_results) + 1]] <- analysis$results
      all_distances[[length(all_distances) + 1]] <- analysis$distances
      all_nn_distances[[length(all_nn_distances) + 1]] <- analysis$nn_distances
      
    } else {
      cat(sprintf("Skipping region %s: insufficient cells (%d)\n", 
                  region_name, nrow(region_data)))
    }
  }
  
  # Return both combined results
  return(list(
    results = do.call(rbind, all_results),
    distances = do.call(rbind, all_distances),
    nn_distances = do.call(rbind, all_nn_distances)
  ))
}

# Run analysis on all slides and create visualizations
run_full_analysis <- function(viz_root, all_subtypes = T, subtype_only = '30_Novel') {
  # Create visualization directory
  viz_dir <- file.path(viz_root, "NN_visualizations")
  dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Get all lasso files
  lasso_files <- list.files(lasso_root, pattern = "\\d+_HQ_Lasso_coordinates\\.csv$", 
                            full.names = TRUE)
  cat(sprintf("Found %d lasso files\n", length(lasso_files)))
  
  # Initialize list for results
  # Initialize lists for results
  all_results <- list()
  all_distances <- list()
  all_nn_distances <- list()
  
  # Process each slide
  for(file in lasso_files) {
    # Get slide ID
    slide_id <- as.numeric(gsub("(\\d+)_HQ_Lasso_coordinates\\.csv", "\\1", 
                                basename(file)))
    cat(sprintf("\nProcessing slide %d...\n", slide_id))
    
    # Load lasso regions
    lasso_list <- load_lasso_regions(file)
    
    # Filter RGC data for this slide
    slide_data <- rgc_df %>% filter(slide == slide_id)
    
    if(nrow(slide_data) > 0) {
      # Create point process plots for each region
      for(i in seq_along(lasso_list)) {
        region_name <- names(lasso_list)[i]
        p <- plot_lasso_region(lasso_list[[i]], slide_data, region_name)
        ggsave(file.path(viz_dir, sprintf("slide%d_%s_points.png", slide_id, region_name)), 
               p, width = 10, height = 10)
      }
      
      if (!all_subtypes) {
        slide_data <- slide_data %>%
          filter(Prediction == subtype_only)
      }
      
      # Run analysis
      analysis <- analyze_all_regions(slide_data, lasso_list, slide_id)
      all_results[[length(all_results) + 1]] <- analysis$results
      all_distances[[length(all_distances) + 1]] <- analysis$distances
      all_nn_distances[[length(all_nn_distances) + 1]] <- analysis$nn_distances
    } else {
      cat(sprintf("Skipping slide %d: no matching RGC data\n", slide_id))
    }
  }
  
  # Combine and save all results
  final_results <- do.call(rbind, all_results)
  final_distances <- do.call(rbind, all_distances)
  final_nn_distances <- do.call(rbind, all_nn_distances)
  
  # Create summary statistics
  summary_results <- create_summary_statistics(final_results)
  
  write_csv(filter(summary_results, Observed), file.path(viz_dir, "summary_spatial_analysis_real.csv"))
  
  write_csv(filter(summary_results, !Observed), file.path(viz_dir, "summary_spatial_analysis_boot.csv"))
  
  write_csv(final_results, file.path(viz_dir, "full_spatial_analysis.csv"))
  write_csv(final_distances, file.path(viz_dir, "neighbor_distances.csv"))
  write_csv(final_nn_distances, file.path(viz_dir, "nearest_neighbor_distances.csv"))
  
  # Create combined region identifier
  final_distances <- final_distances %>%
    mutate(region_id = paste(slide, region, sep="_"))
  
  # Create histogram for each subtype
  subtypes <- unique(final_distances$Subtype)
  for(subtype in subtypes) {
    subtype_data <- final_distances %>% filter(Subtype == subtype)
    
    p <- ggplot(subtype_data, aes(x = distance, fill = region_id)) +
      geom_histogram(position = "stack", bins = 50) +
      scale_fill_viridis_d() +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid = element_blank()
      ) +
      labs(title = paste("Neighbor Distance Distribution:", subtype),
           x = "Distance",
           y = "Count")
    
    ggsave(file.path(viz_dir, sprintf("neighbor_dist_%s.png", subtype)), 
           p, width = 10, height = 6)
    
    subtype_data <- final_nn_distances %>% filter(Subtype == subtype)
    
    p <- ggplot(subtype_data, aes(x = distance)) +
      geom_histogram(bins = 50) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid = element_blank()
      ) +
      labs(title = paste("Nearest Neighbor Distance Distribution:", subtype),
           x = "Distance",
           y = "Count")
    
    ggsave(file.path(viz_dir, sprintf("nn_dist_%s.png", subtype)), 
           p, width = 10, height = 6)
  }
  # Create NNRI visualization
  p_nnri <- final_results %>%
    filter(N > 10) %>%
    ggplot(aes(x = NNRI, y = reorder(Subtype, N), color = Observed)) +
    geom_jitter(size = 0.4, height = 0.2, alpha = 0.5) +
    geom_point(data = filter(final_results, Observed == TRUE & N > 10), 
               color = 'black') +
    theme_minimal() +
    theme(
      legend.position = 'None',
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(color = "black")
    ) +
    labs(title = "Spatial Analysis Results",
         subtitle = "NNRI Distribution by Subtype") +
    facet_grid(slide ~ region)
  
  # Create VDRI visualization
  p_vdri <- final_results %>%
    filter(N > 10) %>%
    ggplot(aes(x = VDRI, y = reorder(Subtype, N), color = Observed)) +
    geom_jitter(size = 0.4, height = 0.2, alpha = 0.5) +
    geom_point(data = filter(final_results, Observed == TRUE & N > 10), 
               color = 'black') +
    theme_minimal() +
    theme(
      legend.position = 'None',
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(color = "black")
    ) +
    labs(title = "Spatial Analysis Results",
         subtitle = "VDRI Distribution by Subtype") +
    facet_grid(slide ~ region)
  
  # Save plots
  ggsave(file.path(viz_dir, "NNRI_analysis.png"), p_nnri, 
         width = 15, height = 20, limitsize = FALSE)
  ggsave(file.path(viz_dir, "VDRI_analysis.png"), p_vdri, 
         width = 15, height = 20, limitsize = FALSE)
  
  return(list(
    results = final_results,
    distances = final_distances,
    nn_distances = final_nn_distances
  ))
}
#################################################################################

####################################################################################
# Function to add significance stars
get_significance_stars <- function(p_value) {
  if(is.na(p_value)) return("NA")
  if(p_value <= 0.001) return("***")
  if(p_value <= 0.01) return("**")
  if(p_value <= 0.05) return("*")
  return("ns")
}

# Function to add significance stars
# Create label function for sample sizes
create_label <- function(row, sig_col) {
  sprintf("%s (n=%d)", 
          row[[sig_col]], 
          row$n_observed)
}

# Create statistical summary and perform tests
create_statistical_summary <- function(results, viz_root, min_cells = 10, min_observations = 3) {
  # Create visualization directory
  viz_dir <- file.path(viz_root, "NN_visualizations")
  dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)
  
  # First filter for minimum cells before aggregation
  filtered_results <- results %>% filter(N >= min_cells)
  
  # Create summary statistics as before
  test_results <- filtered_results %>%
    group_by(Subtype, slide, region, Observed) %>%
    summarize(
      NNRI = mean(NNRI),
      VDRI = mean(VDRI),
      N = mean(N),
      N_region = mean(N_region),
      .groups = 'drop'
    )
  
  # Sort subtypes alphabetically
  subtypes <- sort(unique(test_results$Subtype), decreasing = TRUE)
  stat_results <- list()
  
  for(st in subtypes) {
    subtype_data <- test_results %>% filter(Subtype == st)
    
    # Get observed and null values
    observed <- subtype_data %>% filter(Observed)
    null <- subtype_data %>% filter(!Observed)
    
    # Initialize results with NA
    t_test_result <- list(
      Subtype = st,
      mean_N = mean(observed$N),
      n_observed = nrow(observed),
      n_null = nrow(null),
      mean_observed_nnri = mean(observed$NNRI),
      mean_null_nnri = mean(null$NNRI),
      p_value_nnri = NA,
      mean_observed_vdri = mean(observed$VDRI),
      mean_null_vdri = mean(null$VDRI),
      p_value_vdri = NA
    )
    
    # Only perform tests if we have enough observations
    if(nrow(observed) >= min_observations && nrow(null) >= min_observations) {
      # Try NNRI t-test
      tryCatch({
        t_test_nnri <- t.test(observed$NNRI, null$NNRI)
        t_test_result$p_value_nnri <- t_test_nnri$p.value
      }, error = function(e) {
        warning(sprintf("Could not perform NNRI t-test for subtype %s: %s", st, e$message))
      })
      
      # Try VDRI t-test
      tryCatch({
        t_test_vdri <- t.test(observed$VDRI, null$VDRI)
        t_test_result$p_value_vdri <- t_test_vdri$p.value
      }, error = function(e) {
        warning(sprintf("Could not perform VDRI t-test for subtype %s: %s", st, e$message))
      })
    }
    
    # Add significance stars
    t_test_result$significance_nnri <- get_significance_stars(t_test_result$p_value_nnri)
    t_test_result$significance_vdri <- get_significance_stars(t_test_result$p_value_vdri)
    
    # Store results
    stat_results[[st]] <- as.data.frame(t_test_result)
  }
  
  # Combine statistical results
  stat_summary <- bind_rows(stat_results)
  
  # Create label function for sample sizes
  create_label <- function(row, sig_col) {
    sprintf("%s (n=%d)", 
            row[[sig_col]], 
            row$n_observed)
  }
  
  # Add labels to summary
  stat_summary$label_nnri <- mapply(create_label, split(stat_summary, 1:nrow(stat_summary)), 
                                    MoreArgs = list(sig_col = "significance_nnri"))
  stat_summary$label_vdri <- mapply(create_label, split(stat_summary, 1:nrow(stat_summary)), 
                                    MoreArgs = list(sig_col = "significance_vdri"))
  
  # Define colors
  box_colors <- c("TRUE" = "black", "FALSE" = "grey80")
  
  # Function to create plot
  create_metric_plot <- function(data, stat_data, metric, label_col) {
    ggplot(data %>% mutate(Subtype = factor(Subtype, levels = subtypes)), 
           aes(x = .data[[metric]], y = Subtype)) +
      geom_boxplot(aes(color = Observed, fill = Observed), 
                   outlier.size = 0.4,
                   alpha = 0.5) +
      scale_color_manual(values = box_colors) +
      scale_fill_manual(values = box_colors) +
      scale_x_continuous(expand = expansion(mult = c(0.05, 0.2))) +
      geom_text(data = stat_data,
                aes(x = max(data[[metric]]), 
                    label = .data[[label_col]]),
                hjust = 0,
                size = 3,
                vjust = 0.5) +
      theme_minimal() +
      theme(
        legend.position = 'None',
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(color = "black")
      ) +
      labs(title = "Spatial Analysis Results",
           subtitle = sprintf("%s Distribution by Subtype (min cells = %d, min obs = %d)", 
                              metric, min_cells, min_observations),
           caption = "* p<0.05, ** p<0.01, *** p<0.001")
  }
  
  # Create plots
  p_nnri <- create_metric_plot(test_results, stat_summary, "NNRI", "label_nnri")
  p_vdri <- create_metric_plot(test_results, stat_summary, "VDRI", "label_vdri")
  
  # Save plots
  ggsave(file.path(viz_dir, "statistical_summary_NNRI.png"), p_nnri, 
         width = 12, height = 15)
  ggsave(file.path(viz_dir, "statistical_summary_VDRI.png"), p_vdri, 
         width = 12, height = 15)
  
  # Save statistical summary
  write_csv(stat_summary, file.path(viz_dir, "statistical_summary.csv"))
  
  return(list(
    summary = stat_summary,
    plot_nnri = p_nnri,
    plot_vdri = p_vdri
  ))
}

########################################################################################################
# Function to get text annotation positioned correctly in the plot
get_stats_label <- function(stats, window) {
  x_range <- diff(window$xrange)
  y_range <- diff(window$yrange)
  data.frame(
    x = window$xrange[1] + x_range * 0.05,
    y = window$yrange[2] - y_range * 0.1,
    label = sprintf("NNRI: %.3f\nVDRI: %.3f", stats$nnri, stats$vdri)
  )
}

# Main inspection function
inspect_region_subtype <- function(rgc_df, 
                                   target_slide, target_region, target_subtype,
                                   results_df = NULL,
                                   lasso_root = '/media/sam/New Volume/Xenium_Data/HQ_NearestNeighbor_Zones') {
  # Generate new random seed for this run
  set.seed(as.integer(Sys.time()))
  
  # Load the lasso region
  lasso_file <- list.files(lasso_root, 
                           pattern = paste0(target_slide, "_HQ_Lasso_coordinates\\.csv$"), 
                           full.names = TRUE)
  
  if(length(lasso_file) == 0) {
    stop(sprintf("Could not find lasso file for slide %d", target_slide))
  }
  
  # Load lasso regions
  lasso_list <- load_lasso_regions(lasso_file[1])
  region_idx <- which(names(lasso_list) == target_region)
  
  if(length(region_idx) == 0) {
    stop(sprintf("Could not find region %s", target_region))
  }
  
  # Get the region geometry and create window
  current_region <- lasso_list[[region_idx]]
  region_window <- sf_to_owin(current_region)
  
  # Filter points for this region
  region_data <- filter_points_in_window(
    rgc_df %>% filter(slide == target_slide), 
    region_window
  )
  
  # Get target and background cells
  target_cells <- region_data %>% filter(Prediction == target_subtype)
  other_cells <- region_data %>% filter(Prediction != target_subtype)
  
  # Create target point pattern and calculate observed stats
  target_ppp <- ppp(
    x = target_cells$x,
    y = target_cells$y,
    window = region_window,
    check = TRUE
  )
  obs_stats <- calculate_spatial_stats(target_ppp)
  
  # Get stats annotation for observed data
  obs_stats_label <- data.frame(
    x = region_window$xrange[1] + diff(region_window$xrange) * 0.05,
    y = region_window$yrange[2] - diff(region_window$yrange) * 0.1,
    label = sprintf("NNRI: %.3f\nVDRI: %.3f", obs_stats$nnri, obs_stats$vdri)
  )
  
  # Create Voronoi tessellation for target cells
  voro <- dirichlet(target_ppp)
  voro_tiles <- tiles(voro)
  voro_polygons <- lapply(voro_tiles, function(tile) {
    coords <- cbind(tile$bdry[[1]]$x, tile$bdry[[1]]$y)
    coords <- rbind(coords, coords[1,]) # Close polygon
    st_polygon(list(coords))
  })
  voro_sf <- st_sf(geometry = st_sfc(voro_polygons))
  
  # Perform one bootstrap sample (fresh sample each time)
  sampled_indices <- sample(nrow(region_data), nrow(target_cells))
  sampled_points <- region_data[sampled_indices, ]
  non_sampled_points <- region_data[-sampled_indices, ]
  
  # Calculate bootstrap stats
  boot_ppp <- ppp(
    x = sampled_points$x,
    y = sampled_points$y,
    window = region_window,
    check = TRUE
  )
  boot_stats <- calculate_spatial_stats(boot_ppp)
  
  # Get stats annotation for bootstrap
  boot_stats_label <- data.frame(
    x = region_window$xrange[1] + diff(region_window$xrange) * 0.05,
    y = region_window$yrange[2] - diff(region_window$yrange) * 0.1,
    label = sprintf("NNRI: %.3f\nVDRI: %.3f", boot_stats$nnri, boot_stats$vdri)
  )
  
  # Convert region to sf for plotting
  region_sf <- st_sf(geometry = st_sfc(current_region))
  
  # Create plot comparing observed and bootstrap
  p_boot <- ggplot() +
    geom_sf(data = region_sf, fill = NA, color = "black") +
    geom_point(data = non_sampled_points, 
               aes(x = x, y = y), 
               color = "grey70", size = 0.5, alpha = 0.5) +
    geom_point(data = sampled_points, 
               aes(x = x, y = y), 
               color = "blue", size = 1) +
    geom_text(data = boot_stats_label,
              aes(x = x, y = y, label = label),
              hjust = 0, vjust = 1, size = 3) +
    coord_sf() +
    theme_minimal() +
    theme(
      legend.position = 'None',
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank()
    ) +
    labs(title = sprintf("Bootstrap Sample: %s", target_subtype))
  
  # Create point distribution plot
  p_points <- ggplot() +
    geom_sf(data = region_sf, fill = NA, color = "black") +
    geom_point(data = other_cells, aes(x = x, y = y), 
               color = "grey70", size = 0.5, alpha = 0.5) +
    geom_point(data = target_cells, aes(x = x, y = y), 
               color = "red", size = 1) +
    geom_text(data = obs_stats_label,
              aes(x = x, y = y, label = label),
              hjust = 0, vjust = 1, size = 3) +
    coord_sf() +
    theme_minimal() +
    theme(
      legend.position = 'None',
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank()
    ) +
    labs(title = sprintf("Observed Data: %s", target_subtype))
  
  # Create Voronoi plot
  p_voronoi <- ggplot() +
    geom_sf(data = region_sf, fill = NA, color = "black") +
    geom_sf(data = voro_sf, fill = NA, color = "red", alpha = 0.5) +
    geom_point(data = other_cells, aes(x = x, y = y), 
               color = "grey70", size = 0.5, alpha = 0.5) +
    geom_point(data = target_cells, aes(x = x, y = y), 
               color = "red", size = 1) +
    coord_sf() +
    theme_minimal() +
    theme(
      legend.position = 'None',
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank()
    ) +
    labs(title = sprintf("Voronoi Tessellation: %s", target_subtype))
  
  # Print statistics
  cat(sprintf("\nAnalysis for %s (Slide %d, Region %s):\n", 
              target_subtype, target_slide, target_region))
  cat("\nObserved Statistics:\n")
  cat(sprintf("NNRI: %.3f\n", obs_stats$nnri))
  cat(sprintf("VDRI: %.3f\n", obs_stats$vdri))
  cat("\nBootstrap Sample Statistics:\n")
  cat(sprintf("NNRI: %.3f\n", boot_stats$nnri))
  cat(sprintf("VDRI: %.3f\n", boot_stats$vdri))
  
  # If results_df is provided, print comparison
  if (!is.null(results_df)) {
    stored_results <- results_df %>%
      filter(
        slide == target_slide,
        region == target_region,
        Subtype == target_subtype
      )
    
    if (nrow(stored_results) > 0) {
      cat("\nStored Results Comparison:\n")
      stored_obs <- stored_results %>% filter(Observed)
      stored_null <- stored_results %>% filter(!Observed)
      cat("\nStored Observed Values:\n")
      cat(sprintf("NNRI: %.3f\n", stored_obs$NNRI[1]))
      cat(sprintf("VDRI: %.3f\n", stored_obs$VDRI[1]))
      cat("\nStored Null Distribution Summary:\n")
      cat(sprintf("NNRI mean: %.3f (sd: %.3f)\n", 
                  mean(stored_null$NNRI), sd(stored_null$NNRI)))
      cat(sprintf("VDRI mean: %.3f (sd: %.3f)\n", 
                  mean(stored_null$VDRI), sd(stored_null$VDRI)))
    }
  }
  
  # Return all plots
  return(list(
    points = p_points, 
    voronoi = p_voronoi, 
    bootstrap = p_boot
  ))
}



visualize_edge_correction <- function(results_df, rgc_df, 
                                      target_slide, target_region, target_subtype,
                                      lasso_root = '/media/sam/New Volume/Xenium_Data/HQ_NearestNeighbor_Zones') {
  # Load the lasso region (reuse previous code)
  lasso_file <- list.files(lasso_root, 
                           pattern = paste0(target_slide, "_HQ_Lasso_coordinates\\.csv$"), 
                           full.names = TRUE)
  
  if(length(lasso_file) == 0) {
    stop(sprintf("Could not find lasso file for slide %d", target_slide))
  }
  
  # Load lasso regions
  lasso_list <- load_lasso_regions(lasso_file[1])
  region_idx <- which(names(lasso_list) == target_region)
  
  if(length(region_idx) == 0) {
    stop(sprintf("Could not find region %s", target_region))
  }
  
  # Get the region geometry and create window
  current_region <- lasso_list[[region_idx]]
  region_window <- sf_to_owin(current_region)
  
  # Filter points for this region
  region_data <- filter_points_in_window(
    rgc_df %>% filter(slide == target_slide), 
    region_window
  )
  
  # Get target cells and create ppp
  target_cells <- region_data %>% filter(Prediction == target_subtype)
  target_ppp <- ppp(
    x = target_cells$x,
    y = target_cells$y,
    window = region_window,
    check = TRUE
  )
  
  # Create Voronoi tessellation
  voro <- dirichlet(target_ppp)
  
  # Calculate buffer distance (10% of sqrt of area)
  buffer_distance <- sqrt(area(region_window)) * 0.1
  
  # Identify edge and interior points
  edge_distances <- bdist.points(target_ppp)
  interior_points <- edge_distances > buffer_distance
  edge_points <- !interior_points
  
  # Create data frame for plotting
  plot_data <- data.frame(
    x = target_ppp$x,
    y = target_ppp$y,
    edge_distance = edge_distances,
    is_edge = edge_points,
    nn_dist = nndist(target_ppp, correction="border")
  )
  
  # Convert Voronoi tessellation to sf objects
  voro_tiles <- tiles(voro)
  voro_polygons <- lapply(seq_along(voro_tiles), function(i) {
    tile <- voro_tiles[[i]]
    coords <- cbind(tile$bdry[[1]]$x, tile$bdry[[1]]$y)
    coords <- rbind(coords, coords[1,]) # Close polygon
    st_polygon(list(coords))
  })
  voro_sf <- st_sf(
    geometry = st_sfc(voro_polygons),
    is_edge = edge_points
  )
  
  # Create region sf object
  region_sf <- st_sf(geometry = st_sfc(current_region))
  
  # Create buffer visualization
  buffer_region <- st_buffer(region_sf, dist = -buffer_distance)
  
  # Create edge correction visualization plot
  p <- ggplot() +
    # Plot region boundary
    geom_sf(data = region_sf, fill = NA, color = "black") +
    # Plot buffer zone
    geom_sf(data = buffer_region, fill = NA, color = "blue", linetype = "dashed") +
    # Plot Voronoi cells with edge status
    geom_sf(data = voro_sf, aes(fill = is_edge), 
            alpha = 0.3, color = "red") +
    # Plot points with edge status
    geom_point(data = plot_data, 
               aes(x = x, y = y, color = is_edge, size = nn_dist)) +
    # Add other cells in background
    geom_point(data = region_data %>% filter(Prediction != target_subtype), 
               aes(x = x, y = y), color = "grey70", size = 0.5, alpha = 0.3) +
    # Customize appearance
    scale_fill_manual(values = c("FALSE" = "green", "TRUE" = "red"),
                      labels = c("Interior", "Edge"),
                      name = "Cell Location") +
    scale_color_manual(values = c("FALSE" = "darkgreen", "TRUE" = "darkred"),
                       labels = c("Interior", "Edge"),
                       name = "Cell Location") +
    scale_size_continuous(name = "NN Distance",
                          guide = guide_legend(override.aes = list(color = "black"))) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank()
    ) +
    labs(title = sprintf("Edge Correction Visualization: %s", target_subtype),
         subtitle = paste("Blue dashed line: Buffer zone boundary\n",
                          "Red polygons: Voronoi cells\n",
                          "Point size: Edge-corrected nearest neighbor distance"))
  
  # Print summary statistics
  cat(sprintf("\nEdge Correction Analysis for %s (Slide %d, Region %s):\n", 
              target_subtype, target_slide, target_region))
  cat(sprintf("Total points: %d\n", npoints(target_ppp)))
  cat(sprintf("Edge points: %d (%.1f%%)\n", 
              sum(edge_points), 100*mean(edge_points)))
  cat(sprintf("Interior points: %d (%.1f%%)\n", 
              sum(interior_points), 100*mean(interior_points)))
  cat(sprintf("Buffer distance: %.2f\n", buffer_distance))
  
  return(p)
}
#################################################################################
########################### Analysis#############################################
#################################################################################

re_run_all <- F

cat("\nStarting analysis pipeline\n")
root <- '/media/sam/Data2/baysor_rbpms_consolidated'
viz_root <- "/home/sam/MappingAllMurineRGCs/Data/Figure3_Outputs"
lasso_root <- '/media/sam/New Volume/Xenium_Data/HQ_NearestNeighbor_Zones'

# Load RGC data
cat("Loading RGC data...")
rgc_df <- read_csv(paste0(root,'/OldModel/merged_rgc_prediction_expmat.csv'))
cat(" done\n")


if (re_run_all){
  analysis_output <- run_full_analysis(viz_root, all_subtypes = T, subtype_only = '01_W3D1.1')
  test_results <- analysis_output$results
  distance_data <- analysis_output$distances
} else{
  test_results <- read_csv('/media/sam/Data2/baysor_rbpms_consolidated/NN_visualizations/full_spatial_analysis.csv')
  distance_data <- read_csv('/media/sam/Data2/baysor_rbpms_consolidated/NN_visualizations/neighbor_distances.csv')
}
test_results <- test_results %>%
  group_by(Subtype, slide, region, Observed) %>%
  summarize(NNRI = mean(NNRI),
            VDRI = mean(VDRI),
            N = mean(N),
            N_region = mean(N_region))


test_results %>%
  filter(N > 10) %>%
  ggplot(aes(x = NNRI, y = reorder(Subtype, N), color = Observed)) +
  geom_boxplot() +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(color = "black")
  ) +
  labs(title = "Spatial Analysis Results",
       subtitle = "NNRI Distribution by Subtype") #+
  # facet_grid(slide ~ region)












create_statistical_summary <- function(results, viz_root, min_cells = 10, min_observations = 3) {
  # Create visualization directory
  viz_dir <- file.path(viz_root, "NN_visualizations")
  dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)
  
  # First filter for minimum cells before aggregation
  filtered_results <- results %>% filter(N >= min_cells)
  
  # Create summary statistics
  test_results <- filtered_results %>%
    group_by(Subtype, slide, region, Observed) %>%
    summarize(
      NNRI = mean(NNRI),
      VDRI = mean(VDRI),
      N = mean(N),
      N_region = mean(N_region),
      .groups = 'drop'
    )
  
  # Sort subtypes alphabetically
  subtypes <- sort(unique(test_results$Subtype), decreasing = TRUE)
  stat_results <- list()
  diagnostic_info <- list()
  
  cat("\nBeginning statistical analysis for each subtype:\n")
  cat("------------------------------------------------\n")
  
  for(st in subtypes) {
    cat(sprintf("\nAnalyzing subtype: %s\n", st))
    subtype_data <- test_results %>% filter(Subtype == st)
    
    # Get observed and null values
    observed <- subtype_data %>% filter(Observed)
    null <- subtype_data %>% filter(!Observed)
    
    # Initialize diagnostic information
    diagnostic <- list(
      subtype = st,
      total_samples = nrow(subtype_data),
      observed_samples = nrow(observed),
      null_samples = nrow(null),
      observed_variance_nnri = if(nrow(observed) > 0) var(observed$NNRI) else NA,
      null_variance_nnri = if(nrow(null) > 0) var(null$NNRI) else NA,
      observed_variance_vdri = if(nrow(observed) > 0) var(observed$VDRI) else NA,
      null_variance_vdri = if(nrow(null) > 0) var(null$VDRI) else NA,
      status_nnri = "Not tested yet",
      status_vdri = "Not tested yet"
    )
    
    # Print diagnostic information
    cat(sprintf("  Total samples: %d\n", diagnostic$total_samples))
    cat(sprintf("  Observed samples: %d\n", diagnostic$observed_samples))
    cat(sprintf("  Null samples: %d\n", diagnostic$null_samples))
    cat(sprintf("  Observed NNRI variance: %s\n", 
                ifelse(is.na(diagnostic$observed_variance_nnri), 
                       "NA", sprintf("%.6f", diagnostic$observed_variance_nnri))))
    cat(sprintf("  Null NNRI variance: %s\n", 
                ifelse(is.na(diagnostic$null_variance_nnri), 
                       "NA", sprintf("%.6f", diagnostic$null_variance_nnri))))
    
    # Initialize results with NA
    t_test_result <- list(
      Subtype = st,
      mean_N = mean(observed$N),
      n_observed = nrow(observed),
      n_null = nrow(null),
      mean_observed_nnri = mean(observed$NNRI),
      mean_null_nnri = mean(null$NNRI),
      p_value_nnri = NA,
      mean_observed_vdri = mean(observed$VDRI),
      mean_null_vdri = mean(null$VDRI),
      p_value_vdri = NA
    )
    
    # Function to perform t-test with diagnostics
    perform_t_test <- function(observed_vals, null_vals, metric_name) {
      if(length(observed_vals) < min_observations) {
        return(list(p_value = NA, status = sprintf("Insufficient observed samples (n=%d < %d)", 
                                                  length(observed_vals), min_observations)))
      }
      if(length(null_vals) < min_observations) {
        return(list(p_value = NA, status = sprintf("Insufficient null samples (n=%d < %d)", 
                                                  length(null_vals), min_observations)))
      }
      if(var(observed_vals) == 0) {
        return(list(p_value = NA, status = "Zero variance in observed data"))
      }
      if(var(null_vals) == 0) {
        return(list(p_value = NA, status = "Zero variance in null data"))
      }
      
      tryCatch({
        test_result <- t.test(observed_vals, null_vals)
        return(list(p_value = test_result$p.value, status = "Success"))
      }, error = function(e) {
        return(list(p_value = NA, status = paste("T-test error:", e$message)))
      })
    }
    
    # Perform NNRI test
    nnri_test <- perform_t_test(observed$NNRI, null$NNRI, "NNRI")
    t_test_result$p_value_nnri <- nnri_test$p_value
    diagnostic$status_nnri <- nnri_test$status
    
    # Perform VDRI test
    vdri_test <- perform_t_test(observed$VDRI, null$VDRI, "VDRI")
    t_test_result$p_value_vdri <- vdri_test$p_value
    diagnostic$status_vdri <- vdri_test$status
    
    # Print test results
    cat(sprintf("  NNRI test status: %s\n", diagnostic$status_nnri))
    cat(sprintf("  VDRI test status: %s\n", diagnostic$status_vdri))
    
    # Add significance stars
    t_test_result$significance_nnri <- get_significance_stars(t_test_result$p_value_nnri)
    t_test_result$significance_vdri <- get_significance_stars(t_test_result$p_value_vdri)
    
    # Create labels with sample sizes
    t_test_result$label_nnri <- sprintf("%s (n=%d)", 
                                       t_test_result$significance_nnri, 
                                       t_test_result$n_observed)
    t_test_result$label_vdri <- sprintf("%s (n=%d)", 
                                       t_test_result$significance_vdri, 
                                       t_test_result$n_observed)
    
    # Store results
    stat_results[[st]] <- as.data.frame(t_test_result)
    diagnostic_info[[st]] <- diagnostic
  }
  
  # Combine statistical results
  stat_summary <- bind_rows(stat_results)
  
  # Create plots
  create_metric_plot <- function(data, stat_data, metric, label_col) {
    ggplot(data %>% mutate(Subtype = factor(Subtype, levels = subtypes)), 
           aes(x = .data[[metric]], y = Subtype)) +
      geom_boxplot(aes(color = Observed, fill = Observed), 
                   outlier.size = 0.4,
                   alpha = 0.5) +
      scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey80")) +
      scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "grey80")) +
      scale_x_continuous(expand = expansion(mult = c(0.05, 0.2))) +
      geom_text(data = stat_data,
                aes(x = max(data[[metric]]), 
                    label = .data[[label_col]]),
                hjust = 0,
                size = 3,
                vjust = 0.5) +
      theme_minimal() +
      theme(
        legend.position = 'None',
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(color = "black")
      ) +
      labs(title = "Spatial Analysis Results",
           subtitle = sprintf("%s Distribution by Subtype (min cells = %d, min obs = %d)", 
                            metric, min_cells, min_observations),
           caption = "* p<0.05, ** p<0.01, *** p<0.001")
  }
  
  # Create and save plots
  p_nnri <- create_metric_plot(test_results, stat_summary, "NNRI", "label_nnri")
  p_vdri <- create_metric_plot(test_results, stat_summary, "VDRI", "label_vdri")
  
  ggsave(file.path(viz_dir, "statistical_summary_NNRI.png"), p_nnri, 
         width = 12, height = 15)
  ggsave(file.path(viz_dir, "statistical_summary_VDRI.png"), p_vdri, 
         width = 12, height = 15)
  
  # Save statistical summary
  write_csv(stat_summary, file.path(viz_dir, "statistical_summary.csv"))
  
  # Save diagnostic information
  diagnostic_df <- bind_rows(diagnostic_info)
  write_csv(diagnostic_df, file.path(viz_dir, "statistical_diagnostic.csv"))
  
  return(list(
    summary = stat_summary,
    diagnostics = diagnostic_df,
    plot_nnri = p_nnri,
    plot_vdri = p_vdri
  ))
}



stat_results <- create_statistical_summary(drop_na(test_results), viz_root, min_cells = 5, min_observations = 3)



























###########################################################################################
# Create visualization directory if it doesn't exist
viz_dir <- file.path(viz_root, "NN_visualizations")
dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)

# Create combined region identifier if it doesn't exist
distance_data <- distance_data %>%
  mutate(region_id = paste(slide, region, sep="_"))

# Create histogram for each subtype
subtypes <- unique(distance_data$Subtype)
for(subtype in subtypes) {
  subtype_data <- distance_data %>% filter(Subtype == subtype)
  
  p <- ggplot(subtype_data, aes(x = distance, fill = region_id)) +
    geom_histogram(position = "stack", bins = 50) +
    scale_fill_viridis_d() +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank()
    ) +
    labs(title = paste("Neighbor Distance Distribution:", subtype),
         x = "Distance",
         y = "Count")
  
  ggsave(file.path(viz_dir, sprintf("neighbor_dist_%s.png", subtype)), 
         p, width = 10, height = 6)
  
  # Print progress
  cat(sprintf("Created histogram for subtype: %s\n", subtype))
}


#############################################################################################################
################################## Sanity checking visualizations ###########################################
#############################################################################################################
plot_slide <- 18432
plot_region <- "HQ_Lasso_4"
plot_type<- "01_W3D1.1"#"43_Alpha_ONS"
for (i in 1:2){
  
   plots <- inspect_region_subtype(
    rgc_df = rgc_df,
    target_slide = plot_slide,
    target_region = plot_region,
    target_subtype = plot_type, #"12_ooDS_NT", #   "02_W3D1.2",
    results_df = test_results 
  )
  gridExtra::grid.arrange(plots$points, plots$voronoi, plots$bootstrap, ncol = 3)
  Sys.sleep(2)
}
edge_plot <- visualize_edge_correction(
  results_df = results,
  rgc_df = rgc_df,
  target_slide = plot_slide,  # example slide ID
  target_region = plot_region,  # example region name
  target_subtype = plot_type  # example subtype
)
print(edge_plot)




visualize_edge_correction_publication <- function(results_df, rgc_df, 
                                      target_slide, target_region, target_subtype,
                                      lasso_root = '/media/sam/New Volume/Xenium_Data/HQ_NearestNeighbor_Zones',
                                      plot_label = "A") {
  # Load the lasso region (reuse previous code)
  lasso_file <- list.files(lasso_root, 
                           pattern = paste0(target_slide, "_HQ_Lasso_coordinates\\.csv$"), 
                           full.names = TRUE)
  
  if(length(lasso_file) == 0) {
    stop(sprintf("Could not find lasso file for slide %d", target_slide))
  }
  
  # Load lasso regions
  lasso_list <- load_lasso_regions(lasso_file[1])
  region_idx <- which(names(lasso_list) == target_region)
  
  if(length(region_idx) == 0) {
    stop(sprintf("Could not find region %s", target_region))
  }
  
  # Get the region geometry and create window
  current_region <- lasso_list[[region_idx]]
  region_window <- sf_to_owin(current_region)
  
  # Filter points for this region
  region_data <- filter_points_in_window(
    rgc_df %>% filter(slide == target_slide), 
    region_window
  )
  
  # Get target cells and create ppp
  target_cells <- region_data %>% filter(Prediction == target_subtype)
  target_ppp <- ppp(
    x = target_cells$x,
    y = target_cells$y,
    window = region_window,
    check = TRUE
  )
  
  # Create Voronoi tessellation
  voro <- dirichlet(target_ppp)
  
  # Calculate buffer distance (10% of sqrt of area)
  buffer_distance <- sqrt(area(region_window)) * 0.1
  
  # Identify edge and interior points
  edge_distances <- bdist.points(target_ppp)
  interior_points <- edge_distances > buffer_distance
  edge_points <- !interior_points
  
  # Create data frame for plotting
  plot_data <- data.frame(
    x = target_ppp$x,
    y = target_ppp$y,
    edge_distance = edge_distances,
    is_edge = edge_points,
    nn_dist = nndist(target_ppp, correction="border")
  )
  
  # Convert Voronoi tessellation to sf objects
  voro_tiles <- tiles(voro)
  voro_polygons <- lapply(seq_along(voro_tiles), function(i) {
    tile <- voro_tiles[[i]]
    coords <- cbind(tile$bdry[[1]]$x, tile$bdry[[1]]$y)
    coords <- rbind(coords, coords[1,]) # Close polygon
    st_polygon(list(coords))
  })
  voro_sf <- st_sf(
    geometry = st_sfc(voro_polygons),
    is_edge = edge_points
  )
  
  # Create region sf object
  region_sf <- st_sf(geometry = st_sfc(current_region))
  
  # Create buffer visualization
  buffer_region <- st_buffer(region_sf, dist = -buffer_distance)
  
  # Create edge correction visualization plot
  p <- ggplot() +
    geom_sf(data = region_sf, fill = NA, color = "black") +
    geom_sf(data = buffer_region, fill = NA, color = "blue", linetype = "dashed") +
    geom_sf(data = voro_sf, aes(fill = is_edge), 
            alpha = 0.3, color = "red") +
    geom_point(data = plot_data %>% filter(!is_edge), 
               aes(x = x, y = y, size = nn_dist),
               color = "darkgreen", fill = "darkgreen", shape = 21) +
    geom_point(data = plot_data %>% filter(is_edge), 
               aes(x = x, y = y, size = nn_dist),
               color = "darkgreen", fill = "white", shape = 21) +
    geom_point(data = region_data %>% filter(Prediction != target_subtype), 
               aes(x = x, y = y), color = "grey70", size = 0.5, alpha = 0.3) +
    scale_fill_manual(values = c("FALSE" = "green", "TRUE" = "red"),
                      labels = c("Interior", "Edge"),
                      name = "Cell Location") +
    scale_size_continuous(name = "NN Distance",
                          guide = guide_legend(override.aes = list(color = "black"))) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0, face = "bold", size = 19)
    ) +
    labs(title = plot_label)
  
  return(p)
  
  # Print summary statistics
  cat(sprintf("\nEdge Correction Analysis for %s (Slide %d, Region %s):\n", 
              target_subtype, target_slide, target_region))
  cat(sprintf("Total points: %d\n", npoints(target_ppp)))
  cat(sprintf("Edge points: %d (%.1f%%)\n", 
              sum(edge_points), 100*mean(edge_points)))
  cat(sprintf("Interior points: %d (%.1f%%)\n", 
              sum(interior_points), 100*mean(interior_points)))
  cat(sprintf("Buffer distance: %.2f\n", buffer_distance))
  
  return(p)
}

visualize_edge_correction_comparison <- function(rgc_df, 
                                                 target_slide, target_region, target_subtype,
                                                 lasso_root = '/media/sam/New Volume/Xenium_Data/HQ_NearestNeighbor_Zones') {
  # Load the lasso region
  lasso_file <- list.files(lasso_root, 
                           pattern = paste0(target_slide, "_HQ_Lasso_coordinates\\.csv$"), 
                           full.names = TRUE)
  
  if(length(lasso_file) == 0) {
    stop(sprintf("Could not find lasso file for slide %d", target_slide))
  }
  
  # Load lasso regions
  lasso_list <- load_lasso_regions(lasso_file[1])
  region_idx <- which(names(lasso_list) == target_region)
  
  if(length(region_idx) == 0) {
    stop(sprintf("Could not find region %s", target_region))
  }
  
  # Get the region geometry and create window
  current_region <- lasso_list[[region_idx]]
  region_window <- sf_to_owin(current_region)
  
  # Filter points for this region
  region_data <- filter_points_in_window(
    rgc_df %>% filter(slide == target_slide), 
    region_window
  )
  
  # Create bootstrap data by shuffling labels
  shuffled_data <- region_data
  shuffled_data$Prediction <- sample(shuffled_data$Prediction)
  
  # Plot A: Original data
  p1 <- visualize_edge_correction_publication(
    results_df = NULL,
    rgc_df = rgc_df,
    target_slide = target_slide,
    target_region = target_region,
    target_subtype = target_subtype,
    lasso_root = lasso_root,
    plot_label = "a"
  )
  
  # Plot B: Bootstrapped data
  p2 <- visualize_edge_correction_publication(
    results_df = NULL,
    rgc_df = shuffled_data,
    target_slide = target_slide,
    target_region = target_region,
    target_subtype = target_subtype,
    lasso_root = lasso_root,
    plot_label = "b"
  )
  
  # Return both plots
  return(list(observed = p1, bootstrap = p2))
}

plots <- visualize_edge_correction_comparison(
  rgc_df = rgc_df,
  target_slide = plot_slide,
  target_region = plot_region,
  target_subtype = plot_type
)

# Display plots side by side
(combined_plot <- gridExtra::grid.arrange(plots$observed, plots$bootstrap, ncol = 2))
ggsave(filename ="/home/sam/MappingAllMurineRGCs/Data/Figure3_Outputs/Sup8_localstats.png", combined_plot, width = 12, height = 5)

