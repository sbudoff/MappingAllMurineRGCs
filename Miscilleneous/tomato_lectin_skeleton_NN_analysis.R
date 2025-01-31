library(tidyverse)
library(lme4)
library(cowplot)
library(sf)
library(spatstat)
library(stats)  # For ANOVA and Tukey's test


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

compute_vascular_distances <- function(rgc_data, vasculature_df, slide_id) {
  # Initialize empty list to store results
  results_list <- list()
  
  # Get unique subtypes
  subtypes <- unique(rgc_data$Prediction)
  
  # Split vasculature data into small and large
  small_vasc <- vasculature_df %>% filter(small == TRUE) %>% select(x, y)
  large_vasc <- vasculature_df %>% filter(small == FALSE) %>% select(x, y)
  
  # Create matrices for faster computation
  small_vasc_matrix <- as.matrix(small_vasc)
  large_vasc_matrix <- as.matrix(large_vasc)
  
  # Loop through each subtype
  for(subtype in subtypes) {
    # Filter RGC data for current subtype
    subtype_data <- rgc_data %>% 
      filter(Prediction == subtype) %>%
      select(x, y)
    
    if(nrow(subtype_data) > 0) {
      # Convert to matrix for faster computation
      rgc_matrix <- as.matrix(subtype_data)
      
      # Compute distances to small vasculature
      small_dist_matrix <- fields::rdist(rgc_matrix, small_vasc_matrix)
      small_nn <- apply(small_dist_matrix, 1, min)
      
      # Compute distances to large vasculature
      large_dist_matrix <- fields::rdist(rgc_matrix, large_vasc_matrix)
      large_nn <- apply(large_dist_matrix, 1, min)
      
      # Create results dataframe for this subtype
      results_df <- data.frame(
        slide = slide_id,
        Subtype = subtype,
        X = subtype_data$x,
        Y = subtype_data$y,
        small_vasc_NN = small_nn,
        large_vasc_NN = large_nn
      )
      
      results_list[[subtype]] <- results_df
    }
  }
  
  # Combine all results
  final_results <- do.call(rbind, results_list)
  return(final_results)
}

plot_region_with_vasculature <- function(lasso_region, rgc_data, vasculature_df, region_name) {
  # Convert lasso region to sf object and owin
  lasso_sf <- st_sf(geometry = st_sfc(lasso_region))
  lasso_window <- sf_to_owin(lasso_region)
  
  # Filter RGC points to only those within the lasso region
  filtered_points <- filter_points_in_window(rgc_data, lasso_window)
  
  # Get bounding box of lasso region to filter vasculature
  bbox <- st_bbox(lasso_sf)
  
  # Filter vasculature to bounding box region and split by size
  filtered_vasculature <- vasculature_df %>%
    filter(x >= bbox[1], x <= bbox[3],
           y >= bbox[2], y <= bbox[4])
  
  small_vasc <- filtered_vasculature %>% 
    filter(small == TRUE) %>%
    group_by(group)
  
  large_vasc <- filtered_vasculature %>% 
    filter(small == FALSE) %>%
    group_by(group)
  
  # Create base plot
  p <- ggplot() +
    # Add lasso region outline
    geom_sf(data = lasso_sf, fill = NA, color = "black", linewidth = 0.5) +
    
    # Add vasculature lines
    geom_line(data = small_vasc, aes(x = x, y = y, group = group), 
              color = "red", linewidth = 0.5, alpha = 0.7) +
    geom_line(data = large_vasc, aes(x = x, y = y, group = group), 
              color = "darkred", linewidth = 1.2, alpha = 0.7) +
    
    # Add RGC points
    geom_point(data = filtered_points, 
               aes(x = x, y = y, color = Prediction), 
               size = 0.8, alpha = 0.6) +
    
    # Styling
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "none",  # Remove legend
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10)
    ) +
    labs(title = paste("Region:", region_name),
         subtitle = sprintf("RGC points: %d | Small vessels: %d | Large vessels: %d", 
                            nrow(filtered_points), 
                            n_distinct(small_vasc$group),
                            n_distinct(large_vasc$group))) +
    coord_sf()
  
  return(p)
}


# Function to perform pairwise comparisons
perform_pairwise <- function(vessel_type_data) {
  subtypes <- unique(vessel_type_data$Subtype)
  results <- list()
  
  for(i in 1:(length(subtypes)-1)) {
    for(j in (i+1):length(subtypes)) {
      st1 <- subtypes[i]
      st2 <- subtypes[j]
      
      group1_data <- vessel_type_data %>% filter(Subtype == st1) %>% pull(distance)
      group2_data <- vessel_type_data %>% filter(Subtype == st2) %>% pull(distance)
      
      if(length(group1_data) >= 3 && length(group2_data) >= 3) {
        test <- wilcox.test(group1_data, group2_data)
        
        results[[length(results) + 1]] <- data.frame(
          comparison_type = "between",
          group1 = st1,
          group2 = st2,
          vessel_type = unique(vessel_type_data$vessel_type),
          p.value = test$p.value,
          n_cells = min(length(group1_data), length(group2_data))
        )
      }
    }
  }
  return(do.call(rbind, results))
}

analyze_vascular_distances <- function(final_results) {
  require(stats)
  require(dplyr)
  require(tidyr)
  
  # Reshape data for analysis
  long_data <- final_results %>%
    pivot_longer(
      cols = c(small_vasc_NN, large_vasc_NN),
      names_to = "vessel_type",
      values_to = "distance"
    ) %>%
    mutate(
      vessel_type = factor(vessel_type, 
                           levels = c("small_vasc_NN", "large_vasc_NN"),
                           labels = c("Small", "Large")),
      cell_id = paste(slide, X, Y, sep = "_")
    )
  
  # 1. Within-group comparisons (large vs small for each subtype)
  subtypes <- unique(long_data$Subtype)
  within_group_tests <- list()
  
  for(st in subtypes) {
    subtype_data <- long_data %>% 
      filter(Subtype == st) %>%
      arrange(cell_id, vessel_type)
    
    if(n_distinct(subtype_data$cell_id) >= 3) {
      small_vals <- subtype_data %>% 
        filter(vessel_type == "Small") %>% 
        pull(distance)
      large_vals <- subtype_data %>% 
        filter(vessel_type == "Large") %>% 
        pull(distance)
      
      if(length(small_vals) == length(large_vals) && length(small_vals) > 0) {
        test <- wilcox.test(small_vals, large_vals, paired = TRUE)
        
        within_group_tests[[st]] <- data.frame(
          comparison_type = "within",
          group1 = st,
          group2 = st,
          vessel_comparison = "Large-Small",
          p.value = test$p.value,
          n_cells = length(small_vals)
        )
      }
    }
  }
  
  # 2. Between-group comparisons for each vessel type
  between_group_tests <- list()
  
  
  
  # Perform between-group comparisons for each vessel type
  for(vtype in c("Small", "Large")) {
    vessel_data <- long_data %>% filter(vessel_type == vtype)
    between_group_tests[[vtype]] <- perform_pairwise(vessel_data)
  }
  
  # Combine all results
  within_group_df <- do.call(rbind, within_group_tests)
  between_group_df <- do.call(rbind, between_group_tests)
  
  # Adjust p-values separately for within and between group comparisons
  if(nrow(within_group_df) > 0) {
    within_group_df$adj.p.value <- p.adjust(within_group_df$p.value, method = "BH")
  }
  if(nrow(between_group_df) > 0) {
    between_group_df$adj.p.value <- p.adjust(between_group_df$p.value, method = "BH")
  }
  
  # Calculate summary statistics
  summary_stats <- long_data %>%
    group_by(Subtype, vessel_type) %>%
    summarise(
      mean_dist = mean(distance),
      sd_dist = sd(distance),
      n = n(),
      .groups = "drop"
    )
  
  return(list(
    within_group = within_group_df,
    between_group = between_group_df,
    summary = summary_stats
  ))
}

plot_vascular_distances_with_heatmap <- function(final_results, stat_results) {
  require(ggplot2)
  require(patchwork)
  require(cowplot)
  
  # Reshape data for plotting
  long_data <- final_results %>%
    pivot_longer(
      cols = c(small_vasc_NN, large_vasc_NN),
      names_to = "vessel_type",
      values_to = "distance"
    ) %>%
    mutate(
      vessel_type = factor(vessel_type, 
                           levels = c("small_vasc_NN", "large_vasc_NN"),
                           labels = c("Small", "Large")),
      subtype_num = as.numeric(gsub("([0-9]+).*", "\\1", Subtype))
    ) %>%
    arrange(desc(subtype_num))
  
  # Create ordered factor for consistent ordering across plots
  subtype_levels <- unique(long_data$Subtype)
  long_data$Subtype <- factor(long_data$Subtype, levels = subtype_levels)
  
  # Create boxplot
  boxplot <- ggplot(long_data, aes(x = distance, y = Subtype)) +
    geom_boxplot(aes(fill = vessel_type), 
                 position = position_dodge(width = 0.8),
                 outlier.shape = NA,
                 width = 0.7) +
    scale_fill_manual(values = c("Small" = "red", "Large" = "darkred")) +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "top",
      legend.title = element_text(size = 10),
      axis.text.y = element_text(size = 8)
    ) +
    labs(
      x = "Distance to Vessel",
      y = NULL,
      fill = "Vessel Type"
    )
  
  # Prepare comparison data 
  all_comparisons <- bind_rows(
    stat_results$between_group %>%
      mutate(comparison_type = "between"),
    stat_results$within_group %>%
      mutate(
        group2 = group1,
        comparison_type = "within"
      )
  ) %>%
    mutate(
      log_p = -log10(adj.p.value),
      log_p = pmin(log_p, 3),  # Cap at p=0.001
      significance = case_when(
        adj.p.value < 0.001 ~ "***",
        adj.p.value < 0.01 ~ "**",
        adj.p.value < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      group1 = factor(group1, levels = subtype_levels),
      group2 = factor(group2, levels = subtype_levels)
    )
  
  # Create heatmap
  heatmap <- ggplot(all_comparisons, aes(x = group2, y = group1)) +
    geom_tile(aes(fill = log_p)) +
    geom_text(aes(label = ifelse(group1 == group2, significance, "")), 
              size = 3, color = "black") +
    scale_fill_gradientn(
      colors = c("white", "#fee090", "#fc8d59", "#d73027"),  # ColorBrewer YlOrRd palette
      values = scales::rescale(c(0, 1.301, 2, 3)),
      limits = c(0, 3),
      na.value = "white",
      name = "-log10(adj.p)     ",  # Added spacing after name
      breaks = c(0, 1.301, 2, 3),
      labels = c("ns    ", "0.05    ", "0.01    ", "0.001"),  # Added spacing between values
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        label.hjust = 1,
        barwidth = 10,
        barheight = 0.5,
        nbin = 100
      )
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.text.y = element_blank(),  # Hide y-axis text
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "top",
      legend.title = element_text(size = 10),
      plot.margin = margin(5.5, 5.5, 5.5, 0)
    )
  
  # Combine plots using patchwork
  combined_plot <- boxplot + heatmap +
    plot_layout(
      widths = c(3, 2),
      guides = "collect",  # Collect legends together
      design = "
        AB"
    ) &
    theme(
      legend.position = "top",
      legend.box = "horizontal",
      legend.margin = margin(0, 20, 10, 20),
      legend.spacing.x = unit(0.5, "cm")
    ) &
    plot_annotation(
      title = "Distance to Vasculature by RGC Subtype",
      theme = theme(
        plot.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(hjust = 0.5, size = 12)
      )
    )
  
  return(combined_plot)
}

compute_vessel_ratios <- function(rgc_data, vasculature_df, lasso_region, radius = 90) {
  # Input validation and window bounds setup
  if (nrow(rgc_data) == 0) {
    warning("No RGC data provided")
    return(data.frame())
  }
  
  if (nrow(vasculature_df) == 0) {
    warning("No vasculature data provided")
    return(data.frame())
  }
  
  # Get window bounds from the lasso region coordinates
  coords <- st_coordinates(lasso_region)
  window_bounds <- list(
    xmin = min(coords[, "X"]),
    xmax = max(coords[, "X"]),
    ymin = min(coords[, "Y"]),
    ymax = max(coords[, "Y"])
  )
  
  window_width <- window_bounds$xmax - window_bounds$xmin
  window_height <- window_bounds$ymax - window_bounds$ymin
  
  cat(sprintf("\nWindow dimensions: %.1f x %.1f\n", window_width, window_height))
  cat(sprintf("Analysis radius: %.1f\n", radius))
  cat(sprintf("RGC X range: %.1f to %.1f\n", min(rgc_data$x), max(rgc_data$x)))
  cat(sprintf("RGC Y range: %.1f to %.1f\n", min(rgc_data$y), max(rgc_data$y)))
  
  cat(sprintf("Processing %d cells...\n", nrow(rgc_data)))
  
  # Initialize results list
  results_list <- list()
  valid_cells <- 0
  boundary_cells <- 0
  no_vessel_cells <- 0
  
  for(i in 1:nrow(rgc_data)) {
    cell_x <- rgc_data$x[i]
    cell_y <- rgc_data$y[i]
    
    # Check if cell's radius touches window boundary
    is_boundary <- cell_x - radius <= window_bounds$xmin || 
      cell_x + radius >= window_bounds$xmax ||
      cell_y - radius <= window_bounds$ymin ||
      cell_y + radius >= window_bounds$ymax
    
    if(is_boundary) {
      boundary_cells <- boundary_cells + 1
      next
    }
    
    # Find vessel points within radius
    dist_to_vessels <- sqrt((vasculature_df$x - cell_x)^2 + (vasculature_df$y - cell_y)^2)
    vessels_in_radius <- vasculature_df[dist_to_vessels <= radius, ]
    
    # Skip cells with no vessels within radius
    if(nrow(vessels_in_radius) == 0) {
      no_vessel_cells <- no_vessel_cells + 1
      next
    }
    
    # Count small and large vessels
    small_count <- sum(vessels_in_radius$small)
    large_count <- sum(!vessels_in_radius$small)
    
    # Calculate ratio (small/large)
    ratio <- if(large_count > 0) small_count / large_count else small_count
    
    # Calculate index that ranges from -1 to 1
    # -1: exclusively large vessels
    # 0: equal proportion
    # 1: exclusively small vessels
    total_count <- small_count + large_count
    index <- (small_count - large_count) / total_count
    
    # Store results
    valid_cells <- valid_cells + 1
    results_list[[valid_cells]] <- data.frame(
      X = cell_x,
      Y = cell_y,
      small = small_count,
      large = large_count,
      ratio = ratio,
      index = index,
      Prediction = rgc_data$Prediction[i]
    )
  }
  
  cat(sprintf("Found %d valid cells:\n", valid_cells))
  cat(sprintf("- %d cells excluded due to boundary\n", boundary_cells))
  cat(sprintf("- %d cells excluded due to no vessels within radius\n", no_vessel_cells))
  
  # Combine all results
  if (length(results_list) == 0) {
    warning("No valid cells found in region")
    return(data.frame())
  }
  
  results_df <- do.call(rbind, results_list)
  return(results_df)
}
# compute_vessel_ratios <- function(rgc_data, vasculature_df, lasso_region, radius = 90) {
#   # Input validation and window bounds setup remain the same
#   if (nrow(rgc_data) == 0) {
#     warning("No RGC data provided")
#     return(data.frame())
#   }
#   
#   if (nrow(vasculature_df) == 0) {
#     warning("No vasculature data provided")
#     return(data.frame())
#   }
#   
#   # Get window bounds from the lasso region coordinates
#   coords <- st_coordinates(lasso_region)
#   window_bounds <- list(
#     xmin = min(coords[, "X"]),
#     xmax = max(coords[, "X"]),
#     ymin = min(coords[, "Y"]),
#     ymax = max(coords[, "Y"])
#   )
#   
#   window_width <- window_bounds$xmax - window_bounds$xmin
#   window_height <- window_bounds$ymax - window_bounds$ymin
#   
#   cat(sprintf("\nWindow dimensions: %.1f x %.1f\n", window_width, window_height))
#   cat(sprintf("Analysis radius: %.1f\n", radius))
#   cat(sprintf("RGC X range: %.1f to %.1f\n", min(rgc_data$x), max(rgc_data$x)))
#   cat(sprintf("RGC Y range: %.1f to %.1f\n", min(rgc_data$y), max(rgc_data$y)))
#   
#   cat(sprintf("Processing %d cells...\n", nrow(rgc_data)))
#   
#   # Initialize results list
#   results_list <- list()
#   valid_cells <- 0
#   boundary_cells <- 0
#   
#   for(i in 1:nrow(rgc_data)) {
#     cell_x <- rgc_data$x[i]
#     cell_y <- rgc_data$y[i]
#     
#     # Check if cell's radius touches window boundary
#     is_boundary <- cell_x - radius <= window_bounds$xmin || 
#       cell_x + radius >= window_bounds$xmax ||
#       cell_y - radius <= window_bounds$ymin ||
#       cell_y + radius >= window_bounds$ymax
#     
#     if(is_boundary) {
#       boundary_cells <- boundary_cells + 1
#       next
#     }
#     
#     # Find vessel points within radius
#     dist_to_vessels <- sqrt((vasculature_df$x - cell_x)^2 + (vasculature_df$y - cell_y)^2)
#     vessels_in_radius <- vasculature_df[dist_to_vessels <= radius, ]
#     
#     # Count small and large vessels
#     small_count <- sum(vessels_in_radius$small)
#     large_count <- sum(!vessels_in_radius$small)
#     
#     # Calculate ratio (now small/large instead of large/small)
#     # When large_count is 0, set to maximum value (determined by counts)
#     ratio <- if(large_count > 0) small_count / large_count else small_count
#     
#     # Calculate new index that ranges from -1 to 1
#     # -1: exclusively large vessels
#     # 0: equal proportion of small and large vessels
#     # 1: exclusively small vessels
#     if(small_count == 0 && large_count == 0) {
#       index <- 0  # No vessels nearby
#     } else {
#       total_count <- small_count + large_count
#       index <- (small_count - large_count) / total_count
#     }
#     
#     # Store results
#     valid_cells <- valid_cells + 1
#     results_list[[valid_cells]] <- data.frame(
#       X = cell_x,
#       Y = cell_y,
#       small = small_count,
#       large = large_count,
#       ratio = ratio,
#       index = index,
#       Prediction = rgc_data$Prediction[i]
#     )
#   }
#   
#   cat(sprintf("Found %d valid cells (%d excluded due to boundary)\n", 
#               valid_cells, boundary_cells))
#   
#   # Combine all results
#   if (length(results_list) == 0) {
#     warning("No valid cells found in region")
#     return(data.frame())
#   }
#   
#   results_df <- do.call(rbind, results_list)
#   return(results_df)
# }

analyze_vessel_proximity <- function(rgc_data, vasculature_df, lasso_region, radius = 90) {
  # Get window bounds from the lasso region coordinates
  coords <- st_coordinates(lasso_region)
  window_bounds <- list(
    xmin = min(coords[, "X"]),
    xmax = max(coords[, "X"]),
    ymin = min(coords[, "Y"]),
    ymax = max(coords[, "Y"])
  )
  
  # Initialize results storage
  subtypes <- unique(rgc_data$Prediction)
  results <- data.frame(
    Subtype = subtypes,
    total_cells = 0,
    cells_near_small = 0,
    cells_near_large = 0,
    cells_near_both = 0,
    cells_near_neither = 0
  )
  
  # Process each subtype
  for(st in subtypes) {
    subtype_data <- rgc_data[rgc_data$Prediction == st,]
    total_cells <- nrow(subtype_data)
    near_small <- 0
    near_large <- 0
    near_both <- 0
    near_neither <- 0
    
    # Process each cell
    for(i in 1:nrow(subtype_data)) {
      cell_x <- subtype_data$x[i]
      cell_y <- subtype_data$y[i]
      
      # Skip if cell is too close to boundary
      if(cell_x - radius <= window_bounds$xmin || 
         cell_x + radius >= window_bounds$xmax ||
         cell_y - radius <= window_bounds$ymin || 
         cell_y + radius >= window_bounds$ymax) {
        total_cells <- total_cells - 1
        next
      }
      
      # Find vessels within radius
      dist_to_vessels <- sqrt((vasculature_df$x - cell_x)^2 + (vasculature_df$y - cell_y)^2)
      vessels_in_radius <- vasculature_df[dist_to_vessels <= radius,]
      
      # Count vessel types
      has_small <- any(vessels_in_radius$small)
      has_large <- any(!vessels_in_radius$small)
      
      if(has_small && has_large) near_both <- near_both + 1
      else if(has_small) near_small <- near_small + 1
      else if(has_large) near_large <- near_large + 1
      else near_neither <- near_neither + 1
    }
    
    # Update results
    results$total_cells[results$Subtype == st] <- total_cells
    results$cells_near_small[results$Subtype == st] <- near_small
    results$cells_near_large[results$Subtype == st] <- near_large
    results$cells_near_both[results$Subtype == st] <- near_both
    results$cells_near_neither[results$Subtype == st] <- near_neither
  }
  
  # Calculate percentages
  results <- results %>%
    mutate(
      pct_near_small = 100 * cells_near_small / total_cells,
      pct_near_large = 100 * cells_near_large / total_cells,
      pct_near_both = 100 * cells_near_both / total_cells,
      pct_near_neither = 100 * cells_near_neither / total_cells
    )
  
  return(results)
}

plot_vessel_proximity_bars <- function(proximity_data, slide_id = NULL) {
  # Prepare data for plotting
  plot_data <- proximity_data %>%
    pivot_longer(
      cols = c(pct_near_small, pct_near_large, pct_near_both, pct_near_neither),
      names_to = "category",
      values_to = "percentage"
    ) %>%
    mutate(
      category = factor(category,
                        levels = c("pct_near_both", "pct_near_small", 
                                   "pct_near_large", "pct_near_neither"),
                        labels = c("Both", "Small Only", 
                                   "Large Only", "Neither"))
    )
  
  # Create the plot
  title_text <- if(is.null(slide_id)) {
    "Cell Distribution by Vessel Proximity"
  } else {
    sprintf("Cell Distribution by Vessel Proximity - Slide %s", slide_id)
  }
  
  p <- ggplot(plot_data, aes(x = Subtype, y = percentage, fill = category)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c(
      "Both" = "#4CAF50",      # Green
      "Small Only" = "#2196F3", # Blue
      "Large Only" = "#F44336", # Red
      "Neither" = "#9E9E9E"     # Gray
    )) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "right",
      legend.title = element_blank()
    ) +
    labs(
      title = title_text,
      x = "Cell Type",
      y = "Percentage of Cells"
    )
  
  return(p)
}

plot_vessel_proximity_boxplots <- function(all_slides_data) {
  require(dplyr)
  require(tidyr)
  
  # Calculate the mean values for sorting
  mean_values <- all_slides_data %>%
    group_by(Subtype) %>%
    summarise(
      mean_small = mean(cells_near_small + cells_near_both),
      mean_large = mean(cells_near_large + cells_near_both)
    ) %>%
    arrange(desc(mean_small))  # Sort by small vessel counts
  
  # Create ordered factor levels based on small vessel counts
  subtype_order <- mean_values$Subtype
  
  # Prepare data for plotting
  plot_data <- all_slides_data %>%
    mutate(
      small_total = cells_near_small + cells_near_both,
      large_total = cells_near_large + cells_near_both,
      Subtype = factor(Subtype, levels = subtype_order)  # Apply ordering
    ) %>%
    pivot_longer(
      cols = c(small_total, large_total),
      names_to = "vessel_type",
      values_to = "count"
    ) %>%
    mutate(
      vessel_type = factor(vessel_type,
                           levels = c("small_total", "large_total"),
                           labels = c("Small Vessels", "Large Vessels"))
    )
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = count, y = Subtype)) +
    geom_boxplot(aes(fill = vessel_type), width = 0.7) +
    scale_fill_manual(values = c(
      "Small Vessels" = "#2196F3",  # Blue
      "Large Vessels" = "#F44336"   # Red
    )) +
    facet_wrap(~vessel_type, scales = "free_x", ncol = 2) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "none",
      strip.text = element_text(size = 12, face = "bold"),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 14, hjust = 0.5)
    ) +
    labs(
      title = "Cell Counts by Vessel Proximity Across Slides",
      x = "Number of Cells"
    )
  
  return(p)
}

# Example of how to add summary statistics to the plotting data
add_summary_stats <- function(all_slides_data) {
  summary_stats <- all_slides_data %>%
    group_by(Subtype) %>%
    summarise(
      small_mean = mean(cells_near_small + cells_near_both),
      small_sd = sd(cells_near_small + cells_near_both),
      large_mean = mean(cells_near_large + cells_near_both),
      large_sd = sd(cells_near_large + cells_near_both),
      n_slides = n()
    ) %>%
    arrange(desc(small_mean))  # Sort by small vessel counts
  
  return(summary_stats)
}

analyze_normalized_vessel_proximity <- function(rgc_data, vasculature_df, lasso_region, radius = 90) {
  # Get window bounds from the lasso region coordinates
  coords <- st_coordinates(lasso_region)
  window_bounds <- list(
    xmin = min(coords[, "X"]),
    xmax = max(coords[, "X"]),
    ymin = min(coords[, "Y"]),
    ymax = max(coords[, "Y"])
  )
  
  # Initialize results storage
  subtypes <- unique(rgc_data$Prediction)
  results <- data.frame(
    Subtype = character(),
    total_cells_in_region = integer(),
    cells_near_small = integer(),
    cells_near_large = integer(),
    cells_near_both = integer(),
    proportion_near_small = numeric(),
    proportion_near_large = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Process each subtype
  for(st in subtypes) {
    subtype_data <- rgc_data[rgc_data$Prediction == st,]
    valid_cells <- 0
    near_small <- 0
    near_large <- 0
    near_both <- 0
    
    # Count total valid cells (not too close to boundary)
    total_cells <- sum(!( 
      subtype_data$x - radius <= window_bounds$xmin | 
        subtype_data$x + radius >= window_bounds$xmax |
        subtype_data$y - radius <= window_bounds$ymin | 
        subtype_data$y + radius >= window_bounds$ymax
    ))
    
    if(total_cells > 0) {  # Only process if we have valid cells
      # Process each cell
      for(i in 1:nrow(subtype_data)) {
        cell_x <- subtype_data$x[i]
        cell_y <- subtype_data$y[i]
        
        # Skip if cell is too close to boundary
        if(cell_x - radius <= window_bounds$xmin || 
           cell_x + radius >= window_bounds$xmax ||
           cell_y - radius <= window_bounds$ymin || 
           cell_y + radius >= window_bounds$ymax) {
          next
        }
        
        # Find vessels within radius
        dist_to_vessels <- sqrt((vasculature_df$x - cell_x)^2 + (vasculature_df$y - cell_y)^2)
        vessels_in_radius <- vasculature_df[dist_to_vessels <= radius,]
        
        # Count vessel types
        has_small <- any(vessels_in_radius$small)
        has_large <- any(!vessels_in_radius$small)
        
        if(has_small && has_large) near_both <- near_both + 1
        else if(has_small) near_small <- near_small + 1
        else if(has_large) near_large <- near_large + 1
      }
      
      # Calculate proportions
      total_near_small <- (near_small + near_both) / total_cells
      total_near_large <- (near_large + near_both) / total_cells
      
      # Add to results
      results <- rbind(results, data.frame(
        Subtype = st,
        total_cells_in_region = total_cells,
        cells_near_small = near_small,
        cells_near_large = near_large,
        cells_near_both = near_both,
        proportion_near_small = total_near_small,
        proportion_near_large = total_near_large,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(results)
}

plot_normalized_vessel_proximity <- function(all_regions_data) {
  require(dplyr)
  require(tidyr)
  
  # Calculate mean values for sorting
  mean_values <- all_regions_data %>%
    group_by(Subtype) %>%
    summarise(
      mean_small = mean(proportion_near_small, na.rm = TRUE)
    ) %>%
    arrange(desc(mean_small))
  
  # Create ordered factor levels based on small vessel proportions
  subtype_order <- mean_values$Subtype
  
  # Prepare data for plotting
  plot_data <- all_regions_data %>%
    mutate(Subtype = factor(Subtype, levels = subtype_order)) %>%
    pivot_longer(
      cols = c(proportion_near_small, proportion_near_large),
      names_to = "vessel_type",
      values_to = "proportion"
    ) %>%
    mutate(
      vessel_type = factor(vessel_type,
                           levels = c("proportion_near_small", "proportion_near_large"),
                           labels = c("Small Vessels", "Large Vessels"))
    )
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = proportion * 100, y = Subtype)) +  # Convert to percentage
    geom_boxplot(aes(fill = vessel_type), width = 0.7) +
    scale_fill_manual(values = c(
      "Small Vessels" = "#2196F3",  # Blue
      "Large Vessels" = "#F44336"   # Red
    )) +
    facet_wrap(~vessel_type, scales = "free_x", ncol = 2) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "none",
      strip.text = element_text(size = 12, face = "bold"),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 14, hjust = 0.5)
    ) +
    labs(
      title = "Proportion of Cells Near Vessels by Region",
      x = "Percentage of Cells (%)"
    )
  
  return(p)
}
################################################################################

analyze_vessel_ratios <- function(lasso_list, rgc_df, interpolated_vasculature, slide_id, radius) {
  if (length(lasso_list) == 0) {
    warning("No lasso regions provided")
    return(data.frame())
  }
  
  all_results <- list()
  valid_regions <- 0
  
  # Process each lasso region
  for(i in seq_along(lasso_list)) {
    region_name <- names(lasso_list)[i]
    cat(paste("\nProcessing region:", region_name, "\n"))
    
    # Get current lasso region
    current_lasso <- lasso_list[[i]]
    
    # Create owin object for filtering
    lasso_window <- sf_to_owin(current_lasso)
    
    # Filter RGC data for this slide and region
    rgc_data <- rgc_df %>% 
      filter(slide == slide_id) %>%
      filter_points_in_window(., lasso_window)
    
    # Compute vessel ratios
    ratio_results <- compute_vessel_ratios(
      rgc_data = rgc_data,
      vasculature_df = interpolated_vasculature,
      lasso_region = current_lasso,
      radius = radius
    )
    
    # Only store results if we found valid cells
    if (nrow(ratio_results) > 0) {
      valid_regions <- valid_regions + 1
      ratio_results$region <- region_name
      ratio_results$slide <- slide_id
      all_results[[valid_regions]] <- ratio_results
    }
  }
  
  cat(sprintf("\nProcessed %d regions with valid cells\n", valid_regions))
  
  if (length(all_results) == 0) {
    warning("No valid results found in any region")
    return(data.frame())
  }
  
  final_results <- do.call(rbind, all_results)
  return(final_results)
}


#############################################################################################
################################################################################################
##############################################################################################

cat("\nStarting analysis pipeline\n")
rgc_root <- '/media/sam/Data2/baysor_rbpms_consolidated'
lasso_root <- '/media/sam/New Volume/Xenium_Data/HQ_NearestNeighbor_Zones'

# Load RGC data
cat("Loading RGC data...")
rgc_df <- read_csv(paste0(rgc_root,'/OldModel/merged_rgc_prediction_expmat.csv'))
cat("Done\n")

# Find all vasculature data files
root <- '/media/sam/New Volume/Xenium_Data/IHC_experiments/Keyence/TomatoLectin/'
vasculature_files <- list.files(root, pattern = "retina_\\d+\\.csv$", full.names = TRUE)

# Extract slide IDs from filenames
slide_ids <- str_extract(basename(vasculature_files), "\\d+")
cat(sprintf("\nFound %d slides to process: %s\n", length(slide_ids), paste(slide_ids, collapse = ", ")))

# Create results directory if it doesn't exist
results_path <- file.path(rgc_root, "vascular_analysis")
dir.create(results_path, showWarnings = FALSE)

# Initialize lists to store combined results
all_final_results <- list()
all_summary_stats <- list()
all_vessel_ratio_results <- list()

# Define slides to skip interpolation (empty vector means interpolate all slides)
no_interpolation_slides <- c()  # Add slide IDs here, e.g., c("18300")





################################################################################3
################################################################################3
interpolate_vasculature <- function(vasculature_df, spacing = 0.5, max_segment_length = 500) {
  # Remove any NA values first
  vasculature_df <- na.omit(vasculature_df)
  
  cat(sprintf("\nInterpolating vasculature with %d points and %d unique groups\n", 
              nrow(vasculature_df), length(unique(vasculature_df$group))))
  
  # Initialize list to store interpolated points for each group
  interpolated_list <- list()
  next_group_id <- max(vasculature_df$group) + 1
  
  # Process each vessel group
  for(g in unique(vasculature_df$group)) {
    tryCatch({
      # Get current vessel points
      vessel <- vasculature_df %>% 
        filter(group == g) %>%
        arrange(x, y)
      
      # Get small/large status for this group
      is_small <- unique(vessel$small)
      if(length(is_small) != 1) {
        warning(sprintf("Group %d has inconsistent small/large status. Using first value.", g))
        is_small <- is_small[1]
      }
      
      cat(sprintf("\nProcessing group %d with %d points\n", g, nrow(vessel)))
      
      if(nrow(vessel) > 1) {
        # Convert to matrix for calculations
        points <- as.matrix(vessel[, c("x", "y")])
        
        # Calculate segment lengths
        diffs <- diff(points, 1)
        segment_lengths <- sqrt(rowSums(diffs^2))
        
        cat(sprintf("Segment lengths: min=%.2f, max=%.2f, mean=%.2f\n",
                    min(segment_lengths), max(segment_lengths), mean(segment_lengths)))
        
        # Process each sub-segment between original points
        for(i in 1:(nrow(points)-1)) {
          start_point <- points[i,]
          end_point <- points[i+1,]
          segment_length <- segment_lengths[i]
          
          # If segment is too long, split it into smaller pieces
          if(segment_length > max_segment_length) {
            # Calculate how many pieces to split into
            n_splits <- ceiling(segment_length / max_segment_length)
            sub_length <- segment_length / n_splits
            
            cat(sprintf("Splitting segment %d (length %.2f) into %d pieces\n", 
                        i, segment_length, n_splits))
            
            # Create intermediate points
            for(j in 1:n_splits) {
              # Calculate start and end points for this sub-segment
              sub_start <- start_point + (j-1) * (end_point - start_point) / n_splits
              sub_end <- start_point + j * (end_point - start_point) / n_splits
              
              # Calculate number of interpolation points for this sub-segment
              n_points <- max(2, ceiling(sub_length/spacing))
              
              # Create sequence of points
              t <- seq(0, 1, length.out = n_points)
              new_x <- sub_start[1] + (sub_end[1] - sub_start[1]) * t
              new_y <- sub_start[2] + (sub_end[2] - sub_start[2]) * t
              
              # Create data frame with new group ID
              interpolated_list[[length(interpolated_list) + 1]] <- data.frame(
                x = new_x,
                y = new_y,
                group = next_group_id,
                small = is_small
              )
              next_group_id <- next_group_id + 1
            }
          } else {
            # Process normal-length segment
            n_points <- max(2, ceiling(segment_length/spacing))
            
            # Create sequence of points
            t <- seq(0, 1, length.out = n_points)
            new_x <- start_point[1] + (end_point[1] - start_point[1]) * t
            new_y <- start_point[2] + (end_point[2] - start_point[2]) * t
            
            interpolated_list[[length(interpolated_list) + 1]] <- data.frame(
              x = new_x,
              y = new_y,
              group = g,
              small = is_small
            )
          }
        }
      } else {
        # If vessel only has one point, keep it as is
        interpolated_list[[length(interpolated_list) + 1]] <- vessel
      }
    }, error = function(e) {
      warning(sprintf("Error processing group %d: %s", g, e$message))
    })
  }
  
  # Combine all vessels, checking for empty list
  if(length(interpolated_list) == 0) {
    warning("No valid interpolated points generated")
    return(vasculature_df)  # Return original data if interpolation failed
  }
  
  interpolated_df <- do.call(rbind, interpolated_list)
  
  cat(sprintf("\nInterpolation complete. Original points: %d, Interpolated points: %d\n",
              nrow(vasculature_df), nrow(interpolated_df)))
  cat(sprintf("Original groups: %d, Final groups: %d\n",
              length(unique(vasculature_df$group)), length(unique(interpolated_df$group))))
  
  return(interpolated_df)
}
#################################################################################




radius = 20

all_proximity_results <- list()
all_normalized_proximity_results <- list()
# Process each slide
for(slide_idx in seq_along(slide_ids)) {
  slide_id <- slide_ids[slide_idx]
  cat(sprintf("\n\nProcessing slide %s (%d of %d)\n", slide_id, slide_idx, length(slide_ids)))
  
  # Load Vasculature data
  cat("Loading Vasculature data...")
  path <- file.path(root, sprintf("retina_%s.csv", slide_id))
  
  # Read data with explicit column names
  vasculature_df <- read_csv(path, 
                             col_names = c("x", "y", "group", "small"),
                             col_types = cols(
                               x = col_double(),
                               y = col_double(),
                               group = col_double(),
                               small = col_logical()
                             ))
  
  # Validate data
  if (!all(c("x", "y", "group", "small") %in% names(vasculature_df))) {
    stop(sprintf("Invalid data format in file for slide %s", slide_id))
  }
  
  # Print summary of data
  cat(sprintf("\nLoaded %d points with %d unique vessel groups\n", 
              nrow(vasculature_df), length(unique(vasculature_df$group))))
  
  # Check for any NA values
  na_count <- sum(is.na(vasculature_df))
  if (na_count > 0) {
    warning(sprintf("Found %d NA values in data for slide %s", na_count, slide_id))
    # Remove rows with NA values
    vasculature_df <- na.omit(vasculature_df)
  }
  

  # Check if this slide should skip interpolation
  if(slide_id %in% no_interpolation_slides) {
    cat("Skipping interpolation for slide", slide_id, "\n")
    interpolated_vasculature <- vasculature_df
  } else {
    # Interpolate vasculature before analysis
    cat("Interpolating vasculature points...")
    interpolated_vasculature <- interpolate_vasculature(vasculature_df, spacing = 0.5)
    cat("Done\n")
  }
  
  # Load lasso regions for this slide
  cat(paste("Loading Lasso data for slide", slide_id, "...\n"))
  lasso_file <- file.path(lasso_root, paste0(slide_id, "_HQ_Lasso_coordinates.csv"))
  if (!file.exists(lasso_file)) {
    warning(paste("Lasso file not found for slide", slide_id, "- skipping"))
    next
  }
  lasso_list <- load_lasso_regions(lasso_file)
  cat("Done\n")
  
  # Initialize list to store results for this slide
  slide_results <- list()
  
  # Process each lasso region
  slide_proximity_results <- list()
  slide_normalized_proximity_results <- list()
  
  # Process each lasso region
  for (i in seq_along(lasso_list)) {
    region_name <- names(lasso_list)[i]
    current_lasso <- lasso_list[[i]]
    lasso_window <- sf_to_owin(current_lasso)
    
    # Filter RGC data for this slide and region
    rgc_data <- rgc_df %>% 
      filter(slide == slide_id) %>%
      filter_points_in_window(., lasso_window)
    
    # Add proximity analysis
    proximity_results <- analyze_vessel_proximity(
      rgc_data = rgc_data,
      vasculature_df = interpolated_vasculature,
      lasso_region = current_lasso,
      radius = radius
    )
    
    # Store results
    slide_proximity_results[[region_name]] <- proximity_results
    
    # Add normalized proximity analysis
    normalized_proximity_results <- analyze_normalized_vessel_proximity(
      rgc_data = rgc_data,
      vasculature_df = interpolated_vasculature,
      lasso_region = current_lasso,
      radius = radius
    )
    
    # Add region identifier
    normalized_proximity_results$region <- region_name
    normalized_proximity_results$slide <- slide_id
    
    # Store results
    slide_normalized_proximity_results[[region_name]] <- normalized_proximity_results
    
    # Compute nearest neighbor distances using original vasculature data
    NN_dist_df <- compute_vascular_distances(rgc_data, interpolated_vasculature, slide_id)
    
    # Create and save vessel distribution plot for this region
    plot <- plot_region_with_vasculature(
      lasso_region = current_lasso,
      rgc_data = rgc_data,
      vasculature_df = interpolated_vasculature,  
      region_name = region_name
    )
    
    # Save vessel distribution plot
    plot_path <- file.path(results_path, sprintf("region_plot_%s_%s.pdf", slide_id, region_name))
    ggsave(plot_path, plot, width = 10, height = 8)
    
    # Create and save statistical plots for this region if we have data
    if(nrow(NN_dist_df) > 0) {
      region_stat_results <- analyze_vascular_distances(NN_dist_df)
      p_region_dist <- plot_vascular_distances_with_heatmap(NN_dist_df, region_stat_results)
      ggsave(file.path(results_path, sprintf("distance_analysis_%s_region_%s.pdf", slide_id, region_name)),
             p_region_dist, width = 12, height = 8)
    }
    
    # Store results
    slide_results[[region_name]] <- NN_dist_df
  }
  
  # Combine normalized proximity results for this slide
  all_normalized_proximity_results[[slide_id]] <- do.call(rbind, slide_normalized_proximity_results)
  
  # Combine proximity results for this slide
  combined_proximity <- do.call(rbind, slide_proximity_results) %>%
    group_by(Subtype) %>%
    summarise(
      total_cells = sum(total_cells),
      cells_near_small = sum(cells_near_small),
      cells_near_large = sum(cells_near_large),
      cells_near_both = sum(cells_near_both),
      cells_near_neither = sum(cells_near_neither)
    ) %>%
    mutate(
      pct_near_small = 100 * cells_near_small / total_cells,
      pct_near_large = 100 * cells_near_large / total_cells,
      pct_near_both = 100 * cells_near_both / total_cells,
      pct_near_neither = 100 * cells_near_neither / total_cells,
      slide = slide_id
    )
  
  # Store combined proximity results
  all_proximity_results[[slide_id]] <- combined_proximity
  
  # Create and save proximity plot for this slide
  p_proximity <- plot_vessel_proximity_bars(combined_proximity, slide_id)
  ggsave(file.path(results_path, sprintf("vessel_proximity_plot_slide_%s.pdf", slide_id)),
         p_proximity, width = 12, height = 8)
  
  
  # Combine results for this slide
  final_results <- do.call(rbind, slide_results)
  all_final_results[[slide_id]] <- final_results
  
  # Add summary statistics
  summary_stats <- final_results %>%
    group_by(Subtype) %>%
    summarise(
      n_cells = n(),
      mean_small_dist = mean(small_vasc_NN),
      sd_small_dist = sd(small_vasc_NN),
      mean_large_dist = mean(large_vasc_NN),
      sd_large_dist = sd(large_vasc_NN)
    )
  all_summary_stats[[slide_id]] <- summary_stats
  
  # Run vessel ratio analysis
  vessel_ratio_results <- analyze_vessel_ratios(
    lasso_list = lasso_list,
    rgc_df = rgc_df,
    interpolated_vasculature = interpolated_vasculature,
    slide_id = slide_id,
    radius = radius
  )
  all_vessel_ratio_results[[slide_id]] <- vessel_ratio_results
  
  # Save individual slide results
  write_csv(final_results, 
            file.path(results_path, sprintf("vascular_distances_slide_%s.csv", slide_id)))
  write_csv(summary_stats, 
            file.path(results_path, sprintf("vascular_summary_slide_%s.csv", slide_id)))
  if (nrow(vessel_ratio_results) > 0) {
    write_csv(vessel_ratio_results, 
              file.path(results_path, sprintf("vessel_ratios_slide_%s.csv", slide_id)))
  }
  
  # Create and save statistical plots for this slide
  stat_results <- analyze_vascular_distances(final_results)
  p_dist <- plot_vascular_distances_with_heatmap(final_results, stat_results)
  ggsave(file.path(results_path, sprintf("distance_analysis_slide_%s.pdf", slide_id)),
         p_dist, width = 12, height = 8)
  
  # Create and save vessel preference index plot for this slide
  if (nrow(vessel_ratio_results) > 0) {
    # Calculate statistics for each cell type
    cell_type_stats <- vessel_ratio_results %>%
      group_by(Prediction) %>%
      summarise(
        mean_index = mean(index),
        sd_index = sd(index),
        n = n(),
        se_index = sd_index / sqrt(n)
      ) %>%
      arrange(mean_index)  # Sort by mean index
    
    # Create ordered factor levels
    cell_order <- cell_type_stats$Prediction
    vessel_ratio_results$Prediction <- factor(vessel_ratio_results$Prediction, 
                                              levels = cell_order)
    cell_type_stats$Prediction <- factor(cell_type_stats$Prediction, 
                                         levels = cell_order)
    
    # Create the index plot
    p_index <- ggplot() +
      geom_boxplot(data = vessel_ratio_results,
                   aes(x = Prediction, y = index),
                   width = 0.2, alpha = 0.8, outlier.shape = NA) +
      geom_point(data = cell_type_stats,
                 aes(x = Prediction, y = mean_index),
                 color = "red", size = 2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      scale_fill_manual(values = c("lightgray", "darkred")) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "none"
      ) +
      labs(
        title = sprintf("Vessel Preference by Cell Type - Slide %s", slide_id),
        subtitle = sprintf("1 indicates only small vessels in %.1f um radius, -1 indicates only large vessels", radius),
        x = "Cell Type",
        y = "Vessel Index\n(-1 = large vessels, 1 = small vessels)"
      )
    
    ggsave(file.path(results_path, sprintf("vessel_preference_plot_slide_%s.pdf", slide_id)),
           p_index, width = 10, height = 8)
  }
}
# Combine all normalized proximity results
all_regions_normalized <- do.call(rbind, all_normalized_proximity_results)

# Save combined normalized data
write_csv(all_regions_normalized, 
          file.path(results_path, "vessel_proximity_normalized_all_regions.csv"))

# Create summary statistics for normalized data
normalized_summary <- all_regions_normalized %>%
  group_by(Subtype) %>%
  summarise(
    mean_proportion_small = mean(proportion_near_small, na.rm = TRUE),
    sd_proportion_small = sd(proportion_near_small, na.rm = TRUE),
    mean_proportion_large = mean(proportion_near_large, na.rm = TRUE),
    sd_proportion_large = sd(proportion_near_large, na.rm = TRUE),
    n_regions = n()
  )

# Save normalized summary statistics
write_csv(normalized_summary,
          file.path(results_path, "vessel_proximity_normalized_summary.csv"))

# Create and save normalized boxplot
p_normalized_proximity <- plot_normalized_vessel_proximity(all_regions_normalized)
ggsave(file.path(results_path, "vessel_proximity_normalized_boxplot.pdf"),
       p_normalized_proximity, width = 15, height = 10)

# After combining all proximity results
all_slides_proximity <- do.call(rbind, all_proximity_results)

# Save combined proximity data
write_csv(all_slides_proximity, 
          file.path(results_path, "vessel_proximity_all_slides.csv"))

# Create and save summary statistics
summary_stats <- add_summary_stats(all_slides_proximity)
write_csv(summary_stats,
          file.path(results_path, "vessel_proximity_summary_stats.csv"))

# Create and save boxplot for all slides
p_proximity_boxplot <- plot_vessel_proximity_boxplots(all_slides_proximity)
ggsave(file.path(results_path, "vessel_proximity_boxplot_all_slides.pdf"),
       p_proximity_boxplot, width = 15, height = 10)

# Combine all results
combined_final_results <- do.call(rbind, all_final_results)
combined_summary_stats <- do.call(rbind, all_summary_stats) %>%
  group_by(Subtype) %>%
  summarise(
    total_cells = sum(n_cells),
    mean_small_dist = mean(mean_small_dist, na.rm = TRUE),
    sd_small_dist = mean(sd_small_dist, na.rm = TRUE),
    mean_large_dist = mean(mean_large_dist, na.rm = TRUE),
    sd_large_dist = mean(sd_large_dist, na.rm = TRUE)
  )
combined_vessel_ratios <- do.call(rbind, all_vessel_ratio_results)

# Save combined results
write_csv(combined_final_results, 
          file.path(results_path, "vascular_distances_all_slides.csv"))
write_csv(combined_summary_stats, 
          file.path(results_path, "vascular_summary_all_slides.csv"))
write_csv(combined_vessel_ratios, 
          file.path(results_path, "vessel_ratios_all_slides.csv"))

# Create combined statistical analysis
stat_results_combined <- analyze_vascular_distances(combined_final_results)
p_combined_dist <- plot_vascular_distances_with_heatmap(combined_final_results, stat_results_combined)
ggsave(file.path(results_path, "distance_analysis_all_slides.pdf"),
       p_combined_dist, width = 12, height = 8)

# Create combined vessel preference index plot
if (nrow(combined_vessel_ratios) > 0) {
  # Calculate statistics for each cell type
  combined_cell_type_stats <- combined_vessel_ratios %>%
    group_by(Prediction) %>%
    summarise(
      mean_index = mean(index),
      sd_index = sd(index),
      n = n(),
      se_index = sd_index / sqrt(n)
    ) %>%
    arrange(mean_index)  # Sort by mean index
  
  # Create ordered factor levels
  cell_order <- combined_cell_type_stats$Prediction
  combined_vessel_ratios$Prediction <- factor(combined_vessel_ratios$Prediction, 
                                              levels = cell_order)
  combined_cell_type_stats$Prediction <- factor(combined_cell_type_stats$Prediction, 
                                                levels = cell_order)
  
  # Create the combined index plot
  p_combined_index <- ggplot() +
    geom_boxplot(data = combined_vessel_ratios,
                 aes(x = Prediction, y = index),
                 width = 0.2, alpha = 0.8, outlier.shape = NA) +
    geom_point(data = combined_cell_type_stats,
               aes(x = Prediction, y = mean_index),
               color = "red", size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    scale_fill_manual(values = c("lightgray", "darkred")) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "none"
    ) +
    labs(
      title = "Vessel Preference by Cell Type - All Slides Combined",
      subtitle = sprintf("1 indicates only small vessels in %.1f um radius, -1 indicates only large vessels", radius),
      x = "Cell Type",
      y = "Vessel Index\n(-1 = large vessels, 1 = small vessels)"
    )
  
  ggsave(file.path(results_path, "vessel_preference_plot_all_slides.pdf"),
         p_combined_index, width = 10, height = 8)
}

cat("\nAnalysis complete. Combined results saved to:", results_path, "\n")