library(spatstat)
library(sf)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

# Analysis hyperparameters
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


# Modified load_squares function with debugging
load_squares <- function(filepath) {
  cat(sprintf("\nLoading squares from %s\n", basename(filepath)))
  
  # Read and print first few lines of raw file for debugging
  raw_lines <- readLines(filepath, n = 5)
  cat("\nFirst 5 lines of file:\n")
  print(raw_lines)
  
  # Skip the first two metadata rows, read from row 3
  coords_df <- read.csv(filepath, skip = 2)
  print_debug_info("Loaded coordinates dataframe", coords_df)
  
  # Split into list of squares based on Selection ID
  square_list <- split(coords_df, coords_df$Selection) %>%
    lapply(function(square) {
      if (!identical(square[1, c("X","Y")], square[nrow(square), c("X","Y")])) {
        square <- rbind(square, square[1,])
      }
      coords_matrix <- as.matrix(square[, c("X","Y")])
      st_polygon(list(coords_matrix))
    })
  
  print_debug_info("Created square list", square_list)
  return(square_list)
}


# Function to calculate spatial statistics for a point pattern
calculate_spatial_stats <- function(ppp_obj) {
  if (npoints(ppp_obj) < 3) {
    return(list(nnri = NA, vdri = NA))
  }
  
  # Calculate edge-corrected nearest neighbor distances
  nn <- nndist(ppp_obj, correction="border")  # border correction for nn distances
  
  # Calculate mean nearest neighbor distance for CSR (theoretical)
  lambda <- npoints(ppp_obj) / area(Window(ppp_obj))
  r_mean_theo <- 1 / (2 * sqrt(lambda))
  
  # NNRI with edge correction
  nnri <- mean(nn) / r_mean_theo
  
  # Calculate VDRI using tessellation with edge handling
  # Create Voronoi tessellation
  voro <- dirichlet(ppp_obj)
  
  # Create buffer to identify edge cells
  win <- Window(ppp_obj)
  buffer_distance <- sqrt(area(win)) * 0.1  # 10% of square root of area as buffer
  
  # Identify points that are far enough from the boundary
  inner_points <- bdist.points(ppp_obj) > buffer_distance
  
  if(sum(inner_points) < 3) {
    return(list(nnri = nnri, vdri = NA))
  }
  
  # Get areas only for cells not on the edge
  voro_areas <- tile.areas(voro)[inner_points]
  
  # Calculate VDRI using only interior cells
  vdri <- sd(voro_areas) / mean(voro_areas)
  
  return(list(
    nnri = nnri,
    vdri = vdri,
    n_edge = sum(!inner_points),  # number of edge cells
    n_interior = sum(inner_points)  # number of interior cells
  ))
}

# Function to create point pattern with proper window handling
create_ppp <- function(points_data, window_obj) {
  # Filter points to only those within the window
  temp_ppp <- ppp(
    x = points_data$x,
    y = points_data$y,
    window = window_obj,
    check = FALSE
  )
  
  inside <- inside.owin(temp_ppp$x, temp_ppp$y, window_obj)
  
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

# Modified analyze_group function with cell count threshold
# Modified analyze_group function to include edge statistics
analyze_group <- function(group_data, all_data, square_window, group_name, 
                          slide_id, square_id, n_bootstrap = 100) {
  
  # Pre-filter points to window
  all_data_filtered <- filter_points_in_window(all_data, square_window)
  group_data_filtered <- filter_points_in_window(group_data, square_window)
  
  # Check if we meet the minimum cell threshold
  if (nrow(group_data_filtered) < MIN_CELLS_PER_REGION) {
    cat(sprintf("\n    Skipping group %s: insufficient cells (%d) in square\n", 
                group_name, nrow(group_data_filtered)))
    return(NULL)
  }
  
  # Create observed point pattern
  ppp_obj <- ppp(
    x = group_data_filtered$x,
    y = group_data_filtered$y,
    window = square_window,
    check = TRUE
  )
  
  if(is.null(ppp_obj) || npoints(ppp_obj) < 3) {
    return(NULL)
  }
  
  n_points <- npoints(ppp_obj)
  obs_stats <- calculate_spatial_stats(ppp_obj)
  intensity <- n_points / area(square_window)
  
  # Bootstrap analysis
  bootstrap_nnri <- numeric(n_bootstrap)
  bootstrap_vdri <- numeric(n_bootstrap)
  
  for(i in 1:n_bootstrap) {
    sampled_indices <- sample(nrow(all_data_filtered), n_points)
    sampled_points <- all_data_filtered[sampled_indices, ]
    
    bootstrap_ppp <- ppp(
      x = sampled_points$x,
      y = sampled_points$y,
      window = square_window,
      check = TRUE
    )
    
    if(!is.null(bootstrap_ppp)) {
      boot_stats <- calculate_spatial_stats(bootstrap_ppp)
      bootstrap_nnri[i] <- boot_stats$nnri
      bootstrap_vdri[i] <- boot_stats$vdri
    }
  }
  
  # Calculate summary statistics
  nnri_null_mean <- mean(bootstrap_nnri, na.rm = TRUE)
  vdri_null_mean <- mean(bootstrap_vdri, na.rm = TRUE)
  p_value_nnri <- sum(bootstrap_nnri <= obs_stats$nnri, na.rm = TRUE) / sum(!is.na(bootstrap_nnri))
  p_value_vdri <- sum(bootstrap_vdri <= obs_stats$vdri, na.rm = TRUE) / sum(!is.na(bootstrap_vdri))
  
  return(list(
    Prediction = group_name,
    nnri = obs_stats$nnri,
    vdri = obs_stats$vdri,
    nnri_null = nnri_null_mean,
    vdri_null = vdri_null_mean,
    N = n_points,
    n_edge = obs_stats$n_edge,
    n_interior = obs_stats$n_interior,
    intensity = intensity,
    p_value_nnri = p_value_nnri,
    p_value_vdri = p_value_vdri,
    points_outside = nrow(group_data) - nrow(group_data_filtered)
  ))
}

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

# Modified analyze_slide function for square regions
analyze_slide <- function(rgc_data, square_list, slide_id) {
  cat(sprintf("\nAnalyzing slide %s\n", slide_id))
  cat(sprintf("  Total RGC data points: %d\n", nrow(rgc_data)))
  cat(sprintf("  Number of square regions: %d\n", length(square_list)))
  
  results <- list()
  
  for (square_idx in seq_along(square_list)) {
    cat(sprintf("  Processing square %d/%d\n", square_idx, length(square_list)))
    
    current_square <- square_list[[square_idx]]
    window <- as.owin(current_square)
    
    groups <- unique(rgc_data$Prediction)
    cat(sprintf("    Processing %d cell types\n", length(groups)))
    
    for (group in groups) {
      cat(sprintf("    Analyzing group %s...", group))
      
      group_data <- rgc_data %>% filter(Prediction == group)
      
      group_results <- analyze_group(
        group_data = group_data,
        all_data = rgc_data,
        square_window = window,
        group_name = group,
        slide_id = slide_id,
        square_id = square_idx
      )
      
      if (!is.null(group_results)) {
        group_results$slide <- slide_id
        group_results$square <- square_idx
        results[[length(results) + 1]] <- group_results
      }
      
      cat(" done\n")
    }
  }
  
  return(results)
}

# Main analysis pipeline
cat("\nStarting analysis pipeline\n")
root <- '/media/sam/Data2/baysor_rbpms_consolidated'
squares_root <- '/media/sam/New Volume/Xenium_Data/NearestNeighborSquares'

# Load RGC data
cat("Loading RGC data...")
rgc_df <- read_csv(paste0(root,'/OldModel/all_rgc_prediction_expmat.csv'))
cat(" done\n")

# Get list of square files
cat("Finding square files...\n")
square_files <- list.files(squares_root, pattern = "\\d+_200x200squares\\.csv$", full.names = TRUE)
cat(sprintf("Found %d square files\n", length(square_files)))

# Process each slide
all_results <- list()
for (i in seq_along(square_files)) {
  square_file <- square_files[i]
  cat(sprintf("\nProcessing file %d of %d: %s\n", i, length(square_files), basename(square_file)))
  
  # Extract slide ID - modified to get just the numeric part
  slide_id <- as.numeric(gsub("(\\d+)_200x200squares\\.csv", "\\1", basename(square_file)))
  cat("\nExtracted slide_id:", slide_id, "\n")
  
  # Load squares
  square_list <- load_squares(square_file)
  
  # Filter RGC data for this slide
  cat("\nFiltering RGC data for slide", slide_id, "\n")
  slide_data <- rgc_df %>% filter(slide == slide_id)
  print_debug_info("Filtered slide data", slide_data)
  
  if (nrow(slide_data) == 0) {
    cat("\nWARNING: No matching data found for slide", slide_id, "\n")
    cat("Available slides in RGC data:\n")
    print(head(unique(rgc_df$slide)))
    next
  }
  
  # Analyze
  slide_results <- analyze_slide(slide_data, square_list, slide_id)
  all_results <- c(all_results, slide_results)
}


# Convert results to data frame
results_df <-  bind_rows(all_results) %>%
  group_by(square, slide) %>%
  mutate(N_square = sum(N)) %>%
  ungroup() %>%
  drop_na() %>%
  mutate(delta_nnri = nnri - nnri_null,
         delta_vdri = vdri - vdri_null) 

# Compute summary statistics
spatial_summary <- results_df %>%
  group_by(Prediction) %>%
  summarise(
    n_regions = n(),
    total_cells = sum(N),
    mean_cells_per_region = mean(N),
    
    mean_nnri = mean(nnri, na.rm = TRUE),
    sd_nnri = sd(nnri, na.rm = TRUE),
    se_nnri = sd_nnri / sqrt(n_regions),
    mean_nnri_null = mean(nnri_null, na.rm = TRUE),
    
    mean_vdri = mean(vdri, na.rm = TRUE),
    sd_vdri = sd(vdri, na.rm = TRUE),
    se_vdri = sd_vdri / sqrt(n_regions),
    mean_vdri_null = mean(vdri_null, na.rm = TRUE),
    
    mean_p_value_nnri = mean(p_value_nnri, na.rm = TRUE),
    mean_p_value_vdri = mean(p_value_vdri, na.rm = TRUE)
  ) %>%
  arrange(desc(total_cells))

# Print summary statistics
cat("\nSummary statistics:\n\n")
print(spatial_summary, digits = 3)

# Create visualization
plot_data <- results_df %>%
  pivot_longer(
    cols = c("nnri", "nnri_null", "vdri", "vdri_null"),
    names_to = c("metric", "distribution"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  mutate(
    metric = toupper(metric),
    distribution = ifelse(is.na(distribution), "observed", distribution)
  )

p <- ggplot(plot_data, aes(x = value, y = reorder(Prediction, -N), fill = distribution)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    alpha = 0.7,
    outlier.size = 0.5
  ) +
  facet_wrap(~ metric, scales = "free_x", nrow = 1) +
  scale_fill_manual(
    values = c("observed" = "black", "null" = "grey80"),
    labels = c("Observed", "Null Distribution")
  ) +
  theme_bw() +
  labs(
    x = "Value",
    title = "Spatial Statistics by Cell Type (Square Regions)"
  )

ggsave(paste0(root, "/spatial_statistics_squares.png"), p, width = 10, height = 8, dpi = 300)
# Save the detailed results dataframe
write_csv(results_df, file.path(root, "spatial_stats_results.csv"))

# Save the summary
write_csv(spatial_summary, file.path(root, "spatial_stats_summary.csv"))


################################################################################
results_df %>%
  filter(N_square>=250) %>%
  ggplot(aes(x = N, y = delta_nnri, color = N_square)) +
  geom_point() +
  geom_hline(yintercept=0, color = 'red') +
  facet_wrap(~Prediction, nrow = 7) +
  xlim(0,20)


# Function to reconstruct point pattern from raw data
create_square_ppp <- function(x_min, x_max, y_min, y_max, points_df) {
  # Create square window
  win <- owin(c(x_min, x_max), c(y_min, y_max))
  
  # Create ppp object
  ppp_obj <- ppp(
    x = points_df$x,
    y = points_df$y,
    window = win,
    check = TRUE
  )
  
  return(ppp_obj)
}

# Function to plot Voronoi comparison
# Function to create voronoi polygons from point pattern
create_voronoi_polygons <- function(points_df, window) {
  # Create ppp object
  ppp_obj <- ppp(
    x = points_df$x,
    y = points_df$y,
    window = window,
    check = TRUE
  )
  
  # Create Voronoi tessellation
  voro <- dirichlet(ppp_obj)
  
  # Extract tiles as polygons
  tiles <- tiles(voro)
  
  # Convert tiles to sf polygons
  polygons <- lapply(tiles, function(tile) {
    coords <- cbind(tile$bdry[[1]]$x, tile$bdry[[1]]$y)
    # Close the polygon by repeating first point
    coords <- rbind(coords, coords[1,])
    st_polygon(list(coords))
  })
  
  # Create sf object
  voro_sf <- st_sf(
    geometry = st_sfc(polygons),
    index = 1:length(polygons)
  ) %>%
    st_set_crs(NA)
  
  return(voro_sf)
}

# Modified plot_selected_squares function to properly select based on cell counts
plot_selected_squares <- function(all_results, target_prediction, rgc_df, squares_root, output_path, factor = "N") {
  # Create subdirectory
  dir_name <- paste0(target_prediction, "_", factor)
  dir_path <- file.path(output_path, dir_name)
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  
  # Filter results for target prediction and remove any NA values in factor
  filtered_results <- all_results %>%
    filter(Prediction == target_prediction) %>%
    filter(!is.na(!!sym(factor)))
  
  if(nrow(filtered_results) == 0) {
    stop("No results found for prediction: ", target_prediction)
  }
  
  # Find max, median, and min factor values
  sorted_values <- sort(unique(filtered_results[[factor]]))
  max_value <- max(sorted_values)
  min_value <- min(sorted_values)
  
  # For median, get the true median value
  n_values <- length(sorted_values)
  if(n_values %% 2 == 0) {
    median_indices <- c(n_values/2, n_values/2 + 1)
    median_value <- sorted_values[median_indices[1]]  # Take the lower of the two middle values
  } else {
    median_value <- sorted_values[ceiling(n_values/2)]
  }
  
  # Function to get one random result for a specific factor value
  get_random_result <- function(value) {
    matching_results <- filtered_results[filtered_results[[factor]] == value, ]
    matching_results[sample(nrow(matching_results), 1), ]
  }
  
  # Get results for plotting
  max_result <- get_random_result(max_value)
  median_result <- get_random_result(median_value)
  min_result <- get_random_result(min_value)
  
  # Create plots for each selected result
  plot_cases <- list(
    list(result = max_result, label = "max"),
    list(result = median_result, label = "median"),
    list(result = min_result, label = "min")
  )
  
  for(case in plot_cases) {
    result <- case$result
    output_file <- file.path(dir_path, 
                             paste0(dir_name, "_", case$label, ".png"))
    
    plot_voronoi_comparison(
      rgc_df = rgc_df,
      square_info = result,
      squares_root = squares_root,
      output_path = output_file
    )
  }
  
  return(list(
    max = max_result,
    median = median_result,
    min = min_result
  ))
}

# Function to create voronoi polygons from point pattern
create_voronoi_polygons <- function(points_df, window) {
  # Create ppp object
  ppp_obj <- ppp(
    x = points_df$x,
    y = points_df$y,
    window = window,
    check = TRUE
  )
  
  # Create Voronoi tessellation
  voro <- dirichlet(ppp_obj)
  
  # Extract tiles as polygons
  tiles <- tiles(voro)
  
  # Convert tiles to sf polygons
  polygons <- lapply(tiles, function(tile) {
    coords <- cbind(tile$bdry[[1]]$x, tile$bdry[[1]]$y)
    # Close the polygon by repeating first point
    coords <- rbind(coords, coords[1,])
    st_polygon(list(coords))
  })
  
  # Create sf object
  voro_sf <- st_sf(
    geometry = st_sfc(polygons),
    index = 1:length(polygons)
  ) %>%
    st_set_crs(NA)
  
  return(voro_sf)
}

# Updated plot_voronoi_comparison function
plot_voronoi_comparison <- function(rgc_df, square_info, squares_root, output_path) {
  # Find the correct square file for this slide
  square_file <- list.files(squares_root, 
                            pattern = paste0(square_info$slide, "_200x200squares\\.csv$"), 
                            full.names = TRUE)
  
  if(length(square_file) == 0) {
    stop("Could not find square file for slide ", square_info$slide)
  }
  
  # Load squares for this slide
  square_list <- load_squares(square_file[1])
  
  # Get the specific square geometry
  if(square_info$square > length(square_list)) {
    stop("Square index ", square_info$square, " exceeds number of squares in file")
  }
  current_square <- square_list[[square_info$square]]
  
  # Get square bounds
  square_bbox <- st_bbox(current_square)
  x_min <- square_bbox[["xmin"]]
  x_max <- square_bbox[["xmax"]]
  y_min <- square_bbox[["ymin"]]
  y_max <- square_bbox[["ymax"]]
  
  # Create window for voronoi tessellation
  square_window <- owin(c(x_min, x_max), c(y_min, y_max))
  
  # Filter data for this specific slide and square
  square_data <- rgc_df %>%
    filter(slide == square_info$slide) %>%
    filter(x >= x_min, x <= x_max,
           y >= y_min, y <= y_max)
  
  # Get target cell data
  target_cells <- square_data %>%
    filter(Prediction == square_info$Prediction) %>%
    select(x, y)
  
  # Get non-target cells for null distribution
  other_cells <- square_data %>%
    filter(Prediction != square_info$Prediction) %>%
    select(x, y)
  
  # Get background cells (all other cells in the region)
  background_cells <- square_data %>%
    filter(Prediction != square_info$Prediction) %>%
    select(x, y)
  
  # Sample from non-target cells for null distribution
  set.seed(42)  # for reproducibility
  null_sample <- other_cells %>%
    slice_sample(n = nrow(target_cells), replace = FALSE)
  
  # Create Voronoi tessellations
  target_voronoi <- create_voronoi_polygons(target_cells, square_window)
  null_voronoi <- create_voronoi_polygons(null_sample, square_window)
  
  # Create plot with ggplot
  p <- ggplot() +
    # White background
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      plot.caption = element_text(hjust = 0.5),
      legend.position = "right",
      aspect.ratio = 1
    ) +
    # Plot background cells first
    geom_point(data = background_cells, aes(x = x, y = y), 
               color = "grey50", alpha = 0.2, size = 0.5) +
    # Plot Voronoi tessellations with only borders
    geom_sf(data = null_voronoi, color = "grey50", fill = NA, linewidth = 0.3) +
    geom_sf(data = target_voronoi, color = "red", fill = NA, linewidth = 0.3) +
    # Plot the square boundary
    geom_sf(data = current_square, fill = NA, color = "black") +
    # Plot sampled points
    geom_point(data = null_sample, aes(x = x, y = y, color = "Random sample"), 
               size = 2) +
    geom_point(data = target_cells, aes(x = x, y = y, color = "Observed cells"), 
               size = 2) +
    # Force square aspect ratio
    coord_sf(expand = FALSE) +
    # Colors and legend
    scale_color_manual(values = c("Observed cells" = "red", 
                                  "Random sample" = "grey50")) +
    # Labels
    labs(title = square_info$Prediction,
         subtitle = sprintf("Slide: %d, Square: %d\nN cells: %d, N_square: %d",
                            square_info$slide, square_info$square,
                            square_info$N, square_info$N_square),
         caption = sprintf(
           "NNRI = %.3f (null: %.3f), p = %.3f\nVDRI = %.3f (null: %.3f), p = %.3f",
           square_info$nnri, square_info$nnri_null, square_info$p_value_nnri,
           square_info$vdri, square_info$vdri_null, square_info$p_value_vdri),
         color = "Cells")
  
  # Save the plot
  ggsave(output_path, p, width = 10, height = 10, dpi = 300, bg = "white")
}


plot_selected_squares(results_df, 
                      target_prediction = "02_W3D1.2",
                      rgc_df = rgc_df,
                      squares_root = squares_root,  # Path to directory with square files
                      output_path = root,
                      factor = "N_square")
