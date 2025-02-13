# Load required libraries
library(tiff)
library(smerc)
library(tidyverse)
library(sf)
library(patchwork)
library(spdep)

#=============================================================================
# Configuration Parameters
#=============================================================================
  #=============================================================================
  # File Paths and Directories
  #=============================================================================
  ROOT_DIR <- "/home/sam/FinalRGC_xenium/"               # Base directory for analysis
  TIF_PATH <- "RGC_stack_pixel_norm.tif"        # Input TIF file containing cell type maps
  OUTPUT_DIR <- "GlobalStatistics_PixelNormalized"                 # Main output directory name
  MORANS_DIR <- "Morans"                                 # Directory for Moran's I results
  SCAN_DIR <- "Scan"                                     # Directory for scan statistics results
  VISUAL_SCENE_DIR <- "VisualScene"                      # Directory for visual scene analysis
  
  #=============================================================================
  # Scan Analysis Parameters
  #=============================================================================
  # Grid and Buffer
  GRID_SIZE <- 15          # Number of cells per side in analysis grid
  BUFFER_SIZE <- 0         # Pixels to remove from edge to prevent boundary effects
  
  # Scan Statistics
  ALPHA <- 0.05           # Statistical significance threshold
  N_SIM <- 99           # Number of Monte Carlo simulations for scan statistics
  UBPOP <- 2/3           # Maximum population size for potential clusters
  MIN_CASES <- 2         # Minimum number of cases required in a cluster
  MIN_CLUSTER_SIZE <- 5  # Minimum number of grid cells in a cluster
  NORMALIZE <- FALSE     # Whether to normalize values before analysis
  INVERT <- FALSE       # Whether to invert normalized values
  
  # Bootstrap Analysis
  N_BOOTSTRAP <- 1000    # Number of bootstrap iterations for visual scene analysis
  MIN_COVERAGE <- 25     # Minimum coverage percentage for visualization
  BOOTSTRAP_METHOD = "gaussian" # How the bootstrap is performed, options are "shuffle", "gaussian", "spin"
  
  #=============================================================================
  # Visual Scene Analysis Parameters
  #=============================================================================
  # Mask Generation
  N_ANNULI <- 6         # Number of concentric rings in retinal masks
  N_SECTORS <- 24       # Number of angular sectors in retinal masks
  TEMPORAL_LEFT <- FALSE # Orientation of temporal retinal region
  SLOPE <- 15           # Angle narrowing per annulus for binocular zone
  SECTOR_FILL_0 <- 140  # Start angle for binocular zone
  SECTOR_FILL_1 <- 315  # End angle for binocular zone
  OUTER_FILL <- 2       # Number of outer annuli to include in masks
  INNER_FILL <- 0       # Number of inner annuli to include in masks
  
  #=============================================================================
  # Visualization Parameters
  #=============================================================================
  PLOT_WIDTH <- 8       # Default plot width in inches
  PLOT_HEIGHT <- 8      # Default plot height in inches
  PLOT_DPI <- 300       # Plot resolution
  BOUNDARY_COLOR <- "green3"     # Color for cluster boundary visualization
  OVERLAY_COLOR <- "#FF000040"   # Color for mask overlay visualization (red with alpha)
  MAX_PLOT_SIZE <- 15   # Maximum dimension for any plot output

#=============================================================================
# Directory Setup
#=============================================================================

setup_directories <- function() {
  # Create main output directory
  main_dir <- file.path(ROOT_DIR, OUTPUT_DIR)
  dir.create(main_dir, showWarnings = FALSE)
  
  # Create subdirectories
  dir.create(file.path(main_dir, MORANS_DIR), showWarnings = FALSE)
  dir.create(file.path(main_dir, SCAN_DIR), showWarnings = FALSE)
  dir.create(file.path(main_dir, VISUAL_SCENE_DIR), showWarnings = FALSE)
  
  return(list(
    main = main_dir,
    morans = file.path(main_dir, MORANS_DIR),
    scan = file.path(main_dir, SCAN_DIR),
    visual_scene = file.path(main_dir, VISUAL_SCENE_DIR)
  ))
}

#=============================================================================
# Utility Functions
#=============================================================================
#' Create a study region by processing the input matrix and adding a buffer zone
#' @param matrix Input matrix to process
#' @param buffer_size Number of pixels to remove from edge (default BUFFER_SIZE)
#' @return Matrix with NAs outside study region and buffer zone
create_study_region <- function(matrix, buffer_size = BUFFER_SIZE) {
  # Convert values <= 0 to NA
  matrix[matrix <= 0] <- NA
  
  # Create buffer by identifying edge pixels and their neighbors
  n_rows <- nrow(matrix)
  n_cols <- ncol(matrix)
  
  # Create a copy to store original non-NA positions
  original_mask <- !is.na(matrix)
  
  # For each pixel in the buffer size
  for(i in 1:buffer_size) {
    # Find current edge pixels
    edge_pixels <- which(original_mask & 
                           (!original_mask[c(2:n_rows, n_rows), ] | # check above
                              !original_mask[c(1, 1:(n_rows-1)), ] | # check below
                              !original_mask[, c(2:n_cols, n_cols)] | # check right
                              !original_mask[, c(1, 1:(n_cols-1))]),  # check left
                         arr.ind = TRUE)
    
    # Remove edge pixels from mask
    if(nrow(edge_pixels) > 0) {
      original_mask[edge_pixels] <- FALSE
    }
  }
  
  # Apply buffered mask to matrix
  matrix[!original_mask] <- NA
  
  return(matrix)
}

#' Divide region into rectangular grid cells for analysis
#' @param matrix Input matrix with NAs outside study region
#' @param grid_size Number of cells per side
#' @return Matrix with numbered grid cells
create_grid <- function(matrix, grid_size) {
  # Get dimensions
  n_rows <- nrow(matrix)
  n_cols <- ncol(matrix)
  
  # Calculate cell dimensions
  cell_height <- ceiling(n_rows/grid_size)
  cell_width <- ceiling(n_cols/grid_size)
  
  # Create empty grid
  grid <- matrix(NA, n_rows, n_cols)
  
  # Fill grid with cell numbers
  cell_num <- 1
  for(i in seq(1, n_rows, by=cell_height)) {
    for(j in seq(1, n_cols, by=cell_width)) {
      # Define cell boundaries
      row_end <- min(i + cell_height - 1, n_rows)
      col_end <- min(j + cell_width - 1, n_cols)
      
      # Only number cells that overlap with study region
      if(any(!is.na(matrix[i:row_end, j:col_end]))) {
        grid[i:row_end, j:col_end] <- cell_num
        cell_num <- cell_num + 1
      }
    }
  }
  
  # Only keep grid cells within study region
  grid[is.na(matrix)] <- NA
  
  return(grid)
}

#' Find boundary points of cluster regions
#' @param grid Reference grid matrix
#' @param cluster_cells Vector of cell numbers in cluster
#' @return Matrix marking cluster boundaries
find_cluster_boundaries <- function(grid, cluster_cells) {
  boundary <- matrix(FALSE, nrow(grid), ncol(grid))
  
  # For each cell in cluster
  for(cell in cluster_cells) {
    cell_locs <- which(grid == cell, arr.ind = TRUE)
    
    # Check each point in the cell
    for(i in 1:nrow(cell_locs)) {
      row <- cell_locs[i,1]
      col <- cell_locs[i,2]
      
      # Check neighbors
      neighbors <- grid[max(1,row-1):min(nrow(grid),row+1),
                        max(1,col-1):min(ncol(grid),col+1)]
      
      # If any neighbor is NA or not in cluster, this is a boundary point
      if(any(is.na(neighbors)) || any(!is.na(neighbors) & !(neighbors %in% cluster_cells))) {
        boundary[row, col] <- TRUE
      }
    }
  }
  
  return(boundary)
}

#=============================================================================
# Moran's I Analysis
#=============================================================================
#' Calculate Moran's I statistic for spatial autocorrelation
#' @param matrix Input matrix of values
#' @param grid Reference grid matrix
#' @param nsim Number of Monte Carlo simulations
#' @return List containing Moran's I statistics and test results
calculate_morans_i <- function(matrix, grid, nsim = N_SIM) {
  # Get positions of valid (non-NA) cells and their values
  valid_positions <- which(!is.na(grid) & !is.na(matrix), arr.ind = TRUE)
  values <- matrix[valid_positions]
  
  # Check if we have enough data
  if(length(values) < 2) {
    warning("Not enough valid data points for Moran's I calculation")
    return(list(
      statistic = NA,
      p_value = NA,
      observed = NA,
      expected = NA,
      variance = NA,
      mc_results = NA,
      alternative = "greater",
      method = "Monte-Carlo simulation",
      data_name = "Normalized spatial data"
    ))
  }
  
  # Create neighbors list
  n_valid <- length(values)
  adj_list <- vector("list", n_valid)
  
  # For each valid cell, find its valid neighbors
  for(i in 1:n_valid) {
    current_row <- valid_positions[i, 1]
    current_col <- valid_positions[i, 2]
    
    neighbor_rows <- current_row + c(-1,-1,-1, 0,0, 1,1,1)
    neighbor_cols <- current_col + c(-1, 0, 1,-1,1,-1,0,1)
    
    valid_neighbors <- neighbor_rows >= 1 & 
      neighbor_rows <= nrow(grid) &
      neighbor_cols >= 1 & 
      neighbor_cols <= ncol(grid)
    
    neighbor_rows <- neighbor_rows[valid_neighbors]
    neighbor_cols <- neighbor_cols[valid_neighbors]
    
    neighbors <- integer(0)
    for(j in seq_along(neighbor_rows)) {
      match_idx <- which(valid_positions[,1] == neighbor_rows[j] & 
                           valid_positions[,2] == neighbor_cols[j])
      if(length(match_idx) > 0) {
        neighbors <- c(neighbors, as.integer(match_idx))
      }
    }
    
    adj_list[[i]] <- neighbors
  }
  
  # Set attributes for neighbor list
  class(adj_list) <- "nb"
  attr(adj_list, "region.id") <- as.character(1:n_valid)
  attr(adj_list, "call") <- match.call()
  attr(adj_list, "sym") <- TRUE
  attr(adj_list, "type") <- "queen"
  
  # Create weights and calculate Moran's I
  listw <- nb2listw(adj_list, style="W", zero.policy=TRUE)
  
  moran_result <- tryCatch({
    moran.test(values, listw, alternative="greater",
               zero.policy=TRUE, na.action=na.omit)
  }, error = function(e) {
    warning("Error in moran.test: ", e$message)
    return(NULL)
  })
  
  if(!is.null(moran_result)) {
    mc_result <- tryCatch({
      moran.mc(values, listw, nsim=nsim, zero.policy=TRUE)
    }, error = function(e) {
      warning("Error in moran.mc: ", e$message)
      return(NULL)
    })
  } else {
    mc_result <- NULL
  }
  
  if(!is.null(moran_result) && !is.null(mc_result)) {
    return(list(
      statistic = moran_result$statistic,
      p_value = mc_result$p.value,
      observed = moran_result$estimate[1],
      expected = moran_result$estimate[2],
      variance = moran_result$estimate[3],
      mc_results = mc_result,
      alternative = "greater",
      method = "Monte-Carlo simulation",
      data_name = "Normalized spatial data"
    ))
  } else {
    warning("Could not calculate Moran's I")
    return(list(
      statistic = NA,
      p_value = NA,
      observed = NA,
      expected = NA,
      variance = NA,
      mc_results = NA,
      alternative = "greater",
      method = "Monte-Carlo simulation",
      data_name = "Normalized spatial data"
    ))
  }
}

run_morans_analysis <- function(tif_stack, grid, output_dir) {
  results <- data.frame(
    Layer = integer(),
    MoransI = numeric(),
    Expected = numeric(),
    PValue = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(i in seq_along(tif_stack)) {
    current_layer <- tif_stack[[i]]
    processed_layer <- create_study_region(current_layer)
    
    # Normalize values
    valid_values <- processed_layer[!is.na(processed_layer)]
    if(length(valid_values) > 0) {
      min_val <- min(valid_values)
      max_val <- max(valid_values)
      if(max_val > min_val) {
        processed_layer[!is.na(processed_layer)] <- 
          (processed_layer[!is.na(processed_layer)] - min_val) / (max_val - min_val)
      }
    }
    
    moran_result <- calculate_morans_i(processed_layer, grid)
    
    results <- rbind(results, data.frame(
      Layer = i,
      MoransI = moran_result$observed,
      Expected = moran_result$expected,
      PValue = moran_result$p_value
    ))
  }
  
  # Save results
  write.csv(results, file.path(output_dir, "morans_results.csv"), row.names = FALSE)
  
  return(results)
}

#=============================================================================
# Scan Statistics Analysis
#=============================================================================
#' Prepare data for scan statistics analysis
#' @param matrix Value matrix
#' @param grid Matrix of grid cells
#' @param normalize Boolean for value normalization
#' @param invert Boolean for inverting normalized values
#' @return List with coordinates and values for scan statistics
prepare_scan_data <- function(matrix, grid, normalize = NORMALIZE, invert = INVERT) {
  if(invert) normalize <- TRUE
  
  data <- data.frame(
    cell = as.vector(grid),
    value = as.vector(matrix)
  ) %>%
    filter(!is.na(cell), value > 0)
  
  if(normalize) {
    min_val <- min(data$value)
    max_val <- max(data$value)
    
    if(max_val > min_val) {
      data$value <- (data$value - min_val) / (max_val - min_val)
    }
    
    if(invert) {
      data$value <- 1 - data$value
    }
  }
  
  centroids <- data.frame(cell = unique(data$cell)) %>%
    rowwise() %>%
    mutate(
      indices = list(which(grid == cell, arr.ind = TRUE)),
      x = mean(indices[,2]),
      y = mean(indices[,1])
    ) %>%
    select(cell, x, y)
  
  values <- data %>%
    group_by(cell) %>%
    summarise(
      total_value = sum(value, na.rm=TRUE),
      count = n()
    )
  
  return(list(
    coords = as.matrix(centroids[,c("x", "y")]),
    values = values$total_value,
    pop = values$count
  ))
}

#' Run scan statistics analysis and save results
#' @param tif_stack List of matrices from TIF file
#' @param grid Reference grid matrix
#' @param output_dir Directory to save results
#' @return List containing results, plots, and grid
run_scan_analysis <- function(tif_stack, grid, output_dir) {
  results <- list()
  plots <- list()
  
  for(i in seq_along(tif_stack)) {
    scan_data <- prepare_scan_data(tif_stack[[i]], grid)
    
    scaling_factor <- if(NORMALIZE) 10000 else 1000
    
    scan_results <- scan.test(
      coords = scan_data$coords,
      cases = scan_data$values * scaling_factor,
      pop = rep(scaling_factor, length(scan_data$values)),
      type = "poisson",
      simdist = "poisson",
      nsim = N_SIM,
      alpha = ALPHA,
      ubpop = UBPOP,
      min.cases = MIN_CASES
    )
    
    results[[i]] <- scan_results
    
    # Create visualization
    plot <- visualize_scan_results(tif_stack[[i]], grid, scan_results, i)
    if (!is.null(plot)) {
      plots[[i]] <- plot
      
      ggsave(
        filename = file.path(output_dir, paste0(i, "_scan_results.png")),
        plot = plots[[i]],
        width = 8,
        height = 8,
        dpi = 300
      )
    }
  }
  
  plots <- plots[!sapply(plots, is.null)]
  saveRDS(results, file = file.path(output_dir, "scan_results.rds"))
  
  # Create portable results directory and export results
  portable_dir <- file.path(output_dir, "ScanPortable")
  dir.create(portable_dir, showWarnings = FALSE)
  
  # Process each layer's results for portable format
  for(i in seq_along(results)) {
    p_value_matrix <- matrix(0, nrow=nrow(grid), ncol=ncol(grid))
    
    current_result <- results[[i]]
    if(!is.null(current_result$clusters)) {
      for(j in seq_along(current_result$clusters)) {
        cluster <- current_result$clusters[[j]]
        if(cluster$pvalue <= 0.05) {
          cluster_cells <- cluster$locids
          cluster_positions <- which(grid %in% cluster_cells, arr.ind = TRUE)
          p_value_matrix[cluster_positions] <- cluster$pvalue
        }
      }
    }
    
    filename <- file.path(portable_dir, sprintf("layer_%d_pvalues.csv", i))
    write.csv(p_value_matrix, filename, row.names = FALSE)
  }
  
  # Create readme file
  readme_text <- paste(
    "Scan Statistics Results - P-value Matrices",
    "",
    "File Format:",
    "- Each CSV file corresponds to one layer from the original analysis",
    "- Files are named 'layer_X_pvalues.csv' where X is the layer number",
    "- Matrix values:",
    "  * 0: No significant cluster or outside study region",
    "  * Other values: P-value of significant cluster (p <= 0.05) at that position",
    "",
    "Created:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    sep = "\n"
  )
  writeLines(readme_text, file.path(portable_dir, "README.txt"))
  
  return(list(
    results = results,
    plots = plots,
    grid = grid
  ))
}

#' Visualize scan statistics results
#' @param matrix Original value matrix
#' @param grid Grid cell matrix
#' @param scan_results Scan test results
#' @param layer_name Name/number of current layer
#' @return ggplot object or NULL if no significant clusters
visualize_scan_results <- function(matrix, grid, scan_results, layer_name) {
  plot_data <- expand.grid(
    row = 1:nrow(matrix),
    col = 1:ncol(matrix)
  ) %>%
    mutate(
      value = as.vector(matrix)
    )
  
  plot_data$value[plot_data$value <= 0] <- NA
  
  sig_clusters <- smerc::clusters(scan_results)
  if (is.null(sig_clusters)) {
    message(paste("No significant clusters found for layer", layer_name))
    return(NULL)
  }
  
  cluster_mask <- matrix(NA, nrow(matrix), ncol(matrix))
  for(i in seq_along(sig_clusters)) {
    cluster_cells <- sig_clusters[[i]]
    cells_in_cluster <- which(grid %in% cluster_cells, arr.ind = TRUE)
    cluster_mask[cells_in_cluster] <- matrix[cells_in_cluster]
  }
  
  boundary_points <- data.frame()
  for(i in seq_along(sig_clusters)) {
    cluster_cells <- sig_clusters[[i]]
    boundary <- find_cluster_boundaries(grid, cluster_cells)
    if(any(boundary)) {
      boundary_points <- rbind(
        boundary_points,
        data.frame(
          row = which(boundary, arr.ind = TRUE)[,1],
          col = which(boundary, arr.ind = TRUE)[,2],
          cluster = i
        )
      )
    }
  }
  
  plot_data_clusters <- plot_data %>%
    mutate(value = as.vector(cluster_mask))
  
  num_clusters <- length(sig_clusters)
  display_clusters <- min(5, num_clusters)
  
  p_values_text <- paste(
    sapply(1:display_clusters, 
           function(i) sprintf("Cluster %d: p = %.3f", 
                               i, 
                               scan_results$clusters[[i]]$pvalue)),
    collapse = "\n"
  )
  
  if(num_clusters > 5) {
    p_values_text <- paste(p_values_text, 
                           sprintf("\n(+ %d more clusters)", num_clusters - 5))
  }
  
  p1 <- ggplot(plot_data, aes(x=col, y=row)) +
    geom_tile(aes(fill=value)) +
    scale_fill_viridis_c(na.value="white") +
    coord_equal() +
    theme_minimal() +
    labs(x = "X", y = "Y", 
         subtitle=paste("Full data with clusters outlined\n", p_values_text))
  
  p2 <- ggplot(plot_data_clusters, aes(x=col, y=row)) +
    geom_tile(aes(fill=value)) +
    scale_fill_viridis_c(na.value="white") +
    coord_equal() +
    theme_minimal() +
    labs(x = "X", y = "Y", subtitle="Significant clusters only")
  
  if(nrow(boundary_points) > 0) {
    p1 <- p1 + 
      geom_point(data=boundary_points, color="red", size=0.5) +
      geom_text(data=boundary_points %>% 
                  group_by(cluster) %>% 
                  summarize(row=mean(row), 
                            col=mean(col), 
                            cluster=first(cluster)),
                aes(label=cluster),
                color="black", size=3)
  }
  
  combined_plot <- p1 + p2 + 
    plot_layout(ncol=2) +
    plot_annotation(
      title = paste("Cluster:", layer_name),
      theme = theme(plot.title = element_text(hjust = 0.5))
    )
  
  return(combined_plot)
}

#=============================================================================
# Visual Scene Analysis
#=============================================================================
#' Create anatomical mask for retinal analysis
#' @param grid Reference grid matrix
#' @param n_annuli Number of concentric rings
#' @param n_sectors Number of angular sectors
#' @param temporal_left Boolean for temporal region orientation
#' @param slope Degrees per annulus
#' @param sector_fill_0 Start angle
#' @param sector_fill_1 End angle
#' @param outer_fill Outer annuli to fill
#' @param inner_fill Inner annuli to fill
#' @return Mask matrix
create_retinal_mask <- function(grid, n_annuli = N_ANNULI, n_sectors = N_SECTORS, 
                                temporal_left = TEMPORAL_LEFT, slope = SLOPE,
                                sector_fill_0 = SECTOR_FILL_0, sector_fill_1 = SECTOR_FILL_1,
                                outer_fill = OUTER_FILL, inner_fill = INNER_FILL) {
  n_rows <- nrow(grid)
  n_cols <- ncol(grid)
  center_row <- n_rows/2
  center_col <- n_cols/2
  
  mask <- matrix(0, n_rows, n_cols)
  
  valid_positions <- which(!is.na(grid), arr.ind = TRUE)
  distances <- sqrt((valid_positions[,1] - center_row)^2 + 
                      (valid_positions[,2] - center_col)^2)
  max_radius <- max(distances)
  
  for(i in 1:n_rows) {
    for(j in 1:n_cols) {
      if(!is.na(grid[i,j])) {
        dy <- i - center_row
        dx <- j - center_col
        
        radius <- sqrt(dx^2 + dy^2)
        norm_radius <- radius/max_radius
        annulus_idx <- floor(norm_radius * n_annuli)
        
        angle <- (atan2(dy, dx) * 180/pi + 90) %% 360
        
        if(temporal_left) {
          angle <- (360 - angle) %% 360
        }
        
        if(annulus_idx >= (n_annuli - outer_fill) && annulus_idx < n_annuli) {
          depth_from_edge <- (n_annuli - 1) - annulus_idx
          
          adjusted_fill_0 <- sector_fill_0 + (slope * depth_from_edge)
          adjusted_fill_1 <- sector_fill_1 - (slope * depth_from_edge)
          
          angle_in_range <- FALSE
          if(adjusted_fill_0 <= adjusted_fill_1) {
            angle_in_range <- angle >= adjusted_fill_0 && angle <= adjusted_fill_1
          } else {
            angle_in_range <- angle >= adjusted_fill_0 || angle <= adjusted_fill_1
          }
          
          if(angle_in_range) {
            mask[i,j] <- 1
          }
        } 
        else if(annulus_idx < inner_fill) {
          angle_in_range <- FALSE
          if(sector_fill_0 <= sector_fill_1) {
            angle_in_range <- angle >= sector_fill_0 && angle <= sector_fill_1
          } else {
            angle_in_range <- angle >= sector_fill_0 || angle <= sector_fill_1
          }
          
          if(angle_in_range) {
            mask[i,j] <- 1
          }
        }
      }
    }
  }
  
  return(mask)
}

#' Create standard set of visual scene masks
#' @param grid Reference grid matrix
#' @return List of mask matrices
create_visual_scene_masks <- function(grid) {
  binocular_mask <- create_retinal_mask(grid, n_annuli = 6, n_sectors = 24, 
                                        temporal_left = T, slope = 15,
                                        sector_fill_0 = 140, sector_fill_1 = 315,
                                        outer_fill = 2, inner_fill = 0)
  
  peripheral_mask <- create_retinal_mask(grid, n_annuli = 1, n_sectors = 1, 
                                         temporal_left = T, slope = 0,
                                         sector_fill_0 = 0, sector_fill_1 = 365,
                                         outer_fill = 1, inner_fill = 0) - binocular_mask
  
  visual_sky_mask <- create_retinal_mask(grid, n_annuli = 30, n_sectors = 24, 
                                         temporal_left = T, slope = 3,
                                         sector_fill_0 = 270, sector_fill_1 = 90,
                                         outer_fill = 20, inner_fill = 0)
  
  visual_ground_mask <- create_retinal_mask(grid, n_annuli = 1, n_sectors = 1, 
                                            temporal_left = T, slope = 0,
                                            sector_fill_0 = 0, sector_fill_1 = 365,
                                            outer_fill = 1, inner_fill = 0) - visual_sky_mask
  
  binocular_sky_mask <- visual_sky_mask * binocular_mask
  binocular_ground_mask <- visual_ground_mask * binocular_mask
  peripheral_sky_mask <- visual_sky_mask * peripheral_mask
  peripheral_ground_mask <- visual_ground_mask * peripheral_mask
  
  mask_list <- list(
    binocular = binocular_mask,
    peripheral = peripheral_mask,
    visual_sky = visual_sky_mask,
    visual_ground = visual_ground_mask,
    binocular_sky = binocular_sky_mask,
    binocular_ground = binocular_ground_mask,
    peripheral_sky = peripheral_sky_mask,
    peripheral_ground = peripheral_ground_mask
  )
  
  return(list(
    binocular = binocular_mask,
    peripheral = peripheral_mask,
    visual_sky = visual_sky_mask,
    visual_ground = visual_ground_mask,
    binocular_sky = binocular_sky_mask,
    binocular_ground = binocular_ground_mask,
    peripheral_sky = peripheral_sky_mask,
    peripheral_ground = peripheral_ground_mask
  ))
}

#' Generate bootstrapped version of cluster mask with rotation and flips only
#' @param cluster_mask Original cluster mask matrix
#' @param grid Reference grid matrix
#' @param study_region Boolean matrix of valid study region
#' @return Transformed cluster mask matrix
bootstrap_cluster_mask <- function(cluster_mask, grid, study_region) {
  # Get dimensions
  n_rows <- nrow(grid)
  n_cols <- ncol(grid)
  center_row <- n_rows/2
  center_col <- n_cols/2
  
  # Store original mask for boundary checks
  original_mask <- cluster_mask
  
  # 1. Random rotation
  angle <- runif(1, 0, 360)
  if(angle > 0) {
    theta <- angle * pi / 180
    rotated_mask <- matrix(0, n_rows, n_cols)
    
    # Track if any valid points were placed
    valid_points <- FALSE
    
    for(i in 1:n_rows) {
      for(j in 1:n_cols) {
        if(cluster_mask[i,j] == 1) {
          x <- j - center_col
          y <- i - center_row
          
          x_rot <- x * cos(theta) - y * sin(theta)
          y_rot <- x * sin(theta) + y * cos(theta)
          
          new_col <- round(x_rot + center_col)
          new_row <- round(y_rot + center_row)
          
          if(new_row >= 1 && new_row <= n_rows && 
             new_col >= 1 && new_col <= n_cols && 
             study_region[new_row, new_col]) {
            rotated_mask[new_row, new_col] <- 1
            valid_points <- TRUE
          }
        }
      }
    }
    
    # Only update if we placed valid points
    if(valid_points) {
      cluster_mask <- rotated_mask
    }
  }
  
  # 2. Random flips
  if(runif(1) > 0.5) {
    cluster_mask <- cluster_mask[n_rows:1,] # Vertical flip
  }
  if(runif(1) > 0.5) {
    cluster_mask <- cluster_mask[,n_cols:1] # Horizontal flip
  }
  
  return(cluster_mask)
}

#' Summarize bootstrap distribution statistics
#' @param bootstrap_distributions List of bootstrap distributions per mask
#' @param confidence_levels Vector of confidence levels
#' @return Data frame of summary statistics
summarize_bootstrap_results <- function(bootstrap_distributions, confidence_levels = c(0.95, 0.99)) {
  summary_stats <- data.frame(
    mask = names(bootstrap_distributions),
    mean = sapply(bootstrap_distributions, mean),
    median = sapply(bootstrap_distributions, median)
  )
  
  # Add confidence intervals
  for(level in confidence_levels) {
    alpha <- 1 - level
    ci_lower <- sapply(bootstrap_distributions, quantile, probs = alpha/2)
    ci_upper <- sapply(bootstrap_distributions, quantile, probs = 1 - alpha/2)
    
    summary_stats[[paste0("ci_", level, "_lower")]] <- ci_lower
    summary_stats[[paste0("ci_", level, "_upper")]] <- ci_upper
  }
  
  return(summary_stats)
}

#' Generate bootstrapped version of cluster mask with rotation and flips
#' @param cluster_mask Original cluster mask matrix
#' @param grid Reference grid matrix
#' @param study_region Boolean matrix of valid study region
#' @param method Bootstrap method ("spin", "normal_draw", or "pixel_shuffle")
#' @param stats Optional statistics for normal_draw method
#' @return Transformed cluster mask matrix
bootstrap_cluster_mask <- function(cluster_mask, grid, study_region, 
                                   method = "spin", stats = NULL) {
  if(method == "spin") {
    # Original spin method with rotation and flips
    n_rows <- nrow(grid)
    n_cols <- ncol(grid)
    center_row <- n_rows/2
    center_col <- n_cols/2
    
    # Store original mask
    original_mask <- cluster_mask
    
    # Random rotation
    angle <- runif(1, 0, 360)
    if(angle > 0) {
      theta <- angle * pi / 180
      rotated_mask <- matrix(0, n_rows, n_cols)
      valid_points <- FALSE
      
      for(i in 1:n_rows) {
        for(j in 1:n_cols) {
          if(cluster_mask[i,j] == 1) {
            x <- j - center_col
            y <- i - center_row
            
            x_rot <- x * cos(theta) - y * sin(theta)
            y_rot <- x * sin(theta) + y * cos(theta)
            
            new_col <- round(x_rot + center_col)
            new_row <- round(y_rot + center_row)
            
            if(new_row >= 1 && new_row <= n_rows && 
               new_col >= 1 && new_col <= n_cols && 
               study_region[new_row, new_col]) {
              rotated_mask[new_row, new_col] <- 1
              valid_points <- TRUE
            }
          }
        }
      }
      
      if(valid_points) {
        cluster_mask <- rotated_mask
      }
    }
    
    # Random flips
    if(runif(1) > 0.5) {
      cluster_mask <- cluster_mask[n_rows:1,] # Vertical flip
    }
    if(runif(1) > 0.5) {
      cluster_mask <- cluster_mask[,n_cols:1] # Horizontal flip
    }
    
  } else if(method == "normal_draw") {
    # Draw from normal distribution based on original statistics
    if(is.null(stats)) {
      stop("Stats required for normal_draw method")
    }
    
    # Generate values from normal distribution
    new_values <- rnorm(length(cluster_mask[!is.na(grid)]),
                        mean = stats$mean,
                        sd = stats$sd)
    
    # Create new mask matching original density
    cluster_mask <- matrix(0, nrow(grid), ncol(grid))
    valid_positions <- which(!is.na(grid))
    threshold <- quantile(new_values, 1 - mean(original_mask, na.rm=TRUE))
    cluster_mask[valid_positions] <- as.numeric(new_values > threshold)
    
  } else if(method == "pixel_shuffle") {
    # Randomly shuffle significant pixels while maintaining total number
    n_pixels <- sum(cluster_mask == 1, na.rm = TRUE)
    valid_positions <- which(!is.na(study_region))
    
    # Create new empty mask
    new_mask <- matrix(0, nrow(grid), ncol(grid))
    
    # Randomly select positions for significant pixels
    selected_positions <- sample(valid_positions, n_pixels)
    new_mask[selected_positions] <- 1
    
    cluster_mask <- new_mask
  }
  
  return(cluster_mask)
}

#' Generate set of bootstrap masks using different methods
#' @param scan_results List of scan test results
#' @param grid Reference grid matrix
#' @param study_region Boolean matrix of valid study region
#' @param n_bootstrap Number of bootstrap iterations
#' @param method Bootstrap method ("spin", "gaussian", or "shuffle")
#' @return List of bootstrap mask matrices
generate_bootstrap_masks <- function(scan_results, grid, study_region, 
                                     n_bootstrap = N_BOOTSTRAP, 
                                     method = "spin") {
  message("Collecting significant cluster masks...")
  all_cluster_masks <- list()
  cluster_areas <- numeric()
  mask_counter <- 1
  
  # First collect all significant cluster masks
  for(i in seq_along(scan_results)) {
    current_result <- scan_results[[i]]
    if(!is.null(current_result$clusters)) {
      coverage_matrix <- matrix(0, nrow(grid), ncol(grid))
      for(cluster in current_result$clusters) {
        if(cluster$pvalue <= 0.05) {
          cluster_cells <- cluster$locids
          cells_in_cluster <- which(grid %in% cluster_cells, arr.ind = TRUE)
          coverage_matrix[cells_in_cluster] <- 1
        }
      }
      if(any(coverage_matrix == 1)) {
        all_cluster_masks[[mask_counter]] <- coverage_matrix
        cluster_areas[mask_counter] <- sum(coverage_matrix)
        mask_counter <- mask_counter + 1
      }
    }
  }
  
  if(length(all_cluster_masks) == 0) {
    stop("No significant clusters found for bootstrapping")
  }
  
  # Setup progress bar
  pb <- txtProgressBar(min = 0, max = n_bootstrap, style = 3)
  
  # Generate bootstrap masks based on method
  bootstrap_masks <- list()
  if(method == "shuffle") {
    message("\nUsing pixel shuffle method...")
    for(i in 1:n_bootstrap) {
      # Randomly select a real mask
      selected_mask_idx <- sample(1:length(all_cluster_masks), 1)
      original_mask <- all_cluster_masks[[selected_mask_idx]]
      
      # Get valid positions from study region
      valid_positions <- which(!is.na(study_region))
      
      # Count number of significant pixels in original mask
      n_significant <- sum(original_mask == 1, na.rm = TRUE)
      
      # Create new empty mask
      shuffled_mask <- matrix(0, nrow(grid), ncol(grid))
      
      # Randomly select positions for significant pixels
      new_positions <- sample(valid_positions, n_significant)
      shuffled_mask[new_positions] <- 1
      
      bootstrap_masks[[i]] <- shuffled_mask
      
      setTxtProgressBar(pb, i)
      if(i %% max(1, round(n_bootstrap/10)) == 0) {
        message(sprintf("\nCompleted %d%% of bootstrap iterations", 
                        round(i/n_bootstrap * 100)))
      }
    }
  } else {
    # Calculate frequencies of each cluster area for spin method
    area_table <- table(cluster_areas)
    area_prob <- as.numeric(area_table)/length(cluster_areas)
    unique_areas <- as.numeric(names(area_table))
    
    # Create lookup table for masks by area
    masks_by_area <- lapply(unique_areas, function(area) {
      which(cluster_areas == area)
    })
    names(masks_by_area) <- unique_areas
    
    message(sprintf("\nFound %d significant cluster masks with %d unique areas", 
                    length(all_cluster_masks), length(unique_areas)))
    
    for(i in 1:n_bootstrap) {
      selected_area <- sample(unique_areas, 1, prob = area_prob)
      available_masks <- masks_by_area[[as.character(selected_area)]]
      selected_mask <- all_cluster_masks[[sample(available_masks, 1)]]
      
      bootstrap_masks[[i]] <- bootstrap_cluster_mask(
        selected_mask, grid, study_region, method = method
      )
      
      setTxtProgressBar(pb, i)
      if(i %% max(1, round(n_bootstrap/10)) == 0) {
        message(sprintf("\nCompleted %d%% of bootstrap iterations", 
                        round(i/n_bootstrap * 100)))
      }
    }
  }
  
  close(pb)
  message("\nBootstrap mask generation complete.")
  
  # Verify area distribution
  final_areas <- sapply(bootstrap_masks, sum)
  message("\nArea distribution verification:")
  message("Original distribution:")
  print(table(cluster_areas)/length(cluster_areas))
  message("Bootstrap distribution:")
  print(table(final_areas)/length(final_areas))
  
  return(bootstrap_masks)
}

#' Shuffle a binary mask while preserving cluster size
#' @param mask Original binary mask matrix
#' @param study_region Boolean matrix of valid study region
#' @return Shuffled binary mask matrix
shuffle_mask <- function(mask, study_region) {
  # Get valid positions from study region
  valid_positions <- which(!is.na(study_region))
  
  # Count number of significant pixels in original mask
  n_sig_pixels <- sum(mask == 1, na.rm = TRUE)
  
  # Create new empty mask
  shuffled_mask <- matrix(0, nrow(mask), ncol(mask))
  
  # Randomly select positions for significant pixels
  new_positions <- sample(valid_positions, n_sig_pixels)
  shuffled_mask[new_positions] <- 1
  
  return(shuffled_mask)
}

#' Generate set of bootstrap masks using pixel shuffling
#' @param scan_results List of scan test results
#' @param grid Reference grid matrix
#' @param study_region Boolean matrix of valid study region
#' @param n_bootstrap Number of bootstrap iterations
#' @return List of bootstrap mask matrices
generate_shuffle_masks <- function(scan_results, grid, study_region, n_bootstrap) {
  # Collect all significant cluster masks
  all_cluster_masks <- list()
  mask_counter <- 1
  
  message("Collecting significant cluster masks...")
  for(i in seq_along(scan_results)) {
    current_result <- scan_results[[i]]
    if(!is.null(current_result$clusters)) {
      coverage_matrix <- matrix(0, nrow(grid), ncol(grid))
      has_significant <- FALSE
      
      for(cluster in current_result$clusters) {
        if(cluster$pvalue <= 0.05) {
          cluster_cells <- cluster$locids
          cells_in_cluster <- which(grid %in% cluster_cells, arr.ind = TRUE)
          coverage_matrix[cells_in_cluster] <- 1
          has_significant <- TRUE
        }
      }
      
      if(has_significant) {
        all_cluster_masks[[mask_counter]] <- coverage_matrix
        mask_counter <- mask_counter + 1
      }
    }
  }
  
  if(length(all_cluster_masks) == 0) {
    stop("No significant clusters found for bootstrapping")
  }
  
  message(sprintf("Found %d masks with significant clusters", length(all_cluster_masks)))
  
  # Setup progress bar
  pb <- txtProgressBar(min = 0, max = n_bootstrap, style = 3)
  
  # Generate bootstrap masks by randomly selecting and shuffling
  bootstrap_masks <- list()
  for(i in 1:n_bootstrap) {
    # Randomly select a mask
    selected_mask <- all_cluster_masks[[sample(length(all_cluster_masks), 1)]]
    
    # Shuffle the selected mask
    bootstrap_masks[[i]] <- shuffle_mask(selected_mask, study_region)
    
    setTxtProgressBar(pb, i)
    
    if(i %% max(1, round(n_bootstrap/10)) == 0) {
      message(sprintf("\nCompleted %d%% of bootstrap iterations", round(i/n_bootstrap * 100)))
    }
  }
  
  close(pb)
  message("\nBootstrap mask generation complete.")
  
  return(bootstrap_masks)
}
#' Perform bootstrap analysis of coverage using Gaussian sampling
#' @param coverage_data Numeric vector of coverage percentages
#' @param n_bootstrap Number of bootstrap samples
#' @param conf_levels Vector of confidence levels
#' @return List containing confidence intervals and distribution
bootstrap_coverage_gaussian <- function(coverage_data, n_bootstrap = 10000, 
                                        conf_levels = c(0.95, 0.99)) {
  # Calculate mean and sd of original data
  orig_mean <- mean(coverage_data)
  orig_sd <- sd(coverage_data)
  
  # Generate bootstrap samples from Gaussian
  bootstrap_samples <- matrix(
    rnorm(n_bootstrap * length(coverage_data), 
          mean = orig_mean, 
          sd = orig_sd),
    nrow = n_bootstrap
  )
  
  # Calculate mean for each bootstrap sample
  bootstrap_means <- rowMeans(bootstrap_samples)
  
  # Calculate confidence intervals for each confidence level
  intervals <- list()
  for(level in conf_levels) {
    alpha <- 1 - level
    intervals[[as.character(level)]] <- quantile(bootstrap_means, 
                                                 probs = c(alpha/2, 1 - alpha/2))
  }
  
  return(list(
    intervals = intervals,
    mean = orig_mean,
    sd = orig_sd,
    bootstrap_dist = bootstrap_means
  ))
}

#' Analyze visual scene coverage with choice of bootstrap method
#' @param scan_results List of scan test results
#' @param grid Reference grid matrix
#' @param mask_list List of mask matrices
#' @param n_bootstrap Number of bootstrap iterations
#' @param bootstrap_method Character: "gaussian", "spin", or "shuffle"
#' @param output_dir Output directory path
#' @return List containing results and visualizations
analyze_visual_scene_coverage <- function(scan_results, grid, mask_list, 
                                          n_bootstrap = N_BOOTSTRAP,
                                          bootstrap_method = "gaussian",
                                          output_dir) {
  # Create bootstrap subdirectory
  bootstrap_dir <- file.path(output_dir, "Bootstrap")
  dir.create(bootstrap_dir, showWarnings = FALSE)
  
  # Get total valid area
  total_area <- sum(!is.na(grid))
  study_region <- !is.na(grid)
  
  # Initialize results storage
  results_cols <- c("Layer", "TotalClusters", "CoveredArea", "PercentCoverage")
  for(mask_name in names(mask_list)) {
    results_cols <- c(results_cols, 
                      paste0(mask_name, "_Area"),
                      paste0(mask_name, "_Percent"))
  }
  
  coverage_results <- as.data.frame(matrix(0, 0, length(results_cols)))
  names(coverage_results) <- results_cols
  
  # Process each layer for actual data
  message("Processing actual data layers...")
  for(i in seq_along(scan_results)) {
    current_result <- scan_results[[i]]
    current_row <- list(
      Layer = i,
      TotalClusters = 0,
      CoveredArea = 0,
      PercentCoverage = 0
    )
    
    if(!is.null(current_result$clusters)) {
      coverage_matrix <- matrix(FALSE, nrow(grid), ncol(grid))
      cluster_sizes <- numeric()
      
      for(cluster in current_result$clusters) {
        if(cluster$pvalue <= 0.05) {
          cluster_cells <- cluster$locids
          cells_in_cluster <- which(grid %in% cluster_cells, arr.ind = TRUE)
          coverage_matrix[cells_in_cluster] <- TRUE
          cluster_sizes <- c(cluster_sizes, length(cluster_cells))
        }
      }
      
      covered_area <- sum(coverage_matrix)
      current_row$TotalClusters <- length(cluster_sizes)
      current_row$CoveredArea <- covered_area
      current_row$PercentCoverage <- (covered_area / total_area) * 100
      
      for(mask_name in names(mask_list)) {
        mask <- mask_list[[mask_name]]
        masked_area <- sum(coverage_matrix * mask)
        masked_percent <- if(covered_area > 0) (masked_area/covered_area) * 100 else 0
        
        current_row[[paste0(mask_name, "_Area")]] <- masked_area
        current_row[[paste0(mask_name, "_Percent")]] <- masked_percent
      }
    }
    coverage_results <- rbind(coverage_results, as.data.frame(current_row))
  }
  
  # Initialize bootstrap results storage
  bootstrap_coverage <- list()
  for(mask_name in names(mask_list)) {
    bootstrap_coverage[[mask_name]] <- numeric(n_bootstrap)
  }
  
  bootstrap_summary <- data.frame(
    mask = character(),
    mean = numeric(),
    median = numeric(),
    ci_0.95_lower = numeric(),
    ci_0.95_upper = numeric(),
    ci_0.99_lower = numeric(),
    ci_0.99_upper = numeric(),
    stringsAsFactors = FALSE
  )
  
  message(sprintf("Performing bootstrap analysis using %s method...", bootstrap_method))
  
  if(bootstrap_method %in% c("spin", "shuffle")) {
    # Generate bootstrap masks using appropriate method
    bootstrap_masks <- if(bootstrap_method == "spin") {
      generate_bootstrap_masks(scan_results, grid, study_region, n_bootstrap)
    } else {
      generate_shuffle_masks(scan_results, grid, study_region, n_bootstrap)
    }
    
    # Save bootstrap masks
    dir.create(file.path(bootstrap_dir, "Masks"), showWarnings = FALSE)
    for(i in seq_along(bootstrap_masks)) {
      write.csv(bootstrap_masks[[i]], 
                file = file.path(bootstrap_dir, "Masks", 
                                 sprintf("bootstrap_mask_%04d.csv", i)),
                row.names = FALSE)
    }
    
    # Calculate coverage for each mask and visual scene
    for(mask_name in names(mask_list)) {
      message(sprintf("Processing mask: %s", mask_name))
      
      for(i in seq_along(bootstrap_masks)) {
        masked_area <- sum(bootstrap_masks[[i]] * mask_list[[mask_name]], na.rm = TRUE)
        total_cluster_area <- sum(bootstrap_masks[[i]], na.rm = TRUE)
        bootstrap_coverage[[mask_name]][i] <- if(total_cluster_area > 0) {
          (masked_area / total_cluster_area) * 100
        } else {
          0
        }
      }
      
      # Calculate summary statistics
      ci_95 <- quantile(bootstrap_coverage[[mask_name]], c(0.025, 0.975))
      ci_99 <- quantile(bootstrap_coverage[[mask_name]], c(0.005, 0.995))
      
      bootstrap_summary <- rbind(bootstrap_summary, data.frame(
        mask = mask_name,
        mean = mean(bootstrap_coverage[[mask_name]]),
        median = median(bootstrap_coverage[[mask_name]]),
        ci_0.95_lower = ci_95[1],
        ci_0.95_upper = ci_95[2],
        ci_0.99_lower = ci_99[1],
        ci_0.99_upper = ci_99[2]
      ))
    }
  } else if(bootstrap_method == "gaussian") {
    # Gaussian bootstrap for each mask
    for(mask_name in names(mask_list)) {
      percent_col <- paste0(mask_name, "_Percent")
      coverage_data <- coverage_results[[percent_col]]
      
      if(all(coverage_data == 0)) next
      
      # Generate bootstrap samples
      n_samples <- length(coverage_data)
      bootstrap_samples <- matrix(
        rnorm(n_bootstrap * n_samples, 
              mean = mean(coverage_data), 
              sd = sd(coverage_data)),
        nrow = n_bootstrap
      )
      bootstrap_coverage[[mask_name]] <- rowMeans(bootstrap_samples)
      
      # Calculate summary statistics
      ci_95 <- quantile(bootstrap_coverage[[mask_name]], c(0.025, 0.975))
      ci_99 <- quantile(bootstrap_coverage[[mask_name]], c(0.005, 0.995))
      
      bootstrap_summary <- rbind(bootstrap_summary, data.frame(
        mask = mask_name,
        mean = mean(bootstrap_coverage[[mask_name]]),
        median = median(bootstrap_coverage[[mask_name]]),
        ci_0.95_lower = ci_95[1],
        ci_0.95_upper = ci_95[2],
        ci_0.99_lower = ci_99[1],
        ci_0.99_upper = ci_99[2]
      ))
    }
  }
  
  # Save bootstrap results
  bootstrap_results_df <- data.frame(
    iteration = 1:n_bootstrap
  )
  for(mask_name in names(bootstrap_coverage)) {
    bootstrap_results_df[[paste0(mask_name, "_coverage")]] <- bootstrap_coverage[[mask_name]]
  }
  write.csv(bootstrap_results_df, 
            file = file.path(bootstrap_dir, "bootstrap_coverage_all.csv"),
            row.names = FALSE)
  write.csv(bootstrap_summary,
            file = file.path(bootstrap_dir, "bootstrap_summary_stats.csv"),
            row.names = FALSE)
  
  # Classify results using bootstrap distributions
  classification_results <- data.frame(Layer = coverage_results$Layer)
  for(mask_name in names(mask_list)) {
    percent_col <- paste0(mask_name, "_Percent")
    coverage_data <- coverage_results[[percent_col]]
    
    if(all(coverage_data == 0)) next
    
    # Get confidence intervals from summary
    ci_95 <- c(
      bootstrap_summary$ci_0.95_lower[bootstrap_summary$mask == mask_name],
      bootstrap_summary$ci_0.95_upper[bootstrap_summary$mask == mask_name]
    )
    ci_99 <- c(
      bootstrap_summary$ci_0.99_lower[bootstrap_summary$mask == mask_name],
      bootstrap_summary$ci_0.99_upper[bootstrap_summary$mask == mask_name]
    )
    
    # Classify based on confidence intervals
    classification <- rep("Typical distribution", length(coverage_data))
    classification[coverage_data > ci_99[2]] <- "Very high"
    classification[coverage_data > ci_95[2] & coverage_data <= ci_99[2]] <- "High"
    classification[coverage_data < ci_99[1]] <- "Very low"
    classification[coverage_data < ci_95[1] & coverage_data >= ci_99[1]] <- "Low"
    
    classification_results[[mask_name]] <- classification
  }
  
  # Create plots using standardized format
  plots <- create_coverage_distribution_plots(
    results = coverage_results,
    classification = classification_results,
    bootstrap_summary = bootstrap_summary,
    output_dir = output_dir
  )
  
  return(list(
    results = coverage_results,
    classification = classification_results,
    bootstrap_summary = bootstrap_summary,
    plots = plots
  ))
}
#' Create coverage distribution plots for each mask
#' @param results Coverage results data frame
#' @param classification Classification results data frame
#' @param bootstrap_summary Bootstrap summary statistics dataframe
#' @param output_dir Output directory
#' @return Combined plot object
create_coverage_distribution_plots <- function(results, classification, 
                                               bootstrap_summary, output_dir) {
  plots <- list()
  
  for(mask_name in names(classification)[-1]) { # Skip Layer column
    percent_col <- paste0(mask_name, "_Percent")
    if(all(results[[percent_col]] == 0)) next
    
    # Get confidence intervals from bootstrap summary
    mask_summary <- bootstrap_summary[bootstrap_summary$mask == mask_name, ]
    ci_95 <- c(mask_summary$ci_0.95_lower, mask_summary$ci_0.95_upper)
    ci_99 <- c(mask_summary$ci_0.99_lower, mask_summary$ci_0.99_upper)
    
    # Create plot data
    plot_data <- data.frame(
      Layer = results$Layer,
      Coverage = results[[percent_col]],
      Classification = factor(classification[[mask_name]],
                              levels = c("Very high", "High", "Typical distribution", 
                                         "Low", "Very low"))
    )
    
    p <- ggplot(plot_data, 
                aes(x = reorder(factor(Layer), Coverage), 
                    y = Coverage,
                    fill = Classification)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = ci_95,
                 linetype = "dashed", color = "black") +
      geom_hline(yintercept = ci_99,
                 linetype = "dotted", color = "black") +
      scale_fill_manual(values = c("Very high" = "#DC3545",
                                   "High" = "#E57373",
                                   "Typical distribution" = "#E9ECEF",
                                   "Low" = "#90CAF9",
                                   "Very low" = "#2196F3")) +
      guides(fill = guide_legend(reverse = TRUE)) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        panel.grid.minor = element_blank(),
        plot.margin = margin(10, 10, 10, 10),
        legend.position = "right",
        legend.key.size = unit(0.8, "lines")
      ) +
      labs(title = paste("Cluster Coverage -", gsub("_", " ", mask_name)),
           subtitle = paste("Dashed: 95% CI; Dotted: 99% CI"),
           x = "Layer",
           y = "Percent Coverage",
           fill = "Distribution Type")
    
    plots[[mask_name]] <- p
  }
  
  # Combine all plots
  if(length(plots) > 0) {
    combined_plot <- wrap_plots(plots, ncol = 2) +
      plot_layout(guides = "collect") &
      theme(plot.margin = margin(20, 20, 20, 20))
    
    ggsave(
      file.path(output_dir, "visual_scene_coverage_plot.png"),
      plot = combined_plot,
      width = 15,
      height = ceiling(length(plots)/2) * 6,
      limitsize = FALSE
    )
    
    return(combined_plot)
  } else {
    return(NULL)
  }
}

create_visual_scene_grid_plots <- function(visual_scene_results, scan_results, grid, 
                                           mask_list, tif_stack, output_dir,
                                           min_coverage = 5,
                                           plot_low_significance = FALSE) {
  # Function to create mask visualization
  create_mask_plot <- function(mask, study_grid, title) {
    # Create plot data using the study grid to define valid areas
    plot_data <- expand.grid(
      row = 1:nrow(mask),
      col = 1:ncol(mask)
    ) %>%
      mutate(
        value = as.vector(mask),
        valid = !is.na(as.vector(study_grid))  # Use grid to determine valid areas
      )
    
    # Show mask only in valid areas
    plot_data$value[!plot_data$valid] <- NA
    
    # Create base plot
    p <- ggplot(plot_data, aes(x=col, y=row)) +
      # Plot valid area in light grey
      geom_tile(data = subset(plot_data, valid), 
                fill = "grey90") +
      # Plot mask area
      geom_tile(data = subset(plot_data, !is.na(value) & value == 1), 
                fill = "steelblue", alpha = 0.7) +
      coord_equal() +
      theme_void() +
      labs(title = title) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 10),
        plot.background = element_rect(fill = "white", color = NA)
      )
    
    return(p)
  }
  
  # Rest of the helper functions remain the same
  create_map_plot <- function(matrix, grid, scan_result, mask = NULL) {
    plot_data <- expand.grid(
      row = 1:nrow(matrix),
      col = 1:ncol(matrix)
    ) %>%
      mutate(
        value = as.vector(matrix)
      )
    
    plot_data$value[plot_data$value <= 0] <- NA
    
    sig_clusters <- smerc::clusters(scan_result)
    cluster_mask <- matrix(FALSE, nrow(matrix), ncol(matrix))
    boundary_points <- data.frame()
    
    if(!is.null(sig_clusters)) {
      for(i in seq_along(sig_clusters)) {
        cluster_cells <- sig_clusters[[i]]
        cells_in_cluster <- which(grid %in% cluster_cells, arr.ind = TRUE)
        cluster_mask[cells_in_cluster] <- TRUE
        
        boundary <- find_cluster_boundaries(grid, cluster_cells)
        if(any(boundary)) {
          boundary_points <- rbind(
            boundary_points,
            data.frame(
              row = which(boundary, arr.ind = TRUE)[,1],
              col = which(boundary, arr.ind = TRUE)[,2]
            )
          )
        }
      }
    }
    
    mask_overlay <- if(!is.null(mask)) {
      cluster_mask & mask == 1
    } else {
      NULL
    }
    
    p <- ggplot(plot_data, aes(x=col, y=row)) +
      geom_tile(aes(fill=value)) +
      scale_fill_gradient(low="white", high="grey20", na.value="white") +
      coord_equal() +
      theme_void() +
      theme(legend.position = "none")
    
    if(nrow(boundary_points) > 0) {
      p <- p + geom_point(data=boundary_points, color="green3", size=0.5)
    }
    
    if(!is.null(mask_overlay)) {
      mask_points <- which(mask_overlay, arr.ind = TRUE)
      if(nrow(mask_points) > 0) {
        p <- p + geom_tile(data=data.frame(
          row=mask_points[,1],
          col=mask_points[,2]
        ), fill="#FF000040")
      }
    }
    
    return(p)
  }
  
  for(mask_name in names(mask_list)) {
    # Create data frame with required columns
    if(!(paste0(mask_name, "_Percent") %in% names(visual_scene_results$results))) {
      message(sprintf("Skipping mask %s: no coverage data found", mask_name))
      next
    }
    
    if(!(mask_name %in% names(visual_scene_results$classification))) {
      message(sprintf("Skipping mask %s: no classification data found", mask_name))
      next
    }
    
    percent_col <- paste0(mask_name, "_Percent")
    if(!(percent_col %in% names(visual_scene_results$results))) {
      message(sprintf("Skipping mask %s: no percent coverage data found", mask_name))
      next
    }
    
    mask_data <- tryCatch({
      data.frame(
        Layer = visual_scene_results$results$Layer,
        Coverage = visual_scene_results$results[[percent_col]],
        Classification = visual_scene_results$classification[[mask_name]]
      )
    }, error = function(e) {
      message(sprintf("Error creating data frame for mask %s: %s", mask_name, e$message))
      return(NULL)
    })
    
    if(is.null(mask_data)) next
    
    # Process high significance layers
    high_sig_layers <- mask_data$Layer[mask_data$Classification %in% c("Very high", "High") & 
                                         mask_data$Coverage >= min_coverage]
    
    if(length(high_sig_layers) > 0) {
      layer_order <- high_sig_layers[order(mask_data$Coverage[match(high_sig_layers, mask_data$Layer)], 
                                           decreasing = TRUE)]
      
      # Create mask visualization with the full study grid
      mask_plot <- create_mask_plot(
        mask_list[[mask_name]], 
        grid,
        paste("Mask:", gsub("_", " ", mask_name))
      )
      
      sig_plots <- lapply(layer_order, function(layer) {
        p <- create_map_plot(matrix = tif_stack[[layer]],
                             grid = grid,
                             scan_result = scan_results[[layer]],
                             mask = mask_list[[mask_name]])
        
        p <- p + 
          labs(title = sprintf("C%d", layer)) +
          theme(plot.title = element_text(hjust = 0.5, size = 10),
                plot.subtitle = element_text(hjust = 0.5, size = 8)) +
          labs(subtitle = sprintf("%.1f%%", mask_data$Coverage[mask_data$Layer == layer]))
        
        return(p)
      })
      
      n_plots <- length(sig_plots)
      n_cols <- ceiling(sqrt(n_plots))
      n_rows <- ceiling(n_plots/n_cols)
      
      combined_plot <- wrap_plots(
        c(list(mask_plot), sig_plots),
        ncol = n_cols + 1,
        widths = c(1, rep(1, n_cols)),
        heights = rep(1, n_rows)
      )
      
      ggsave(file.path(output_dir, sprintf("%s_high_significance.png", mask_name)),
             plot = combined_plot,
             width = min(3 * (n_cols + 1), 15),
             height = min(3 * n_rows, 15),
             dpi = 300)
    }
    
    # Process low significance layers only if requested
    if(plot_low_significance) {
      low_sig_layers <- mask_data$Layer[mask_data$Classification %in% c("Very low", "Low") &
                                          mask_data$Coverage >= min_coverage]
      
      if(length(low_sig_layers) > 0) {
        # ... [Rest of the low significance plotting code remains the same]
        layer_order <- low_sig_layers[order(mask_data$Coverage[match(low_sig_layers, mask_data$Layer)], 
                                            decreasing = TRUE)]
        
        mask_plot <- create_mask_plot(
          mask_list[[mask_name]], 
          grid,
          paste("Mask:", gsub("_", " ", mask_name))
        )
        
        sig_plots <- lapply(layer_order, function(layer) {
          p <- create_map_plot(matrix = tif_stack[[layer]],
                               grid = grid,
                               scan_result = scan_results[[layer]],
                               mask = mask_list[[mask_name]])
          
          p <- p + 
            labs(title = sprintf("C%d", layer)) +
            theme(plot.title = element_text(hjust = 0.5, size = 10),
                  plot.subtitle = element_text(hjust = 0.5, size = 8)) +
            labs(subtitle = sprintf("%.1f%%", mask_data$Coverage[mask_data$Layer == layer]))
          
          return(p)
        })
        
        n_plots <- length(sig_plots)
        n_cols <- ceiling(sqrt(n_plots))
        n_rows <- ceiling(n_plots/n_cols)
        
        combined_plot <- wrap_plots(
          c(list(mask_plot), sig_plots),
          ncol = n_cols + 1,
          widths = c(1, rep(1, n_cols)),
          heights = rep(1, n_rows)
        )
        
        ggsave(file.path(output_dir, sprintf("%s_low_significance.png", mask_name)),
               plot = combined_plot,
               width = min(3 * (n_cols + 1), 15),
               height = min(3 * n_rows, 15),
               dpi = 300)
      }
    }
  }
}
#=============================================================================
# Main Analysis Flow
#=============================================================================

# Setup directories
main_dir <- file.path(ROOT_DIR, OUTPUT_DIR)
dir.create(main_dir, showWarnings = FALSE)

dirs <- list(
  main = main_dir,
  morans = file.path(main_dir, MORANS_DIR),
  scan = file.path(main_dir, SCAN_DIR),
  visual_scene = file.path(main_dir, VISUAL_SCENE_DIR)
)

# Create subdirectories
lapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE)

# Read TIF stack
tif_stack <- readTIFF(file.path(ROOT_DIR, TIF_PATH), all=TRUE)
tif_stack <- lapply(tif_stack, FUN = function(mat) {
  mat[is.nan(mat)] <- -1
  return(mat)
})

tif_stack <- lapply(tif_stack, FUN = function(mat) {
  mat[mat == -1] <- 0
  return(mat)
})

# Create study region and grid
study_region <- create_study_region(tif_stack[[1]])
grid <- create_grid(study_region, GRID_SIZE)

# Run Moran's I analysis
message("Running Moran's I analysis...")
morans_results <- run_morans_analysis(tif_stack, grid, dirs$morans)

# Run scan statistics analysis
message("Running scan statistics analysis...")
scan_results <- run_scan_analysis(tif_stack, grid, dirs$scan)

# Create and analyze visual scene masks
message("Running visual scene analysis...")
masks <- create_visual_scene_masks(grid)
saveRDS(masks, file.path(dirs$visual_scene, "visual_scene_masks.rds"))
# Analyze visual scene coverage with new bootstrap approach
visual_scene_results <- analyze_visual_scene_coverage(
  scan_results$results, 
  grid, 
  masks, 
  n_bootstrap = N_BOOTSTRAP,
  bootstrap_method = BOOTSTRAP_METHOD,
  output_dir = dirs$visual_scene
)

# Save visual scene results
write.csv(
  visual_scene_results$results,
  file = file.path(dirs$visual_scene, "visual_scene_coverage.csv"),
  row.names = FALSE
)

write.csv(
  visual_scene_results$classification,
  file = file.path(dirs$visual_scene, "visual_scene_classification.csv"),
  row.names = FALSE
)


create_visual_scene_grid_plots(
  visual_scene_results = visual_scene_results,
  scan_results = scan_results$results,  
  grid = grid,
  mask_list = masks,
  tif_stack = tif_stack,
  output_dir = dirs$visual_scene,
  min_coverage = MIN_COVERAGE
)
