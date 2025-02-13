# This script allows for mask exploration so parameters can be quickly visualized

# Load required libraries
library(tiff)
library(smerc)
library(tidyverse)
library(sf)
library(patchwork)


root_dir <- "/home/sam/FinalRGC_xenium/"
tif_path <- "M_RGCtypes_norm_circle_012025.tif"
tif_stack <- readTIFF(file.path(root_dir, tif_path), all=TRUE)

#' Create circular mask and identify study region with buffer
#' @param matrix Input matrix
#' @param buffer_size Number of pixels to remove from edge (default 2)
#' @return Binary matrix marking study region (NA for values <= 0 and buffer zone)
create_study_region <- function(matrix, buffer_size = 2) {
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
#' Divide region into rectangular grid cells
#' @param matrix Input matrix with NAs outside study region
#' @param grid_size Number of cells per side (e.g., 8 creates an 8x8 grid)
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

# Create study region and grid (only needs to be done once)
study_region <- create_study_region(tif_stack[[1]], buffer_size = 3)
grid <- create_grid(study_region, grid_size=71)


#' Create anatomical masks for retinal analysis with sloping annuli
#' The default behavior produces the isi-targetting region (binocular zone) characterized by Johnson 2021
#' @param grid Reference grid matrix
#' @param n_annuli Number of concentric rings (default 6)
#' @param n_sectors Number of angular sectors (default 24)
#' @param sector_fill_0 Start angle for sector fill (default 140)
#' @param sector_fill_1 End angle for sector fill (default 315)
#' @param temporal_left Boolean indicating if temporal region is on left (default TRUE)
#' @param outer_fill Number of outer annuli to fill (default 2)
#' @param inner_fill Number of inner annuli to fill (default 0)
#' @param slope Degrees to narrow angle range per annulus (default 15)
#' @return List containing mask matrix and visualization plot
create_retinal_mask <- function(grid, n_annuli = 6, n_sectors = 24, 
                                temporal_left = F, slope = 15,
                                sector_fill_0 = 140, sector_fill_1 = 315,
                                outer_fill = 2, inner_fill = 0) {
  require(ggplot2)
  require(tidyverse)
  
  # Get dimensions and center of grid
  n_rows <- nrow(grid)
  n_cols <- ncol(grid)
  center_row <- n_rows/2
  center_col <- n_cols/2
  
  # Create empty mask
  mask <- matrix(0, n_rows, n_cols)
  
  # Calculate maximum radius (to edge of valid grid points)
  valid_positions <- which(!is.na(grid), arr.ind = TRUE)
  distances <- sqrt((valid_positions[,1] - center_row)^2 + 
                      (valid_positions[,2] - center_col)^2)
  max_radius <- max(distances)
  
  # Create polar coordinates for each point
  for(i in 1:n_rows) {
    for(j in 1:n_cols) {
      # Only process if grid point is valid
      if(!is.na(grid[i,j])) {
        # Calculate relative position from center
        dy <- i - center_row
        dx <- j - center_col
        
        # Convert to polar coordinates
        radius <- sqrt(dx^2 + dy^2)
        
        # Calculate normalized radius and annulus index (0 to n_annuli-1)
        norm_radius <- radius/max_radius
        annulus_idx <- floor(norm_radius * n_annuli)
        
        # Calculate angle (0 degrees at top, clockwise)
        angle <- (atan2(dy, dx) * 180/pi + 90) %% 360
        
        # Adjust angle based on temporal orientation
        if(temporal_left) {
          angle <- (360 - angle) %% 360
        }
        
        # For outer annuli, adjust the angle range based on depth from edge
        if(annulus_idx >= (n_annuli - outer_fill) && annulus_idx < n_annuli) {
          # Calculate how many annuli in from the edge
          depth_from_edge <- (n_annuli - 1) - annulus_idx
          
          # Adjust angle range based on depth
          adjusted_fill_0 <- sector_fill_0 + (slope * depth_from_edge)
          adjusted_fill_1 <- sector_fill_1 - (slope * depth_from_edge)
          
          # Check if angle is within adjusted range
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
        # For inner annuli, use original range
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
  
  # Calculate coverage statistics
  total_valid_points <- sum(!is.na(grid))
  masked_points <- sum(mask == 1, na.rm = TRUE)
  coverage_percent <- (masked_points / total_valid_points) * 100
  
  # Create visualization data
  plot_data <- expand.grid(
    row = 1:n_rows,
    col = 1:n_cols
  ) %>%
    mutate(
      grid_value = as.vector(grid),
      mask_value = as.vector(mask)
    )
  
  # Generate radial lines data (every 30 degrees)
  angles <- seq(0, 330, by = 30)
  radial_lines <- data.frame()
  for(angle in angles) {
    rad <- (angle - 90) * pi / 180
    x <- center_col + cos(rad) * max_radius
    y <- center_row + sin(rad) * max_radius
    radial_lines <- rbind(radial_lines,
                          data.frame(
                            x1 = center_col,
                            y1 = center_row,
                            x2 = x,
                            y2 = y,
                            angle = angle
                          ))
  }
  
  # Generate degree labels
  label_radius <- max_radius * 1.1
  degree_labels <- data.frame(
    x = center_col + cos((angles - 90) * pi / 180) * label_radius,
    y = center_row + sin((angles - 90) * pi / 180) * label_radius,
    label = as.character(angles)
  )
  
  # Create enhanced visualization
  base_plot <- function(data, fill_aes) {
    ggplot(data, aes(x = col, y = row)) +
      geom_tile(fill_aes) +
      geom_segment(data = radial_lines,
                   aes(x = x1, y = y1, xend = x2, yend = y2),
                   color = "grey50", linetype = "dashed", size = 0.25) +
      geom_text(data = degree_labels,
                aes(x = x, y = y, label = label),
                size = 3, color = "grey30") +
      coord_equal() +
      scale_y_reverse() +
      theme_minimal() +
      theme(panel.grid = element_blank())
  }
  
  p1 <- base_plot(plot_data, aes(fill = factor(grid_value))) +
    scale_fill_viridis_d(na.value = "white", guide = "none") +
    labs(title = "Original Grid")
  
  p2 <- base_plot(plot_data, aes(fill = factor(mask_value))) +
    scale_fill_manual(
      values = c("0" = "white", "1" = "red"),
      na.value = "grey80",
      guide = "none"
    ) +
    labs(title = "Generated Mask")
  
  # Create description of angle ranges
  range_desc <- sprintf("Outer edge: %d°-%d°", sector_fill_0, sector_fill_1)
  if(outer_fill > 1) {
    for(i in 1:(outer_fill-1)) {
      adj_0 <- sector_fill_0 + (slope * i)
      adj_1 <- sector_fill_1 - (slope * i)
      range_desc <- paste0(range_desc, sprintf("\nLevel %d in: %d°-%d°", i+1, adj_0, adj_1))
    }
  }
  
  # Combine plots
  combined_plot <- p1 + p2 + 
    plot_layout(ncol = 2) +
    plot_annotation(
      title = "Retinal Mask Generation Results",
      subtitle = paste(range_desc, sprintf("\nCoverage: %.1f%% of valid area", coverage_percent)),
      theme = theme(plot.title = element_text(hjust = 0.5),
                    plot.subtitle = element_text(hjust = 0.5))
    )
  
  # Print plot
  print(combined_plot)
  
  return(list(
    mask = mask,
    plot = combined_plot,
    coverage = coverage_percent,
    range_description = range_desc
  ))
}


w3_top_mask <- create_retinal_mask(grid, n_annuli = 8, n_sectors = 24,
                               temporal_left = F, slope = 6,
                               sector_fill_0 = 210, sector_fill_1 = 150,
                               outer_fill = 7, inner_fill = 1)$mask

w3_bottom_mask <- create_retinal_mask(grid, n_annuli = 8, n_sectors = 24,
                                   temporal_left = F, slope = 0,
                                   sector_fill_0 = 120, sector_fill_1 = 240,
                                   outer_fill = 2, inner_fill = 0)$mask


ac_mask <- create_retinal_mask(grid, n_annuli = 8, n_sectors = 24,
                                      temporal_left = F, slope = 15,
                                      sector_fill_0 = 300, sector_fill_1 = 130,
                                      outer_fill = 4, inner_fill = 0)$mask


# # Create biologically interesting masks
binocular_mask <- create_retinal_mask(grid, n_annuli = 8, n_sectors = 24,
                              temporal_left = T, slope = 6,
                              sector_fill_0 = 170, sector_fill_1 = 345,
                              outer_fill = 6, inner_fill = 0)$mask
# 
# peripheral_mask <- create_retinal_mask(grid, n_annuli = 1, n_sectors = 1, 
#                                       temporal_left = T, slope = 0,
#                                       sector_fill_0 = 0, sector_fill_1 = 365,
#                                     outer_fill = 1, inner_fill = 0)$mask - binocular_mask
# 
visual_ground_mask <- create_retinal_mask(grid, n_annuli = 30, n_sectors = 24,
                                          temporal_left = T, slope = 3,
                                          sector_fill_0 = 270, sector_fill_1 = 90,
                                          outer_fill = 20, inner_fill = 0)$mask


visual_ground_mask <- create_retinal_mask(grid, n_annuli = 30, n_sectors = 24,
                                          temporal_left = T, slope = 3,
                                          sector_fill_0 = 270, sector_fill_1 = 90,
                                          outer_fill = 20, inner_fill = 0)$mask
visual_ground_mask_f <- create_retinal_mask(grid, n_annuli = 30, n_sectors = 24,
                                          temporal_left = T, slope = 3,
                                          sector_fill_0 = 90, sector_fill_1 = 270,
                                          outer_fill = 20, inner_fill = 0)$mask







sopsin_mask <- create_retinal_mask(grid, n_annuli = 3, n_sectors = 24,
                                          temporal_left = T, slope = 0,
                                          sector_fill_0 = 90, sector_fill_1 = 250,
                                          outer_fill = 2, inner_fill = 1)
# 
# 
# visual_sky_mask <- create_retinal_mask(grid, n_annuli = 1, n_sectors = 1, 
#                                        temporal_left = T, slope = 0,
#                                        sector_fill_0 = 0, sector_fill_1 = 365,
#                                        outer_fill = 1, inner_fill = 0)$mask - visual_ground_mask
# 
# binocular_sky_mask <- visual_sky_mask * binocular_mask
# binocular_ground_mask <- visual_ground_mask * binocular_mask
# peripheral_sky_mask <- visual_sky_mask * peripheral_mask
# peripheral_ground_mask <- visual_ground_mask * peripheral_mask
# 
# masks <- c(binocular_mask, peripheral_mask, 
#            visual_sky_mask, visual_ground_mask,
#            binocular_sky_mask, binocular_ground_mask, 
#            peripheral_sky_mask, peripheral_ground_mask
#            )



# Create biologically interesting masks
binocular_mask <- create_retinal_mask(grid, n_annuli = 6, n_sectors = 24, 
                                      temporal_left = T, slope = 15,
                                      sector_fill_0 = 140, sector_fill_1 = 315,
                                      outer_fill = 2, inner_fill = 0)$mask

peripheral_mask <- create_retinal_mask(grid, n_annuli = 1, n_sectors = 1, 
                                       temporal_left = T, slope = 0,
                                       sector_fill_0 = 0, sector_fill_1 = 365,
                                       outer_fill = 1, inner_fill = 0)$mask - binocular_mask

visual_sky_mask <- create_retinal_mask(grid, n_annuli = 30, n_sectors = 24, 
                                          temporal_left = T, slope = 3,
                                          sector_fill_0 = 270, sector_fill_1 = 90,
                                          outer_fill = 20, inner_fill = 0)$mask

visual_ground_mask <- create_retinal_mask(grid, n_annuli = 1, n_sectors = 1, 
                                       temporal_left = T, slope = 0,
                                       sector_fill_0 = 0, sector_fill_1 = 365,
                                       outer_fill = 1, inner_fill = 0)$mask - visual_sky_mask

binocular_sky_mask <- visual_sky_mask * binocular_mask
binocular_ground_mask <- visual_ground_mask * binocular_mask
peripheral_sky_mask <- visual_sky_mask * peripheral_mask
peripheral_ground_mask <- visual_ground_mask * peripheral_mask
