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
TIF_PATH <- "RGC_stack_genes_norm_smooth.tif"        # Input TIF file containing cell type maps
OUTPUT_DIR <- "GlobalStatistics_Genes"                 # Main output directory name
MORANS_DIR <- "Morans"                                 # Directory for Moran's I results
SCAN_DIR <- "Scan"                                     # Directory for scan statistics results
VISUAL_SCENE_DIR <- "VisualScene"                      # Directory for visual scene analysis

#=============================================================================
# Visual Scene Analysis
#=============================================================================


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
        dy <- center_row - i
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
  ipsi_mask <- create_retinal_mask(grid, n_annuli = 6, n_sectors = 24, 
                                   temporal_left = T, slope = 15,
                                   sector_fill_0 = 140, sector_fill_1 = 315,
                                   outer_fill = 2, inner_fill = 0)
  
  #44% holmgren tilt
  binocular_mask <- create_retinal_mask(grid, n_annuli = 8, n_sectors = 24,
                                        temporal_left = T, slope = 4,
                                        sector_fill_0 = 170, sector_fill_1 = 345,
                                        outer_fill = 7, inner_fill = 0)
  
  peripheral_mask <- create_retinal_mask(grid, n_annuli = 1, n_sectors = 1, 
                                         temporal_left = T, slope = 0,
                                         sector_fill_0 = 0, sector_fill_1 = 365,
                                         outer_fill = 1, inner_fill = 0) - binocular_mask
  
  visual_sky_mask <- create_retinal_mask(grid, n_annuli = 3, n_sectors = 24,
                                         temporal_left = T, slope = 0,
                                         sector_fill_0 = 90, sector_fill_1 = 250,
                                         outer_fill = 2, inner_fill = 1)
  
  visual_ground_mask <-  create_retinal_mask(grid, n_annuli = 1, n_sectors = 1, 
                                             temporal_left = T, slope = 0,
                                             sector_fill_0 = 0, sector_fill_1 = 365,
                                             outer_fill = 1, inner_fill = 0) - visual_sky_mask
  
  visual_floor_mask <- create_retinal_mask(grid, n_annuli = 30, n_sectors = 24, 
                                           temporal_left = T, slope = 3,
                                           sector_fill_0 = 270, sector_fill_1 = 90,
                                           outer_fill = 20, inner_fill = 0)
  
  binocular_sky_mask <- visual_sky_mask * binocular_mask
  binocular_ground_mask <- visual_ground_mask * binocular_mask
  binocular_floor_mask <- visual_floor_mask * binocular_mask
  peripheral_sky_mask <- visual_sky_mask * peripheral_mask
  peripheral_ground_mask <- visual_ground_mask * peripheral_mask
  peripheral_floor_mask <- visual_floor_mask * peripheral_mask
  
  return(list(
    ipsi = ipsi_mask,
    binocular = binocular_mask,
    peripheral = peripheral_mask,
    visual_sky = visual_sky_mask,
    visual_ground = visual_ground_mask,
    visual_floor = visual_floor_mask,
    binocular_sky = binocular_sky_mask,
    binocular_ground = binocular_ground_mask,
    peripheral_sky = peripheral_sky_mask,
    peripheral_ground = peripheral_ground_mask,
    peripheral_floor = peripheral_floor_mask,
    binocular_floor = binocular_floor_mask
  ))
}



#' Create a grid matrix from x,y coordinates in rgc df
#' @param x x values observed vector
#' @param y y values observed vector
#' @param grid_size default 100 to match reference grid
#' @return Grid matrix
create_coordinate_grid <- function(x, y, grid_size = 100) {
  # Find the boundaries of the circular region
  max_radius <- max(sqrt(x^2 + y^2))
  
  # Create a square grid that encompasses the circular region
  # Add a small buffer to ensure we capture the full circle
  buffer <- max_radius * 0.1
  grid_extent <- max_radius + buffer
  
  # Create sequence of points for the grid
  seq_len <- seq(-grid_extent, grid_extent, length.out = grid_size)
  
  # Create matrices for x and y coordinates
  x_mat <- matrix(rep(seq_len, grid_size), grid_size, grid_size)
  y_mat <- matrix(rep(seq_len, each = grid_size), grid_size, grid_size)
  
  # Create the grid matrix
  # Points inside the circle get sequential numbers, outside get NA
  grid <- matrix(NA, grid_size, grid_size)
  cell_num <- 1
  
  for(i in 1:grid_size) {
    for(j in 1:grid_size) {
      # Check if point is inside circle
      if(sqrt(x_mat[i,j]^2 + y_mat[i,j]^2) <= max_radius) {
        grid[i,j] <- cell_num
        cell_num <- cell_num + 1
      }
    }
  }
  
  return(grid)
}

#' Function to check if a point is inside a mask region
#' @param x x values observed vector
#' @param y y values observed vector
#' @param mask boolean mask to use for comparison
#' @param grid_size default 100 to match reference grid
check_point_in_mask <- function(x, y, mask, grid_size = 100) {
  # Convert x,y to grid coordinates
  max_radius <- max(abs(c(x, y)))
  buffer <- max_radius * 0.1
  grid_extent <- max_radius + buffer
  
  # Convert x,y to matrix indices
  i <- round((y + grid_extent) * (grid_size - 1)/(2 * grid_extent)) + 1
  j <- round((x + grid_extent) * (grid_size - 1)/(2 * grid_extent)) + 1
  
  # Ensure indices are within bounds
  i <- pmax(1, pmin(grid_size, i))
  j <- pmax(1, pmin(grid_size, j))
  
  # Return mask value at those coordinates
  return(mask[cbind(i, j)] == 1)
}

# Function to convert x,y coordinates to grid indices
coords_to_indices <- function(x, y, grid_size = 100) {
  # Scale coordinates to grid indices
  # Assuming x,y are in [-1.2, 1.2] range based on your data
  scale_factor <- (grid_size - 1) / 2.4  # 2.4 is the range (-1.2 to 1.2)
  
  i <- round((y + 1.2) * scale_factor) + 1
  j <- round((x + 1.2) * scale_factor) + 1
  
  # Ensure indices are within bounds
  i <- pmax(1, pmin(grid_size, i))
  j <- pmax(1, pmin(grid_size, j))
  
  return(list(i = i, j = j))
}

# Function to apply masks to dataframe
apply_masks <- function(df, masks) {
  # Get indices for all points
  indices <- coords_to_indices(df$x, df$y)
  
  # Add mask columns using case_when
  for(mask_name in names(masks)) {
    df[[mask_name]] <- as.numeric(masks[[mask_name]][cbind(indices$i, indices$j)] == 1)
  }
  
  return(df)
}

# Visualization function
plot_masks <- function(df, mask_cols) {
  # Create long format data
  df_long <- df %>%
    pivot_longer(cols = all_of(mask_cols),
                 names_to = "mask",
                 values_to = "value")
  
  # Create plot
  ggplot(df_long, aes(x = x, y = y, color = factor(value))) +
    geom_point(alpha = 0.1) +
    facet_wrap(~mask) +
    scale_color_manual(values = c("0" = "gray80", "1" = "blue")) +
    theme_minimal() +
    labs(color = "In Mask")
}


#=============================================================================
# Main Analysis Flow
#=============================================================================

# Identify index of valid maps
rgc_path <- file.path(ROOT_DIR, "rgc_expMat_with_Studyregions_transformedCoords.csv")

rgc_metadata <- c("index", "cell", "x", "y", "z", "volume", "x_range",        
                  "y_range", "z_range", "rect_vol", "elongation", "flatness",  
                  "slide", "slice", "dapi_max", "dapi_min", "dapi_mean",     
                  "dapi_sd", "sample", "retina", "Class", "dapi_max_norm", 
                  "dapi_min_norm", "dapi_mean_norm", "dapi_sd_norm",    
                  "volume_hull", "nn_dist", "nn_id", "compactness", "sphericity",
                  "study_region", "X_dv", "Y_dv")

rgc_df <- read_csv(rgc_path) %>%
  select(-rgc_metadata) %>%
  rename(x = X_circ,
         y = Y_circ)

# rgc_df %>%
#   summarize(max(x),
#             max(y),
#             min(x),
#             min(y))
# `max(x)` `max(y)` `min(x)` `min(y)`
# <dbl>    <dbl>    <dbl>    <dbl>
#   1     1.14     1.15    -1.17    -1.16

# Derive Grid from x, y and create scene masks
grid <- create_coordinate_grid(rgc_df$x, rgc_df$y)
masks <- create_visual_scene_masks(grid)

mask_cols <- names(masks)

# check each RGC against each mask
rgc_df <- apply_masks(rgc_df, masks)

#Visualize _masks
print(plot_masks(rgc_df, mask_cols))

#=============================================================================
# Generate Volcano Plot DFs and visualizations
#=============================================================================
# Create volcanoes directory
dir.create(file.path(ROOT_DIR, "volcanoes"), showWarnings = FALSE)

# Function to calculate differential expression statistics
calculate_de_stats <- function(data, genes) {
  # Calculate mean expression for each condition
  condition_means <- data %>%
    group_by(target) %>%
    summarize(across(all_of(genes), ~ mean(.x, na.rm = TRUE))) %>%
    pivot_longer(cols = -target, 
                 names_to = "gene", 
                 values_to = "mean_expr")
  
  # Calculate fold changes
  fc_data <- condition_means %>%
    pivot_wider(names_from = target, 
                values_from = mean_expr) %>%
    mutate(
      log2FC = log2((`1` + 1) / (`0` + 1))
    )
  
  # Calculate chi-square test for each gene
  chi_square_results <- map_dfr(genes, function(gene) {
    contingency <- table(data[[gene]] > 0, data$target)
    if(nrow(contingency) > 1 && ncol(contingency) > 1) {
      test <- chisq.test(contingency)
      data.frame(
        gene = gene,
        p_value = test$p.value
      )
    } else {
      data.frame(
        gene = gene,
        p_value = 1
      )
    }
  })
  
  # Combine results
  results <- fc_data %>%
    left_join(chi_square_results, by = "gene") %>%
    mutate(
      adj_p_value = p.adjust(p_value, method = "BH"),
      neg_log10_p = -log10(adj_p_value)
    )
  
  return(results)
}

# Function to create volcano plot
create_volcano_plot <- function(de_results, mask_name, cell_type) {
  # Identify top genes to label
  top_genes <- de_results %>%
    filter(abs(log2FC) > 0.5 & adj_p_value < 0.05) %>%
    top_n(10, wt = abs(log2FC))
  
  # Create plot
  p <- ggplot(de_results, aes(x = log2FC, y = neg_log10_p)) +
    geom_point(aes(color = adj_p_value < 0.05 & abs(log2FC) > 1),
               alpha = 0.6) +
    geom_text_repel(
      data = top_genes,
      aes(label = gene),
      max.overlaps = 15
    ) +
    scale_color_manual(values = c("grey", "red")) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    labs(
      title = paste("Volcano Plot -", mask_name, "\nCell Type:", cell_type),
      x = "log2 Fold Change",
      y = "-log10 Adjusted P-value",
      color = "Significant"
    ) +
    theme(legend.position = "bottom")
  
  return(p)
}

# Main analysis loop
volcanoes <- list()
all_results <- data.frame()  # To store all results for consolidated CSV

for (mask in mask_cols) {
  # Create directory for this mask
  mask_dir <- file.path(ROOT_DIR, "volcanoes", mask)
  dir.create(mask_dir, showWarnings = FALSE)
  
  # Initialize results storage for this mask
  volcanoes[[mask]] <- list()
  
  # Get data for this mask
  rgc_df_i <- rgc_df %>%
    mutate(target = !!sym(mask)) %>%
    select(-all_of(mask_cols), -x, -y)
  
  gene_cols <- names(rgc_df_i)[!names(rgc_df_i) %in% c("Prediction", "target")]
  
  # Get unique cell types
  cell_types <- unique(rgc_df_i$Prediction)
  
  # For each cell type
  for (cell_type in cell_types) {
    # Filter data for this cell type
    cell_type_data <- rgc_df_i %>%
      filter(Prediction == cell_type)
    
    # Skip if not enough data
    if (nrow(cell_type_data) < 10) next
    
    # Calculate statistics
    results <- calculate_de_stats(cell_type_data, gene_cols)
    results$cell_type <- cell_type
    results$mask <- mask
    
    # Add to consolidated results
    all_results <- bind_rows(all_results, results)
    
    # Create and save volcano plot
    p <- create_volcano_plot(results, mask, cell_type)
    
    # Save plot
    ggsave(
      filename = file.path(mask_dir, paste0(cell_type, "_volcano.png")),
      plot = p,
      width = 10,
      height = 8,
      dpi = 300
    )
    
    # Store results
    volcanoes[[mask]][[cell_type]] <- list(
      data = results,
      plot = p
    )
  }
}

# Save consolidated results
write_csv(all_results, file.path(ROOT_DIR, "volcanoes", "all_volcano_results.csv"))

# Print summary of results
cat("\nAnalysis Summary:\n")
for (mask in mask_cols) {
  cat("\nMask:", mask, "\n")
  mask_results <- all_results %>% filter(mask == mask)
  for (cell_type in unique(mask_results$cell_type)) {
    sig_genes <- mask_results %>%
      filter(cell_type == !!cell_type,
             adj_p_value < 0.05 & abs(log2FC) > 0.5) %>%
      nrow()
    cat("  Cell Type:", cell_type, "- Significant DE genes:", sig_genes, "\n")
  }
}

hits <- all_results %>%
  filter(adj_p_value < 0.05 & abs(log2FC) > 1 )
write_csv(hits, file.path(ROOT_DIR, "volcanoes", "signifigant_volcano_results.csv"))


# Create custom p-value border color function
get_pvalue_color <- function(p) {
  case_when(
    p >= 0.05 ~ "white",
    p >= 0.01 ~ "white",
    p >= 0.001 ~ "grey",
    TRUE ~ "black"
  )
}

# Add border color column
hits <- hits %>%
  mutate(border_color = get_pvalue_color(adj_p_value))

# Define mask subsets
subset1_masks <- c("binocular", "ipsi", "visual_ground", "visual_floor")
subset2_masks <- c("binocular_floor", "binocular_ground", "binocular_sky", 
                   "binocular", "peripheral_ground", "peripheral_sky")

# Create base plotting function
create_heatmap <- function(data, title) {
  ggplot(data, aes(x = cell_type, y = gene)) +
    geom_point(aes(fill = log2FC, color = border_color), 
               shape = 21, size = 1) +
    scale_fill_gradient2(
      low = "blue", 
      mid = "white", 
      high = "red", 
      midpoint = 0,
      name = "Log2 Fold Change"
    ) +
    scale_color_identity() +
    facet_wrap(~mask, scales = "free_y", ncol = 2) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 4)  # Smaller gene labels
    ) +
    labs(
      x = "Cell Type",
      y = "Gene",
      title = title
    )
}

# Create and save first subset plot
p1 <- create_heatmap(
  hits %>% filter(mask %in% subset1_masks),
  "Differential Expression Patterns - Basic Regions"
)

ggsave(
  filename = file.path(ROOT_DIR, "volcanoes", "differential_expression_heatmap_subset1.png"),
  plot = p1,
  width = 12,
  height = 15,
  dpi = 300
)

# Create and save second subset plot
p2 <- create_heatmap(
  hits %>% filter(mask %in% subset2_masks),
  "Differential Expression Patterns - Specialized Regions"
)

ggsave(
  filename = file.path(ROOT_DIR, "volcanoes", "differential_expression_heatmap_subset2.png"),
  plot = p2,
  width = 12,
  height = 15,
  dpi = 300
)