library(dplyr)
library(spatstat)
library(sf)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(patchwork)

# Copy over required functions from previous script
filter_points_in_window <- function(points_data, window_obj) {
  temp_ppp <- ppp(
    x = points_data$x,
    y = points_data$y,
    window = window_obj,
    check = FALSE
  )
  inside <- inside.owin(temp_ppp$x, temp_ppp$y, window_obj)
  points_data[inside, ]
}

load_lasso_regions <- function(filepath) {
  cat(sprintf("\nLoading lasso regions from %s\n", basename(filepath)))
  coords_df <- read.csv(filepath, skip = 2)
  coords_df <- coords_df[grep("^HQ_Lasso_\\d+$", coords_df$Selection), ]
  
  lasso_list <- split(coords_df, coords_df$Selection) %>%
    lapply(function(region) {
      if (!identical(region[1, c("X","Y")], region[nrow(region), c("X","Y")])) {
        region <- rbind(region, region[1,])
      }
      coords_matrix <- as.matrix(region[, c("X","Y")])
      st_polygon(list(coords_matrix))
    })
  
  return(lasso_list)
}

sf_to_owin <- function(sf_poly) {
  coords <- st_coordinates(sf_poly)
  owin(poly = list(x = coords[, "X"], y = coords[, "Y"]))
}


root <- '/media/sam/Data2/baysor_rbpms_consolidated'
out_root <- '/home/sam/FinalRGC_xenium'
lasso_root <- '/media/sam/New Volume/Xenium_Data/HQ_NearestNeighbor_Zones'

# Create output directory
exclusion_dir <- file.path(out_root, "exclusion_zones")
dir.create(exclusion_dir, recursive = TRUE, showWarnings = FALSE)

# Output file path
output_file <- file.path(exclusion_dir, "cross_subtype_distances.csv")

# Load RGC data
cat("Loading RGC data...\n")
rgc_df <- read_csv(file.path(root, 'OldModel/merged_rgc_prediction_expmat.csv')) %>%
  tibble::rowid_to_column("index") %>%
  mutate(x = round(x,1),
         y = round(y,1))

alignment_df <- data.frame(X_dv=double(), Y_dv=double(), 
                             X_circ=double(), Y_circ=double(), 
                             index=double(), 
                             x=double(), y=double())

j = 1
for (i in seq(1,35,7)){

  temp_df <-  read_csv(file.path(root, 'aligned_XY_coordinates.csv'), 
                          col_names = F, col_select = c(i:(i+6))) %>%
    drop_na()
  
  colnames(temp_df) <- c("X_dv", "Y_dv", "X_circ", "Y_circ", "index", "x", "y") 
  
  temp_df <- temp_df %>%
    mutate(retina = j,
           x = round(x,1),
           y = round(y,1))
  alignment_df <- rbind(alignment_df, temp_df)
  j = j + 1
}

alignment_df <- alignment_df %>%
  mutate(index = index+1,
         )

(summary <- alignment_df %>%
  group_by(retina) %>%
  summarize(n()))

rgc_df_final <- left_join(rgc_df, alignment_df, by = c("index", "retina", "x", 'y') ) %>%
  drop_na()

# Get all lasso files
lasso_files <- list.files(lasso_root,
                          pattern = "\\d+_HQ_Lasso_coordinates\\.csv$",
                          full.names = TRUE)

#######################################################################

library(dplyr)
library(sf)
library(spatstat)

# Get all lasso files
lasso_files <- list.files(lasso_root, 
                          pattern = "\\d+_HQ_Lasso_coordinates\\.csv$", 
                          full.names = TRUE)

# Process each slide
for(file in lasso_files) {
  # Get slide ID
  slide_id <- as.numeric(gsub("(\\d+)_HQ_Lasso_coordinates\\.csv", "\\1", basename(file)))
  cat(sprintf("\nProcessing slide %d...\n", slide_id))
  
  # Load lasso regions
  lasso_list <- load_lasso_regions(file)
  
  # Filter RGC data for this slide
  slide_data <- rgc_df_final %>% filter(slide == slide_id)
  
  if(nrow(slide_data) > 0) {
    # Process each region
    for(i in seq_along(lasso_list)) {
      region_name <- names(lasso_list)[i]
      region_id <- paste(slide_id, region_name, sep="_")
      
      # Convert lasso region to sf object and owin
      lasso_sf <- st_sf(geometry = st_sfc(lasso_list[[i]]))
      lasso_window <- sf_to_owin(lasso_list[[i]])
      
      # Filter points in this region
      points_in_region <- filter_points_in_window(slide_data, lasso_window)
      
      if(nrow(points_in_region) > 0) {
        # Mark these points as belonging to this region
        rgc_df_final$study_region[rgc_df_final$index %in% points_in_region$index] <- region_id
      }
    }
  }
}

# Save the updated dataframe
write_csv(rgc_df_final, file.path(exclusion_dir, "rgc_df_with_regions_new.csv"))
###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
# Function to create plot for a single region
create_region_plot <- function(data, region_name, lasso_list, show_legend = FALSE) {
  # Get region data
  region_data <- data %>%
    filter(study_region == region_name)
  
  # Get slide ID and lasso name
  slide_id <- as.numeric(strsplit(region_name, "_")[[1]][1])
  lasso_name <- paste(strsplit(region_name, "_")[[1]][-1], collapse="_")
  
  # Get lasso geometry
  lasso_geom <- lasso_list[[which(names(lasso_list) == lasso_name)]]
  lasso_sf <- st_sf(geometry = st_sfc(lasso_geom))
  
  # Get lasso boundaries for cropping
  lasso_coords <- st_coordinates(lasso_sf)
  x_min <- min(lasso_coords[,"X"])
  x_max <- max(lasso_coords[,"X"])
  y_min <- min(lasso_coords[,"Y"])
  y_max <- max(lasso_coords[,"Y"])
  
  # Add padding (100 units) to boundaries
  x_padding <- 100
  y_padding <- 100
  x_min <- x_min - x_padding
  x_max <- x_max + x_padding
  y_min <- y_min - y_padding
  y_max <- y_max + y_padding
  
  # Get ALL points from the same slide within the expanded boundaries
  slide_data <- data %>%
    filter(slide == slide_id,
           x >= x_min, x <= x_max,
           y >= y_min, y <= y_max) %>%
    mutate(
      point_type = case_when(
        study_region == region_name ~ "in_region",
        !is.na(study_region) ~ "other_region",
        TRUE ~ "no_region"
      )
    )
  
  # Create plot
  p <- ggplot() +
    # Add points with appropriate sizing and alpha
    geom_point(
      data = slide_data,
      aes(x = x, y = y, color = Prediction),
      size = case_when(
        slide_data$point_type == "in_region" ~ 0.5,
        slide_data$point_type == "other_region" ~ 0.5,
        TRUE ~ 0.5
      ),
      alpha = case_when(
        slide_data$point_type == "in_region" ~ 0.1,
        slide_data$point_type == "other_region" ~ 0.3,
        TRUE ~ 0.3
      )
    ) +
    # Add lasso boundary
    geom_sf(data = lasso_sf, fill = NA, color = "black", linewidth = 0.5) +
    coord_sf(
      xlim = c(x_min, x_max),
      ylim = c(y_min, y_max)
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = if(show_legend) "right" else "none",
      plot.title = element_text(size = 10)
    ) +
    labs(
      title = region_name,
      x = "X (µm)",
      y = "Y (µm)"
    )
  
  return(p)
}
# Function to create plot for a transformed region (dv or circ)
create_transformed_region_plot <- function(data, region_name, lasso_list, transform = "dv", show_legend = FALSE) {
  # Get region data and slide ID
  slide_id <- as.numeric(strsplit(region_name, "_")[[1]][1])
  lasso_name <- paste(strsplit(region_name, "_")[[1]][-1], collapse="_")
  
  # Get lasso geometry to determine boundaries in original space
  lasso_geom <- lasso_list[[which(names(lasso_list) == lasso_name)]]
  lasso_sf <- st_sf(geometry = st_sfc(lasso_geom))
  
  # Get lasso boundaries for filtering
  lasso_coords <- st_coordinates(lasso_sf)
  x_min <- min(lasso_coords[,"X"])
  x_max <- max(lasso_coords[,"X"])
  y_min <- min(lasso_coords[,"Y"])
  y_max <- max(lasso_coords[,"Y"])
  
  # Add padding
  x_padding <- 100
  y_padding <- 100
  x_min <- x_min - x_padding
  x_max <- x_max + x_padding
  y_min <- y_min - y_padding
  y_max <- y_max + y_padding
  
  # Filter points based on original xy coordinates
  slide_data <- data %>%
    filter(slide == slide_id,
           x >= x_min, x <= x_max,
           y >= y_min, y <= y_max)
  
  # Determine which coordinates to use based on transform
  if(transform == "dv") {
    slide_data <- slide_data %>%
      mutate(
        plot_x = X_dv,
        plot_y = Y_dv
      )
  } else if(transform == "circ") {
    slide_data <- slide_data %>%
      mutate(
        plot_x = X_circ,
        plot_y = Y_circ
      )
  }
  
  # Add point type and color information
  slide_data <- slide_data %>%
    mutate(
      point_type = case_when(
        study_region == region_name ~ "in_region",
        !is.na(study_region) ~ "other_region",
        TRUE ~ "no_region"
      ),
      point_color = case_when(
        study_region == region_name ~ as.character(Prediction),
        TRUE ~ "grey80"
      )
    )
  
  # Create plot
  p <- ggplot() +
    geom_point(
      data = slide_data,
      aes(x = plot_x, y = plot_y, color = point_color),
      size = case_when(
        slide_data$point_type == "in_region" ~ 0.5,
        slide_data$point_type == "other_region" ~ 0.5,
        TRUE ~ 0.5
      ),
      alpha = case_when(
        slide_data$point_type == "in_region" ~ 0.8,
        slide_data$point_type == "other_region" ~ 0.3,
        TRUE ~ 0.3
      )
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = if(show_legend) "right" else "none",
      plot.title = element_text(size = 10)
    ) +
    labs(
      title = region_name,
      x = if(transform == "dv") "X_dv" else "X_circ",
      y = if(transform == "dv") "Y_dv" else "Y_circ",
      color = "Cell Type"
    )
  
  return(p)
}

# # Updated plot_rgc_regions function
# plot_rgc_regions <- function(data, lasso_files, transform = "none") {
#   # Validate transform parameter
#   if(!transform %in% c("none", "dv", "circ")) {
#     stop("transform must be one of: 'none', 'dv', 'circ'")
#   }
#   
#   # Process each study region
#   study_regions <- unique(data$study_region[!is.na(data$study_region)])
#   
#   # Create plot list
#   plot_list <- list()
#   for(i in seq_along(study_regions)) {
#     region <- study_regions[i]
#     slide_id <- as.numeric(strsplit(region, "_")[[1]][1])
#     lasso_file <- lasso_files[grep(paste0(slide_id, "_HQ_Lasso"), lasso_files)]
#     lasso_list <- load_lasso_regions(lasso_file)
#     
#     # Choose plotting function based on transform
#     if(transform == "none") {
#       plot_list[[region]] <- create_region_plot(data, region, lasso_list, 
#                                                 show_legend = (i == length(study_regions)))
#     } else {
#       plot_list[[region]] <- create_transformed_region_plot(data, region, lasso_list, 
#                                                             transform = transform,
#                                                             show_legend = (i == length(study_regions)))
#     }
#   }
#   
#   # Set title based on transform
#   title <- switch(transform,
#                   "none" = "Local Study Regions on Xenium Slide",
#                   "dv" = "Local Study Regions Along Dorsal/Ventral Axis",
#                   "circ" = "Local Study Regions in Circular Coordinates")
#   
#   # Create combined plot
#   combined_plot <- wrap_plots(plot_list, ncol = 3) +
#     plot_layout(guides = "collect") +
#     plot_annotation(
#       title = title,
#       theme = theme(
#         plot.title = element_text(size = 16, hjust = 0.5)
#       )
#     )
#   
#   return(combined_plot)
# }
# 
# # Original xy coordinates
# (original_plot <- plot_rgc_regions(rgc_df_final, lasso_files, transform = "none"))
# 
# # Dorsal-ventral coordinates
# (dv_plot <- plot_rgc_regions(rgc_df_final, lasso_files, transform = "dv"))
# 
# # Circular coordinates
# (circ_plot <- plot_rgc_regions(rgc_df_final, lasso_files, transform = "circ"))
# 

#######################################################################################
################################################################################################################
##################################################################################################################
# Improved function to compute rotation angle between two sets of points
compute_rotation_angle <- function(x1, y1, x2, y2) {
  # Center the points
  x1_c <- x1 - mean(x1)
  y1_c <- y1 - mean(y1)
  x2_c <- x2 - mean(x2)
  y2_c <- y2 - mean(y2)
  
  # Compute principal direction using PCA for both coordinate systems
  pca1 <- prcomp(cbind(x1_c, y1_c))
  pca2 <- prcomp(cbind(x2_c, y2_c))
  
  # Get the angles of the first principal components
  angle1 <- atan2(pca1$rotation[2,1], pca1$rotation[1,1])
  angle2 <- atan2(pca2$rotation[2,1], pca2$rotation[1,1])
  
  # Compute rotation needed to align the principal directions
  angle_diff <- (angle2 - angle1)
  
  # Normalize to [-pi, pi]
  angle_diff <- (angle_diff + pi) %% (2 * pi) - pi
  
  return(angle_diff * (180/pi))
}

# Function to find modal rotation for a group of points


# Improved rotation calculation using principal components
find_modal_rotation <- function(points_df) {
  # Center both coordinate sets
  xy_centered <- scale(points_df[, c("x", "y")], scale = FALSE)
  dv_centered <- scale(points_df[, c("X_dv", "Y_dv")], scale = FALSE)
  
  # Compute principal components for both coordinate systems
  xy_pca <- prcomp(xy_centered, scale = FALSE)
  dv_pca <- prcomp(dv_centered, scale = FALSE)
  
  # Get the first principal component direction for both
  xy_direction <- xy_pca$rotation[, 1]
  dv_direction <- dv_pca$rotation[, 1]
  
  # Compute angle between principal directions
  angle <- atan2(
    xy_direction[1] * dv_direction[2] - xy_direction[2] * dv_direction[1],
    xy_direction[1] * dv_direction[1] + xy_direction[2] * dv_direction[2]
  )
  
  # Convert to degrees
  angle_deg <- angle * (180/pi)
  
  # Check if we need to flip (if the determinant is negative)
  det_sign <- sign(det(solve(xy_pca$rotation) %*% dv_pca$rotation))
  if(det_sign < 0) {
    angle_deg <- angle_deg + 180
  }
  
  return(angle_deg)
}

# Modified apply_rotation function to handle reflection if needed
apply_rotation <- function(x, y, angle_deg) {
  angle_rad <- angle_deg * (pi/180)
  
  # Create rotation matrix
  R <- matrix(c(cos(angle_rad), -sin(angle_rad),
                sin(angle_rad), cos(angle_rad)), 2, 2)
  
  # Apply transformation
  coords <- matrix(c(x, y), ncol = 2)
  rotated <- coords %*% R
  
  return(list(x_rot = rotated[,1], y_rot = rotated[,2]))
}

# Function to compute and add rotated coordinates
add_rotated_coordinates <- function(data) {
  # Compute rotations for each retina+slide+study_region group
  rotation_map <- data %>%
    filter(!is.na(study_region)) %>%
    group_by(retina, slide, study_region) %>%
    group_modify(~{
      modal_rot <- find_modal_rotation(.x)
      tibble(rotation = modal_rot)
    }) %>%
    ungroup()
  
  # Add rotated coordinates to the data
  data_with_rot <- data %>%
    left_join(rotation_map, by = c("retina", "slide", "study_region")) %>%
    rowwise() %>%
    mutate(
      rot_coords = list(apply_rotation(x, y, rotation)),
      X_rot = rot_coords$x_rot,
      Y_rot = rot_coords$y_rot
    ) %>%
    select(-rot_coords)
  
  return(list(
    data = data_with_rot,
    rotation_map = rotation_map
  ))
}


# Modified create_transformed_region_plot to handle rotated coordinates
create_transformed_region_plot <- function(data, region_name, lasso_list, 
                                           transform = "dv", rotation_map = NULL,
                                           show_legend = FALSE) {
  # Get region data and slide ID
  slide_id <- as.numeric(strsplit(region_name, "_")[[1]][1])
  lasso_name <- paste(strsplit(region_name, "_")[[1]][-1], collapse="_")
  
  # Get lasso geometry and boundaries
  lasso_geom <- lasso_list[[which(names(lasso_list) == lasso_name)]]
  lasso_sf <- st_sf(geometry = st_sfc(lasso_geom))
  
  # Get original boundaries for filtering
  lasso_coords <- st_coordinates(lasso_sf)
  x_min <- min(lasso_coords[,"X"]) - 100
  x_max <- max(lasso_coords[,"X"]) + 100
  y_min <- min(lasso_coords[,"Y"]) - 100
  y_max <- max(lasso_coords[,"Y"]) + 100
  
  # Filter points
  slide_data <- data %>%
    filter(slide == slide_id,
           x >= x_min, x <= x_max,
           y >= y_min, y <= y_max)
  
  # Set coordinates based on transform type
  slide_data <- slide_data %>%
    mutate(
      plot_x = case_when(
        transform == "dv" ~ X_dv,
        transform == "circ" ~ X_circ,
        transform == "rot" ~ X_rot,
        TRUE ~ x
      ),
      plot_y = case_when(
        transform == "dv" ~ Y_dv,
        transform == "circ" ~ Y_circ,
        transform == "rot" ~ Y_rot,
        TRUE ~ y
      ),
      point_type = case_when(
        study_region == region_name ~ "in_region",
        !is.na(study_region) ~ "other_region",
        TRUE ~ "no_region"
      ),
      point_color = case_when(
        study_region == region_name ~ as.character(Prediction),
        TRUE ~ "grey80"
      )
    )
  
  # Create base plot
  p <- ggplot() +
    geom_point(
      data = slide_data,
      aes(x = plot_x, y = plot_y, color = point_color),
      size = case_when(
        slide_data$point_type == "in_region" ~ 0.5,
        slide_data$point_type == "other_region" ~ 0.5,
        TRUE ~ 0.5
      ),
      alpha = case_when(
        slide_data$point_type == "in_region" ~ 0.8,
        slide_data$point_type == "other_region" ~ 0.3,
        TRUE ~ 0.3
      )
    )
  
  # Add rotated lasso if using rot transform
  if(transform == "rot" && !is.null(rotation_map)) {
    # Get rotation angle for this region
    region_data <- slide_data %>% 
      filter(study_region == region_name) %>% 
      select(retina, slide) %>% 
      distinct()
    
    if(nrow(region_data) > 0) {
      rotation_angle <- rotation_map %>%
        filter(retina == region_data$retina[1],
               slide == region_data$slide[1]) %>%
        pull(rotation)
      
      if(length(rotation_angle) > 0) {
        rotated_lasso <- rotate_lasso(lasso_sf, rotation_angle)
        p <- p + geom_sf(data = st_sf(geometry = st_sfc(rotated_lasso)),
                         fill = NA, color = "black", linewidth = 0.5)
      }
    }
  }
  
  # Add theme and labels
  p <- p + theme_minimal() +
    theme(
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = if(show_legend) "right" else "none",
      plot.title = element_text(size = 10)
    ) +
    labs(
      title = region_name,
      x = switch(transform,
                 "dv" = "X_dv",
                 "circ" = "X_circ",
                 "rot" = "X_rot",
                 "X"),
      y = switch(transform,
                 "dv" = "Y_dv",
                 "circ" = "Y_circ",
                 "rot" = "Y_rot",
                 "Y"),
      color = "Cell Type"
    )
  
  return(p)
}

# Updated plot_rgc_regions function
plot_rgc_regions <- function(data, lasso_files, transform = "none", rotation_map = NULL) {
  # Validate transform parameter
  if(!transform %in% c("none", "dv", "circ", "rot")) {
    stop("transform must be one of: 'none', 'dv', 'circ', 'rot'")
  }
  
  # Process each study region
  study_regions <- unique(data$study_region[!is.na(data$study_region)])
  
  # Create plot list
  plot_list <- list()
  for(i in seq_along(study_regions)) {
    region <- study_regions[i]
    slide_id <- as.numeric(strsplit(region, "_")[[1]][1])
    lasso_file <- lasso_files[grep(paste0(slide_id, "_HQ_Lasso"), lasso_files)]
    lasso_list <- load_lasso_regions(lasso_file)
    
    # Choose plotting function based on transform
    if(transform == "none") {
      plot_list[[region]] <- create_region_plot(data, region, lasso_list, 
                                                show_legend = (i == length(study_regions)))
    } else {
      plot_list[[region]] <- create_transformed_region_plot(
        data, region, lasso_list, transform = transform,
        rotation_map = rotation_map,
        show_legend = (i == length(study_regions))
      )
    }
  }
  
  # Set title based on transform
  title <- switch(transform,
                  "none" = "Local Study Regions on Xenium Slide",
                  "dv" = "Local Study Regions Along Dorsal/Ventral Axis",
                  "circ" = "Local Study Regions in Circular Coordinates",
                  "rot" = "Local Study Regions with Modal Rotation")
  
  # Create combined plot
  combined_plot <- wrap_plots(plot_list, ncol = 3) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = title,
      theme = theme(
        plot.title = element_text(size = 16, hjust = 0.5)
      )
    )
  
  return(combined_plot)
}

# Example usage:
# First compute rotated coordinates and get rotation map
result <- add_rotated_coordinates(rgc_df_final)
rgc_df_final <- result$data
rotation_map <- result$rotation_map

# Then create plots
# circ_plot <- plot_rgc_regions(rgc_df_final, lasso_files, transform = "circ")
(original_plot <- plot_rgc_regions(rgc_df_final, lasso_files, transform = "none"))
(rot_plot <- plot_rgc_regions(rgc_df_final, lasso_files, transform = "rot",
                            rotation_map = rotation_map))
(dv_plot <- plot_rgc_regions(rgc_df_final, lasso_files, transform = "dv"))