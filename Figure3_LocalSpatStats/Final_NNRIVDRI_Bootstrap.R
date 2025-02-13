library(tidyverse)
library(patchwork)
library(spatstat)
library(sf)
################################################################################################
######################## FUNCTIONS FOR LOADING & JOING LASSOS ##################################
################################################################################################
load_lasso_regions <- function(lasso_root, study_region_levels) {
  # Get all lasso files
  lasso_files <- list.files(lasso_root, pattern = "\\d+_HQ_Lasso_coordinates\\.csv$", 
                            full.names = TRUE)
  
  regions_unnested <- list()
  k <- 1
  
  for(file in lasso_files) {
    # Get slide ID from filename
    slide_id <- as.numeric(gsub("(\\d+)_HQ_Lasso_coordinates\\.csv", "\\1", basename(file)))
    
    # Read coords, skipping first two rows
    coords_df <- read.csv(file, skip = 2)
    
    # Filter for HQ_Lasso regions
    coords_df <- coords_df[grep("^HQ_Lasso_\\d+$", coords_df$Selection), ]
    
    # Process each selection
    region_list <- split(coords_df, coords_df$Selection)
    
    for(region_name in names(region_list)) {
      region <- region_list[[region_name]]
      
      # Create study_region identifier matching rgc_df
      study_region <- paste0(slide_id, "_", region_name)
      
      # Only process if this study_region exists in our data
      if(study_region %in% study_region_levels) {
        # Get region_ID based on factor level - exactly matching rgc_df transformation
        region_ID <- as.numeric(factor(study_region, levels = study_region_levels))
        
        # Ensure polygon is closed
        if (!identical(region[1, c("X","Y")], region[nrow(region), c("X","Y")])) {
          region <- rbind(region, region[1,])
        }
        
        # Transform X coordinates using region_ID
        region$X <- region$X + region_ID * 10000
        
        # Store transformed coordinates
        regions_unnested[[k]] <- list(
          geometry = st_polygon(list(as.matrix(region[, c("X","Y")]))),
          study_region = study_region,
          region_ID = region_ID
        )
        k <- k + 1
      }
    }
  }
  
  return(regions_unnested)
}

transform_and_plot_regions <- function(regions_unnested) {
  all_coords <- data.frame()
  
  for(i in 1:length(regions_unnested)) {
    region <- regions_unnested[[i]]
    coords <- st_coordinates(region$geometry)
    
    all_coords <- rbind(all_coords, data.frame(
      x = coords[,"X"],
      y = coords[,"Y"],
      study_region = region$study_region
    ))
  }
  
  # Create plot matching rgc_df plot
  p <- ggplot(all_coords, aes(x = x, y = y, color = study_region, group = study_region)) +
    geom_path() +
    theme_minimal() +
    ggtitle("Splayed Out Study Regions") +
    theme(legend.position = "right")
  
  return(list(
    coordinates = all_coords,
    plot = p
  ))
}

add_corridors_to_regions <- function(regions_unnested, corridor_width = 1) {
  region_ids <- sort(unique(sapply(regions_unnested, function(x) x$region_ID)))
  connected_regions <- list()
  corridors <- list()
  
  # Function to find closest points between two polygons
  find_closest_points <- function(poly1_coords, poly2_coords) {
    min_dist <- Inf
    p1_best <- NULL
    p2_best <- NULL
    
    # Compare all points from poly1 to all points from poly2
    for(i in 1:nrow(poly1_coords)) {
      p1 <- poly1_coords[i, c("X","Y")]
      for(j in 1:nrow(poly2_coords)) {
        p2 <- poly2_coords[j, c("X","Y")]
        dist <- sqrt(sum((p1 - p2)^2))
        if(dist < min_dist) {
          min_dist <- dist
          p1_best <- p1
          p2_best <- p2
        }
      }
    }
    return(list(p1 = p1_best, p2 = p2_best))
  }
  
  # Function to create corridor between two points
  create_corridor <- function(start_point, end_point, width = corridor_width) {
    # Vector from start to end
    vec <- end_point - start_point
    # Normalized perpendicular vector
    perp <- c(-vec[2], vec[1])
    perp <- perp / sqrt(sum(perp^2)) * width/2
    
    # Create corridor as thin rectangle
    corridor_points <- matrix(c(
      start_point[1] + perp[1], start_point[2] + perp[2],
      start_point[1] - perp[1], start_point[2] - perp[2],
      end_point[1] - perp[1], end_point[2] - perp[2],
      end_point[1] + perp[1], end_point[2] + perp[2],
      start_point[1] + perp[1], start_point[2] + perp[2]
    ), ncol = 2, byrow = TRUE)
    
    st_polygon(list(corridor_points))
  }
  
  # For each region (except the last)
  for(i in 1:(length(region_ids) - 1)) {
    current_region <- regions_unnested[[i]]
    next_region <- regions_unnested[[i + 1]]
    
    # Get boundary coordinates
    current_coords <- st_coordinates(current_region$geometry)
    next_coords <- st_coordinates(next_region$geometry)
    
    # Find closest pair of points
    closest_points <- find_closest_points(current_coords, next_coords)
    
    # Create corridor
    corridor <- create_corridor(closest_points$p1, closest_points$p2)
    corridors[[i]] <- list(
      geometry = corridor,
      is_corridor = TRUE,
      start_point = closest_points$p1,
      end_point = closest_points$p2
    )
  }
  
  # Combine all geometries
  all_geometries <- c(
    lapply(regions_unnested, function(x) list(geometry = x$geometry, is_corridor = FALSE)),
    corridors
  )
  
  return(all_geometries)
}

# Modified plot function to show connection points if desired
plot_connected_regions <- function(connected_regions, show_connection_points = TRUE) {
  plot_data <- data.frame()
  connection_points <- data.frame()
  
  for(i in seq_along(connected_regions)) {
    coords <- st_coordinates(connected_regions[[i]]$geometry)
    plot_data <- rbind(plot_data, data.frame(
      x = coords[,"X"],
      y = coords[,"Y"],
      group = i,
      is_corridor = connected_regions[[i]]$is_corridor
    ))
    
    # Store connection points if available
    if(show_connection_points && !is.null(connected_regions[[i]]$start_point)) {
      connection_points <- rbind(connection_points,
                                 data.frame(
                                   x = c(connected_regions[[i]]$start_point[1], 
                                         connected_regions[[i]]$end_point[1]),
                                   y = c(connected_regions[[i]]$start_point[2], 
                                         connected_regions[[i]]$end_point[2])
                                 )
      )
    }
  }
  
  p <- ggplot(plot_data, aes(x = x, y = y, group = group)) +
    geom_path(aes(color = is_corridor, size = is_corridor)) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
    scale_size_manual(values = c("TRUE" = 0.5, "FALSE" = 1))
  
  if(show_connection_points && nrow(connection_points) > 0) {
    p <- p + geom_point(data = connection_points, aes(x = x, y = y),
                        inherit.aes = FALSE, color = "green", size = 2)
  }
  
  p <- p + theme_minimal() +
    ggtitle("Connected Study Regions") +
    theme(legend.position = "right")
  
  return(p)
}

create_connected_window <- function(connected_regions) {
  # First combine all geometries (regions and corridors) into a single polygon
  # Extract all coordinates
  all_polys <- lapply(connected_regions, function(x) st_coordinates(x$geometry))
  
  # Convert to owin format
  # We need to separate the polygons into a list of x and y coordinates
  poly_list <- lapply(all_polys, function(coords) {
    list(x = coords[,"X"], y = coords[,"Y"])
  })
  
  # Create owin object using union.owin
  window_obj <- NULL
  for(poly in poly_list) {
    poly_owin <- owin(poly = poly)
    if(is.null(window_obj)) {
      window_obj <- poly_owin
    } else {
      window_obj <- union.owin(window_obj, poly_owin)
    }
  }
  
  return(window_obj)
}
################################################################################################
######################## FUNCTIONS FOR SPATIAL STATISTICS ######################################
################################################################################################

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
  nn_mean <- mean(nn)
  nn_sd <- sd(nn)
  nnri <- nn_mean / nn_sd
  
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
  area_mean <- mean(voro_areas)
  area_sd <- sd(voro_areas)
  vdri <- area_mean / area_sd
  
  return(list(
    nnri = nnri,
    vdri = vdri,
    n_edge = sum(!inner_points),  # number of edge cells
    n_interior = sum(inner_points),  # number of interior cells
    nn = nn, # nearest neighbor distances
    voro_areas = voro_areas, # all tesselation areas
    nn_mean = nn_mean,
    nn_sd = nn_sd,
    area_mean = area_mean,
    area_sd = area_sd
  ))
}

################################################################################################
######################## LOAD DATASETS #########################################################
################################################################################################

set.seed(18)
n_bootstrap = 1000
data_root <- '/home/sam/MappingAllMurineRGCs/Data/'
out_root <- file.path(data_root, 'Figure3_Outputs')
lasso_root <- '/media/sam/New Volume/Xenium_Data/HQ_NearestNeighbor_Zones'
NN_dir <- file.path(out_root, "NNRI_VDRI_Final_Analysis")
rgc_path <- file.path(data_root, "rgc_expMat_with_Studyregions_transformedCoords.csv")

rgc_df <- read_csv(rgc_path) %>%
  drop_na() %>%
  mutate(study_region = factor(study_region),
         subtype = factor(Prediction),
         region_ID = as.numeric(study_region),
         x = x + region_ID*10000)  %>%
  select(x,y,study_region, subtype) %>%
  group_by(study_region, subtype) %>% # Filter out all under sampled subtypes to prevent outlier driven results
  mutate(n = n()) %>%
  filter(n > 2) %>%
  select(-n) %>%
  ungroup()
  


################################################################################################
######################## LOAD AND JOIN STUDY WINDOWS ###########################################
################################################################################################

# Get the study region levels from rgc_df
study_region_levels <- levels(rgc_df$study_region)
# Load the lasso regions
regions_unnested <- load_lasso_regions(lasso_root, study_region_levels)
# Place the lasso regions onto the common "slide"
transformed_data <- transform_and_plot_regions(regions_unnested)

# Connect the lasso-regions
connected_regions <- add_corridors_to_regions(regions_unnested)
# Create the window object
window_obj <- create_connected_window(connected_regions)


################################################################################################
######################## VERIFY ALL DATA LOADING IS CORRECT ####################################
################################################################################################
ggplot(rgc_df, aes(x,y, color=study_region)) +
  geom_point() +
  theme_minimal() +
  ggtitle("Splayed Out Study Regions")


# Visualize and demonstrate the bootstrapping approach
p_boot <- ggplot(rgc_df, aes(x,y, color=subtype))+
  geom_point(size = 1) +
  theme_void() +
  ggtitle("Original") +
  theme(legend.position = "none")+ 
  xlim(17850,18500)+
  ylim(18000,19000)
for (i in 1:3) {
  rgc_shuffled <- rgc_df %>%
    group_by(study_region) %>%
    mutate(subtype = sample(subtype)) %>%
    ungroup()
  
  p <- ggplot(rgc_shuffled, aes(x,y, color=subtype)) +
    geom_point(size = 1) +
    theme_void() +
    ggtitle("Bootstrap") +
    theme(legend.position = "none")+
    xlim(17850,18500)+
    ylim(18000,19000)
  
  p_boot <- p_boot + p
}
print(p_boot)
#Verify all steps worked with plotting
print(transformed_data$plot)
p <- plot_connected_regions(connected_regions, show_connection_points = TRUE)
print(p)

print(p +
        xlim(17500,38000)+
        ylim(17000,20500))

# Test creating a ppp object with some sample points
test_ppp <- ppp(
  x = rgc_df$x,
  y = rgc_df$y,
  window = window_obj,
  check = TRUE
)
plot(test_ppp, main="Points in Connected Study Region")

################################################################################################
######################## BOOTSTRAP FUNCTIONS ###################################################
################################################################################################

precompute_spatial_data <- function(rgc_df, window_obj) {
  # Create initial ppp object
  full_ppp <- ppp(
    x = rgc_df$x,
    y = rgc_df$y,
    window = window_obj,
    check = TRUE
  )
  
  cat("Computing full distance matrix...\n")
  # Precompute ALL pairwise distances
  dist_mat <- crossdist(full_ppp$x, full_ppp$y, full_ppp$x, full_ppp$y)
  diag(dist_mat) <- Inf  # Set diagonal to Inf to exclude self-distances
  
  # Store study region information for each point
  region_indices <- split(1:nrow(rgc_df), rgc_df$study_region)
  
  cat("Precomputation complete.\n")
  return(list(
    ppp = full_ppp,
    dist_mat = dist_mat,
    region_indices = region_indices,
    window = window_obj
  ))
}

calculate_stats_from_precomputed <- function(precomp_data, point_indices, minimum_distance = 10) {
  if(length(point_indices) < 3) {
    warning("Fewer than 3 points provided")
    return(NULL)
  }
  
  # Get distances for this subset of points
  subset_dist_mat <- precomp_data$dist_mat[point_indices, point_indices]
  nn_dists <- apply(subset_dist_mat, 1, min)
  
  # Filter distances
  valid_nn <- nn_dists[nn_dists >= minimum_distance]
  
  if(length(valid_nn) < 3) {
    warning("Fewer than 3 valid distances after minimum distance filtering")
    return(NULL)
  }
  
  # Calculate NNRI
  nn_mean <- mean(valid_nn)
  nn_sd <- sd(valid_nn)
  nnri <- nn_mean / nn_sd
  
  # Create temporary ppp object for Voronoi calculations
  temp_ppp <- ppp(
    x = precomp_data$ppp$x[point_indices],
    y = precomp_data$ppp$y[point_indices],
    window = precomp_data$window,
    check = TRUE
  )
  
  # Calculate Voronoi tessellation
  voro <- dirichlet(temp_ppp)
  
  # Identify edge points
  edge_dists <- bdist.points(temp_ppp)
  buffer_distance <- sqrt(area.owin(precomp_data$window)) * 0.05
  inner_points <- edge_dists > buffer_distance
  
  if(sum(inner_points) < 3) {
    return(list(
      nnri = nnri,
      vdri = NA,
      n_points = length(point_indices),
      n_valid_nn = length(valid_nn),
      nn_mean = nn_mean,
      nn_sd = nn_sd,
      nn_distances = valid_nn,
      voro_areas = NA
    ))
  }
  
  # Calculate areas only for interior cells
  voro_areas <- tile.areas(voro)[inner_points]
  area_mean <- mean(voro_areas)
  area_sd <- sd(voro_areas)
  vdri <- area_mean / area_sd
  
  # print(paste("Points before min dist:", length(point_indices)))
  # print(paste("Points after min dist:", length(valid_nn)))
  # print(paste("Interior points:", sum(inner_points)))
  
  return(list(
    nnri = nnri,
    vdri = vdri,
    n_points = length(point_indices),
    n_valid_nn = length(valid_nn),
    n_edge = sum(!inner_points),
    n_interior = sum(inner_points),
    nn_mean = nn_mean,
    nn_sd = nn_sd,
    area_mean = area_mean,
    area_sd = area_sd,
    nn_distances = valid_nn,
    voro_areas = voro_areas
  ))
}

run_bootstrap_analysis <- function(rgc_df, precomp_data, output_dir, n_bootstrap = 1000, minimum_distance = 10) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Function to safely access list elements
  safe_get <- function(x, field, default = NA) {
    if (is.null(x) || is.null(x[[field]])) {
      return(default)
    }
    return(x[[field]])
  }
  
  # Function to load existing data
  # Function to load existing data
  load_existing_data <- function(subtype) {
    stats_file <- file.path(output_dir, sprintf("temp_%s_bootstrap_stats.csv", subtype))
    nn_file <- file.path(output_dir, sprintf("temp_%s_bootstrap_nn_distances.csv", subtype))
    areas_file <- file.path(output_dir, sprintf("temp_%s_bootstrap_voro_areas.csv", subtype))
    
    existing_data <- list(
      stats = if(file.exists(stats_file)) suppressMessages(read_csv(stats_file, show_col_types = FALSE)) else data.frame(),
      nn_distances = if(file.exists(nn_file)) suppressMessages(read_csv(nn_file, show_col_types = FALSE)) else data.frame(),
      voro_areas = if(file.exists(areas_file)) suppressMessages(read_csv(areas_file, show_col_types = FALSE)) else data.frame()
    )
    
    # Get maximum bootstrap_id
    max_bootstrap <- if(nrow(existing_data$stats) > 0) max(existing_data$stats$bootstrap_id) else 0
    
    return(list(data = existing_data, max_bootstrap = max_bootstrap))
  }
  
  # Initialize overall results storage
  observed_stats <- data.frame()
  
  subtypes <- levels(rgc_df$subtype)
  
  for(st in subtypes) {
    cat(sprintf("\nAnalyzing subtype %s\n", st))
    
    # Load any existing data for this subtype
    existing <- load_existing_data(st)
    start_bootstrap <- existing$max_bootstrap + 1
    
    # Get observed indices for this subtype
    obs_indices <- which(rgc_df$subtype == st)
    
    # Calculate observed statistics if not already done
    obs_file <- file.path(output_dir, sprintf("temp_%s_observed.rds", st))
    
    if(file.exists(obs_file)) {
      cat("Loading existing observed statistics...\n")
      obs_stats <- readRDS(obs_file)
    } else {
      obs_stats <- calculate_stats_from_precomputed(precomp_data, obs_indices, minimum_distance)
      if(!is.null(obs_stats)) {
        saveRDS(obs_stats, obs_file)
      }
    }
    
    # Create a single row data frame with safe value extraction
    new_stats_row <- data.frame(
      subtype = st,
      nnri = safe_get(obs_stats, "nnri"),
      vdri = safe_get(obs_stats, "vdri"),
      n_points = safe_get(obs_stats, "n_points", length(obs_indices)),
      n_valid_nn = safe_get(obs_stats, "n_valid_nn", 0),
      n_edge = safe_get(obs_stats, "n_edge"),
      n_interior = safe_get(obs_stats, "n_interior"),
      nn_mean = safe_get(obs_stats, "nn_mean"),
      nn_sd = safe_get(obs_stats, "nn_sd"),
      area_mean = safe_get(obs_stats, "area_mean"),
      area_sd = safe_get(obs_stats, "area_sd"),
      status = if(is.null(obs_stats)) "insufficient_data" else "analyzed"
    )
    
    # Add new row to observed_stats
    observed_stats <- rbind(observed_stats, new_stats_row)
    
    # Skip bootstrapping if we don't have valid observed statistics
    if(is.null(obs_stats) || new_stats_row$status == "insufficient_data") {
      cat(sprintf("Skipping bootstrap for subtype %s: insufficient valid observations\n", st))
      next
    }
    
    # Rest of your bootstrap code remains the same...
    # Only run remaining bootstrap iterations
    if(start_bootstrap <= n_bootstrap) {
      remaining_bootstraps <- n_bootstrap - start_bootstrap + 1
      if(remaining_bootstraps <= 0) {
        cat("All requested bootstrap iterations already completed for this subtype.\n")
        next
      }
      
      cat(sprintf("Running %d additional bootstrap iterations (%d to %d)...\n", 
                  remaining_bootstraps, start_bootstrap, n_bootstrap))
      
      # Only create progress bar if more than one iteration
      use_progress_bar <- remaining_bootstraps > 1
      if(use_progress_bar) {
        pb <- txtProgressBar(min = 0, max = remaining_bootstraps, style = 3)
      }
      
      for(i in start_bootstrap:n_bootstrap) {
        # Set seed based on bootstrap iteration
        set.seed(18 + i)
        
        # Shuffle labels within each study region
        boot_indices <- unlist(lapply(precomp_data$region_indices, function(idx) {
          n_select <- sum(obs_indices %in% idx)
          if(n_select > 0) sample(idx, n_select) else integer(0)
        }))
        
        # Calculate bootstrap statistics
        boot_stats <- calculate_stats_from_precomputed(precomp_data, boot_indices, minimum_distance)
        
        if(!is.null(boot_stats)) {
          # Create new rows for this iteration using safe value extraction
          new_stats <- data.frame(
            subtype = st,
            bootstrap_id = i,
            nnri = safe_get(boot_stats, "nnri"),
            vdri = safe_get(boot_stats, "vdri"),
            n_points = safe_get(boot_stats, "n_points", length(boot_indices)),
            n_valid_nn = safe_get(boot_stats, "n_valid_nn", 0),
            n_edge = safe_get(boot_stats, "n_edge"),
            n_interior = safe_get(boot_stats, "n_interior"),
            nn_mean = safe_get(boot_stats, "nn_mean"),
            nn_sd = safe_get(boot_stats, "nn_sd"),
            area_mean = safe_get(boot_stats, "area_mean"),
            area_sd = safe_get(boot_stats, "area_sd")
          )
          
          new_nn_distances <- if(!is.null(boot_stats$nn_distances)) {
            data.frame(
              subtype = st,
              bootstrap_id = i,
              distance = boot_stats$nn_distances
            )
          } else {
            data.frame()
          }
          
          new_voro_areas <- if(!is.null(boot_stats$voro_areas) && !all(is.na(boot_stats$voro_areas))) {
            data.frame(
              subtype = st,
              bootstrap_id = i,
              area = boot_stats$voro_areas
            )
          } else {
            data.frame()
          }
          
          # Append to existing data
          existing$data$stats <- rbind(existing$data$stats, new_stats)
          if(nrow(new_nn_distances) > 0) {
            existing$data$nn_distances <- rbind(existing$data$nn_distances, new_nn_distances)
          }
          if(nrow(new_voro_areas) > 0) {
            existing$data$voro_areas <- rbind(existing$data$voro_areas, new_voro_areas)
          }
          
          # Save intermediate results every 10 iterations
          if(i %% 10 == 0 || i == n_bootstrap) {
            write_csv(existing$data$stats, 
                      file.path(output_dir, sprintf("temp_%s_bootstrap_stats.csv", st)))
            if(nrow(existing$data$nn_distances) > 0) {
              write_csv(existing$data$nn_distances, 
                        file.path(output_dir, sprintf("temp_%s_bootstrap_nn_distances.csv", st)))
            }
            if(nrow(existing$data$voro_areas) > 0) {
              write_csv(existing$data$voro_areas, 
                        file.path(output_dir, sprintf("temp_%s_bootstrap_voro_areas.csv", st)))
            }
          }
        }
        
        if(use_progress_bar) {
          setTxtProgressBar(pb, i - start_bootstrap + 1)
        }
      }
      if(use_progress_bar) {
        close(pb)
      }
    } else {
      cat("All requested bootstrap iterations already completed for this subtype.\n")
    }
  }
  
  # Combine all results at the end
  all_bootstrap_stats <- data.frame()
  all_bootstrap_nn <- data.frame()
  all_bootstrap_areas <- data.frame()
  
  for(st in subtypes) {
    existing <- load_existing_data(st)
    all_bootstrap_stats <- rbind(all_bootstrap_stats, existing$data$stats)
    all_bootstrap_nn <- rbind(all_bootstrap_nn, existing$data$nn_distances)
    all_bootstrap_areas <- rbind(all_bootstrap_areas, existing$data$voro_areas)
  }
  
  # Save final combined results
  write_csv(observed_stats, file.path(output_dir, "observed_statistics.csv"))
  write_csv(all_bootstrap_stats, file.path(output_dir, "bootstrap_statistics.csv"))
  write_csv(all_bootstrap_nn, file.path(output_dir, "bootstrap_nn_distances.csv"))
  write_csv(all_bootstrap_areas, file.path(output_dir, "bootstrap_voro_areas.csv"))
  
  # Save observed nn distances and voro areas
  observed_nn_distances <- data.frame()
  observed_voro_areas <- data.frame()
  
  for(i in seq_along(subtypes)) {
    st <- subtypes[i]
    obs_file <- file.path(output_dir, sprintf("temp_%s_observed.rds", st))
    if(file.exists(obs_file)) {
      obs_stats <- readRDS(obs_file)
      if(!is.null(obs_stats)) {
        # Add nn distances
        if(!is.null(obs_stats$nn_distances)) {
          observed_nn_distances <- rbind(
            observed_nn_distances,
            data.frame(
              subtype = st,
              distance = obs_stats$nn_distances
            )
          )
        }
        
        # Add voro areas
        if(!is.null(obs_stats$voro_areas) && !all(is.na(obs_stats$voro_areas))) {
          observed_voro_areas <- rbind(
            observed_voro_areas,
            data.frame(
              subtype = st,
              area = obs_stats$voro_areas
            )
          )
        }
      }
    }
  }
  
  # Save observed metrics
  write_csv(observed_nn_distances, file.path(output_dir, "observed_nn_distances.csv"))
  write_csv(observed_voro_areas, file.path(output_dir, "observed_voro_areas.csv"))
  
  # Create comprehensive statistical summary
  create_statistical_summary <- function(stats_df, nn_df, voro_df) {
    # Helper functions for safe statistics
    safe_stat <- function(x, fn) {
      if(all(is.na(x))) return(NA)
      fn(x, na.rm = TRUE)
    }
    
    safe_quantile <- function(x, prob) {
      if(all(is.na(x))) return(NA)
      quantile(x, prob, na.rm = TRUE)
    }
    
    # Calculate quantiles and CIs for NNRI and VDRI
    metrics_summary <- stats_df %>%
      group_by(subtype) %>%
      summarise(
        # NNRI statistics
        nnri_n = sum(!is.na(nnri)),
        nnri_mean = safe_stat(nnri, mean),
        nnri_median = safe_stat(nnri, median),
        nnri_min = safe_stat(nnri, min),
        nnri_max = safe_stat(nnri, max),
        nnri_ci_95_lower = safe_quantile(nnri, 0.025),
        nnri_ci_95_upper = safe_quantile(nnri, 0.975),
        nnri_ci_99_lower = safe_quantile(nnri, 0.005),
        nnri_ci_99_upper = safe_quantile(nnri, 0.995),
        
        # VDRI statistics
        vdri_n = sum(!is.na(vdri)),
        vdri_mean = safe_stat(vdri, mean),
        vdri_median = safe_stat(vdri, median),
        vdri_min = safe_stat(vdri, min),
        vdri_max = safe_stat(vdri, max),
        vdri_ci_95_lower = safe_quantile(vdri, 0.025),
        vdri_ci_95_upper = safe_quantile(vdri, 0.975),
        vdri_ci_99_lower = safe_quantile(vdri, 0.005),
        vdri_ci_99_upper = safe_quantile(vdri, 0.995)
      )
    
    # Calculate NN distance statistics
    nn_summary <- nn_df %>%
      group_by(subtype) %>%
      summarise(
        nn_dist_n = n(),
        nn_dist_mean = mean(distance, na.rm = TRUE),
        nn_dist_sd = sd(distance, na.rm = TRUE)
      )
    
    # Calculate Voronoi area statistics
    voro_summary <- voro_df %>%
      group_by(subtype) %>%
      summarise(
        voro_area_n = n(),
        voro_area_mean = mean(area, na.rm = TRUE),
        voro_area_sd = sd(area, na.rm = TRUE)
      )
    
    # Combine all summaries
    full_summary <- metrics_summary %>%
      left_join(nn_summary, by = "subtype") %>%
      left_join(voro_summary, by = "subtype")
    
    return(full_summary)
  }
  
  # Create and save statistical summaries
  observed_summary <- create_statistical_summary(
    observed_stats,
    observed_nn_distances,
    observed_voro_areas
  )
  
  bootstrap_summary <- create_statistical_summary(
    all_bootstrap_stats,
    all_bootstrap_nn,
    all_bootstrap_areas
  )
  
  # Save summaries
  write_csv(observed_summary, file.path(output_dir, "observed_comprehensive_summary.csv"))
  write_csv(bootstrap_summary, file.path(output_dir, "bootstrap_comprehensive_summary.csv"))
  
  return(list(
    observed_stats = observed_stats,
    bootstrap_stats = all_bootstrap_stats,
    bootstrap_nn_distances = all_bootstrap_nn,
    bootstrap_voro_areas = all_bootstrap_areas,
    observed_nn_distances = observed_nn_distances,
    observed_voro_areas = observed_voro_areas,
    observed_summary = observed_summary,
    bootstrap_summary = bootstrap_summary
  ))
}

########################################################################################333
create_spatial_statistics_plot <- function(observed_stats, bootstrap_stats, output_dir) {
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  
  # Helper function to calculate significance
  get_significance_stars <- function(p_value) {
    if(is.na(p_value)) return("NA")
    if(p_value <= 0.001) return("***")
    if(p_value <= 0.01) return("**")
    if(p_value <= 0.05) return("*")
    return("ns")
  }
  
  # Helper function for safe mean calculation
  safe_mean <- function(x) {
    if(all(is.na(x))) return(NA)
    mean(x, na.rm = TRUE)
  }
  
  # Helper function for safe median calculation
  safe_median <- function(x) {
    if(all(is.na(x))) return(NA)
    median(x, na.rm = TRUE)
  }
  
  # Calculate statistics and p-values for each subtype
  stat_summary <- observed_stats %>%
    group_by(subtype) %>%
    summarize(
      obs_nnri = safe_mean(nnri),
      obs_vdri = safe_mean(vdri)
    ) %>%
    left_join(
      bootstrap_stats %>%
        group_by(subtype) %>%
        summarize(
          n_bootstrap = n(),
          p_value_nnri = mean(nnri >= safe_mean(observed_stats$nnri[observed_stats$subtype == first(subtype)]), 
                              na.rm = TRUE),
          p_value_vdri = mean(vdri >= safe_mean(observed_stats$vdri[observed_stats$subtype == first(subtype)]), 
                              na.rm = TRUE)
        ),
      by = "subtype"
    ) %>%
    mutate(
      sig_nnri = sapply(p_value_nnri, get_significance_stars),
      sig_vdri = sapply(p_value_vdri, get_significance_stars),
      label_nnri = sprintf("%s (n=%d)", sig_nnri, n_bootstrap),
      label_vdri = sprintf("%s (n=%d)", sig_vdri, n_bootstrap)
    )
  
  # Combine observed and bootstrap data
  plot_data <- bind_rows(
    observed_stats %>% mutate(type = "Observed"),
    bootstrap_stats %>% mutate(type = "Bootstrap")
  )
  
  # Function to create plot for each metric
  create_metric_plot <- function(data, metric, metric_name, summary_data, obs_metric) {
    # Remove NA values for reordering
    valid_data <- data[!is.na(data[[metric]]), ]
    if(nrow(valid_data) == 0) {
      return(ggplot() + 
               theme_void() + 
               ggtitle(paste("No valid", metric_name, "data available")))
    }
    
    # Calculate order based on median of valid values
    order_data <- valid_data %>%
      group_by(subtype) %>%
      summarize(median_val = median(.data[[metric]], na.rm = TRUE)) %>%
      arrange(median_val)
    
    # Set factor levels based on calculated order
    data$subtype <- factor(data$subtype, levels = order_data$subtype)
    
    p <- ggplot(data, aes(x = .data[[metric]], y = subtype)) +
      geom_boxplot(aes(fill = type), outlier.size = 0.4) +
      scale_fill_manual(values = c("Observed" = "black", "Bootstrap" = "grey80"))
    
    # Only add text annotations if we have valid statistics
    if(nrow(summary_data) > 0) {
      max_val <- max(data[[metric]], na.rm = TRUE)
      p <- p + geom_text(
        data = summary_data,
        aes(x = max_val * 1.1,
            y = subtype,
            label = .data[[paste0("label_", tolower(metric))]]),
        hjust = 0,
        size = 3
      )
    }
    
    p <- p + theme_minimal() +
      theme(
        legend.position = "none",
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid = element_blank(),
        axis.title.y = element_blank()
      ) +
      labs(title = paste(metric_name, "Distribution by Subtype"),
           x = metric_name)
    
    return(p)
  }
  
  # Create individual plots
  p_nnri <- create_metric_plot(plot_data, "nnri", "NNRI", stat_summary, "obs_nnri")
  p_vdri <- create_metric_plot(plot_data, "vdri", "VDRI", stat_summary, "obs_vdri")
  
  # Combine plots
  combined_plot <- p_nnri + p_vdri +
    plot_layout(guides = "collect") +
    plot_annotation(
      caption = "* p<0.05, ** p<0.01, *** p<0.001",
      theme = theme(
        plot.caption = element_text(hjust = 1, size = 8)
      )
    )
  
  # Save plot
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(output_dir, "spatial_statistics_summary.png"),
         combined_plot,
         width = 12,
         height = 15)
  
  # Save statistical summary
  write_csv(stat_summary, file.path(output_dir, "statistical_summary.csv"))
  
  return(list(
    plot = combined_plot,
    statistics = stat_summary
  ))
}
# plot_stats <- create_spatial_statistics_plot(
#   bs_results$observed_stats,
#   bs_results$bootstrap_stats,
#   bs_dir
# )  
#################################################################################
visualize_comprehensive_stats <- function(output_dir) {
  library(tidyverse)
  library(patchwork)
  
  # Load the comprehensive summaries
  observed <- read_csv(file.path(output_dir, "observed_comprehensive_summary.csv"))
  bootstrap <- read_csv(file.path(output_dir, "bootstrap_statistics.csv"))
  
  # Calculate p-values for each subtype
  calculate_p_values <- function(subtype_obs, subtype_boot) {
    # For NNRI
    p_nnri <- if(is.na(subtype_obs$nnri_mean) || all(is.na(subtype_boot$nnri))) {
      NA
    } else {
      mean(subtype_boot$nnri >= subtype_obs$nnri_mean, na.rm = TRUE)
    }
    
    # For VDRI
    p_vdri <- if(is.na(subtype_obs$vdri_mean) || all(is.na(subtype_boot$vdri))) {
      NA
    } else {
      mean(subtype_boot$vdri >= subtype_obs$vdri_mean, na.rm = TRUE)
    }
    
    return(c(p_nnri, p_vdri))
  }
  
  # Add p-values to observed statistics
  observed <- observed %>%
    group_by(subtype) %>%
    mutate(
      p_value_nnri = calculate_p_values(
        observed[observed$subtype == subtype,],
        bootstrap[bootstrap$subtype == subtype,]
      )[1],
      p_value_vdri = calculate_p_values(
        observed[observed$subtype == subtype,],
        bootstrap[bootstrap$subtype == subtype,]
      )[2]
    ) %>%
    ungroup()
  
  # Save updated observed statistics with p-values
  write_csv(observed, file.path(output_dir, "observed_comprehensive_summary_with_pvalues.csv"))
  
  bootstrap <- read_csv(file.path(output_dir, "bootstrap_comprehensive_summary.csv"))
  
  # Function to create plot for each metric (NNRI or VDRI)
  create_metric_plot <- function(observed_data, bootstrap_data, metric_prefix, title) {
    # Get column names for this metric
    mean_col <- paste0(metric_prefix, "_mean")
    median_col <- paste0(metric_prefix, "_median")
    ci95_lower <- paste0(metric_prefix, "_ci_95_lower")
    ci95_upper <- paste0(metric_prefix, "_ci_95_upper")
    ci99_lower <- paste0(metric_prefix, "_ci_99_lower")
    ci99_upper <- paste0(metric_prefix, "_ci_99_upper")
    p_value_col <- paste0("p_value_", tolower(metric_prefix))
    
    # Order subtypes by bootstrap median
    subtype_order <- bootstrap_data %>%
      filter(!is.na(.data[[median_col]])) %>%
      arrange(.data[[median_col]]) %>%
      pull(subtype)
    
    # Set factor levels for subtypes
    bootstrap_data <- bootstrap_data %>%
      mutate(subtype = factor(subtype, levels = subtype_order))
    
    observed_data <- observed_data %>%
      mutate(subtype = factor(subtype, levels = subtype_order))
    
    # Calculate actual n (non-NA values) for each subtype
    n_values <- bootstrap_data %>%
      group_by(subtype) %>%
      summarise(
        n = sum(!is.na(.data[[mean_col]]))
      )
    
    # Create the plot
    p <- ggplot() +
      # Bootstrap statistics
      geom_errorbar(data = bootstrap_data,
                    aes(x = subtype, 
                        ymin = .data[[ci99_lower]], 
                        ymax = .data[[ci99_upper]]),
                    width = 0.5,
                    color = "lightblue",
                    size = 0.5,
                    alpha = 0.5) +
      geom_errorbar(data = bootstrap_data,
                    aes(x = subtype, 
                        ymin = .data[[ci95_lower]], 
                        ymax = .data[[ci95_upper]]),
                    width = 0.3,
                    color = "blue",
                    size = 0.7,
                    alpha = 0.5) +
      # Bootstrap mean/median
      geom_point(data = bootstrap_data,
                 aes(x = subtype, y = .data[[mean_col]]),
                 color = "blue",
                 size = 2) +
      geom_point(data = bootstrap_data,
                 aes(x = subtype, y = .data[[median_col]]),
                 color = "darkblue",
                 size = 2) +
      # Observed statistics with p-value significance
      geom_point(data = observed_data,
                 aes(x = subtype, y = .data[[mean_col]]),
                 color = "red",
                 size = 3,
                 shape = 18) +
      geom_point(data = observed_data,
                 aes(x = subtype, y = .data[[median_col]]),
                 color = "darkred",
                 size = 3,
                 shape = 18) +
      # Add sample sizes and p-values
      geom_text(data = left_join(n_values, observed_data),
                aes(x = subtype, 
                    y = .data[[ci99_upper]],
                    label = sprintf("n=%d\np=%.3f", 
                                    n,
                                    .data[[p_value_col]])),
                vjust = -0.5,
                size = 2.5) +
      # Customize theme
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(hjust = 0.5)
      ) +
      labs(title = title,
           x = "Subtype",
           y = title) +
      coord_flip()
    
    return(p)
  }
  
  # Create plots for NNRI and VDRI
  p_nnri <- create_metric_plot(observed, bootstrap, "nnri", "NNRI Distribution")
  p_vdri <- create_metric_plot(observed, bootstrap, "vdri", "VDRI Distribution")
  
  # Combine plots
  combined_plot <- p_nnri + p_vdri +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Spatial Statistics Summary",
      subtitle = paste("Red diamonds: Observed values (mean/median)",
                       "Blue points = mean, Dark blue points = median", 
                       "Light blue bars = 99% CI, Blue bars = 95% CI",
                       "p-values: proportion of bootstrap values >= observed",
                       sep = "\n"),
      theme = theme(
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 10)
      )
    )
  
  # Save plot
  ggsave(file.path(output_dir, "comprehensive_statistics_summary.png"),
         combined_plot,
         width = 15,
         height = 12)
  
  return(list(
    plot = combined_plot,
    observed_with_pvalues = observed,
    bootstrap = bootstrap
  ))
}
################################################################################################
######################## EXCECUTE BOOTSTRAP############## ######################################
################################################################################################

# Create the new directory for plots
bs_dir <- file.path(out_root, "NNRI-VDRI_BootStrap")
if (!dir.exists(bs_dir)) {
  dir.create(bs_dir)
}

precomp_data <- precompute_spatial_data(rgc_df, window_obj)
bs_results <- run_bootstrap_analysis(rgc_df, precomp_data, bs_dir, n_bootstrap = n_bootstrap )

vis_results <- visualize_comprehensive_stats(bs_dir)

################################################################################################
######################## Distributions ############## ######################################
################################################################################################
library(tidyverse)
library(patchwork)

create_distribution_plots <- function(nn_dir) {
  # Create distributions directory
  dist_dir <- file.path(nn_dir, "distributions")
  dir.create(dist_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Read all data files
  observed_nn <- read_csv(file.path(nn_dir, "observed_nn_distances.csv"))
  observed_voro <- read_csv(file.path(nn_dir, "observed_voro_areas.csv"))
  bootstrap_nn <- read_csv(file.path(nn_dir, "bootstrap_nn_distances.csv"))
  bootstrap_voro <- read_csv(file.path(nn_dir, "bootstrap_voro_areas.csv"))
  
  # Get unique subtypes
  subtypes <- unique(observed_nn$subtype)
  
  # Create plots for each subtype
  for(st in subtypes) {
    # Filter data for current subtype
    obs_nn_sub <- observed_nn %>% filter(subtype == st)
    obs_voro_sub <- observed_voro %>% filter(subtype == st)
    boot_nn_sub <- bootstrap_nn %>% filter(subtype == st)
    boot_voro_sub <- bootstrap_voro %>% filter(subtype == st)
    
    # Calculate statistics for subtitle
    stats <- list(
      obs_nn_n = nrow(obs_nn_sub),
      obs_nn_na = sum(is.na(obs_nn_sub$distance)),
      obs_voro_n = nrow(obs_voro_sub),
      obs_voro_na = sum(is.na(obs_voro_sub$area)),
      boot_nn_n = nrow(boot_nn_sub),
      boot_nn_na = sum(is.na(boot_nn_sub$distance)),
      boot_voro_n = nrow(boot_voro_sub),
      boot_voro_na = sum(is.na(boot_voro_sub$area))
    )
    
    # Create subplot for NN distances - Observed
    p1 <- ggplot(obs_nn_sub, aes(x = distance)) +
      geom_histogram(fill = "blue", alpha = 0.5, bins = 30) +
      theme_minimal() +
      labs(title = "Observed NN Distances",
           subtitle = sprintf("n = %d, NA = %d", stats$obs_nn_n, stats$obs_nn_na)) +
      theme(plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8))
    
    # Create subplot for NN distances - Bootstrap
    p2 <- ggplot(boot_nn_sub, aes(x = distance)) +
      geom_histogram(fill = "red", alpha = 0.5, bins = 30) +
      theme_minimal() +
      labs(title = "Bootstrap NN Distances",
           subtitle = sprintf("n = %d, NA = %d", stats$boot_nn_n, stats$boot_nn_na)) +
      theme(plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8))
    
    # Create subplot for Voronoi areas - Observed
    p3 <- ggplot(obs_voro_sub, aes(x = area)) +
      geom_histogram(fill = "blue", alpha = 0.5, bins = 30) +
      theme_minimal() +
      labs(title = "Observed Tesselation Areas",
           subtitle = sprintf("n = %d, NA = %d", stats$obs_voro_n, stats$obs_voro_na)) +
      theme(plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8))
    
    # Create subplot for Voronoi areas - Bootstrap
    p4 <- ggplot(boot_voro_sub, aes(x = area)) +
      geom_histogram(fill = "red", alpha = 0.5, bins = 30) +
      theme_minimal() +
      labs(title = "Bootstrap Tesselation Areas",
           subtitle = sprintf("n = %d, NA = %d", stats$boot_voro_n, stats$boot_voro_na)) +
      theme(plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8))
    
    # Combine plots
    combined_plot <- (p1 + p2) / (p3 + p4) +
      plot_annotation(
        title = sprintf("Distribution Plots for %s", st),
        theme = theme(plot.title = element_text(size = 12, face = "bold"))
      )
    
    # Save the plot
    ggsave(
      filename = file.path(dist_dir, sprintf("%s_distributions.png", st)),
      plot = combined_plot,
      width = 12,
      height = 10,
      dpi = 300
    )
  }
  
  # Return the directory path where plots were saved
  return(dist_dir)
}
dist_plots_dir <- create_distribution_plots(bs_dir)
