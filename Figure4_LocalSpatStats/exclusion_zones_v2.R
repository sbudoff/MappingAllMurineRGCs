library(tidyverse)
library(patchwork)
library(sf)

################################################################################################
######################## MAIN FUNCTIONS ######################################################
################################################################################################
d_i_computer <- function(xy_df, R=100, delta_r=10, min_dist =10) {
  N <- xy_df %>%
    filter(!edge_point) %>%
    nrow()
  
  breaks <- seq(0, R, by = delta_r)
  annuli <- seq_len(length(breaks)-1)
  delta_A_i <- pi * (delta_r^2) * (2 * annuli - 1)/1e6
  
  dist_mat <- xy_df %>%
    select(x, y) %>%
    as.matrix() %>%
    dist() %>%
    as.matrix()%>%
    as.numeric() 
  
  centers <- abs(xy_df$edge_point-1)
  
  center_mat <- dist_mat*centers
  
  hist_counts <- center_mat[center_mat > min_dist & center_mat <= 100] %>%
    hist(breaks = breaks,
         include.lowest = FALSE,
         plot = FALSE)
  
  out_df <- data.frame(
    bin = paste(hist_counts$breaks[-length(hist_counts$breaks)], 
                hist_counts$breaks[-1], 
                sep="-"),
    annuli = annuli,
    N = N,
    delta_A_i = delta_A_i,
    n_i = hist_counts$counts,
    d_i = hist_counts$counts/(delta_A_i*N)
  )
  
  return(out_df)
}

Rodieck1991_DensityRecovery <- function(rgc_df_filtered, centered_df, r, subtype_ND, A, delta_r=10, 
                                        subtype_target=F) {
  
  # Get all unique subtypes
  subtypes <- unique(rgc_df_filtered$Prediction)
  
  # Print initial message
  cat("Starting Rodieck DRP analysis...\n")
  # Create progress bar
  pb <- progress::progress_bar$new(
    format = "Processing subtype :what [:bar] :percent eta: :eta",
    total = length(subtypes),
    clear = FALSE,
    width = 60
  )
  
  # Handle subtype_target logic
  if (subtype_target == "all") {
    target_subtypes <- subtypes
  } else if (subtype_target == F) {
    target_subtypes <- NULL  # Will be set per iteration
  } else {
    target_subtypes <- subtype_target
  }
  
  # delta A_i = area of annulus i = pi*delta_r^2 (2i-1)
  radii <- seq(0,r, delta_r) # set of radii for each annulus
  
  # Rodieck performed all his analysis on single subtypes, 
  #thus to get valid recovery values I must do the same
  # with subtype_target I will build in flexibility for heterologous analysis
  rodieck_df <- data.frame(
    subtype = character(),
    N = numeric(),        # total number of points
    D = numeric(),        # overall density
    
    annulus = numeric(),  # annulus distance
    n_i = numeric(),      # number of points in annulus i
    d_i = numeric(),      # distance for annulus i
    lambda_i = numeric(),  # density for annulus i
    delta_V_i = numeric(),  # "dip" from mean annulus i
    V_e = numeric(), # average magnitude of the dead space
    r_e = numeric() # effective radius of a cell
  )
  
  for (subtype in subtypes) {
    
    # Update progress bar
    pb$tick(tokens = list(what = subtype))
    
    # Print detailed progress
    cat(sprintf("\nProcessing subtype %s...\n", subtype))
    
    subtype_ND_i <- subtype_ND %>%
      filter(Prediction==subtype)
    
    # Determine which subtypes to compare against
    comparison_subtypes <- if (is.null(target_subtypes)) {
      c(subtype)  # Compare only against same subtype
    } else {
      target_subtypes  # Compare against specified subtypes or all subtypes
    }
    
    
    cat(sprintf("  Calculating distances for comparison with %d subtypes...\n", 
                length(comparison_subtypes)))
    
    # Use the efficient distance calculator to get n_i and d_i values
    distance_data <- d_i_computer(
      rgc_df_filtered %>% 
        filter(Prediction %in% comparison_subtypes) %>%
        select(x, y, edge_point),
      R = r,
      delta_r = delta_r
    )
    
    cat("  Computing density recovery profile...\n")
    
    # Extract values from the distance calculator output
    n_i <- distance_data$n_i
    d_i <- distance_data$d_i
    delta_A_i <- distance_data$delta_A_i
    
    # N = number of points in A
    N <- unique(distance_data$N)
    
    # D = density of points in the region, equal to N/A;
    D <- mean(d_i)
    
    # lambda_i = expected number of points in annulus i, assuming a Poisson distribution
    lambda_i <- N*D*delta_A_i
    
    # Delta_V_i = the "volume" of the dip at bin position i scaled by the annular area
    delta_V_i = 1e6*delta_A_i*(D-d_i)#(lambda_i - n_i)/N
    
    # V_e = the total volume of the dip, obtained by adding the volumes of the annuli until the density value for an annulus exceeds the mean density
    # ie until the calculated volume for the annulus becomes negative
    V_e <- 0
    for (i in seq_along(distance_data$annuli)) {
      if (n_i[i] < lambda_i[i]) {  # Continue until ni exceeds λi
        V_e <- V_e + delta_V_i[i]
      } else {
        break
      }
    }
    
    #r_e = the effective radius. The total decrease in cells within the dip is 
    #then made equivalent to one composed of a step change from zero to the mean density at the effective radius.
    r_e = sqrt(V_e/(pi*D))
    
    cat(sprintf("  Completed analysis for subtype %s (r_e = %.2f)\n", subtype, r_e))
    
    temp_df <- data.frame(subtype = subtype, 
                          N = N,
                          D = D,
                          annulus = distance_data$annuli, 
                          n_i = n_i, 
                          d_i = d_i,
                          lambda_i = lambda_i,
                          delta_V_i = delta_V_i,
                          V_e = V_e,
                          r_e = r_e)
    
    rodieck_df <- rbind(rodieck_df, temp_df) 
  }
  
  rodieck_df <- rodieck_df %>%
    mutate(delta_r = delta_r) 
  
  cat("\nRodieck analysis completed successfully!\n")
  
  return(rodieck_df)
}

################################################################################################
######################## HELPER FUNCTIONS ######################################################
################################################################################################


load_shuffle_results <- function(delta_r, shuffle_dir) {
  # Get list of all CSV files matching our pattern
  shuffle_files <- list.files(
    path = shuffle_dir,
    pattern = paste0("Effective_radius_calculations_binWidth", delta_r, "seed_.*\\.csv"),
    full.names = TRUE
  )
  
  # Check if we found any files
  if(length(shuffle_files) == 0) {
    stop("No shuffle files found matching the specified pattern")
  }
  
  # Read and combine all files
  shuffle_results <- shuffle_files %>%
    map_df(~{
      read_csv(.x, show_col_types = FALSE) %>%
        # Extract r_e for each subtype and seed
        group_by(subtype, random_seed) %>%
        slice(1) %>%  # Take first row since r_e should be same within groups
        select(subtype, random_seed, r_e) %>%
        ungroup()
    })
  
  return(shuffle_results)
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

identify_edge_points <- function(rgc_df, slide_id, lasso_id, r = 100) {
  # Create the full study region identifier
  study_region_id <- paste0(slide_id, "_", lasso_id)
  
  # Filter the dataframe for the specified region
  region_df <- rgc_df %>%
    filter(study_region == study_region_id)
  
  if(nrow(region_df) == 0) {
    warning(sprintf("No points found in region %s", study_region_id))
    return(NULL)
  }
  
  # Get the lasso polygon for this region
  lasso_polygon <- st_sf(geometry = st_sfc(lasso_list[[lasso_id]]))
  
  # Convert points to sf object
  points_sf <- st_as_sf(region_df, coords = c("x", "y"))
  
  # Get the boundary of the polygon
  lasso_boundary <- st_cast(lasso_polygon, "LINESTRING")
  
  # Calculate distances from points to the boundary
  distances <- st_distance(points_sf, lasso_boundary)
  
  # Add edge_point column to the original filtered dataframe
  region_df %>%
    mutate(edge_point = as.vector(distances) <= r)
}

# Function to find minimum exclusion distance for a sample
find_min_distance <- function(data) {
  # Sample with replacement from the centered arrangements
  sampled_centers <- sample(unique(data$center_cell), 
                            size = length(unique(data$center_cell)), 
                            replace = TRUE)
  
  # Get all points from these sampled arrangements
  sample_data <- data %>%
    filter(center_cell %in% sampled_centers)
  
  # Create distance bins
  max_r <- max(sample_data$r)
  bins <- seq(0, max_r, by = bin_width)
  
  # Find first bin with points
  for(bin_start in bins) {
    bin_end <- bin_start + bin_width
    points_in_bin <- sample_data %>%
      filter(r >= bin_start & r < bin_end) %>%
      nrow()
    
    if(points_in_bin > 0) {
      return(bin_start)
    }
  }
  return(NA)
}
compute_soma_exclusion <- function(rgcs_centered, subtype, n_bootstrap = 1000, bin_width = 1, conf_level = 0.95) {
  # Filter data for the specific subtype centers
  subtype_data <- rgcs_centered %>%
    filter(center == subtype) %>%
    # Calculate radial distance for each point
    mutate(r = sqrt(x^2 + y^2))
  
  # Perform bootstrap
  bootstrap_results <- replicate(n_bootstrap, 
                                 find_min_distance(subtype_data))
  
  # Calculate statistics
  mean_exclusion <- mean(bootstrap_results, na.rm = TRUE)
  ci_lower <- quantile(bootstrap_results, (1 - conf_level)/2, na.rm = TRUE)
  ci_upper <- quantile(bootstrap_results, 1 - (1 - conf_level)/2, na.rm = TRUE)
  
  # Return results as a data frame
  result <- data.frame(
    subtype = subtype,
    mean_soma_exclusion = mean_exclusion,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    n_arrangements = length(unique(subtype_data$center_cell)),
    n_bootstrap = n_bootstrap,
    bin_width = bin_width
  )
  
  return(result)
}

# Function to analyze all subtypes and create visualization
analyze_all_soma_exclusions <- function(rgcs_centered, bin_width = 1, n_bootstrap = 1000) {
  subtypes <- unique(rgcs_centered$center)
  
  # Calculate for all subtypes
  results <- do.call(rbind, lapply(subtypes, function(subtype) {
    compute_soma_exclusion(rgcs_centered, subtype, 
                           bin_width = bin_width, 
                           n_bootstrap = n_bootstrap)
  }))
  
  # Create visualization
  p <- ggplot(results, aes(x = reorder(subtype, mean_soma_exclusion), 
                           y = mean_soma_exclusion)) +
    geom_point() +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
    coord_flip() +
    theme_minimal() +
    labs(x = "Cell Type",
         y = "Soma Exclusion Zone (μm)",
         title = "Minimum Cell-Cell Distance by Type",
         subtitle = sprintf("Based on %d bootstrap samples\nError bars show 95%% CI", n_bootstrap)) +
    theme(axis.text.y = element_text(size = 8))
  
  return(list(results = results, plot = p))
}

################################################################################################
######################## PLOTTING FUNCTIONS ####################################################
################################################################################################
plot_exclusion_analysis <- function(rgcs_centered, rodieck_df, subtype, is_paired = FALSE, bootstrap_df = NULL) {
  circle_rad <- ceiling(max(rgcs_centered$x, rgcs_centered$y))
  
  if (is_paired) {
    # Split pair into components
    subtypes <- strsplit(subtype, "-")[[1]]
    subtype1 <- subtypes[1]
    subtype2 <- subtypes[2]
    
    # Filter data for the paired subtypes
    spatial_data <- rgcs_centered %>%
      filter(center == subtype) %>%
      mutate(
        is_subtype1 = Prediction == subtype1,
        is_subtype2 = Prediction == subtype2
      )
  } else {
    spatial_data <- rgcs_centered %>%
      filter(center == subtype) %>%
      mutate(is_target = Prediction == subtype)
  }
  
  # Calculate normalized density values
  density_data <- rodieck_df %>%
    filter(subtype == !!subtype) %>%
    mutate(
      density = d_i,
      bin_start = (annulus - 1) * delta_r,
      bin_end = annulus * delta_r,
      bin_center = (bin_start + bin_end) / 2
    )
  
  r_e_value <- unique(density_data$r_e)
  D_value <- unique(density_data$D)
  N_val <- sum(density_data$n_i)
  
  # Create data frame for annuli circles
  annuli_radii <- rodieck_df %>%
    mutate(bin_end = annulus * delta_r) %>%
    pull(bin_end) %>%
    unique()
  circles_df <- data.frame(
    r = annuli_radii,
    x0 = 0,
    y0 = 0
  )
  
  # Spatial distribution plot
  p1 <- ggplot() +
    ggforce::geom_circle(
      data = circles_df,
      aes(x0 = x0, y0 = y0, r = r),
      color = "grey80",
      linetype = "dashed",
      size = 0.2
    )
  
  if (is_paired) {
    # Add points for paired data
    p1 <- p1 +
      # Other cells (neither subtype)
      geom_point(
        data = filter(spatial_data, Prediction != subtype1 & Prediction != subtype2),
        aes(x = x, y = y),
        color = "grey50",
        alpha = 0.3,
        size = 0.5
      ) +
      # Subtype 1
      geom_point(
        data = filter(spatial_data, is_subtype1),
        aes(x = x, y = y),
        color = "red",
        alpha = 0.7,
        size = 1
      ) +
      # Subtype 2
      geom_point(
        data = filter(spatial_data, is_subtype2),
        aes(x = x, y = y),
        color = "blue",
        alpha = 0.7,
        size = 1
      )
  } else {
    # Original single subtype plotting
    p1 <- p1 +
      geom_point(
        data = filter(spatial_data, !is_target),
        aes(x = x, y = y),
        color = "grey50",
        alpha = 0.3,
        size = 0.5
      ) +
      geom_point(
        data = filter(spatial_data, is_target),
        aes(x = x, y = y),
        color = "red",
        alpha = 0.7,
        size = 1
      )
  }
  
  p1 <- p1 +
    coord_fixed() +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 12),
      axis.text = element_text(size = 8)
    ) +
    labs(
      title = paste("Spatial Distribution -", subtype),
      x = "Distance (μm)",
      y = "Distance (μm)"
    ) +
    lims(
      x = c(-circle_rad, circle_rad),
      y = c(-circle_rad, circle_rad)
    )
  
  # Density recovery profile plot
  p2 <- ggplot(density_data) +
    geom_rect(aes(
      xmin = bin_start,
      xmax = bin_end,
      ymin = 0,
      ymax = density
    ), fill = "grey50") +
    geom_hline(yintercept = D_value, linetype = "dashed", color = "blue", alpha = 0.5) +
    geom_vline(xintercept = r_e_value, color = "red", linetype = "dashed") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12),
      axis.text = element_text(size = 8)
    ) +
    labs(
      title = "Density Recovery Profile",
      subtitle = sprintf("Effective Radius: %.1f μm, Observations: %.0f", r_e_value, N_val),
      x = "Radius (μm)",
      y = "Density (cells/mm²)"
    ) +
    ylim(0, max(max(density_data$density) * 1.1, D_value*1.1))
  
  # If bootstrap data is provided, create p3
  if (!is.null(bootstrap_df)) {
    # Filter bootstrap data for the specific subtype
    bootstrap_density_data <- bootstrap_df %>%
      filter(subtype == !!subtype) %>%
      mutate(
        density = d_i,
        bin_start = (annulus - 1) * delta_r,
        bin_end = annulus * delta_r,
        bin_center = (bin_start + bin_end) / 2
      )
    
    bootstrap_r_e <- unique(bootstrap_density_data$r_e)
    bootstrap_D <- unique(bootstrap_density_data$D)
    bootstrap_N <- sum(bootstrap_density_data$n_i)
    
    p3 <- ggplot(bootstrap_density_data) +
      geom_rect(aes(
        xmin = bin_start,
        xmax = bin_end,
        ymin = 0,
        ymax = density
      ), fill = "grey50") +
      geom_hline(yintercept = bootstrap_D, linetype = "dashed", color = "blue", alpha = 0.5) +
      geom_vline(xintercept = bootstrap_r_e, color = "red", linetype = "dashed") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 12),
        axis.text = element_text(size = 8)
      ) +
      labs(
        title = "Bootstrap Density Recovery",
        subtitle = sprintf("Bootstrap R_e: %.1f μm, Observations: %.0f", bootstrap_r_e, bootstrap_N),
        x = "Radius (μm)",
        y = "Density (cells/mm²)"
      ) +
      ylim(0, max(max(bootstrap_density_data$density) * 1.1, bootstrap_D*1.1))
    
    # Combine all three plots
    combined_plot <- p1 + p2 + p3 +
      plot_layout(widths = c(1, 1, 1)) &
      theme(plot.margin = margin(5, 5, 5, 5))
  } else {
    # Original two-plot layout
    combined_plot <- p1 + p2 +
      plot_layout(widths = c(1, 1)) &
      theme(plot.margin = margin(5, 5, 5, 5))
  }
  
  return(combined_plot)
}

plot_edge_analysis <- function(df, lasso_id, r) {
  # Get the lasso polygon
  lasso_polygon <- st_sf(geometry = st_sfc(lasso_list[[lasso_id]]))
  
  # Calculate summary statistics
  summary <- df %>%
    group_by(edge_point, Prediction) %>%
    summarize(N = n(), .groups = "drop") %>%
    group_by(Prediction) %>%
    mutate(P = N/sum(N))
  
  # Spatial plot
  p1 <- ggplot() +
    # Plot non-edge points (larger, filled circles)
    geom_point(data = filter(df, !edge_point),
               aes(x = x, y = y, color = Prediction), 
               size = 1,
               shape = 16) +
    # Plot edge points (hollow circles)
    geom_point(data = filter(df, edge_point),
               aes(x = x, y = y, color = Prediction), 
               size = 1,
               shape = 1) +
    geom_sf(data = lasso_polygon, 
            fill = NA, 
            color = "black", 
            linewidth = 0.5) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text = element_text(size = 8),
          plot.title = element_text(size = 10)) +
    labs(title = paste("Edge Points Detection\nDistance threshold:", r, "units"))
  
  # Prepare data for stacked bar plot
  counts_data <- summary %>%
    arrange(desc(Prediction))
  
  # Stacked bar plot with custom fill and outline
  p2 <- ggplot(counts_data, aes(x = N, y = Prediction)) +
    # Create the base layer with all data
    geom_col(aes(fill = edge_point), width = 0.7) +
    scale_fill_manual(values = c("TRUE" = "white", "FALSE" = "black")) +
    # Add black borders to the stacked bars
    geom_col(data = counts_data, 
             aes(group = edge_point),
             fill = NA, 
             color = "black",
             width = 0.7) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8),
          axis.title.y = element_blank(),
          legend.position = "none") +
    labs(x = "Count",
         title = "Cell Counts")
  
  # Proportion plot
  p3 <- ggplot(filter(counts_data, !edge_point), 
               aes(x = P, y = Prediction)) +
    geom_col(fill = "grey50") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none") +
    labs(x = "Centered Cell Proportion",
         title = "Centered Cells Proportion") +
    scale_x_continuous(limits = c(0, 1))
  
  # Combine plots
  combined_plot <- p1 + p2 + p3 +
    plot_layout(ncol = 3, widths = c(1.5, 1, 1)) &
    theme(plot.margin = margin(5, 5, 5, 5))
  
  return(combined_plot)
}

plot_centered_subtypes <- function(data) {
  # Get unique subtypes and sort them by numeric prefix
  subtypes <- unique(data$center)
  # Extract numeric prefix and use it for sorting
  sorted_subtypes <- subtypes[order(as.numeric(gsub("^(\\d+).*", "\\1", subtypes)))]
  
  # Create plot list
  plot_list <- lapply(sorted_subtypes, function(subtype) {
    # Filter data for this subtype's centers
    subtype_data <- data %>%
      filter(center == subtype) %>%
      mutate(is_same = Prediction == center)
    
    # Plot with different layers for matching and non-matching points
    ggplot() +
      # Non-matching points (grey, smaller, more transparent)
      geom_point(data = filter(subtype_data, !is_same),
                 aes(x = x, y = y),
                 color = "grey50",
                 alpha = 0.2,
                 size = 0.1) +
      # Matching points (red, larger, less transparent)
      geom_point(data = filter(subtype_data, is_same),
                 aes(x = x, y = y),
                 color = "red",
                 alpha = 0.8,
                 size = 0.3) +
      coord_fixed() +  # Ensure aspect ratio is 1:1
      theme_minimal() +
      theme(legend.position = "none",
            plot.title = element_text(size = 10),
            axis.text = element_text(size = 8)) +
      labs(title = subtype)
  })
  
  # Name the plots
  names(plot_list) <- sorted_subtypes
  
  # Combine plots using patchwork
  combined_plot <- wrap_plots(plot_list, ncol = 9) +
    plot_annotation(
      title = "Cell Type Distribution Around Center Points",
      theme = theme(plot.title = element_text(size = 14, hjust = 0.5))
    )
  
  return(combined_plot)
}
################################################################################################
######################## LOAD DATASETS #########################################################
################################################################################################


out_root <- '/home/sam/FinalRGC_xenium'
lasso_root <- '/media/sam/New Volume/Xenium_Data/HQ_NearestNeighbor_Zones'
exclusion_dir <- file.path(out_root, "exclusion_zones")
rgc_path <- file.path(out_root, "rgc_expMat_with_Studyregions_transformedCoords.csv")

rgc_df <- read_csv(rgc_path) %>%
  drop_na() %>%
  mutate(study_region = factor(study_region))

lasso_files <- list.files(lasso_root,
                          pattern = "\\d+_HQ_Lasso_coordinates\\.csv$",
                          full.names = TRUE)

################################################################################################
######################## IDENTIFY VALID CELLS ##################################################
################################################################################################

# Create the new directory for plots
plot_dir <- file.path(exclusion_dir, "Edge_plots")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}


# Radius of study region
r <- 100

regions <- unique(rgc_df$study_region)
csv_filename <- file.path(plot_dir, paste0("rgc_exp_mat_edge_cells_", r,"um.csv"))

if (!file.exists(csv_filename)) {
  cat("Identifying edge points")
  start = T
  for (region in regions) {
    slide_id <- as.numeric(strsplit(region, "_")[[1]][1])  # This will get the first part before any "_"
    lasso_id <- as.numeric(tail(strsplit(region, "_")[[1]], 1))  # This gets the last number
    cat(sprintf("\nProcessing slide %d...\n", slide_id))
    
    # Load lasso regions
    file <- lasso_files[grep(paste0(slide_id, "_HQ_Lasso"), lasso_files)]
    lasso_list <- load_lasso_regions(file)
  
    lasso_id <- paste0("HQ_Lasso_", lasso_id)
    temp_df <- identify_edge_points(rgc_df, slide_id, lasso_id, r = r)
    p<-plot_edge_analysis(temp_df, lasso_id, r)
  
    # Create filename
    filename <- file.path(plot_dir, paste0(slide_id, "_", lasso_id,"_edge_cells.png"))
    # Save plot with ggsave
    ggsave(
      filename = filename,
      plot = p,
      width = 10,
      height = 5.5,
      dpi = 300,
      units = "in"
    )
    
    if (start){
      rgc_df_filtered <- temp_df
      start = F
    } else {
      rgc_df_filtered <- rbind(rgc_df_filtered, temp_df)
    }
  }
  write_csv(rgc_df_filtered, file=csv_filename)
  cat("Edges computed, and saved")
} else {
  rgc_df_filtered <- read_csv(csv_filename)
  cat("Edges previously computed, and reloaded")
}

summary <- rgc_df_filtered %>%
  group_by(edge_point, Prediction) %>%
  summarize(N = n(), .groups = "drop") %>%
  group_by(Prediction) %>%
  mutate(P = N/sum(N))

print(rgc_df_filtered %>%
      group_by(edge_point) %>%
      summarize(N = n(), .groups = "drop") %>%
      mutate(P = N/sum(N)))

start = T  
subtypes <- unique(rgc_df$Prediction)
for (subtype in subtypes) {
  rgc_s_centers <- rgc_df_filtered %>%
    filter(Prediction==subtype, !edge_point) %>%
    select(cell, x, y, Prediction, study_region) 
  
  regions <- unique(rgc_s_centers$study_region)
  
  
  for (i in 1:nrow(rgc_s_centers)) {
    rgc <- rgc_s_centers[i,]
    
    rgc_i <- rgc_df_filtered %>%
      filter(edge_point, study_region == rgc$study_region) %>%
      select(cell, x, y, Prediction, study_region) %>%
      mutate(x = x - rgc$x,
             y = y- rgc$y,
             center = rgc$Prediction,
             center_cell = rgc$cell)%>%
      filter(sqrt(x^2 + y^2) <= r)
      
    
    if (start){
      rgcs_centered <- rgc_i
      start = F
    } else {
      rgcs_centered <- rbind(rgcs_centered, rgc_i)
    }
  }
}


(p <- plot_centered_subtypes(drop_na(rgcs_centered)))

write.csv(rgcs_centered, file = file.path(exclusion_dir, "centered_RGCs.csv"))
# Create filename
filename <- file.path(exclusion_dir, "All_RGC_Autocorrelation_plots.png")
# Save plot with ggsave
ggsave(
  filename = filename,
  plot = p,
  width = 12,
  height = 12,
  dpi = 600,
  units = "in"
)

################################################################################################
######################## COMPUTE RODIECK PARAMS ################################################
################################################################################################


# To get to mm^2 divide by 1million

# A = area of the region containing the reference points
areas_by_region <- rgc_df_filtered %>%
  # filter(!edge_point) %>%
  select(study_region, x, y) %>%
  group_by(study_region) %>%
  group_modify(~{
    # Create concave hull
    hull <- concaveman::concaveman(cbind(.x$x, .x$y))
    # Convert to polygon and calculate area
    area <- st_area(st_polygon(list(hull)))
    data.frame(area = area)
  })
A <- sum(areas_by_region$area)
A <- A / 1e6 # to get from um^2 to mm^2

































delta_r <- 5

# Create the new directory for plots
plot_dir <- file.path(exclusion_dir, "Homotypical_plots")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}
csv_filename <- file.path(plot_dir, paste0("Homotypical_effective_radius_calculations_binWidth", delta_r, ".csv"))

if (!file.exists(csv_filename)) {
  rodieck_df <- Rodieck1991_DensityRecovery(rgc_df_filtered, rgcs_centered, r, subtype_ND, A, delta_r)
  
  # Loop over subtypes and create plots
  for (subtype in subtypes) {
    p <- plot_exclusion_analysis(rgcs_centered, rodieck_df, subtype)
    
    # Create filename
    filename <- file.path(plot_dir, paste0(subtype, "_effective_radius.png"))
    
    # Save plot with ggsave
    ggsave(
      filename = filename,
      plot = p,
      width = 12,
      height = 6,
      dpi = 300,
      units = "in"
    )
  }
  
  write_csv(rodieck_df, file = csv_filename)
} else {
  rodieck_df <- read.csv(csv_filename)
}
################################################################################################
######################## SHUFFLED DENSITY PROFILE VIA BOOTSTRAP ################################
################################################################################################

BOOTSTRAP_SHUFFLES <- 1000
shuffle_dir <- file.path(exclusion_dir, "EffectiveRadius_SHUFFLES")
if (!dir.exists(shuffle_dir)) {
  dir.create(shuffle_dir)
}

# Initialize progress bar
pb <- txtProgressBar(min = 0, max = BOOTSTRAP_SHUFFLES, style = 3)

for (s in 1:BOOTSTRAP_SHUFFLES) {
  seed <- 18 + s
  csv_filename <- file.path(shuffle_dir, 
                            paste0("Effective_radius_calculations_binWidth", 
                                   delta_r, "seed_", seed, ".csv"))
  
  if (!file.exists(csv_filename)) {
    set.seed(seed)
    
    # Shuffle positions while maintaining x-y pairs
    rgc_df_shuffled <- rgc_df_filtered %>%
      group_by(study_region) %>%
      mutate(
        shuffle_idx = sample(n()),
        x = x[shuffle_idx],
        y = y[shuffle_idx]
      ) %>%
      select(-shuffle_idx) %>%
      ungroup()
    
    rodieck_df_shuffled <- Rodieck1991_DensityRecovery(
      rgc_df_shuffled, 
      rgcs_centered, 
      r,
      subtype_ND, 
      A, 
      delta_r
    ) %>%
      mutate(random_seed = seed)
    
    write_csv(rodieck_df_shuffled, file = csv_filename)
  } else {
    rodieck_df_shuffled<- read.csv(csv_filename)
  }
  
  # Update progress bar
  setTxtProgressBar(pb, s)
}

# Close progress bar
close(pb)


r_e_sum <- rodieck_df %>%
  select(subtype, r_e) %>%
  distinct() %>%
  mutate(random_seed = 0)

shuffled_recoveries <- load_shuffle_results(delta_r, shuffle_dir)
shuffled_recoveries <- rbind(shuffled_recoveries,
                             r_e_sum)

analyze_and_plot_bootstraps <- function(combined_data) {
  # Calculate statistics for each subtype
  bootstrap_stats <- combined_data %>%
    group_by(subtype) %>%
    reframe(
      mean_r_e = mean(r_e[random_seed != 0], na.rm = TRUE),
      min_r_e = min(r_e[random_seed != 0], na.rm = TRUE),
      max_r_e = max(r_e[random_seed != 0], na.rm = TRUE),
      ci_lower99 = quantile(r_e[random_seed != 0], 0.005, na.rm = TRUE),
      ci_upper99 = quantile(r_e[random_seed != 0], 0.995, na.rm = TRUE),
      ci_lower = quantile(r_e[random_seed != 0], 0.025, na.rm = TRUE),
      ci_upper = quantile(r_e[random_seed != 0], 0.975, na.rm = TRUE),
      real_r_e = r_e[random_seed == 0][1],  # Direct indexing of the real value
      p_value = if(any(random_seed == 0)) {
        real_val <- r_e[random_seed == 0][1]
        boot_mean <- mean(r_e[random_seed != 0], na.rm = TRUE)
        mean(abs(r_e[random_seed != 0] - boot_mean) >= abs(real_val - boot_mean), na.rm = TRUE)
      } else {
        NA_real_
      }
    ) %>%
    # Order subtypes by dd number
    arrange(desc(subtype))
  
  # Create the plot
  p <- ggplot(combined_data, aes(y = reorder(subtype, desc(subtype)))) +
    # Add bootstrap points
    geom_jitter(data = . %>% filter(random_seed != 0, !is.na(r_e)),
               aes(x = r_e),
               color = "grey80",
               alpha = 0.3,
               size = 0.5) +
    # Add real values
    geom_point(data = . %>% filter(random_seed == 0, !is.na(r_e)),
               aes(x = r_e),
               color = "red",
               size = 2) +
    # Add confidence intervals
    geom_segment(data = bootstrap_stats %>% filter(!is.na(ci_lower), !is.na(ci_upper)),
                 aes(x = ci_lower, xend = ci_upper,
                     yend = subtype),
                 color = "grey50",
                 size = 1,
                 alpha = 0.8) +
    # Customize appearance
    theme_minimal() +
    labs(x = "Effective Radius",
         y = "Cell Type") +
    theme(axis.text.y = element_text(size = 8),
          panel.grid.minor = element_blank())
  
  # Return both plot and statistics
  return(list(
    plot = p,
    statistics = bootstrap_stats
  ))
}

# Usage example:
results <- analyze_and_plot_bootstraps(shuffled_recoveries)

bootstrap_results <- results$statistics


plot_dir <- file.path(exclusion_dir, "EffectiveRadius_bootstrap_plots")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

write_csv(bootstrap_results, file.path(plot_dir,'summary_bootstrap_effective_radius.csv'))
write_csv(shuffled_recoveries, file.path(plot_dir,'values_bootstrap_effective_radius.csv'))

# Loop over subtypes and create plots
for (subtype in subtypes) {
  p <- plot_exclusion_analysis(rgcs_centered, rodieck_df, subtype, bootstrap_df = rodieck_df_shuffled)
  
  # Create filename
  filename <- file.path(plot_dir, paste0(subtype, "_effective_radius_bootstrap.png"))
  
  # Save plot with ggsave
  ggsave(
    filename = filename,
    plot = p,
    width = 12,
    height = 6,
    dpi = 300,
    units = "in"
  )
}




# View plot
print(results$plot)

filename <- file.path(plot_dir, paste0(subtype, "Effective_radius_Aggregate_Summary.png"))
# Save plot with ggsave
ggsave(
  filename = filename,
  plot = results$plot,
  width = 8,
  height = 12,
  dpi = 300,
  units = "in"
)

################################################################################################
######################## COMPUTE SOMA EXCLUSION ZONE ###########################################
################################################################################################

plot_dir <- file.path(exclusion_dir, "Soma_EffectiveRadius_plots")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}
csv_filename <- file.path(plot_dir, paste0("Soma_effective_radius_calculations_binWidth", delta_r, ".csv"))

if (!file.exists(csv_filename)) {
  rodieck_df <- Rodieck1991_DensityRecovery(rgc_df_filtered, rgcs_centered, r, 
                                            subtype_ND, A, delta_r,
                                            subtype_target = "all"
                                            )
  
  # Loop over subtypes and create plots
  for (subtype in subtypes) {
    p <- plot_exclusion_analysis(rgcs_centered, rodieck_df, subtype)
    
    # Create filename
    filename <- file.path(plot_dir, paste0(subtype, "_effective_radius.png"))
    
    # Save plot with ggsave
    ggsave(
      filename = filename,
      plot = p,
      width = 12,
      height = 6,
      dpi = 300,
      units = "in"
    )
  }
  
  write_csv(rodieck_df, file = csv_filename)
} else {
  rodieck_df <- read.csv(csv_filename)
}



################################################################################################
######################## COMPUTE SOMA EXCLUSION ZONE VIA BOOTSTRAP #############################
################################################################################################
csv_path <- file.path(exclusion_dir, "soma_exclusion_analysis.csv")
if (!file.exists(csv_path)) {
  soma_analysis <- analyze_all_soma_exclusions(rgcs_centered, 
                                               bin_width = 1,  # in μm
                                               n_bootstrap = 1000)
  print(soma_analysis$plot)
  soma_exc_df<-soma_analysis$results
  # Save results
  write_csv(soma_exc_df, 
            file = csv_path)
  # Save plot
  ggsave(file.path(exclusion_dir, "soma_exclusion_analysis.png"),
         plot = soma_analysis$plot,
         width = 8, height = 12, dpi = 300)
} else {
  soma_exc_df <- read_csv(csv_path)
}



analyze_spacing_relationships <- function(rgcs_centered, rodieck_df, soma_results,
                                          rgc_df, n_bootstrap=1000) {
  # Get effective radius from Rodieck analysis
  recovery_results <- rodieck_df %>%
    group_by(subtype) %>%
    summarize(
      effective_radius = first(r_e),
      density = first(D)
    )
  
  # Calculate volume statistics for each subtype
  volume_stats <- rgc_df %>%
    select(Prediction, volume) %>%
    group_by(Prediction) %>%
    summarize(
      mean_volume = mean(volume, na.rm = TRUE),
      sd_volume = sd(volume, na.rm = TRUE),
      n = n(),
      se_volume = sd_volume / sqrt(n)
    )
  
  # Combine all results
  combined_results <- soma_results %>%
    left_join(recovery_results, by = c("subtype" = "subtype")) %>%
    left_join(volume_stats, by = c("subtype" = "Prediction")) %>%
    mutate(
      spacing_ratio = effective_radius / mean_soma_exclusion,
      free_space = effective_radius - mean_soma_exclusion
    )
  
  # Perform linear regressions
  lm_fit_spacing <- lm(effective_radius ~ mean_soma_exclusion, data = combined_results)
  r_squared_spacing <- summary(lm_fit_spacing)$r.squared
  
  lm_fit_vol_exclusion <- lm(mean_soma_exclusion ~ mean_volume, data = combined_results)
  r_squared_vol_exclusion <- summary(lm_fit_vol_exclusion)$r.squared
  
  lm_fit_vol_effective <- lm(effective_radius ~ mean_volume, data = combined_results)
  r_squared_vol_effective <- summary(lm_fit_vol_effective)$r.squared
  
  # 1. Original scatter plot comparing exclusion zone and effective radius
  p1 <- ggplot(combined_results, 
               aes(x = mean_soma_exclusion, y = effective_radius)) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    geom_point(aes(size = density), alpha = 0.6) +
    geom_errorbar(aes(ymin = effective_radius, ymax = effective_radius), 
                  width = 0.5, color = "grey50", alpha = 0.5) +
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), 
                  width = 0.5, color = "grey50", alpha = 0.5) +
    geom_smooth(method = "lm", color = "blue", alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", alpha = 0.5) +
    ggrepel::geom_text_repel(aes(label = subtype), size = 3) +
    annotate("text", x = min(combined_results$mean_soma_exclusion), 
             y = max(combined_results$effective_radius),
             label = sprintf("R² = %.3f", r_squared_spacing),
             hjust = 0, vjust = 1) +
    labs(x = "Soma Exclusion Zone (μm)",
         y = "Effective Radius (μm)",
         size = "Cell Density",
         title = "Physical vs. Functional Cell Spacing",
         subtitle = "Dashed line indicates 1:1 ratio")
  
  # 2. Volume vs Exclusion Zone
  p2 <- ggplot(combined_results, 
               aes(x = mean_volume, y = mean_soma_exclusion)) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    geom_point(alpha = 0.6) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                  width = 0.5, color = "grey50", alpha = 0.5) +
    geom_errorbarh(aes(xmin = mean_volume - se_volume, 
                       xmax = mean_volume + se_volume), 
                   height = 0.5, color = "grey50", alpha = 0.5) +
    geom_smooth(method = "lm", color = "blue", alpha = 0.2) +
    ggrepel::geom_text_repel(aes(label = subtype), size = 3) +
    annotate("text", x = min(combined_results$mean_volume), 
             y = max(combined_results$mean_soma_exclusion),
             label = sprintf("R² = %.3f", r_squared_vol_exclusion),
             hjust = 0, vjust = 1) +
    labs(x = "Cell Volume (μm³)",
         y = "Soma Exclusion Zone (μm)",
         title = "Cell Volume vs Physical Spacing")
  
  # 3. Volume vs Effective Radius
  p3 <- ggplot(combined_results, 
               aes(x = mean_volume, y = effective_radius)) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    geom_point(alpha = 0.6) +
    geom_errorbarh(aes(xmin = mean_volume - se_volume, 
                       xmax = mean_volume + se_volume), 
                   height = 0.5, color = "grey50", alpha = 0.5) +
    geom_smooth(method = "lm", color = "blue", alpha = 0.2) +
    ggrepel::geom_text_repel(aes(label = subtype), size = 3) +
    annotate("text", x = min(combined_results$mean_volume), 
             y = max(combined_results$effective_radius),
             label = sprintf("R² = %.3f", r_squared_vol_effective),
             hjust = 0, vjust = 1) +
    labs(x = "Cell Volume (μm³)",
         y = "Effective Radius (μm)",
         title = "Cell Volume vs Functional Spacing")
  
  # 4. Bar plot of spacing ratios
  p4 <- ggplot(combined_results, 
               aes(x = reorder(subtype, spacing_ratio), y = spacing_ratio)) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    coord_flip() +
    labs(x = "Cell Type",
         y = "Spacing Ratio (Effective Radius / Soma Exclusion)",
         title = "Relative Spacing Efficiency",
         subtitle = "Higher values indicate greater spacing relative to cell size")
  
  # 5. Free space analysis
  p5 <- ggplot(combined_results, 
               aes(x = reorder(subtype, free_space), y = free_space)) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    geom_col(fill = "darkgreen", alpha = 0.7) +
    coord_flip() +
    labs(x = "Cell Type",
         y = "Free Space (μm)",
         title = "Available Space Between Cells",
         subtitle = "Difference between effective radius and soma exclusion")
  
  # Combine plots using patchwork
  combined_plot <- (p1 + p2 + p3) / (p4 + p5) +
    plot_layout(heights = c(1.2, 1)) +
    plot_annotation(
      title = "Comprehensive Cell Size and Spacing Analysis",
      subtitle = sprintf("Based on %d bootstrap samples", n_bootstrap)
    ) &
    theme(plot.margin = margin(5, 5, 5, 5))
  
  # Return results and plots
  return(list(
    results = combined_results,
    spacing_plot = p1,
    volume_exclusion_plot = p2,
    volume_effective_plot = p3,
    spacing_ratio_plot = p4,
    free_space_plot = p5,
    combined_plot = combined_plot,
    models = list(
      spacing = lm_fit_spacing,
      vol_exclusion = lm_fit_vol_exclusion,
      vol_effective = lm_fit_vol_effective
    )
  ))
}


spacing_analysis <- analyze_spacing_relationships(
  rgcs_centered,
  rodieck_df,
  soma_exc_df,
  rgc_df,
  n_bootstrap = 1000
)

# Save results
write_csv(spacing_analysis$results,
         file.path(exclusion_dir, "spacing_analysis_results.csv"))

# Save combined plot
ggsave(file.path(exclusion_dir, "spacing_analysis_combined.png"),
      plot = spacing_analysis$combined_plot,
      width = 15, height = 13, dpi = 300)


################################################################################################
######################## COMPUTE HETEROTYPICAL PARAMS ##########################################
################################################################################################


library(future)
library(furrr)

# Create the new directory for plots - FIXED DIRECTORY
plot_dir <- file.path(exclusion_dir, "Heterotypical_plots")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

csv_filename1 <- file.path(plot_dir, "Heterotypical_rgc_pairs_autocorrelation.csv")
csv_filename2 <- file.path(plot_dir, paste0("Heterotypical_effective_radius_calculations_binWidth", delta_r, ".csv"))

if (!file.exists(csv_filename1)) {
  cat("Computing pair autocorrelations\n")
  
  # Set up parallel processing
  plan(multisession, workers = 18)
  
  # Generate all subtype pairs
  subtypes <- unique(rgc_df$Prediction)
  subtype_pairs <- expand.grid(subtype1 = subtypes, subtype2 = subtypes) %>%
    filter(subtype1 != subtype2) %>%
    as_tibble()
  
  # Function to process a single pair of subtypes
  process_pair <- function(pair, rgc_df_filtered, r) {
    subtype1 <- pair$subtype1
    subtype2 <- pair$subtype2
    pair_name <- paste0(subtype1, "-", subtype2)
    
    # Get centers from first subtype only
    rgc_centers <- rgc_df_filtered %>%
      filter(Prediction == subtype1 & !edge_point) %>%
      select(cell, x, y, Prediction, study_region)
    
    # Process each center
    results <- map_dfr(seq_len(nrow(rgc_centers)), function(i) {
      rgc <- rgc_centers[i,]
      
      # Get edge points of second subtype only
      edge_points <- rgc_df_filtered %>%
        filter(edge_point & Prediction == subtype2,
               study_region == rgc$study_region) %>%
        select(cell, x, y, Prediction, study_region) %>%
        mutate(
          x = x - rgc$x,
          y = y - rgc$y,
          center = pair_name,
          center_cell = rgc$cell,
          # Keep original Prediction to match existing Rodieck function
          Prediction = subtype2
        ) %>%
        filter(sqrt(x^2 + y^2) <= r)
      
      return(edge_points)
    })
    
    return(results)
  }
  
  # Process all pairs in parallel with progress indicator
  cat(sprintf("Processing %d subtype pairs in parallel...\n", nrow(subtype_pairs)))
  rgcs_centered_hetero <- future_map_dfr(
    split(subtype_pairs, 1:nrow(subtype_pairs)),
    ~process_pair(., rgc_df_filtered, r),
    .progress = TRUE
  )
  
  # Verify the results
  cat("\nVerification of results:\n")
  cat("Number of rows:", nrow(rgcs_centered_hetero), "\n")
  cat("Range of x coordinates:", range(rgcs_centered_hetero$x), "\n")
  cat("Range of y coordinates:", range(rgcs_centered_hetero$y), "\n")
  cat("Number of unique centers:", length(unique(rgcs_centered_hetero$center_cell)), "\n")
  
  # Save results
  write_csv(rgcs_centered_hetero, file = csv_filename1)
  cat("Completed and saved pair autocorrelations\n")
  
} else {
  rgcs_centered_hetero <- read_csv(csv_filename1)
  cat("Reloaded pre-computed pair autocorrelations\n")
}

if (!file.exists(csv_filename2)) {
  rodieck_df <- Rodieck1991_DensityRecovery(rgcs_centered_hetero, r, delta_r=delta_r)
  
  # Get unique subtype pairs from the data
  subtype_pairs <- unique(rgcs_centered_hetero$center)
  
  # Loop over subtypes and create plots
  for (pair in subtype_pairs) {
    p <- plot_exclusion_analysis(rgcs_centered_hetero, rodieck_df, pair, is_paired = TRUE)
    
    # Create filename using the pair identifier
    filename <- file.path(plot_dir, paste0(pair, "_effective_radius.png"))
    ggsave(
      filename = filename,
      plot = p,
      width = 12,
      height = 6,
      dpi = 300,
      units = "in"
    )
  }
  
  write_csv(rodieck_df, file = csv_filename2)
} else {
  rgcs_centered_hetero <- read_csv(csv_filename2)
  cat("Reloaded pre-computed statistics for the pairs\n")
}

##################################################################################
# normalized_density_hists <- read_csv(paste0(exclusion_dir, "/normalized_exclusion_histogram.csv"), col_names = F)
# colnames(normalized_density_hists) <- c("subtype", "dist", "d_i")
# normalized_density_hists$subtype <- factor(normalized_density_hists$subtype, 
#                                            levels = 1:length(subtypes),
#                                            labels = subtypes)
# normalized_density_hists <- normalized_density_hists %>%
#   group_by(subtype) %>%
#   mutate(annuli = row_number() - 1) %>%
#   ungroup()
# 
# Normalized_DensityRecovery <- function(normalized_df, A,
#                                        subtype_target=F) {
#   
#   delta_r <- pull(normalized_df[2,2] - normalized_df[1,2])
#   r <- max(normalized_density_hists$dist) + delta_r
#   
#   # Get all unique subtypes
#   subtypes <- unique(normalized_df$subtype)
#   
#   # Handle subtype_target logic
#   if (subtype_target == "all") {
#     target_subtypes <- subtypes
#   } else if (subtype_target == F) {
#     target_subtypes <- NULL  # Will be set per iteration
#   } else {
#     target_subtypes <- subtype_target
#   }
#   
#   # delta A_i = area of annulus i = pi*delta_r^2 (2i-1)
#   radii <- seq(0,r, delta_r) # set of radii for each annulus
#   
#   # Rodieck performed all his analysis on single subtypes, 
#   #thus to get valid recovery values I must do the same
#   # with subtype_target I will build in flexibility for heterologous analysis
#   rodieck_df <- data.frame(
#     subtype = character(),
#     N = numeric(),        # total number of points
#     D = numeric(),        # overall density
#     
#     annulus = numeric(),  # annulus distance
#     n_i = numeric(),      # number of points in annulus i
#     d_i = numeric(),      # distance for annulus i
#     lambda_i = numeric(),  # density for annulus i
#     delta_V_i = numeric(),  # "dip" from mean annulus i
#     V_e = numeric(), # average magnitude of the dead space
#     r_e = numeric() # effective radius of a cell
#   )
#   
#   for (subtype in subtypes) {
#     
#     # Print detailed progress
#     cat(sprintf("\nProcessing subtype %s...\n", subtype))
#     
#     # Determine which subtypes to compare against
#     comparison_subtypes <- if (is.null(target_subtypes)) {
#       c(subtype)  # Compare only against same subtype
#     } else {
#       target_subtypes  # Compare against specified subtypes or all subtypes
#     }
#     
#     distance_data  <- normalized_df %>%
#       filter(subtype==subtype)
#     
#     # Extract values from the distance calculator output
#     d_i <- distance_data$d_i
#     d_i <- distance_data$n_i
#     delta_A_i <- pi * (delta_r^2) * (2 * distance_data$annuli - 1)/1e6
#     
#     
#     # N = number of points in A
#     N <- mean(n_i/(d_i*delta_A_i), na.rm=T)
#     
#     # D = density of points in the region, equal to N/A;
#     D <- mean(d_i)
#     
#     # lambda_i = expected number of points in annulus i, assuming a Poisson distribution
#     lambda_i <- N*D*delta_A_i
#     
#     # Delta_V_i = the "volume" of the dip at bin position i scaled by the annular area
#     delta_V_i = 1e6*delta_A_i*(D-d_i)#(lambda_i - n_i)/N
#     
#     # V_e = the total volume of the dip, obtained by adding the volumes of the annuli until the density value for an annulus exceeds the mean density
#     # ie until the calculated volume for the annulus becomes negative
#     V_e <- 0
#     for (i in seq_along(distance_data$annuli)) {
#       if (n_i[i] < lambda_i[i]) {  # Continue until ni exceeds λi
#         V_e <- V_e + delta_V_i[i]
#       } else {
#         break
#       }
#     }
#     
#     #r_e = the effective radius. The total decrease in cells within the dip is 
#     #then made equivalent to one composed of a step change from zero to the mean density at the effective radius.
#     r_e = sqrt(V_e/(pi*D))
#     
#     cat(sprintf("  Completed analysis for subtype %s (r_e = %.2f)\n", subtype, r_e))
#     
#     temp_df <- data.frame(subtype = subtype, 
#                           N = N,
#                           D = D,
#                           annulus = distance_data$annuli, 
#                           n_i = n_i, 
#                           d_i = d_i,
#                           lambda_i = lambda_i,
#                           delta_V_i = delta_V_i,
#                           V_e = V_e,
#                           r_e = r_e)
#     
#     rodieck_df <- rbind(rodieck_df, temp_df) 
#   }
#   
#   rodieck_df <- rodieck_df %>%
#     mutate(delta_r = delta_r) 
#   
#   cat("\nRodieck analysis completed successfully!\n")
#   
#   return(rodieck_df)
# }
