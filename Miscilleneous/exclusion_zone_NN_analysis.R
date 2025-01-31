library(spatstat)
library(sf)
library(dplyr)
library(readr)
library(tidyr)
library(parallel)
library(doParallel)
library(ggplot2)

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

calculate_subtype_distances <- function(target_data, test_data, max_distance = 500) {
  # Calculate distance matrix between target and test cells
  dist_mat <- as.matrix(dist(rbind(target_data[, c("x", "y")], 
                                   test_data[, c("x", "y")])))
  
  # Subset the matrix to get only target-to-test distances
  n_target <- nrow(target_data)
  n_test <- nrow(test_data)
  target_test_distances <- dist_mat[1:n_target, (n_target + 1):(n_target + n_test)]
  
  # Convert to vector and filter distances
  distances <- as.vector(target_test_distances)
  distances <- distances[distances <= max_distance]
  
  return(distances)
}

# Main analysis function
analyze_cross_subtype_distances <- function(rgc_df, output_file) {
  # Check if output already exists
  if (file.exists(output_file)) {
    cat("Loading existing results from", output_file, "\n")
    return(read_csv(output_file))
  }
  
  # Get all lasso files
  lasso_files <- list.files(lasso_root, 
                            pattern = "\\d+_HQ_Lasso_coordinates\\.csv$", 
                            full.names = TRUE)
  
  all_distances <- list()
  
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
      # Process each region
      for(i in seq_along(lasso_list)) {
        region_name <- names(lasso_list)[i]
        cat(sprintf("Processing region %s...\n", region_name))
        
        # Convert region to window and filter points
        region_window <- sf_to_owin(lasso_list[[i]])
        region_data <- filter_points_in_window(slide_data, region_window)
        
        # Get all subtypes in this region
        subtypes <- unique(region_data$Prediction)
        
        # Loop through each target subtype
        for(target_subtype in subtypes) {
          cat(sprintf("Target subtype: %s\n", target_subtype))
          
          # Get target cells
          target_cells <- region_data
          
          # Only proceed if we have enough target cells
          if(nrow(target_cells) >= 3) {
            # Loop through test subtypes
            for(test_subtype in subtypes) {
              if(test_subtype != target_subtype) {
                # Get test cells
                test_cells <- region_data %>% 
                  filter(Prediction == test_subtype)
                
                # Only proceed if we have test cells
                if(nrow(test_cells) >= 1) {
                  # Calculate distances
                  distances <- calculate_subtype_distances(target_cells, test_cells)
                  
                  # Store results
                  distance_df <- data.frame(
                    distance = distances,
                    target_subtype = target_subtype,
                    test_subtype = test_subtype,
                    slide_region = paste(slide_id, region_name, sep="_")
                  )
                  
                  all_distances[[length(all_distances) + 1]] <- distance_df
                }
              }
            }
          }
        }
      }
    }
  }
  
  # Combine all results
  final_distances <- do.call(rbind, all_distances)
  
  # Save results
  write_csv(final_distances, output_file)
  
  return(final_distances)
}

# Set up directories and run analysis
# Set up directories and run analysis
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
  select(Prediction, slide, x, y)

# Run analysis
cross_subtype_distances <- analyze_cross_subtype_distances(rgc_df, output_file)

# Print summary
cat("\nAnalysis complete. Results saved to:", output_file, "\n")
cat("Total number of distance measurements:", nrow(cross_subtype_distances), "\n")
cat("Number of unique target subtypes:", length(unique(cross_subtype_distances$target_subtype)), "\n")
cat("Number of unique study regions:", length(unique(cross_subtype_distances$slide_region)), "\n")

##########################################################################################
##########################################################################################
##########################################################################################

# Updated RodieckDRP function to handle homotypic and all-cell comparisons
RodieckDRP <- function(distances_df, reference_subtype, target_subtype = NULL,
                       max_radius = 500, n_annuli = 100) {
  
  # Add debug print
  cat(sprintf("\nDRP Debug - Processing: %s vs %s\n", reference_subtype, target_subtype))
  cat(sprintf("Input distances rows: %d\n", nrow(distances_df)))
  
  # Handle special cases for target_subtype
  is_homotypic <- is.null(target_subtype)  # Only compare to same type
  is_all_cells <- !is.null(target_subtype) && target_subtype == "all"  # Compare to all cells
  
  # Filter distances for reference-target pair
  if (is_homotypic) {
    pair_distances <- distances_df %>%
      filter(target_subtype == reference_subtype,
             test_subtype == reference_subtype)
  } else if (is_all_cells) {
    pair_distances <- distances_df %>%
      filter(target_subtype == reference_subtype)
  } else {
    pair_distances <- distances_df %>%
      filter(target_subtype == reference_subtype,
             test_subtype == target_subtype)
  }
  
  # Debug print
  cat(sprintf("Filtered distances rows: %d\n", nrow(pair_distances)))
  
  # Early return if not enough data
  if(nrow(pair_distances) < 10) {  # Minimum threshold for meaningful analysis
    cat("Warning: Not enough data points for analysis\n")
    return(list(
      bins = NULL,
      n_i = NULL,
      d_i = NULL,
      lambda_i = NULL,
      delta_V_i = NULL,
      delta_A_i = NULL,
      D = 0,
      V_e = 0,
      effective_radius = 0,
      bin_centers = NULL
    ))
  }
  
  # Calculate bin width
  bin_width <- max_radius / n_annuli
  
  # Create bins
  bins <- seq(0, max_radius, by = bin_width)
  
  # Initialize results vectors
  n_i <- numeric(n_annuli)  # counts per annulus
  d_i <- numeric(n_annuli)  # density per annulus
  lambda_i <- numeric(n_annuli)  # expected counts
  delta_V_i <- numeric(n_annuli)  # volume of dip per annulus
  delta_A_i <- numeric(n_annuli)  # area of each annulus
  
  # Calculate areas of annuli using equation (1) from paper
  for(i in 1:n_annuli) {
    # Calculate area using eqn (1): ΔA_i = πΔr^2(2i - 1)
    delta_A_i[i] <- pi * bin_width^2 * (2*i - 1)
    
    # Count points in this annulus
    r_inner <- bins[i]
    r_outer <- bins[i + 1]
    n_i[i] <- sum(pair_distances$distance >= r_inner & 
                    pair_distances$distance < r_outer)
  }
  
  # Calculate total study area and mean density (D)
  total_area <- pi * max_radius^2
  N <- nrow(pair_distances)  # total number of points
  D <- N / total_area  # mean density
  
  # Calculate densities and expected counts for each annulus
  cat("\nDensity calculations:\n")
  for(i in 1:n_annuli) {
    # Calculate density for this annulus
    d_i[i] <- n_i[i] / delta_A_i[i]
    
    # Calculate expected count (lambda_i) eqn (2)
    lambda_i[i] <- N * D * delta_A_i[i]
    
    # Calculate volume of dip eqn (6)
    delta_V_i[i] <- (lambda_i[i] - n_i[i]) / N
    
    # Debug print for this bin
    cat(sprintf("Bin %d (%.1f-%.1f μm): count=%d, density=%.2f, mean=%.2f\n", 
                i, bins[i], bins[i+1], n_i[i], d_i[i], D))
  }
  
  # Calculate total volume of dip (V_e) until density exceeds mean
  # Following eqn (7) from paper
  # Modified effective radius calculation
  V_e <- 0
  effective_radius <- NA  # Start with NA
  below_mean <- FALSE   # Track if we've gone below mean
  
  # First pass to calculate V_e
  for(i in 1:n_annuli) {
    # First wait until density drops below mean
    if(!below_mean && d_i[i] <= D) {
      below_mean <- TRUE
      cat(sprintf("Density drops below mean at bin %d (%.1f μm): %.2f < %.2f\n", 
                  i, bins[i], d_i[i], D))
    }
    
    # Then look for recovery
    if(below_mean && d_i[i] >= D) {
      effective_radius <- bins[i]
      cat(sprintf("Density recovers at bin %d (%.1f μm): %.2f >= %.2f\n", 
                  i, effective_radius, d_i[i], D))
      break
    }
    
    # Accumulate volume only after dropping below mean
    if(below_mean) {
      V_e <- V_e + delta_V_i[i]
    }
  }
  
  if(is.na(effective_radius)) {
    effective_radius <- max_radius
    cat("Warning: Density never recovers to mean within max_radius\n")
    cat(sprintf("Min density = %.2f at %.1f μm\n", min(d_i), bins[which.min(d_i)]))
  }
  
  # Return results
  list(
    bins = bins,
    n_i = n_i,
    d_i = d_i,
    lambda_i = lambda_i,
    delta_V_i = delta_V_i,
    delta_A_i = delta_A_i,
    D = D,
    V_e = V_e,
    effective_radius = effective_radius,
    bin_centers = bins[-1] - bin_width/2
  )
}

# Function to plot DRP results
plot_DRP <- function(drp_results, reference_subtype, target_subtype, save_path = NULL) {
  
  # If save_path is provided, ensure it uses the correct base directory
  if(!is.null(save_path)) {
    rel_path <- basename(save_path)
    subtype_dir <- basename(dirname(save_path))
    save_path <- file.path(out_root, "DensityRecoveryProfiles", 
                           subtype_dir, rel_path)  # Remove extra exclusion_zones
  }
  
  # Check for NA values and handle them
  if(is.na(drp_results$effective_radius)) {
    warning("NA effective radius for ", reference_subtype, " vs ", target_subtype)
    return(NULL)
  }
  
  plot_data <- data.frame(
    distance = drp_results$bin_centers,
    density = drp_results$d_i,
    mean_density = rep(drp_results$D, length(drp_results$d_i))
  )
  
  # Remove any NA values
  plot_data <- plot_data %>% 
    filter(!is.na(distance), !is.na(density))
  
  # Only create plot if we have valid data
  if(nrow(plot_data) == 0) {
    warning("No valid data for plotting DRP of ", reference_subtype, " vs ", target_subtype)
    return(NULL)
  }
  
  p <- ggplot(plot_data) +
    geom_line(aes(x = distance, y = density)) +
    geom_hline(yintercept = drp_results$D, linetype = "dashed", color = "blue") +
    geom_vline(xintercept = drp_results$effective_radius, 
               linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(title = sprintf("Density Recovery Profile\n%s vs %s", 
                         reference_subtype, target_subtype),
         x = "Distance (μm)",
         y = "Density",
         subtitle = sprintf("Effective Radius: %.2f μm", 
                            drp_results$effective_radius))
  
  if(!is.null(save_path)) {
    tryCatch({
      ggsave(save_path, p, width = 8, height = 6)
    }, error = function(e) {
      warning("Failed to save plot for ", reference_subtype, " vs ", target_subtype, ": ", e$message)
    })
  }
  
  return(p)
}
# Main analysis function
analyze_exclusion_zones <- function(cross_subtype_distances, out_root) {
  # Create directories with explicit checks
  drp_dir <- file.path(out_root, "DensityRecoveryProfiles")
  intermediate_dir <- file.path(out_root, "intermediates")
  
  for(dir in c(out_root, drp_dir, intermediate_dir)) {
    if(!dir.exists(dir)) {
      cat(sprintf("Creating directory: %s\n", dir))
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  # Check if results already exist
  checkpoint_file <- file.path(intermediate_dir, "analysis_progress.rds")
  
  if(file.exists(checkpoint_file)) {
    cat("Loading existing results...\n")
    saved_results <- readRDS(checkpoint_file)
    return(saved_results)
  }
  
  # Get unique subtypes
  subtypes <- unique(cross_subtype_distances$target_subtype)
  
  # Initialize results storage
  soma_sizes <- data.frame()
  homotypic_radii <- data.frame()
  heterotypic_comparisons <- data.frame()
  
  # 1. Calculate soma sizes (vs all cells)
  cat("\nCalculating soma sizes (correlation with all cells)...\n")
  for(subtype in subtypes) {
    cat(sprintf("Processing subtype: %s\n", subtype))
    
    # Create subtype directory
    subtype_dir <- file.path(drp_dir, subtype)
    dir.create(subtype_dir, showWarnings = FALSE)
    
    # Calculate DRP vs all cells
    soma_drp <- RodieckDRP(cross_subtype_distances, subtype, "all")
    soma_size <- soma_drp$effective_radius
    
    # Save soma DRP plot
    plot_DRP(soma_drp, subtype, "all cells", 
             file.path(subtype_dir, "soma_size.png"))
    
    # Store soma size
    soma_sizes <- rbind(soma_sizes, data.frame(
      subtype = subtype,
      soma_size = soma_size
    ))
    
    # Save intermediate results
    saveRDS(list(
      progress = "soma_sizes",
      current_subtype = subtype,
      soma_sizes = soma_sizes
    ), checkpoint_file)
  }
  
  # Save soma sizes
  write_csv(soma_sizes, file.path(out_root, "soma_sizes.csv"))
  
  # 2. Calculate homotypic exclusion zones
  cat("\nCalculating homotypic exclusion zones...\n")
  for(subtype in subtypes) {
    # Calculate DRP for homotypic interactions
    homo_drp <- RodieckDRP(cross_subtype_distances, subtype, NULL)
    homo_radius <- homo_drp$effective_radius
    
    # Save homotypic plot
    plot_DRP(homo_drp, subtype, sprintf("%s (homotypic)", subtype),
             file.path(file.path(drp_dir, subtype), "homotypic.png"))
    
    # Store homotypic radius with relative measure
    soma_size <- soma_sizes$soma_size[soma_sizes$subtype == subtype]
    homotypic_radii <- rbind(homotypic_radii, data.frame(
      subtype = subtype,
      homotypic_radius = homo_radius,
      relative_to_soma = homo_radius / soma_size
    ))
    
    # Save intermediate results
    saveRDS(list(
      progress = "homotypic",
      current_subtype = subtype,
      soma_sizes = soma_sizes,
      homotypic_radii = homotypic_radii
    ), checkpoint_file)
  }
  
  # Save homotypic radii
  write_csv(homotypic_radii, file.path(out_root, "homotypic_radii.csv"))
  
  # 3. Calculate heterotypic exclusion zones
  cat("\nCalculating heterotypic exclusion zones...\n")
  total_pairs <- length(subtypes) * (length(subtypes) - 1)
  pair_counter <- 0
  
  for(target_subtype in subtypes) {
    other_subtypes <- setdiff(subtypes, target_subtype)
    
    for(test_subtype in other_subtypes) {
      pair_counter <- pair_counter + 1
      cat(sprintf("\rProcessing pair %d/%d (%.1f%%): %s vs %s", 
                  pair_counter, total_pairs,
                  100 * pair_counter/total_pairs,
                  target_subtype, test_subtype))
      
      # Calculate DRP for heterotypic interaction
      hetero_drp <- RodieckDRP(cross_subtype_distances, target_subtype, test_subtype)
      hetero_radius <- hetero_drp$effective_radius
      
      # Save heterotypic plot
      plot_DRP(hetero_drp, target_subtype, test_subtype,
               file.path(file.path(drp_dir, target_subtype), 
                         sprintf("vs_%s.png", test_subtype)))
      
      # Get relevant reference values
      soma_size <- soma_sizes$soma_size[soma_sizes$subtype == target_subtype]
      homo_radius <- homotypic_radii$homotypic_radius[
        homotypic_radii$subtype == target_subtype]
      
      # Store comparison
      heterotypic_comparisons <- rbind(heterotypic_comparisons, data.frame(
        target_subtype = target_subtype,
        test_subtype = test_subtype,
        heterotypic_radius = hetero_radius,
        relative_to_soma = hetero_radius / soma_size,
        relative_to_homo = hetero_radius / homo_radius
      ))
      
      # Save intermediate results periodically
      if(pair_counter %% 10 == 0) {
        saveRDS(list(
          progress = "heterotypic",
          current_pair = pair_counter,
          soma_sizes = soma_sizes,
          homotypic_radii = homotypic_radii,
          heterotypic_comparisons = heterotypic_comparisons
        ), checkpoint_file)
      }
    }
  }
  
  # 4. Identify split candidates
  cat("\nIdentifying split candidates...\n")
  split_candidates <- heterotypic_comparisons %>%
    filter(relative_to_homo >= 1) %>%
    arrange(target_subtype, desc(relative_to_homo))
  
  # Save all results
  write_csv(heterotypic_comparisons, file.path(out_root, "heterotypic_comparisons.csv"))
  write_csv(split_candidates, file.path(out_root, "split_candidates.csv"))
  
  # Create summary visualization
  p_homotypic <- ggplot(homotypic_radii, 
                        aes(x = reorder(subtype, relative_to_soma), 
                            y = relative_to_soma)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Homotypic Radii Relative to Soma Size",
         x = "Subtype",
         y = "Ratio (Homotypic Radius / Soma Size)")
  
  ggsave(file.path(out_root, "homotypic_ratios.png"), p_homotypic, 
         width = 12, height = 6)
  
  # Prepare final results
  results <- list(
    soma_sizes = soma_sizes,
    homotypic_radii = homotypic_radii,
    heterotypic_comparisons = heterotypic_comparisons,
    split_candidates = split_candidates
  )
  
  # Save final results
  saveRDS(results, checkpoint_file)
  
  return(results)
}

# Main script execution
out_root <- file.path(root, "exclusion_zones")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# Run analysis
exclusion_results <- analyze_exclusion_zones(cross_subtype_distances, out_root)

# Print summary of findings
cat("\nAnalysis Summary:\n")
cat(sprintf("Total subtypes analyzed: %d\n", 
            nrow(exclusion_results$soma_sizes)))
cat(sprintf("Number of split candidates identified: %d\n", 
            nrow(exclusion_results$split_candidates)))

if(nrow(exclusion_results$split_candidates) > 0) {
  cat("\nPotential split candidates:\n")
  print(exclusion_results$split_candidates %>%
          select(target_subtype, test_subtype, relative_to_homo) %>%
          arrange(desc(relative_to_homo)))
}

###############################################################################
###############################################################################
###############################################################################
# Statistical testing function for spatial relationships
# bootstrap_spatial_test <- function(reference_subtype, target_subtype, 
#                                    cross_subtype_distances, n_boot = 1000) {
#   # Get homotypic distances
#   homo_distances <- cross_subtype_distances %>%
#     filter(target_subtype == reference_subtype,
#            test_subtype == reference_subtype)
#   
#   # Get heterotypic distances
#   hetero_distances <- cross_subtype_distances %>%
#     filter(target_subtype == reference_subtype,
#            test_subtype == target_subtype)
#   
#   # Calculate observed ratio using full dataset
#   homo_drp <- RodieckDRP(homo_distances, reference_subtype, NULL)
#   hetero_drp <- RodieckDRP(hetero_distances, reference_subtype, target_subtype)
#   observed_ratio <- hetero_drp$effective_radius / homo_drp$effective_radius
#   
#   # Storage for bootstrap results
#   boot_ratios <- numeric(n_boot)
#   
#   # Perform bootstrap
#   for(i in 1:n_boot) {
#     # Resample distances with replacement
#     boot_homo <- homo_distances[sample(1:nrow(homo_distances), replace = TRUE), ]
#     boot_hetero <- hetero_distances[sample(1:nrow(hetero_distances), replace = TRUE), ]
#     
#     # Calculate DRPs for bootstrap sample
#     boot_homo_drp <- RodieckDRP(boot_homo, reference_subtype, NULL)
#     boot_hetero_drp <- RodieckDRP(boot_hetero, reference_subtype, target_subtype)
#     
#     # Store ratio
#     boot_ratios[i] <- boot_hetero_drp$effective_radius / boot_homo_drp$effective_radius
#   }
#   
#   # Calculate confidence intervals
#   ci <- quantile(boot_ratios, c(0.025, 0.975))
#   
#   # Calculate p-value (two-sided test)
#   # Test if ratio significantly differs from 1
#   p_value <- mean(abs(boot_ratios - 1) >= abs(observed_ratio - 1))
#   
#   return(list(
#     p_value = p_value,
#     boot_ratios = boot_ratios,
#     observed_ratio = observed_ratio,
#     ci_lower = ci[1],
#     ci_upper = ci[2]
#   ))
# }

bootstrap_spatial_test <- function(reference_subtype, target_subtype, 
                                   cross_subtype_distances, n_boot = 100) {
  # Get homotypic distances
  homo_distances <- cross_subtype_distances %>%
    filter(target_subtype == reference_subtype,
           test_subtype == reference_subtype)
  
  # Get heterotypic distances
  hetero_distances <- cross_subtype_distances %>%
    filter(target_subtype == reference_subtype,
           test_subtype == target_subtype)
  
  # Calculate observed ratio using full dataset
  homo_drp <- RodieckDRP(homo_distances, reference_subtype, NULL)
  hetero_drp <- RodieckDRP(hetero_distances, reference_subtype, target_subtype)
  observed_ratio <- hetero_drp$effective_radius / homo_drp$effective_radius
  
  # Setup parallel processing
  n_cores <-2
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Perform parallel bootstrap
  boot_ratios <- parLapply(cl, 1:n_boot, function(i) {
    # Resample distances with replacement
    boot_homo <- homo_distances[sample(1:nrow(homo_distances), replace = TRUE), ]
    boot_hetero <- hetero_distances[sample(1:nrow(hetero_distances), replace = TRUE), ]
    
    # Calculate DRPs for bootstrap sample
    boot_homo_drp <- RodieckDRP(boot_homo, reference_subtype, NULL)
    boot_hetero_drp <- RodieckDRP(boot_hetero, reference_subtype, target_subtype)
    
    # Return ratio
    boot_hetero_drp$effective_radius / boot_homo_drp$effective_radius
  })
  
  stopCluster(cl)
  
  # Convert list to vector
  boot_ratios <- unlist(boot_ratios)
  
  # Calculate confidence intervals
  ci <- quantile(boot_ratios, c(0.025, 0.975))
  
  # Calculate p-value (two-sided test)
  p_value <- mean(abs(boot_ratios - 1) >= abs(observed_ratio - 1))
  
  return(list(
    p_value = p_value,
    boot_ratios = boot_ratios,
    observed_ratio = observed_ratio,
    ci_lower = ci[1],
    ci_upper = ci[2]
  ))
}

adjust_pvalues <- function(ratio_df) {
  # Group by target subtype to handle comparisons within each group
  ratio_df %>%
    group_by(target_subtype) %>%
    mutate(
      # Add FDR correction (Benjamini-Hochberg)
      p_adj = p.adjust(p_value, method = "BH"),
      # Add Bonferroni correction for more stringent control
      p_adj_bonf = p.adjust(p_value, method = "bonferroni"),
      # Add significance levels
      significance = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01 ~ "**", 
        p_adj < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    ) %>%
    ungroup()
}

# Function to create all visualizations
create_exclusion_visualizations <- function(results, out_root) {
  # Create subdirectory for intermediate results
  intermediate_dir <- file.path(out_root, "intermediates")
  dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create subdirectory for DRP plots
  drp_dir <- file.path(out_root, "DensityRecoveryProfiles")
  dir.create(drp_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Get subtypes
  subtypes <- unique(c(results$heterotypic_comparisons$target_subtype,
                       results$heterotypic_comparisons$test_subtype))
  
  # Create ratio matrix with confidence intervals
  ratio_df <- data.frame()
  total_pairs <- length(subtypes) * (length(subtypes) - 1)
  pair_counter <- 0
  
  for(target in subtypes) {
    # Create checkpoint file path
    checkpoint_file <- file.path(intermediate_dir, sprintf("%s_ratios.rds", target))
    
    # Check if checkpoint exists
    if(file.exists(checkpoint_file)) {
      cat(sprintf("\nLoading existing results for target %s\n", target))
      target_ratios <- readRDS(checkpoint_file)
      ratio_df <- rbind(ratio_df, target_ratios)
      next
    }
    
    target_ratios <- data.frame()
    
    for(test in subtypes) {
      if(target == test) {
        # Homotypic case
        target_ratios <- rbind(target_ratios, data.frame(
          target_subtype = target,
          test_subtype = test,
          ratio = 1,
          ci_lower = 1,
          ci_upper = 1
        ))
      } else {
        pair_counter <- pair_counter + 1
        cat(sprintf("\rProcessing pair %d/%d (%.1f%%): %s vs %s", 
                    pair_counter, total_pairs,
                    100 * pair_counter/total_pairs,
                    target, test))
        
        # Heterotypic case
        boot_results <- bootstrap_spatial_test(target, test, cross_subtype_distances)
        target_ratios <- rbind(target_ratios, data.frame(
          target_subtype = target,
          test_subtype = test,
          ratio = boot_results$observed_ratio,
          ci_lower = boot_results$ci_lower,
          ci_upper = boot_results$ci_upper,
          p_value = boot_results$p_value
        ))
      }
    }
    
    # Save checkpoint
    saveRDS(target_ratios, checkpoint_file)
    ratio_df <- rbind(ratio_df, target_ratios)
  }

  
  # Create ratio heatmap
  p_heatmap <- ggplot(ratio_df, aes(x = test_subtype, y = target_subtype, fill = ratio)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 1, limits = c(0, max(ratio_df$ratio))) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8)) +
    labs(title = "Relative Exclusion Zone Ratios",
         subtitle = "Color shows ratio of heterotypic to homotypic effective radius",
         x = "Test Subtype",
         y = "Target Subtype",
         fill = "Ratio")
  
  # Apply multiple testing correction
  ratio_df <- ratio_df %>%
    adjust_pvalues()
  
  # Update significance heatmap to use adjusted p-values
  p_significance <- ggplot(ratio_df, 
                           aes(x = test_subtype, y = target_subtype, 
                               fill = -log10(p_adj))) +
    geom_tile() +
    geom_text(aes(label = significance), size = 2) +
    scale_fill_gradient(low = "white", high = "red") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8)) +
    labs(title = "Statistical Significance of Exclusion Zone Ratios",
         subtitle = "Stars indicate FDR-adjusted significance level",
         x = "Test Subtype",
         y = "Target Subtype",
         fill = "-log10(p_adj)")
  
  # Save plots
  ggsave(file.path(drp_dir, "size_comparison.png"), p_sizes, width = 10, height = 8)
  ggsave(file.path(drp_dir, "ratio_heatmap.png"), p_heatmap, width = 12, height = 10)
  ggsave(file.path(drp_dir, "significance_heatmap.png"), p_significance, width = 12, height = 10)
  
  # Save detailed results
  write_csv(ratio_df, file.path(drp_dir, "spatial_relationships.csv"))
  
  return(list(
    size_comparison = p_sizes,
    ratio_heatmap = p_heatmap,
    significance_heatmap = p_significance,
    ratio_data = ratio_df
  ))
}

# Main analysis wrapper function
# analyze_exclusion_zones <- function(cross_subtype_distances, out_root) {
#   # Previous analysis code remains the same...
#   
#   # Add visualization and statistical testing
#   viz_results <- create_exclusion_visualizations(results, out_root)
#   
#   # Merge all results
#   results$visualizations <- viz_results
#   results$ratio_statistics <- viz_results$ratio_data
#   
#   return(results)
# }

# Run the analysis
results <- analyze_exclusion_zones(cross_subtype_distances, exclusion_dir)
