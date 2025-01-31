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
  TIF_PATH <- "M_RGCtypes_norm_circle_012025.tif"        # Input TIF file containing cell type maps
  OUTPUT_DIR <- "GlobalStatistics_Densities"                 # Main output directory name
  MORANS_DIR <- "Morans"                                 # Directory for Moran's I results
  SCAN_DIR <- "Scan"                                     # Directory for scan statistics results
  VISUAL_SCENE_DIR <- "VisualScene"                      # Directory for visual scene analysis
  
  #=============================================================================
  # Scan Analysis Parameters
  #=============================================================================
  # Toggle off for density analysis
  GENE_MAP_ANALYSIS <- FALSE
  # Gene threshold
  GENE_THRESHOLD<-0.75
  # Grid and Buffer
  GRID_SIZE <- 71          # Number of cells per side in analysis grid
  BUFFER_SIZE <- 3         # Pixels to remove from edge to prevent boundary effects
  
  # Scan Statistics
  ALPHA <- 0.05           # Statistical significance threshold
  N_SIM <- 1000           # Number of Monte Carlo simulations for scan statistics
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

if (GENE_MAP_ANALYSIS) {
  #' Compile comprehensive statistics table from all analyses with gene/subtype information
  #' @param base_dir Base output directory containing all results
  #' @param valid_rgc_genes Dataframe containing valid genes before Moran's I filtering
  #' @param clustered_rgc_genes Dataframe containing genes that passed Moran's I filtering
  #' @return List containing full, raw embedding, and normalized embedding tables
  compile_statistics <- function(base_dir, valid_rgc_genes, clustered_rgc_genes) {
    # Construct full paths using ROOT_DIR and OUTPUT_DIR
    full_base_dir <- file.path(ROOT_DIR, base_dir)
    
    # Read individual result files with verbose error checking
    morans_path <- file.path(full_base_dir, "Morans", "morans_results.csv")
    f1_path <- file.path(full_base_dir, "VisualScene", "F1", "f1_scores.csv")
    enrich_path <- file.path(full_base_dir, "VisualScene", "Enrichment", "enrichment_statistics.csv")
    
    # Verbose error checking for file existence
    for(path in c(morans_path, f1_path, enrich_path)) {
      if(!file.exists(path)) {
        stop(sprintf("Cannot find required file: %s", path))
      }
    }
    
    # Prepare valid_rgc_genes and clustered_rgc_genes
    valid_rgc_genes <- valid_rgc_genes %>%
      rename(subtype = Prediction, gene = column) %>%
      mutate(original_index = row_number())
    
    clustered_rgc_genes <- clustered_rgc_genes %>%
      rename(subtype = Prediction, gene = column) %>%
      mutate(original_index = row_number())
    
    # Read and process Moran's results
    morans_results <- read.csv(morans_path) %>%
      mutate(original_index = Layer) %>%
      select(-Layer) %>%
      left_join(valid_rgc_genes, by = "original_index") %>%
      select(-original_index)
    
    # Read and process F1 results
    f1_results <- read.csv(f1_path) %>%
      mutate(original_index = Layer) %>%
      select(-Layer) %>%
      left_join(clustered_rgc_genes, by = "original_index") %>%
      select(-original_index)
    
    # Read and process enrichment statistics
    enrich_stats <- read.csv(enrich_path) %>%
      mutate(original_index = Layer) %>%
      select(-Layer) %>%
      left_join(clustered_rgc_genes, by = "original_index") %>%
      select(-original_index)
    
    # Handle the log transformation with protection against zeros
    enrich_wide <- enrich_stats %>%
      mutate(
        neg_log_p = ifelse(P_adj == 0, 
                           -log10(.Machine$double.eps),
                           -log10(P_adj)),
        neglog_col_name = paste0(Mask, "_negative10Friedmanpvalue"),
        padj_col_name = paste0(Mask, "_Padj"),
        p_col_name = paste0(Mask, "_P"),
        odds_col_name = paste0(Mask, "_OddsRatio")
      ) %>%
      select(subtype, gene, neglog_col_name, padj_col_name, odds_col_name, p_col_name,
             neg_log_p, P_adj, P_value, Odds_ratio) %>%
      pivot_longer(
        cols = c(neg_log_p, P_adj, P_value, Odds_ratio),
        names_to = "stat_type",
        values_to = "value"
      ) %>%
      mutate(
        col_name = case_when(
          stat_type == "neg_log_p" ~ neglog_col_name,
          stat_type == "P_adj" ~ padj_col_name,
          stat_type == "P_value" ~ p_col_name,
          stat_type == "Odds_ratio" ~ odds_col_name
        )
      ) %>%
      select(subtype, gene, col_name, value) %>%
      tidyr::pivot_wider(
        names_from = col_name,
        values_from = value
      )
    
    # Get largest cluster centroids from scan results
    scan_dir <- file.path(full_base_dir, "Scan", "ScanPortable")
    scan_centroids <- data.frame(original_index = integer(), 
                                 ScanX = numeric(), 
                                 ScanY = numeric())
    
    # Process each layer's scan results
    layer_files <- list.files(scan_dir, pattern = "layer_.*_pvalues.csv")
    for(file in layer_files) {
      layer_num <- as.integer(gsub("layer_|_pvalues.csv", "", file))
      p_values <- as.matrix(read.csv(file.path(scan_dir, file)))
      
      if(!all(p_values == -1)) {
        sig_points <- which(p_values != -1, arr.ind = TRUE)
        
        if(nrow(sig_points) > 0) {
          tryCatch({
            # Create adjacency matrix
            adj_matrix <- matrix(FALSE, nrow(p_values), ncol(p_values))
            adj_matrix[p_values != -1] <- TRUE
            
            # Find connected components
            clusters <- igraph::components(
              igraph::graph.adjacency(
                adj_matrix,
                mode = "undirected"
              )
            )
            
            if(length(clusters$csize) > 0) {
              largest_cluster <- which.max(clusters$csize)
              cluster_points <- which(clusters$membership == largest_cluster)
              
              # Make sure cluster_points are within bounds of sig_points
              valid_points <- cluster_points[cluster_points <= nrow(sig_points)]
              
              if(length(valid_points) > 0) {
                cluster_coords <- sig_points[valid_points, , drop = FALSE]
                
                centroid_x <- mean(cluster_coords[, 2])
                centroid_y <- mean(cluster_coords[, 1])
                
                scan_centroids <- rbind(scan_centroids,
                                        data.frame(original_index = layer_num,
                                                   ScanX = centroid_x,
                                                   ScanY = centroid_y))
              }
            }
          }, error = function(e) {
            warning(sprintf("Error processing scan results for layer %d: %s", 
                            layer_num, e$message))
          })
        }
      }
    }
    
    # Add gene/subtype information to scan_centroids
    scan_centroids <- scan_centroids %>%
      left_join(clustered_rgc_genes, by = "original_index") %>%
      select(-original_index)
    
    # Merge all statistics using gene and subtype as keys
    full_stats <- valid_rgc_genes %>%
      select(subtype, gene) %>%
      left_join(morans_results, by = c("subtype", "gene")) %>%
      left_join(scan_centroids, by = c("subtype", "gene")) %>%
      left_join(f1_results, by = c("subtype", "gene")) %>%
      left_join(enrich_wide, by = c("subtype", "gene"))
    
    # Create embedding_raw subset - only using rows that passed Moran's I test
    mask_cols <- grep("_F1$|_negative10Friedmanpvalue$", names(full_stats), value = TRUE)
    
    # Use only clustered genes (ones that passed Moran's I)
    embedding_raw <- full_stats %>%
      semi_join(clustered_rgc_genes, by = c("subtype", "gene")) %>%
      select(subtype, gene, MoransI, ScanX, ScanY,
             all_of(mask_cols)) %>%
      # Remove any remaining rows with NA values
      drop_na()
    
    # Create normalized version
    embedding_normalized <- embedding_raw %>%
      select(-subtype, -gene) %>%
      mutate(across(everything(), function(x) {
        if(max(x) == min(x)) return(x)  # No need for na.rm since we dropped NAs
        (x - min(x)) / (max(x) - min(x))
      }))
    
    # Create version with NAs replaced by -1000
    full_stats_igor <- full_stats %>%
      mutate(across(everything(), ~replace(., is.na(.), -1000)))
    
    # Save all versions
    write.csv(full_stats,
              file.path(full_base_dir, "comprehensive_statistics.csv"),
              row.names = FALSE)
    
    write.csv(full_stats_igor,
              file.path(full_base_dir, "comprehensive_statistics_igor.csv"),
              row.names = FALSE)
    
    write.csv(embedding_raw,
              file.path(full_base_dir, "embedding_raw.csv"),
              row.names = FALSE)
    
    write.csv(embedding_normalized,
              file.path(full_base_dir, "embedding_normalized.csv"),
              row.names = FALSE)
    
    message("Files saved to:")
    message(sprintf("  - %s", file.path(full_base_dir, "comprehensive_statistics.csv")))
    message(sprintf("  - %s", file.path(full_base_dir, "comprehensive_statistics_igor.csv")))
    message(sprintf("  - %s", file.path(full_base_dir, "embedding_raw.csv")))
    message(sprintf("  - %s", file.path(full_base_dir, "embedding_normalized.csv")))
    
    return(list(
      full_stats = full_stats,
      embedding_raw = embedding_raw,
      embedding_normalized = embedding_normalized
    ))
  }
} else {
  #' Compile comprehensive statistics table from all analyses
  #' @param base_dir Base output directory containing all results
  #' @return List containing full, raw embedding, and normalized embedding tables
  compile_statistics <- function(base_dir) {
    # Construct full paths using ROOT_DIR and OUTPUT_DIR
    full_base_dir <- file.path(ROOT_DIR, base_dir)
    
    # Read individual result files with verbose error checking
    morans_path <- file.path(full_base_dir, "Morans", "morans_results.csv")
    f1_path <- file.path(full_base_dir, "VisualScene", "F1", "f1_scores.csv")
    enrich_path <- file.path(full_base_dir, "VisualScene", "Enrichment", "enrichment_statistics.csv")
    
    # Verbose error checking for file existence
    for(path in c(morans_path, f1_path, enrich_path)) {
      if(!file.exists(path)) {
        stop(sprintf("Cannot find required file: %s", path))
      }
    }
    
    # Read files with error handling
    morans_results <- read.csv(morans_path)
    f1_results <- read.csv(f1_path)
    enrich_stats <- read.csv(enrich_path)
    
    # Handle the log transformation with protection against zeros
    enrich_wide <- enrich_stats %>%
      mutate(
        # Add small epsilon to prevent log(0), and handle Inf values
        neg_log_p = ifelse(P_adj == 0, 
                           -log10(.Machine$double.eps), # Use smallest possible number instead of 0
                           -log10(P_adj)),
        neglog_col_name = paste0(Mask, "_negative10Friedmanpvalue"),
        padj_col_name = paste0(Mask, "_Padj"),
        p_col_name = paste0(Mask, "_P"),
        odds_col_name = paste0(Mask, "_OddsRatio")
      ) %>%
      select(Layer, neglog_col_name, padj_col_name, odds_col_name, p_col_name,
             neg_log_p, P_adj, P_value, Odds_ratio) %>%
      pivot_longer(
        cols = c(neg_log_p, P_adj, P_value, Odds_ratio),
        names_to = "stat_type",
        values_to = "value"
      ) %>%
      mutate(
        col_name = case_when(
          stat_type == "neg_log_p" ~ neglog_col_name,
          stat_type == "P_adj" ~ padj_col_name,
          stat_type == "P_value" ~ p_col_name,
          stat_type == "Odds_ratio" ~ odds_col_name
        )
      ) %>%
      select(Layer, col_name, value) %>%
      tidyr::pivot_wider(
        names_from = col_name,
        values_from = value
      )
    
    # Get largest cluster centroids from scan results
    scan_dir <- file.path(full_base_dir, "Scan", "ScanPortable")
    scan_centroids <- data.frame(Layer = integer(), ScanX = numeric(), ScanY = numeric())
    
    # Process each layer's scan results
    layer_files <- list.files(scan_dir, pattern = "layer_.*_pvalues.csv")
    for(file in layer_files) {
      layer_num <- as.integer(gsub("layer_|_pvalues.csv", "", file))
      p_values <- as.matrix(read.csv(file.path(scan_dir, file)))
      
      # Find largest significant cluster
      if(all(p_values == -1)) next
      
      # Create matrix of significant points
      sig_points <- which(p_values != -1, arr.ind = TRUE)
      
      if(nrow(sig_points) > 0) {
        # Find connected components
        clusters <- igraph::components(
          igraph::graph.adjacency(
            p_values != -1,
            mode = "undirected"
          )
        )
        
        # Get largest cluster
        largest_cluster <- which.max(clusters$csize)
        cluster_points <- which(clusters$membership == largest_cluster)
        
        # Get coordinates of points in largest cluster
        cluster_coords <- sig_points[cluster_points, , drop = FALSE]
        
        # Calculate centroid (note: columns are now correct)
        centroid_x <- mean(cluster_coords[, 2])  # column coordinates
        centroid_y <- mean(cluster_coords[, 1])  # row coordinates
        
        scan_centroids <- rbind(scan_centroids,
                                data.frame(Layer = layer_num,
                                           ScanX = centroid_x,
                                           ScanY = centroid_y))
      }
    }
    
    
    # Merge all statistics
    full_stats <- morans_results %>%
      select(Layer, MoransI) %>%
      left_join(scan_centroids, by = "Layer") %>%
      left_join(f1_results, by = "Layer") %>%
      left_join(enrich_wide, by = "Layer")
    
    # Create embedding_raw subset
    mask_cols <- grep("_F1$|_negative10Friedmanpvalue$", names(full_stats), value = TRUE)
    
    # Use all GMM columns for each Gaussian
    embedding_raw <- full_stats %>%
      select(Layer, MoransI, ScanX, ScanY,
             all_of(mask_cols))
    
    # Create normalized version
    embedding_normalized <- embedding_raw %>%
      select(-Layer) %>%
      mutate(across(everything(), function(x) {
        if(all(is.na(x))) return(x)
        if(max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) return(x)
        (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
      }))
    
    # Save all versions
    write.csv(full_stats,
              file.path(full_base_dir, "comprehensive_statistics.csv"),
              row.names = FALSE)
    
    write.csv(embedding_raw,
              file.path(full_base_dir, "embedding_raw.csv"),
              row.names = FALSE)
    
    write.csv(embedding_normalized,
              file.path(full_base_dir, "embedding_normalized.csv"),
              row.names = FALSE)
    
    message("Files saved to:")
    message(sprintf("  - %s", file.path(full_base_dir, "comprehensive_statistics.csv")))
    message(sprintf("  - %s", file.path(full_base_dir, "embedding_raw.csv")))
    message(sprintf("  - %s", file.path(full_base_dir, "embedding_normalized.csv")))
    
    return(list(
      full_stats = full_stats,
      embedding_raw = embedding_raw,
      embedding_normalized = embedding_normalized
    ))
  }
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
  
  significant_indices <- integer()  # Track indices of significant layers
  
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
    
    # Track significant layers (p < 0.05)
    if (!is.na(moran_result$p_value) && moran_result$p_value < 0.05) {
      significant_indices <- c(significant_indices, i)
    }
  }
  
  # Save results
  write.csv(results, file.path(output_dir, "morans_results.csv"), row.names = FALSE)
  
  # Also save the list of significant indices
  write.csv(data.frame(SignificantLayers = significant_indices),
            file.path(output_dir, "significant_layers.csv"),
            row.names = FALSE)
  
  return(list(
    results = results,
    significant_indices = significant_indices
  ))
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
    p_value_matrix <- matrix(-1, nrow=nrow(grid), ncol=ncol(grid))
    
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
    "  * -1: No significant cluster or outside study region",
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
  
  ac_mask <- create_retinal_mask(grid, n_annuli = 8, n_sectors = 24,
                                 temporal_left = F, slope = 15,
                                 sector_fill_0 = 300, sector_fill_1 = 130,
                                 outer_fill = 4, inner_fill = 0)
  
  # Not tilted much
  # binocular_mask <- create_retinal_mask(grid, n_annuli = 7, n_sectors = 24,
  #                                       temporal_left = T, slope = 3,
  #                                       sector_fill_0 = 140, sector_fill_1 = 315,
  #                                       outer_fill = 4, inner_fill = 0)
  
  #44% holmgren tilt
  binocular_mask <- create_retinal_mask(grid, n_annuli = 8, n_sectors = 24,
                                        temporal_left = T, slope = 4,
                                        sector_fill_0 = 170, sector_fill_1 = 345,
                                        outer_fill = 7, inner_fill = 0)
  #40% holmgren tilt
  # binocular_mask <- create_retinal_mask(grid, n_annuli = 8, n_sectors = 24,
  #                                       temporal_left = T, slope = 6,
  #                                       sector_fill_0 = 170, sector_fill_1 = 345,
  #                                       outer_fill = 6, inner_fill = 0)
  
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
    ac = ac_mask,
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

#' Calculate confusion matrix and F1 score for a cluster mask against a visual scene mask
#' @param cluster_mask Binary matrix of significant clusters
#' @param scene_mask Binary matrix of visual scene region
#' @param study_region Matrix defining valid study region
#' @return List containing confusion matrix and F1 metrics
compute_f1_metrics <- function(cluster_mask, scene_mask, study_region) {
  # Only consider pixels within study region
  valid_pixels <- which(!is.na(study_region))
  
  # Extract valid pixels from both masks
  cluster_pixels <- cluster_mask[valid_pixels]
  scene_pixels <- scene_mask[valid_pixels]
  
  # Compute confusion matrix elements
  true_positive <- sum(cluster_pixels == 1 & scene_pixels == 1)
  true_negative <- sum(cluster_pixels == 0 & scene_pixels == 0)
  false_positive <- sum(cluster_pixels == 1 & scene_pixels == 0)
  false_negative <- sum(cluster_pixels == 0 & scene_pixels == 1)
  
  # Calculate precision, recall, and F1
  precision <- if(true_positive + false_positive > 0) {
    true_positive / (true_positive + false_positive)
  } else {
    0
  }
  
  recall <- if(true_positive + false_negative > 0) {
    true_positive / (true_positive + false_negative)
  } else {
    0
  }
  
  f1_score <- if(precision + recall > 0) {
    2 * (precision * recall) / (precision + recall)
  } else {
    0
  }
  
  return(list(
    confusion_matrix = matrix(
      c(true_positive, false_positive,
        false_negative, true_negative),
      nrow = 2,
      dimnames = list(
        Predicted = c("Positive", "Negative"),
        Actual = c("Positive", "Negative")
      )
    ),
    metrics = list(
      precision = precision,
      recall = recall,
      f1_score = f1_score
    )
  ))
}

#' Run F1 analysis for all maps and masks
#' @param scan_results List of scan test results
#' @param grid Reference grid matrix
#' @param mask_list List of visual scene masks
#' @param output_dir Base output directory
#' @return List containing results and plots
analyze_f1_scores <- function(scan_results, grid, mask_list, output_dir) {
  # Create directory structure
  f1_dir <- file.path(output_dir, "F1")
  dir.create(f1_dir, showWarnings = FALSE)
  confusion_dir <- file.path(f1_dir, "ConfusionMatrix")
  dir.create(confusion_dir, showWarnings = FALSE)
  
  # Create subdirectories for each mask
  mask_dirs <- lapply(names(mask_list), function(mask_name) {
    mask_dir <- file.path(confusion_dir, mask_name)
    dir.create(mask_dir, showWarnings = FALSE)
    return(mask_dir)
  })
  names(mask_dirs) <- names(mask_list)
  
  # Initialize results dataframe
  results_cols <- c("Layer")
  for(mask_name in names(mask_list)) {
    results_cols <- c(results_cols,
                      paste0(mask_name, "_Precision"),
                      paste0(mask_name, "_Recall"),
                      paste0(mask_name, "_F1"))
  }
  
  f1_results <- as.data.frame(matrix(0, 0, length(results_cols)))
  names(f1_results) <- results_cols
  
  # Process each layer
  for(i in seq_along(scan_results)) {
    current_result <- scan_results[[i]]
    current_row <- list(Layer = i)
    
    if(!is.null(current_result$clusters)) {
      # Create cluster mask
      cluster_mask <- matrix(0, nrow(grid), ncol(grid))
      for(cluster in current_result$clusters) {
        if(cluster$pvalue <= 0.05) {
          cluster_cells <- cluster$locids
          cells_in_cluster <- which(grid %in% cluster_cells, arr.ind = TRUE)
          cluster_mask[cells_in_cluster] <- 1
        }
      }
      
      # Process each mask
      for(mask_name in names(mask_list)) {
        mask <- mask_list[[mask_name]]
        
        # Compute F1 metrics
        f1_metrics <- compute_f1_metrics(cluster_mask, mask, grid)
        
        # Save confusion matrix
        write.csv(f1_metrics$confusion_matrix,
                  file = file.path(mask_dirs[[mask_name]],
                                   sprintf("layer_%d_confusion.csv", i)))
        
        # Add metrics to results
        current_row[[paste0(mask_name, "_Precision")]] <- f1_metrics$metrics$precision
        current_row[[paste0(mask_name, "_Recall")]] <- f1_metrics$metrics$recall
        current_row[[paste0(mask_name, "_F1")]] <- f1_metrics$metrics$f1_score
      }
    }
    
    f1_results <- rbind(f1_results, as.data.frame(current_row))
  }
  
  # Save complete results
  write.csv(f1_results, file = file.path(f1_dir, "f1_scores.csv"), row.names = FALSE)
  
  # Create F1 score plots for each mask
  plots <- list()
  for(mask_name in names(mask_list)) {
    f1_col <- paste0(mask_name, "_F1")
    
    plot_data <- data.frame(
      Layer = f1_results$Layer,
      F1_Score = f1_results[[f1_col]]
    ) %>%
      arrange(desc(F1_Score))
    
    p <- ggplot(plot_data, aes(x = reorder(factor(Layer), F1_Score),
                               y = F1_Score)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_blank()
      ) +
      labs(title = paste("F1 Scores -", gsub("_", " ", mask_name)),
           x = "Layer",
           y = "F1 Score")
    
    plots[[mask_name]] <- p
    
    # Save individual plot
    ggsave(
      file.path(f1_dir, paste0(mask_name, "_f1_scores.png")),
      plot = p,
      width = 8,
      height = 6,
      dpi = 300
    )
  }
  
  # Create combined plot
  if(length(plots) > 0) {
    combined_plot <- wrap_plots(plots, ncol = 2) +
      plot_layout(guides = "collect") &
      theme(plot.margin = margin(20, 20, 20, 20))
    
    ggsave(
      file.path(f1_dir, "all_f1_scores.png"),
      plot = combined_plot,
      width = 15,
      height = ceiling(length(plots)/2) * 6,
      limitsize = FALSE
    )
  }
  
  return(list(
    results = f1_results,
    plots = plots
  ))
}

#' Analyze visual scene enrichment using Fisher's exact test
#' @param f1_results Results from analyze_f1_scores
#' @param output_dir Base output directory
#' @return List containing statistical results and plots
analyze_visual_scene_enrichment <- function(f1_results, output_dir, mask_list) {
  require(tidyverse)
  require(patchwork)
  require(ggrepel)  # Add this package for repelling labels
  # Create output directory
  enrich_dir <- file.path(output_dir, "Enrichment")
  dir.create(enrich_dir, showWarnings = FALSE)
  
  # Get mask names correctly from the input data
  mask_names <- names(mask_list)  # Use all masks from the provided mask_list
  
  if(length(mask_names) == 0) {
    stop("No masks found in mask_list")
  }
  
  # Initialize results storage
  enrichment_results <- data.frame()
  
  # Process each layer and mask
  for(mask_name in mask_names) {
    # Get confusion matrix files for this mask
    conf_dir <- file.path(output_dir, "F1", "ConfusionMatrix", mask_name)
    conf_files <- list.files(conf_dir, pattern = "layer_.*_confusion\\.csv$")
    
    for(file in conf_files) {
      # Extract layer number
      layer <- as.numeric(gsub("layer_|_confusion\\.csv", "", file))
      
      # Read confusion matrix
      conf_matrix <- as.matrix(read.csv(file.path(conf_dir, file), row.names = 1))
      
      # Perform Fisher's exact test
      fisher_result <- fisher.test(conf_matrix, alternative = "greater")
      
      # Store results
      enrichment_results <- rbind(enrichment_results, data.frame(
        Layer = layer,
        Mask = mask_name,
        P_value = fisher_result$p.value,
        Odds_ratio = as.numeric(fisher_result$estimate),
        F1_score = f1_results$results[[paste0(mask_name, "_F1")]][layer]
      ))
    }
  }
  
  # Adjust p-values within each mask group
  enrichment_results <- enrichment_results %>%
    group_by(Mask) %>%
    mutate(P_adj = p.adjust(P_value, method = "BH")) %>%
    ungroup()
  
  # Save results
  write.csv(enrichment_results, 
            file = file.path(enrich_dir, "enrichment_statistics.csv"), 
            row.names = FALSE)
  
  # Create volcano plots
  plots <- list()
  for(mask_name in unique(enrichment_results$Mask)) {
    mask_data <- enrichment_results %>%
      filter(Mask == mask_name)
    
    p <- ggplot(mask_data, aes(x = F1_score, y = -log10(P_adj))) +
      geom_point(aes(color = P_adj < 0.05, size = log2(Odds_ratio))) +
      scale_color_manual(values = c("grey50", "#DC3545"),
                         labels = c("Not Significant", "Significant")) +
      scale_size_continuous(name = "log2(Odds Ratio)") +
      theme_bw() +
      theme(
        legend.position = "right",
        plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 10)
      ) +
      labs(
        title = paste("Enrichment Analysis -", gsub("_", " ", mask_name)),
        x = "F1 Score",
        y = "-log10(adjusted p-value)",
        color = "Significance"
      ) +
      # Add labels for significant points without repel
      geom_text(
        data = subset(mask_data, P_adj < 0.05),
        aes(label = Layer),
        size = 3,
        nudge_y = 0.1
      )
    
    plots[[mask_name]] <- p
    
    # Save individual plot
    ggsave(
      file.path(enrich_dir, paste0(mask_name, "_volcano.png")),
      plot = p,
      width = 8,
      height = 6,
      dpi = 300
    )
  }
  
  # Create combined plot
  if(length(plots) > 0) {
    combined_plot <- wrap_plots(plots, ncol = 2) +
      plot_layout(guides = "collect") &
      theme(plot.margin = margin(20, 20, 20, 20))
    
    ggsave(
      file.path(enrich_dir, "all_volcano_plots.png"),
      plot = combined_plot,
      width = 15,
      height = ceiling(length(plots)/2) * 6,
      limitsize = FALSE
    )
  }
  
  return(list(
    results = enrichment_results,
    plots = plots
  ))
}

#=============================================================================
# Clustering
#=============================================================================

#' Create UMAP visualization from normalized statistics
#' @param emdedding_data Df containing the normalized embeddings loaded from CSV
#' @param base_dir Base output directory containing embedding files
#' @param perplexity Perplexity parameter for UMAP (default 5)
#' @param seed Random seed for reproducibility
#' @return UMAP plot object
create_umap_visualization <- function(embedding_data, base_dir, perplexity=5, seed=18) {
  # Set random seed for reproducibility
  set.seed(seed)
  
  # Construct full path
  full_base_dir <- file.path(ROOT_DIR, base_dir)
  
  # Remove any columns with all NA values
  embedding_data <- embedding_data[, colSums(is.na(embedding_data)) < nrow(embedding_data)]
  
  # Handle any remaining NA values by imputing with column means
  for(col in names(embedding_data)) {
    if(any(is.na(embedding_data[[col]]))) {
      embedding_data[[col]][is.na(embedding_data[[col]])] <- mean(embedding_data[[col]], na.rm=TRUE)
    }
  }
  
  # Create labels for points
  point_labels <- paste0("T", seq_len(nrow(embedding_data)))
  
  # Compute UMAP
  require(umap)
  umap_config <- umap.defaults
  umap_config$n_neighbors <- min(perplexity, nrow(embedding_data) - 1)
  umap_config$random_state <- seed
  
  umap_result <- umap(embedding_data, config=umap_config)
  
  # Create plot data
  plot_data <- data.frame(
    UMAP1 = umap_result$layout[,1],
    UMAP2 = umap_result$layout[,2],
    Label = point_labels
  )
  
  # Create plot with ggplot2
  p <- ggplot(plot_data, aes(x=UMAP1, y=UMAP2, label=Label)) +
    geom_point(size=3, alpha=0.7) +
    geom_text_repel(
      size=3,
      box.padding = 0.5,
      point.padding = 0.2,
      force = 2
    ) +
    theme_minimal() +
    labs(
      title = "UMAP Visualization of RGC Subtypes",
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size=14),
      axis.title = element_text(size=12),
      axis.text = element_text(size=10)
    )
  
  # Save plot
  ggsave(
    filename = file.path(full_base_dir, "umap_visualization.png"),
    plot = p,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  # Save UMAP coordinates
  write.csv(
    plot_data,
    file = file.path(full_base_dir, "umap_coordinates.csv"),
    row.names = FALSE
  )
  
  message("Files saved:")
  message(sprintf("  - Plot: %s", file.path(full_base_dir, "umap_visualization.png")))
  message(sprintf("  - Coordinates: %s", file.path(full_base_dir, "umap_coordinates.csv")))
  
  return(list(
    plot = p,
    coordinates = plot_data,
    umap_object = umap_result
  ))
}

#' Perform hierarchical clustering on UMAP coordinates and update statistics
#' @param base_dir Base output directory
#' @param max_k Maximum number of clusters to try
#' @param seed Random seed for reproducibility
#' @return List containing cluster results and updated statistics
perform_hierarchical_clustering <- function(base_dir, max_k=15, seed=18) {
  require(cluster)
  require(dendextend)
  
  set.seed(seed)
  
  # Construct full path
  full_base_dir <- file.path(ROOT_DIR, base_dir)
  
  # Read UMAP coordinates and full statistics
  umap_coords <- read.csv(file.path(full_base_dir, "umap_coordinates.csv"))
  full_stats <- read.csv(file.path(full_base_dir, "comprehensive_statistics.csv")) %>%
    drop_na()
  
  # Perform hierarchical clustering
  dist_matrix <- dist(umap_coords[, c("UMAP1", "UMAP2")])
  hc <- hclust(dist_matrix, method="ward.D2")
  
  # Try different numbers of clusters and calculate silhouette scores
  sil_scores <- numeric(max_k - 1)
  for(k in 2:max_k) {
    clusters <- cutree(hc, k=k)
    sil <- silhouette(clusters, dist_matrix)
    sil_scores[k-1] <- mean(sil[,3])
  }
  
  # Find optimal k
  optimal_k <- which.max(sil_scores) + 1
  message(sprintf("Optimal number of clusters: %d", optimal_k))
  
  # Get cluster assignments using optimal k
  clusters <- cutree(hc, k=optimal_k)
  
  # Create dendrogram visualization
  dend <- as.dendrogram(hc)
  dend <- color_branches(dend, k=optimal_k)
  
  # Plot silhouette analysis
  sil_plot <- ggplot(data.frame(k = 2:max_k, score = sil_scores), aes(x=k, y=score)) +
    geom_line() +
    geom_point() +
    geom_point(data=data.frame(k=optimal_k, score=sil_scores[optimal_k-1]),
               color="red", size=3) +
    theme_minimal() +
    labs(title="Silhouette Score Analysis",
         x="Number of Clusters (k)",
         y="Average Silhouette Score")
  
  # Plot dendrogram
  png(file.path(full_base_dir, "dendrogram.png"),
      width=1200, height=800, res=150)
  plot(dend, main="Hierarchical Clustering Dendrogram",
       ylab="Height", leaflab="none")
  dev.off()
  
  # Save silhouette plot
  ggsave(file.path(full_base_dir, "silhouette_analysis.png"),
         plot=sil_plot, width=8, height=6, dpi=300)
  
  # Create UMAP plot with cluster colors
  umap_cluster_plot <- ggplot(cbind(umap_coords, Cluster=as.factor(clusters)),
                              aes(x=UMAP1, y=UMAP2, color=Cluster, label=Label)) +
    geom_point(size=3, alpha=0.7) +
    geom_text_repel(size=3, box.padding=0.5, show.legend=FALSE) +
    theme_minimal() +
    labs(title="UMAP Visualization with Cluster Assignments") +
    theme(plot.title = element_text(hjust=0.5))
  
  # Save UMAP cluster plot
  ggsave(file.path(full_base_dir, "umap_clusters.png"),
         plot=umap_cluster_plot, width=10, height=8, dpi=300)
  
  # Add cluster assignments to full statistics
  full_stats$Cluster <- clusters
  
  # Save updated statistics
  write.csv(full_stats,
            file.path(full_base_dir, "comprehensive_statistics.csv"),
            row.names=FALSE)
  
  # Save cluster assignments separately
  cluster_assignments <- data.frame(
    Label = umap_coords$Label,
    Cluster = clusters
  )
  write.csv(cluster_assignments,
            file.path(full_base_dir, "cluster_assignments.csv"),
            row.names=FALSE)
  
  message("Files saved:")
  message(sprintf("  - Dendrogram: %s", file.path(full_base_dir, "dendrogram.png")))
  message(sprintf("  - Silhouette analysis: %s", file.path(full_base_dir, "silhouette_analysis.png")))
  message(sprintf("  - UMAP clusters: %s", file.path(full_base_dir, "umap_clusters.png")))
  message(sprintf("  - Updated statistics: %s", file.path(full_base_dir, "comprehensive_statistics.csv")))
  message(sprintf("  - Cluster assignments: %s", file.path(full_base_dir, "cluster_assignments.csv")))
  
  return(list(
    optimal_k = optimal_k,
    silhouette_scores = sil_scores,
    cluster_assignments = clusters,
    dendrogram = dend,
    plots = list(
      silhouette = sil_plot,
      umap = umap_cluster_plot
    )
  ))
}

#' Create visualization of clustered maps with improved organization
#' @param tif_stack List of matrices containing the map data
#' @param grid Reference grid defining valid regions
#' @param cluster_assignments Data frame containing Label and Cluster columns
#' @param output_dir Directory to save the plot
#' @param plot_width Width of the plot in inches (default 15)
#' @param map_size Size of individual maps in inches (default 2)
#' @return ggplot object
visualize_clustered_maps <- function(tif_stack, grid, cluster_assignments, output_dir,
                                     plot_width = 15, map_size = 2) {
  require(tidyverse)
  require(gridExtra)
  require(cowplot)
  require(grid)
  
  # Extract layer numbers from labels and ensure proper ordering
  cluster_assignments$Layer <- as.numeric(gsub("T", "", cluster_assignments$Label))
  
  # Sort by cluster and then by layer
  cluster_assignments <- cluster_assignments %>%
    arrange(Cluster, Layer)
  
  # Function to create masked and normalized heatmap
  create_heatmap <- function(matrix, grid, title) {
    # Create data frame for plotting
    plot_data <- expand.grid(
      row = 1:nrow(matrix),
      col = 1:ncol(matrix)
    ) %>%
      mutate(
        value = as.vector(matrix)
      )
    
    # Mask values outside valid region
    plot_data$value[is.na(as.vector(grid))] <- NA
    
    # Normalize values within valid region
    valid_values <- plot_data$value[!is.na(plot_data$value) & plot_data$value > 0]
    if(length(valid_values) > 0) {
      min_val <- min(valid_values)
      max_val <- max(valid_values)
      plot_data$value[!is.na(plot_data$value) & plot_data$value > 0] <- 
        (valid_values - min_val) / (max_val - min_val)
    }
    
    # Create plot
    p <- ggplot(plot_data, aes(x = col, y = row, fill = value)) +
      geom_tile() +
      scale_fill_viridis_c(na.value = "white") +
      coord_equal() +
      theme_void() +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(hjust = 0.5, size = 8),
        legend.position = "none"
      ) +
      labs(title = title)
    
    return(p)
  }
  
  # Create a list to hold all plots by cluster
  plots_by_cluster <- list()
  
  # Process each cluster separately
  for(current_cluster in sort(unique(cluster_assignments$Cluster))) {
    # Get maps for this cluster
    cluster_maps <- cluster_assignments %>%
      filter(Cluster == current_cluster)
    
    # Create plots for this cluster
    cluster_plots <- list()
    
    for(i in 1:nrow(cluster_maps)) {
      layer <- cluster_maps$Layer[i]
      map <- tif_stack[[layer]]
      
      p <- create_heatmap(
        matrix = map,
        grid = grid,
        title = sprintf("T%d", layer)
      )
      
      cluster_plots[[i]] <- p
    }
    
    # Create a cluster label
    cluster_label <- textGrob(
      sprintf("Cluster %d", current_cluster),
      x = unit(0, "npc"),
      y = unit(0.5, "npc"),
      just = "left",
      gp = gpar(fontface = "bold", fontsize = 12)
    )
    
    # Store plots and label for this cluster
    plots_by_cluster[[current_cluster]] <- list(
      label = cluster_label,
      plots = cluster_plots
    )
  }
  
  # Create the final layout
  plot_grobs <- list()
  current_grob <- 1
  
  # Calculate maximum maps per row for consistent layout
  max_maps_per_row <- max(table(cluster_assignments$Cluster))
  
  # Create layout matrix
  n_clusters <- length(unique(cluster_assignments$Cluster))
  n_rows <- n_clusters * 2  # One row for label, one for plots per cluster
  layout_matrix <- matrix(NA, nrow = n_rows, ncol = max_maps_per_row + 1)
  
  # Fill layout matrix
  for(cluster in sort(unique(cluster_assignments$Cluster))) {
    cluster_idx <- which(sort(unique(cluster_assignments$Cluster)) == cluster)
    label_row <- (cluster_idx - 1) * 2 + 1
    plots_row <- label_row + 1
    
    # Add label to first column
    plot_grobs[[current_grob]] <- plots_by_cluster[[cluster]]$label
    layout_matrix[label_row, 1] <- current_grob
    current_grob <- current_grob + 1
    
    # Add plots to subsequent columns
    n_plots <- length(plots_by_cluster[[cluster]]$plots)
    for(i in 1:n_plots) {
      plot_grobs[[current_grob]] <- plots_by_cluster[[cluster]]$plots[[i]]
      layout_matrix[plots_row, i + 1] <- current_grob
      current_grob <- current_grob + 1
    }
  }
  
  # Calculate plot dimensions
  plot_height <- n_clusters * map_size * 2  # Double height to account for labels
  
  # Create final combined plot
  combined_plot <- arrangeGrob(
    grobs = plot_grobs,
    layout_matrix = layout_matrix,
    widths = c(1.5, rep(1, max_maps_per_row)),
    heights = rep(c(0.3, 1), n_clusters)
  )
  
  # Save plot
  ggsave(
    filename = file.path(output_dir, "clustered_maps.png"),
    plot = combined_plot,
    width = plot_width,
    height = plot_height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  
  return(combined_plot)
}

#==============================================================================
# UMAP Projections and Clustering
# Calculate F1 score between two masks
#==============================================================================
calculate_mask_f1 <- function(mask1, mask2, grid) {
  # Only consider valid grid positions
  valid_positions <- !is.na(grid)
  
  # Convert to binary vectors considering only valid positions
  mask1_vec <- as.vector(mask1)[valid_positions]
  mask2_vec <- as.vector(mask2)[valid_positions]
  
  # Calculate confusion matrix components
  true_positive <- sum(mask1_vec == 1 & mask2_vec == 1)
  false_positive <- sum(mask1_vec == 1 & mask2_vec == 0)
  false_negative <- sum(mask1_vec == 0 & mask2_vec == 1)
  
  # Calculate precision and recall
  precision <- if(true_positive + false_positive > 0) {
    true_positive / (true_positive + false_positive)
  } else {
    0
  }
  
  recall <- if(true_positive + false_negative > 0) {
    true_positive / (true_positive + false_negative)
  } else {
    0
  }
  
  # Calculate F1 score
  f1_score <- if(precision + recall > 0) {
    2 * (precision * recall) / (precision + recall)
  } else {
    0
  }
  
  return(f1_score)
}

# Calculate centroid for a mask
calculate_mask_centroid <- function(mask, grid) {
  # Get positions of valid mask pixels
  valid_positions <- which(mask == 1 & !is.na(grid), arr.ind = TRUE)
  
  if(nrow(valid_positions) == 0) {
    return(list(x = NA, y = NA))
  }
  
  # Calculate centroid
  centroid_x <- mean(valid_positions[, 2])  # column coordinates
  centroid_y <- mean(valid_positions[, 1])  # row coordinates
  
  return(list(x = centroid_x, y = centroid_y))
}

project_masks_to_umap <- function(mask_features, umap_model) {
  # Print column names for debugging
  print("Columns in mask_features:")
  print(names(mask_features))
  
  print("Columns in embedding_data:")
  print(names(embedding_data))
  
  # Create feature matrix with only the common columns used in the UMAP
  common_cols <- intersect(names(mask_features), names(embedding_data))
  print("Common columns:")
  print(common_cols)
  
  if(length(common_cols) == 0) {
    stop("No common columns found between mask_features and embedding_data")
  }
  
  feature_matrix <- mask_features[, common_cols, drop = FALSE]
  
  # Check for missing columns and fill with 0s if necessary
  missing_cols <- setdiff(names(embedding_data), names(mask_features))
  if(length(missing_cols) > 0) {
    print("Adding missing columns with 0 values:")
    print(missing_cols)
    for(col in missing_cols) {
      feature_matrix[[col]] <- 0
    }
  }
  
  # Ensure columns are in the same order as embedding_data
  feature_matrix <- feature_matrix[, names(embedding_data)]
  
  # Project using UMAP model
  mask_coords <- predict(umap_model, as.matrix(feature_matrix))
  
  # Create results dataframe
  mask_projections <- data.frame(
    Label = mask_features$Label,
    UMAP1 = mask_coords[,1],
    UMAP2 = mask_coords[,2]
  )
  
  return(mask_projections)
}


# Updated mask feature creation
create_mask_features <- function(masks, grid) {
  # Print available masks
  print("Available masks:")
  print(names(masks))
  
  # Define all masks we want to include
  mask_names <- c(
    "ac", "ipsi", "binocular", "peripheral",
    "visual_sky", "visual_ground", "visual_floor",
    "binocular_sky", "binocular_ground", "binocular_floor",
    "peripheral_sky", "peripheral_ground", "peripheral_floor"
  )
  
  # Check which masks are available
  available_masks <- intersect(mask_names, names(masks))
  print("Masks being processed:")
  print(available_masks)
  
  # Initialize results dataframe
  mask_features <- data.frame()
  
  # Calculate all centroids first
  centroids <- lapply(available_masks, function(name) {
    centroid <- calculate_mask_centroid(masks[[name]], grid)
    print(sprintf("Centroid for %s: x=%.2f, y=%.2f", name, centroid$x, centroid$y))
    return(centroid)
  })
  names(centroids) <- available_masks
  
  # Get ranges for normalization
  x_coords <- sapply(centroids, function(c) c$x)
  y_coords <- sapply(centroids, function(c) c$y)
  x_range <- range(x_coords, na.rm = TRUE)
  y_range <- range(y_coords, na.rm = TRUE)
  
  # Pre-calculate all F1 scores
  f1_matrix <- matrix(0, length(available_masks), length(available_masks))
  rownames(f1_matrix) <- available_masks
  colnames(f1_matrix) <- available_masks
  
  for(i in seq_along(available_masks)) {
    for(j in seq_along(available_masks)) {
      f1_matrix[i,j] <- calculate_mask_f1(
        masks[[available_masks[i]]], 
        masks[[available_masks[j]]], 
        grid
      )
    }
  }
  
  print("F1 Score Matrix between masks:")
  print(f1_matrix)
  
  for(mask_name in available_masks) {
    centroid <- centroids[[mask_name]]
    
    # Normalize coordinates to [0,1]
    norm_x <- (centroid$x - x_range[1]) / (x_range[2] - x_range[1])
    norm_y <- (centroid$y - y_range[1]) / (y_range[2] - y_range[1])
    
    # Get F1 scores for this mask against all others
    f1_scores <- setNames(
      f1_matrix[mask_name,],
      paste0(available_masks, "_F1")
    )
    
    # Create row for this mask
    mask_row <- data.frame(
      Label = mask_name,
      ScanX = norm_x,
      ScanY = norm_y
    )
    
    # Add F1 scores
    for(score_name in names(f1_scores)) {
      mask_row[[score_name]] <- f1_scores[score_name]
    }
    
    mask_features <- rbind(mask_features, mask_row)
  }
  
  print("Final mask features:")
  print(mask_features)
  
  return(mask_features)
}

# Also update the visualization function to make masks more visible
create_enhanced_umap_visualization <- function(umap_coords, cluster_assignments, mask_projections,
                                               repel_layers = 2, repel_masks=10) {
  require(deldir)
  require(tidyverse)
  
  # Combine cluster assignments with UMAP coordinates
  plot_data <- merge(umap_coords, cluster_assignments, by = "Label")
  
  # Create a fine grid for the background
  grid_points <- 100
  x_range <- range(plot_data$UMAP1)
  y_range <- range(plot_data$UMAP2)
  margin <- 0.1 * c(diff(x_range), diff(y_range))
  x_seq <- seq(x_range[1] - margin[1], x_range[2] + margin[1], length.out = grid_points)
  y_seq <- seq(y_range[1] - margin[2], y_range[2] + margin[2], length.out = grid_points)
  grid_df <- expand.grid(UMAP1 = x_seq, UMAP2 = y_seq)
  
  # Find nearest point for grid cells
  find_nearest_point <- function(x, y, points) {
    dists <- sqrt((points$UMAP1 - x)^2 + (points$UMAP2 - y)^2)
    return(points$Cluster[which.min(dists)])
  }
  
  # Assign cluster to each grid point
  grid_df$Cluster <- vapply(1:nrow(grid_df), function(i) {
    find_nearest_point(grid_df$UMAP1[i], grid_df$UMAP2[i], plot_data)
  }, FUN.VALUE = numeric(1))
  
  # Split mask labels into groups based on their type
  mask_projections <- mask_projections %>%
    mutate(
      mask_type = case_when(
        grepl("visual", Label) ~ "visual",
        grepl("binocular", Label) ~ "binocular",
        grepl("peripheral", Label) ~ "peripheral",
        TRUE ~ "other"
      ),
      # Add some jitter to help with label placement
      UMAP1 = UMAP1 + runif(n(), -0.1, 0.1),
      UMAP2 = UMAP2 + runif(n(), -0.1, 0.1)
    )
  
  # Create the base plot
  p <- ggplot() +
    # Add tessellated background
    geom_tile(data = grid_df, 
              aes(x = UMAP1, y = UMAP2, fill = factor(Cluster)),
              alpha = 0.2) +  
    # Add data points
    geom_point(data = plot_data,
               aes(x = UMAP1, y = UMAP2, color = factor(Cluster)),
               size = 3) +
    # Add data point labels with minimal overlap
    geom_text_repel(data = plot_data,
                    aes(x = UMAP1, y = UMAP2, label = Label),
                    size = 3,
                    max.overlaps = 20,
                    min.segment.length = 0,
                    box.padding = 0.3,
                    point.padding = 0.1,
                    force = repel_layers,
                    color = "black") +
    # Add mask points
    geom_point(data = mask_projections,
               aes(x = UMAP1, y = UMAP2),
               size = 6,
               shape = 18,
               color = "red",
               alpha = 0.6)
  
  # Add mask labels in layers with different repel parameters
  # Visual masks
  p <- p + geom_text_repel(
    data = filter(mask_projections, mask_type == "visual"),
    aes(x = UMAP1, y = UMAP2, label = Label),
    size = 4,
    color = "red",
    fontface = "bold",
    box.padding = 1,
    point.padding = 0.5,
    force = repel_masks,
    max.overlaps = Inf,
    min.segment.length = 0,
    direction = "both",
    nudge_x = 1
  )
  
  # Binocular masks
  p <- p + geom_text_repel(
    data = filter(mask_projections, mask_type == "binocular"),
    aes(x = UMAP1, y = UMAP2, label = Label),
    size = 4,
    color = "red",
    fontface = "bold",
    box.padding = 1,
    point.padding = 0.5,
    force = repel_masks,
    max.overlaps = Inf,
    min.segment.length = 0,
    direction = "both",
    nudge_y = 1
  )
  
  # Peripheral masks
  p <- p + geom_text_repel(
    data = filter(mask_projections, mask_type == "peripheral"),
    aes(x = UMAP1, y = UMAP2, label = Label),
    size = 4,
    color = "red",
    fontface = "bold",
    box.padding = 1,
    point.padding = 0.5,
    force = repel_masks,
    max.overlaps = Inf,
    min.segment.length = 0,
    direction = "both",
    nudge_x = -1
  )
  
  # Other masks
  p <- p + geom_text_repel(
    data = filter(mask_projections, mask_type == "other"),
    aes(x = UMAP1, y = UMAP2, label = Label),
    size = 4,
    color = "red",
    fontface = "bold",
    box.padding = 1,
    point.padding = 0.5,
    force = repel_masks,
    max.overlaps = Inf,
    min.segment.length = 0,
    direction = "both"
  )
  
  # Customize appearance
  p <- p +
    scale_fill_discrete(name = "Cluster") +
    scale_color_discrete(name = "Cluster") +
    theme_minimal() +
    labs(title = "UMAP Visualization with Anatomical Regions") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.position = "right",
      panel.grid = element_blank(),
      plot.margin = margin(1, 1, 1, 1, "cm")
    )
  
  return(p)
}

#=============================================================================
# Main Analysis Flow
#=============================================================================
if (GENE_MAP_ANALYSIS) {
  # Identify index of valid maps
  rgc_path <- file.path(ROOT_DIR, "rgc_expMat_with_Studyregions_transformedCoords.csv")
  
  rgc_metadata <- c("index", "cell", "x", "y", "z", "volume", "x_range",
                    "y_range", "z_range", "rect_vol", "elongation", "flatness",
                    "slide", "slice", "dapi_max", "dapi_min", "dapi_mean",
                    "dapi_sd", "sample", "retina", "Class", "dapi_max_norm",
                    "dapi_min_norm", "dapi_mean_norm", "dapi_sd_norm",
                    "volume_hull", "nn_dist", "nn_id", "compactness", "sphericity",
                    "X_circ", "Y_circ", "study_region", "X_dv", "Y_dv")
  
  rgc_df <- read_csv(rgc_path)
  
  rgc_gene_summary <- rgc_df %>%
    select(-rgc_metadata) %>%
    group_by(Prediction) %>%
    summarise_all(~mean(.))  %>% # get the mean of each column
    ungroup() %>%
    select(order(colnames(.))) %>%
    arrange(Prediction) %>%
    pivot_longer(-Prediction, names_to = "column", values_to = "value") %>%
    mutate(valid = value >= GENE_THRESHOLD)
  
  valid_maps <- rgc_gene_summary %>%
    pull(valid)
  
  valid_rgc_genes <- rgc_gene_summary %>%
    filter(valid) %>%
    select(-valid)
}

#################################################################################

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

if (GENE_MAP_ANALYSIS) {
  # Filter tif_stack to only include significant layers
  valid_tif_stack <- tif_stack[valid_maps]
  
  message(sprintf("Found %d layers with more than 1 gene on average out of %d total layers",
                  length(valid_tif_stack), length(tif_stack)))
  
  if (length(valid_tif_stack) == 0) {
    stop("No maps showed more than 1 gene on average. Analysis cannot continue.")
  }
  
  tif_stack <- valid_tif_stack
}


# Create study region and grid
study_region <- create_study_region(tif_stack[[1]])
grid <- create_grid(study_region, GRID_SIZE)


stats_path <- file.path(dirs$visual_scene, "spatial_stats_results.RData")

if (!file.exists( stats_path)) {
  # Run Moran's I analysis
  message("Running Moran's I analysis...")
  morans_results <- run_morans_analysis(tif_stack, grid, dirs$morans)
  
  if (GENE_MAP_ANALYSIS) {
    # Filter tif_stack to only include significant layers
    significant_tif_stack <- tif_stack[morans_results$significant_indices]
    clustered_rgc_genes <- valid_rgc_genes %>%
      mutate(index = row_number()) %>%
      filter(index %in% morans_results$significant_indices) %>%
      select(-index)

    message(sprintf("Found %d layers with significant spatial autocorrelation out of %d total layers",
                    length(significant_tif_stack), length(tif_stack)))

    if (length(significant_tif_stack) == 0) {
      stop("No layers showed significant spatial autocorrelation. Analysis cannot continue.")
    }

    tif_stack <- significant_tif_stack
  }
  
  # Run scan statistics analysis
  message("Running scan statistics analysis...")
  scan_results <- run_scan_analysis(tif_stack, grid, dirs$scan)
  
  if (GENE_MAP_ANALYSIS) {
    save(scan_results, morans_results, clustered_rgc_genes, file = stats_path)
    } else {
    save(scan_results, morans_results, file = stats_path)
    }

} else {
  message("Loading pre-computed Scan and Moran's data...")
  load(stats_path)
}


# Create and analyze visual scene masks
message("Running visual scene analysis...")
masks <- create_visual_scene_masks(grid)
saveRDS(masks, file.path(dirs$visual_scene, "visual_scene_masks.rds"))


# Create VisualRegionMasks directory
masks_dir <- file.path(dirs$visual_scene, "VisualRegionMasks")
dir.create(masks_dir, showWarnings = FALSE)

# Save each mask as a CSV
for(mask_name in names(masks)) {
  write.csv(
    masks[[mask_name]], 
    file = file.path(masks_dir, paste0(mask_name, "_mask.csv")), 
    row.names = FALSE
  )
}


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

# create_visual_scene_grid_plots(
#   visual_scene_results = visual_scene_results,
#   scan_results = scan_results$results,  
#   grid = grid,
#   mask_list = masks,
#   tif_stack = tif_stack,
#   output_dir = dirs$visual_scene,
#   min_coverage = MIN_COVERAGE
# )

# After your existing visual scene analysis
message("Running F1 score analysis...")
f1_results <- analyze_f1_scores(
  scan_results$results,
  grid,
  masks,
  dirs$visual_scene
)

# Usage example:
message("Running enrichment analysis...")
enrichment_results <- analyze_visual_scene_enrichment(
  f1_results,
  dirs$visual_scene,
  masks
)

# Compile all stats
message("Compile Statistical Summary Table...")
if (GENE_MAP_ANALYSIS) {
  stats_results <- compile_statistics(OUTPUT_DIR, valid_rgc_genes, clustered_rgc_genes)
} else {
  stats_results <- compile_statistics(OUTPUT_DIR)
}

# Compute UMAP
message("Computing UMAP and Clustering...")
# Read normalized embedding data
embedding_path <- file.path(file.path(ROOT_DIR, OUTPUT_DIR), "embedding_normalized.csv")
if(!file.exists(embedding_path)) {
  stop(sprintf("Cannot find embedding file: %s", embedding_path))
}
embedding_data <- read.csv(embedding_path)%>%
  select(
         ScanX,
         ScanY,
         ac_F1,
         ipsi_F1,
         binocular_F1,
         peripheral_F1, 
         visual_sky_F1,                             
         visual_ground_F1, 
         visual_floor_F1,
         binocular_sky_F1,
         binocular_ground_F1,
         peripheral_sky_F1,
         peripheral_ground_F1
         ) 

# Generate UMAPs
umap_results <- create_umap_visualization(embedding_data, OUTPUT_DIR)

# Cluster the UMAP
message("Perform Hierarchical Clustering...")
cluster_results <- perform_hierarchical_clustering(OUTPUT_DIR)
print(cluster_results$plots$umap)

# Load cluster assignments
if (GENE_MAP_ANALYSIS) {
  cluster_assignments <- read.csv(file.path(ROOT_DIR, OUTPUT_DIR, "cluster_assignments.csv")) %>%
    cbind(clustered_rgc_genes) %>%
    select(-Label)
  write.csv(clustered_rgc_genes, file = file.path(ROOT_DIR, OUTPUT_DIR, "cluster_assignments.csv"))
} else {
  cluster_assignments <- read.csv(file.path(ROOT_DIR, OUTPUT_DIR, "cluster_assignments.csv"))
}

# Create visualization of clusters
plot <- visualize_clustered_maps(
  tif_stack = tif_stack,
  grid = grid,
  cluster_assignments = cluster_assignments,
  output_dir = file.path(ROOT_DIR, OUTPUT_DIR)
)

#######################################################################################333


message("Computing mask features and projecting onto UMAP...")

# Create mask features
mask_features <- create_mask_features(masks, grid)

# Print features for verification
print("Mask Features:")
print(mask_features)

# Create mask features with debugging info
mask_features <- create_mask_features(masks, grid)

# Project masks onto UMAP space
mask_projections <- project_masks_to_umap(mask_features, umap_results$umap_object)

# Create enhanced visualization
enhanced_umap <- create_enhanced_umap_visualization(
  umap_coords = umap_results$coordinates,
  cluster_assignments = cluster_assignments,
  mask_projections = mask_projections,
  repel_layers = 1,
  repel_masks = 10
)

print(enhanced_umap)

# Save results
write.csv(
  mask_features,
  file = file.path(ROOT_DIR, OUTPUT_DIR, "mask_features.csv"),
  row.names = FALSE
)

write.csv(
  mask_projections,
  file = file.path(ROOT_DIR, OUTPUT_DIR, "mask_projections.csv"),
  row.names = FALSE
)


#####################################################################################
# Define path to scan results
DENSITY_MAPS_PATH <- '/home/sam/FinalRGC_xenium/GlobalStatistics_Densities/Scan/ScanPortable'

# Function to load and process a single CSV file
process_density_csv <- function(file_path) {
  # Read CSV
  mask_data <- try(read.csv(file_path, header = TRUE))
  
  if(inherits(mask_data, "try-error")) {
    warning(sprintf("Error reading file: %s", file_path))
    return(NULL)
  }
  
  # Convert to matrix
  mask_matrix <- as.matrix(mask_data)
  
  # Convert to boolean mask: -1 -> 0, everything else -> 1
  boolean_mask <- matrix(0, nrow=nrow(mask_matrix), ncol=ncol(mask_matrix))
  boolean_mask[mask_matrix != -1] <- 1
  
  return(boolean_mask)
}

# List all CSV files
csv_files <- list.files(DENSITY_MAPS_PATH, pattern = "layer_.*_pvalues\\.csv$", full.names = TRUE)

if(length(csv_files) == 0) {
  stop("No layer_*_pvalues.csv files found in specified directory")
}

# Sort files to ensure consistent ordering
csv_files <- sort(csv_files)

# Initialize density_masks list
density_masks <- list()

# Process each file
for(file_path in csv_files) {
  # Extract layer number from filename
  layer_num <- as.numeric(gsub(".*layer_([0-9]+)_pvalues\\.csv$", "\\1", file_path))
  
  if(is.na(layer_num)) {
    warning(sprintf("Could not extract layer number from filename: %s", file_path))
    next
  }
  
  # Process file
  mask <- process_density_csv(file_path)
  
  if(!is.null(mask)) {
    density_masks[[layer_num]] <- mask
  }
}

# Verify results
n_masks <- length(density_masks)
message(sprintf("Loaded %d density masks", n_masks))

# Basic validation checks
if(n_masks > 0) {
  # Check dimensions of first mask
  mask_dims <- dim(density_masks[[1]])
  message(sprintf("Mask dimensions: %d x %d", mask_dims[1], mask_dims[2]))
  
  # Check that all masks have the same dimensions
  dims_match <- all(sapply(density_masks, function(mask) {
    all(dim(mask) == mask_dims)
  }))
  
  if(!dims_match) {
    warning("Not all masks have the same dimensions!")
  }
  
  # Print summary of mask values
  value_summary <- sapply(density_masks, function(mask) {
    c(sum(mask == 1), sum(mask == 0))
  })
  
  message("Summary of mask values (ones/zeros):")
  print(t(value_summary))
}

# Save the density_masks object for future use
saveRDS(density_masks, file = file.path(dirname(DENSITY_MAPS_PATH), "density_masks.rds"))
message(sprintf("Density masks saved to %s", 
                file.path(dirname(DENSITY_MAPS_PATH), "density_masks.rds")))




# First, properly name the density masks
named_density_masks <- list()
for(i in seq_along(density_masks)) {
  mask_name <- paste0("density_", i)
  named_density_masks[[mask_name]] <- density_masks[[i]]
}

# Create density overlaps directory
density_dir <- file.path(main_dir, "density_overlaps")
dir.create(density_dir, showWarnings = FALSE)

# Create subdirectory structure to match original analysis
density_dirs <- list(
  main = density_dir,
  f1 = file.path(density_dir, "F1"),
  enrichment = file.path(density_dir, "Enrichment")
)

# Create all subdirectories
lapply(density_dirs, dir.create, showWarnings = FALSE)

# Run F1 score analysis for density masks
message("Running F1 score analysis for density masks...")
density_f1_results <- analyze_f1_scores(
  scan_results$results,
  grid,
  named_density_masks,  # Using named density masks
  density_dirs$main
)

# Verify F1 results
message("F1 Analysis complete. Checking results...")
message("Number of F1 scores calculated: ", nrow(density_f1_results$results))
print(head(density_f1_results$results))

# Run enrichment analysis for density masks
message("Running enrichment analysis for density masks...")
density_enrichment_results <- analyze_visual_scene_enrichment(
  density_f1_results,
  density_dirs$main,
  named_density_masks  # Using named density masks
)

# Save results
saveRDS(density_f1_results, file.path(density_dirs$main, "density_f1_results.rds"))
saveRDS(density_enrichment_results, file.path(density_dirs$main, "density_enrichment_results.rds"))
saveRDS(named_density_masks, file.path(density_dirs$main, "named_density_masks.rds"))

# Print summary
message("\nAnalysis Complete!")
message(sprintf("Results saved in: %s", density_dir))
message("Number of masks analyzed: ", length(named_density_masks))
message("Number of F1 scores calculated: ", nrow(density_f1_results$results))
message("Number of enrichment tests performed: ", nrow(density_enrichment_results$results))







# Function to update comprehensive statistics with density analysis results
update_comprehensive_statistics <- function(comp_stats_path, 
                                            density_f1_results, 
                                            density_enrichment_results,
                                            is_gene_analysis = FALSE) {
  
  # Read existing comprehensive statistics
  comp_stats <- read.csv(comp_stats_path)
  
  # Extract F1 scores
  f1_scores <- density_f1_results$results
  
  # Process enrichment statistics
  enrichment_stats <- density_enrichment_results$results %>%
    mutate(
      neg_log_p = ifelse(P_adj == 0, 
                         -log10(.Machine$double.eps),
                         -log10(P_adj)),
      neglog_col_name = paste0(Mask, "_negative10Friedmanpvalue"),
      padj_col_name = paste0(Mask, "_Padj"),
      p_col_name = paste0(Mask, "_P"),
      odds_col_name = paste0(Mask, "_OddsRatio")
    ) %>%
    select(Layer, Mask, neg_log_p, P_adj, P_value, Odds_ratio) %>%
    pivot_longer(
      cols = c(neg_log_p, P_adj, P_value, Odds_ratio),
      names_to = "stat_type",
      values_to = "value"
    ) %>%
    mutate(
      col_name = case_when(
        stat_type == "neg_log_p" ~ paste0(Mask, "_negative10Friedmanpvalue"),
        stat_type == "P_adj" ~ paste0(Mask, "_Padj"),
        stat_type == "P_value" ~ paste0(Mask, "_P"),
        stat_type == "Odds_ratio" ~ paste0(Mask, "_OddsRatio")
      )
    ) %>%
    select(Layer, col_name, value) %>%
    pivot_wider(
      names_from = col_name,
      values_from = value
    )
  
  # Initialize new columns for F1 scores
  for(mask_name in names(named_density_masks)) {
    f1_col <- paste0(mask_name, "_F1")
    comp_stats[[f1_col]] <- NA
  }
  
  # Initialize new columns for enrichment statistics
  enrichment_cols <- setdiff(names(enrichment_stats), "Layer")
  for(col in enrichment_cols) {
    comp_stats[[col]] <- NA
  }
  
  # Update values based on whether it's gene analysis or layer analysis
  if(is_gene_analysis) {
    # For gene analysis, we need to match by both subtype and gene
    for(i in seq_len(nrow(f1_scores))) {
      layer <- f1_scores$Layer[i]
      # Find matching row in comp_stats
      matching_row <- which(comp_stats$Layer == layer)
      
      if(length(matching_row) == 1) {
        # Update F1 scores
        for(mask_name in names(named_density_masks)) {
          f1_col <- paste0(mask_name, "_F1")
          comp_stats[matching_row, f1_col] <- f1_scores[i, paste0(mask_name, "_F1")]
        }
        
        # Update enrichment statistics
        if(layer %in% enrichment_stats$Layer) {
          enrich_row <- which(enrichment_stats$Layer == layer)
          for(col in enrichment_cols) {
            comp_stats[matching_row, col] <- enrichment_stats[enrich_row, col]
          }
        }
      }
    }
  } else {
    # For layer analysis, we can match directly by Layer
    for(i in seq_len(nrow(f1_scores))) {
      layer <- f1_scores$Layer[i]
      matching_row <- which(comp_stats$Layer == layer)
      
      if(length(matching_row) == 1) {
        # Update F1 scores
        for(mask_name in names(named_density_masks)) {
          f1_col <- paste0(mask_name, "_F1")
          comp_stats[matching_row, f1_col] <- f1_scores[i, paste0(mask_name, "_F1")]
        }
        
        # Update enrichment statistics
        if(layer %in% enrichment_stats$Layer) {
          enrich_row <- which(enrichment_stats$Layer == layer)
          for(col in enrichment_cols) {
            comp_stats[matching_row, col] <- enrichment_stats[enrich_row, col]
          }
        }
      }
    }
  }
  
  # Save updated comprehensive statistics
  write.csv(comp_stats, comp_stats_path, row.names = FALSE)
  
  # Return updated dataframe
  return(comp_stats)
}

# Use the function
comp_stats_path <- file.path(main_dir, "comprehensive_statistics.csv")
updated_stats <- update_comprehensive_statistics(
  comp_stats_path = comp_stats_path,
  density_f1_results = density_f1_results,
  density_enrichment_results = density_enrichment_results,
  is_gene_analysis = GENE_MAP_ANALYSIS
)

message("Comprehensive statistics updated successfully!")
message("New columns added:")
new_cols <- c(
  paste0("density_", seq_along(named_density_masks), "_F1"),
  names(enrichment_stats)[!names(enrichment_stats) %in% c("Layer")]
)
message(paste(new_cols, collapse = "\n"))