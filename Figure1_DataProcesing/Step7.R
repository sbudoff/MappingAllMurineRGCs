library(tidyverse)
library(dplyr)
library(ggthemes)
library(patchwork)

################################################################################
########################## HELPER FUNCTIONS ####################################
################################################################################

nn_computer <- function(df) {
  if ("nn_dist" %in% colnames(df)){
    df <- select(df, -c(nn_dist, nn_id)) # Remove Old NN columns if present
  }
  
  nn_list <- list()
  for (r in unique(df$retina)) {
    for (s in unique(df$slice)) {
      for (p in unique(df$Prediction)) {
        df_i <- df %>%
          filter(retina == r & slice == s & Prediction == p)
        
        if (nrow(df_i) > 1) {
          df_i <- df_i %>%
            select(cell,x,y,z) %>%
            column_to_rownames("cell")
          dists <- as.matrix(dist(df_i))
          diag(dists) <- NA
          nn_dist <- apply(dists,1,min,na.rm=TRUE)
          nn_id <- apply(dists, 1, function(x) rownames(df_i)[which.min(x)])# find the cell id of the nearest neighbor
          nn_list[[paste0(r,s,p)]] <- data.frame(nn_dist) %>%
            rownames_to_column("cell") %>%
            mutate(nn_id = nn_id)
        }
      }
    }
  }
  nn_df <- bind_rows(nn_list)

  df <- left_join(df, nn_df, by="cell")
  return(df)
}

NNRI <- function(df) {
  # Input validation
  required_cols <- c("x", "y")
  if (!all(required_cols %in% colnames(df))) {
    stop("DataFrame must contain columns: x, y")
  }
  
  dists <- df %>%
    select(x,y) %>%
    dist() %>%
    as.matrix()
  diag(dists) <- NA
  nn_dist <- apply(dists,1,min,na.rm=TRUE)
  nn_mu <- mean(nn_dist)
  nn_sigma <- sd(nn_dist)
  nnri <- nn_mu/nn_sigma
  return(nnri)
}

VDRI <- function(df) {
  # Input validation
  required_cols <- c("x", "y")
  if (!all(required_cols %in% colnames(df))) {
    stop("DataFrame must contain columns: x, y")
  }
  
  veroni <- deldir::deldir(df$x, df$y)
  veroni <- veroni[['summary']]
  veroni_area <- veroni$dir.area
  
  vd_mu <- mean(veroni_area)
  vd_sigma <- sd(veroni_area)
  vdri <- vd_mu/vd_sigma
  return(vdri)
}
################################################################################
########################## DATA LOADING ########################################
################################################################################
root <- '/media/sam/Data2/baysor_rbpms_consolidated'
out_path <- paste0(root,'/all_retinas_prediction_expmat.csv')
rgc_out_path <- paste0(root,'/all_rgc_prediction_expmat.csv')
f_rgc_out_path <- paste0(root,'/filtered_rgc_prediction_expmat.csv')
# Load the gene order for cuttlenet inference
gene_order <- read.csv('/media/sam/Data2/baysor_analysis/CuttleNet/input_order.csv', header = F) %>% pull()

# Find all processed segmentation expression matrices
slides <- list.files(root)
expmat_list <- list()
keys <- c()
paths <- c()
dapi_paths <- c()
for (slide in slides) {
  path <- paste0(root,"/",slide)
  for (p in list.files(path)){
    p_i <- paste0(root,"/",slide,"/",p)
    for (file in list.files(p_i)){
      if (file == 'segmentation.csv') {
        
        # Construct the expected output path
        sample_root <- substr(p_i,1, nchar(p_i))
        sample_id <- substr(p_i,50, nchar(p_i)-7)
        expected_output <- paste0(sample_root, '/', sample_id, 'expression_matrix.csv')
        expected_output_dapi <- paste0(sample_root, '/dapi_statistics.csv')
        # Only add to paths if output doesn't exist
        if (file.exists(expected_output)) {
          paths <- c(paths, expected_output)
          dapi_paths <- c(dapi_paths, paste0(p_i, '/dapi_statistics.csv'))
          print(paste("expression matrix found for", p))
          if (!file.exists(out_path)) {
            expmat_list[[p]] <- read_csv(expected_output) %>%
              mutate(slide = substr(p,1,5),
                     slice = substr(p,7,11)) %>%
              left_join(select(read_csv(expected_output_dapi),-1), by = 'cell')
            print(paste(p, "processed"))
          }
          keys <- c(keys, p)
        } 
      }
    }
  }
}
# Bind all retinas together or load the pre-processed file
if (!file.exists(out_path)) {
  retinas_df <- bind_rows(expmat_list)
  
  retinas_df <- retinas_df %>%
    mutate(slide = factor(slide),
           sample = slice,
           retina = factor(substr(slice,2,2)),
           slice = factor(substr(slice,4,5)),
           Class = case_when(
             substr(Prediction, 1, 2) == "AC" ~ "AC",
             grepl("^[[:digit:]]{2}", Prediction) ~ "RGC",
             substr(Prediction, 2, 3) == "BC" ~ "BC", 
             grepl( "Photoreceptors", Prediction, fixed=T) ~ "Ph",
             TRUE ~ "Other"
           ),
           Class = factor(Class),
           Prediction = factor(Prediction)) %>%
    group_by(slice, Class) %>%
    mutate(
      dapi_max_norm = (dapi_max -  min(dapi_min)) / (max(dapi_max) - min(dapi_min)),
      dapi_min_norm = (dapi_min -  min(dapi_min)) / (max(dapi_max) - min(dapi_min)),
      dapi_mean_norm = (dapi_mean -  min(dapi_min)) / (max(dapi_max) - min(dapi_min)),
      dapi_sd_norm = (dapi_sd  -  min(dapi_min)) / (max(dapi_max) - min(dapi_min))
    ) %>%
    ungroup()
  
  write.csv(retinas_df,out_path, row.names = FALSE)
} else {
  retinas_df <- read_csv(out_path)
}
rm(expmat_list)

retinas_df <- retinas_df %>%
  mutate(slide = factor(slide),
         sample= factor(sample),
         retina = factor(retina),
         slice = factor(slice),
         Class = factor(Class),
         Prediction = factor(Prediction))

rgc_df <- retinas_df %>%
  filter(Class == 'RGC') %>%
  mutate(volume_hull = volume,
         volume =(4/3)*pi*(0.5*x_range)*(0.5*y_range)*(0.5*z_range)) 
%>%
  nn_computer()

if (!file.exists(rgc_out_path)) {
  write.csv(rgc_df,rgc_out_path, row.names = FALSE)
}

# Load the Bae 2018 volume background information
# Get all unique predictions
all_subtypes <- unique(rgc_df$Prediction)

background <- read_csv('/home/sam/RGC_Background/RGC Subtype Unification - Aligned Names.csv') %>%
  select(Tran2019_Clusters, Bae2018_25th, Bae2018_75th) %>%
  drop_na() %>%
  mutate(subtype = sapply(Tran2019_Clusters, function(x) {
    # Convert number to two-digit string format
    pattern <- sprintf("%02d", x)
    # Find matching prediction from rgc_df
    matching_pred <- grep(paste0("^", pattern, "_"), all_subtypes, value = TRUE)
    if(length(matching_pred) > 0) return(matching_pred[1])
    return(NA_character_)
  })) %>%
  select(-Tran2019_Clusters)

# Find missing subtypes from background
missing_subtypes <- setdiff(all_subtypes, background$subtype)

# Create new rows for missing subtypes
if(length(missing_subtypes) > 0) {
  new_rows <- data.frame(
    subtype = missing_subtypes,
    Bae2018_25th = min(background$Bae2018_25th, na.rm = TRUE),
    Bae2018_75th = max(background$Bae2018_75th, na.rm = TRUE)
  )
  
  # Bind the new rows to the original dataset
  background <- bind_rows(background, new_rows)
}

background <- background %>%
  group_by(subtype) %>%
  summarize(Bae2018_25th = min(Bae2018_25th),
            Bae2018_75th = max(Bae2018_75th)) %>%
  ungroup()
################################################################################
########################## Basic Statistics ####################################
################################################################################


retinas_stats <- retinas_df %>%
  group_by(Class, Prediction, retina, slice, sample) %>%
  summarize(volume = mean(volume),
            elongation = mean(elongation),
            flatness = mean(flatness),
            sphericity = mean(sphericity),
            compactness = mean(compactness),
            N = n()) %>%
  ungroup() %>%
  group_by(Class, retina, slice) %>%
  mutate(Proportion = N/sum(N)) %>%
  ungroup()
  
retinas_stats %>%
  filter(Class == "RGC" | Class == "AC") %>%
  ggplot(aes(x = Class, y = volume, color = sample)) +
  geom_jitter() +
  ggtitle("Volume of Predicted Subtypes by Class") +
  theme_tufte()

retinas_stats %>%
  filter(Class == "RGC" | Class == "AC") %>%
  mutate(sample = factor(paste(slice, retina))) %>%
  group_by(Class, Prediction) %>%
  mutate(mean_volume = mean(volume)) %>%
  ungroup() %>%
  ggplot(aes(y = reorder(Prediction, -mean_volume), x = volume)) +
  geom_violin() +
  geom_jitter(aes(color = sample), width = 0.1, alpha = 0.5) +
  ggtitle("Volume of Predicted Subtypes by Descending Mean Volume") +
  labs(x = "Proportion", y = "Prediction") +
  theme_tufte() +
  facet_wrap(~Class)



# Step 1: Calculate mean proportions for unique Prediction-Class combinations
mean_proportions <- retinas_stats %>%
  filter(Class %in% c("RGC", "AC") & volume > 500) %>%
  group_by(Class, Prediction) %>%
  summarize(mean_Proportion = mean(Proportion), .groups = 'drop') %>%
  arrange(Class, desc(mean_Proportion))

# Step 2: Set Prediction as an ordered factor within each Class based on mean_Proportion
retinas_stats %>%
  filter(Class %in% c("RGC", "AC") & volume > 500)   %>%
  mutate(sample = factor(paste(slice, retina))) %>%
  left_join(mean_proportions, by = c("Class", "Prediction")) %>%
  group_by(Class) %>%
  mutate(Prediction_ordered = factor(Prediction, 
                                     levels = unique(Prediction[order(-mean_Proportion)]))) %>%
  ungroup() %>%
  # Step 3: Plot with ordered Prediction
  ggplot(aes(x = Proportion, y = Prediction_ordered, color = sample, group = sample)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_line(alpha = 0.7) +
  ggtitle("Proportion of Predicted Subtypes by Descending Mean Proportion") +
  labs(x = "Proportion", y = "Prediction") +
  facet_grid(rows = vars(Class), 
             cols = vars(sample), 
             scales = "free_y",   # This allows different y-scales for each row
             space = "free_y") +  # This adjusts the spacing based on number of levels
  theme_tufte()

################################################################################
########################## RGC NN Stats ########################################
################################################################################
rgc_df <- retinas_df %>%
  filter(Class == 'RGC' & volume > 500 & dapi_mean_norm > 0.05) %>%
  nn_computer()

if (!file.exists(rgc_out_path)) {
  write.csv(rgc_df,f_rgc_out_path, row.names = FALSE)
}
median_nn_dist <- rgc_df %>%
  group_by(Prediction) %>%
  summarize(median_nn_dist = median(nn_dist, na.rm=T), .groups = 'drop') %>%
  arrange(desc(median_nn_dist))

rgc_df %>%
  left_join(median_nn_dist, by = "Prediction") %>%
  mutate(Prediction_ordered = factor(Prediction, 
                                     levels = unique(Prediction[order(-median_nn_dist)]))) %>%
  ggplot(aes(y = Prediction_ordered, x = nn_dist)) +
  geom_violin() +
  geom_jitter(aes(color = sample), height = 0.2, alpha = 0.05, size=1) +
  ggtitle("NN Dist of Predicted Subtypes by Descending Median Distance") +
  labs(x = "NN Distance", y = "Prediction") +
  scale_x_log10() + 
  theme_tufte() 

rgc_df  %>%
  ggplot(aes(x= Prediction, fill = sample)) +
    geom_bar(position = "dodge")+ 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


retinas_df %>%
  group_by(sample, Prediction) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(prop = count/sum(count)) %>%
  ungroup() %>%
ggplot(aes(x= Prediction, y = prop, fill = sample)) +
  geom_col(position = "identity")+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

################################################################################
############################# Retina Mosaics ###################################
################################################################################
# https://pmc.ncbi.nlm.nih.gov/articles/PMC7368823/#FN1
#Two spatial statistics associated with these mosaics, 
  #the “regularity index” derived from an analysis of nearest neighbor distances, and 
  #the “effective radius” determined from the density recovery profile (informally termed the exclusion zone),
#have each come to be regarded as definitive evidence that a population of labeled neurons is indeed “regular” 

# By expressing the variation in those nearest neighbor distances relative to their average 
  #(specifically, the ratio of the mean nearest neighbor distance to the standard deviation), 
  #a summary statistic, the “regularity index”, seemed to yield a simple quantitation of this property, 
  #increasing towards infinity as a mosaic approached a perfect lattice.

# For each of the real fields and their random simulations, a regularity index can be computed, 
  #being the mean nearest neighbor distance divided by the standard deviation (NNRI). 

find_sample_regions <- function(df, max_diameter = 1000, n_cells = 10, n_regions = 30, max_rolls = 1000) {
  # Remove the incorrect assignment of df
  # df <- updated_rgc_df  # This line was overwriting the input parameter
  
  # Input validation
  required_cols <- c("cell", "Prediction", "slide", "x", "y")
  if (!all(required_cols %in% colnames(df))) {
    stop("DataFrame must contain columns: cell, Prediction, slide, x, y")
  }
  
  # Initialize list to store results
  sample_regions <- list()
  list_index <- 1
  
  slides <- unique(df$slide)
  n_slides <- length(slides)
  samples_per_slide <- ceiling(n_regions/n_slides)
  
  subtypes <- unique(df$Prediction)
  
  # Process each eligible pair
  pb <- progress_bar$new(
    total = length(slides)*length(subtypes),
    format = "Extracting sample regions [:bar] :percent eta: :eta"
  )
  
  # Process each slide separately
  for (current_slide in slides) {
    # Create a copy of slide data that we'll update as we remove points
    slide_data <- df[df$slide == current_slide, ]
    regions_found <- 0
    
    for (subtype in subtypes) {
      pb$tick()
      remaining_data <- slide_data[slide_data$Prediction == subtype, ]
      
      # Skip if not enough cells of this subtype
      if (nrow(remaining_data) < n_cells) {
        next
      }
      
      # Compute distance matrix
      d <- remaining_data %>%
        select(x, y) %>%
        dist() %>%
        as.matrix()
      
      # Initialize variables
      samples_found <- 0
      roll <- 0
      
      while (samples_found < samples_per_slide && roll < max_rolls) {
        roll <- roll + 1
        
        # Skip if not enough remaining data
        if (nrow(remaining_data) < n_cells) {
          break
        }
        
        # Randomly pick a cell and find the distance to its [n_cells]th nearest neighbor within max_diameter
        rand_idx <- sample(nrow(remaining_data), 1)
        rand_cell <- remaining_data[rand_idx, ]
        
        nn_dists <- sort(d[rand_idx, ])
        if (length(nn_dists) >= n_cells) {
          if (nn_dists[n_cells] > max_diameter) {
            next # Skip this cell if no neighbors within max_diameter
          }
        } else {
          next # Skip this cell if too few neighbors within max_diameter
        }
        
        # Get the radius to include exactly n_cells neighbors
        dist_n_cells <- nn_dists[n_cells]
        
        # Calculate the center of the circular area based on the selected random cell
        center_x <- rand_cell$x
        center_y <- rand_cell$y
        
        # Find all cells (regardless of subtype) within the circle in the slide_data
        points_in_circle <- slide_data %>%
          filter(sqrt((x - center_x)^2 + (y - center_y)^2) <= dist_n_cells)
        
        # Verify we have enough cells of the target subtype
        n_subtype_cells <- nrow(points_in_circle[points_in_circle$Prediction == subtype, ])
        if (n_subtype_cells < n_cells) {
          next
        }
        
        # Store the sample region
        sample_regions[[list_index]] <- list()
        sample_regions[[list_index]][["df_i"]] <- points_in_circle
        sample_regions[[list_index]][["x"]] <- center_x
        sample_regions[[list_index]][["y"]] <- center_y
        sample_regions[[list_index]][["r"]] <- dist_n_cells
        sample_regions[[list_index]][["subtype"]] <- subtype
        sample_regions[[list_index]][["N"]] <- n_subtype_cells
        
        # Update remaining data and counters
        regions_found <- regions_found + 1
        samples_found <- samples_found + 1
        list_index <- list_index + 1
        remaining_data <- remaining_data[!(remaining_data$cell %in% points_in_circle$cell), ]
        
        # Update distance matrix for remaining points
        if (nrow(remaining_data) >= 2) {
          d <- remaining_data %>%
            select(x, y) %>%
            dist() %>%
            as.matrix()
        }
      }
    }
  }
  
  return(sample_regions)
}

runif_circle <- function(n, center_x = 0, center_y = 0, diameter = 200) {
  # Input validation
  if (n <= 0 || !is.numeric(n)) {
    stop("n must be a positive number")
  }
  if (!is.numeric(center_x) || !is.numeric(center_y) || !is.numeric(diameter)) {
    stop("center_x, center_y, and diameter must be numeric")
  }
  if (diameter <= 0) {
    stop("diameter must be positive")
  }
  
  radius <- diameter / 2
  points <- data.frame(x = numeric(0), y = numeric(0))
  
  # Generate points using rejection sampling
  while(nrow(points) < n) {
    # Generate more points than needed to account for rejection
    n_extra <- ceiling((n - nrow(points)) * 1.3)  # 30% extra for efficiency
    
    # Generate random points in a square
    candidate_x <- runif(n_extra, -radius, radius)
    candidate_y <- runif(n_extra, -radius, radius)
    
    # Keep only points that fall within the circle
    in_circle <- sqrt(candidate_x^2 + candidate_y^2) <= radius
    
    # Add valid points to our collection
    new_points <- data.frame(
      x = candidate_x[in_circle] + center_x,
      y = candidate_y[in_circle] + center_y
    )
    
    points <- rbind(points, new_points)
  }
  
  # Trim to exact number of points needed
  points <- points[1:n, ]
  
  return(points)
}

# sample_regions <- find_sample_regions(rgc_df)
# sample_regions[[1]] %>%
#   # filter(Prediction == '01_W3D1.1')%>%
#   ggplot(aes(x=x, y = y, color = Prediction))+
#   geom_point(size=1, alpha =0.7) +
#   theme_tufte() 



RI_bootstrapper <- function(df_i, N, simulations = 100){
  NNRIs <- c()
  VDRIs <- c()
  for (i in 1:simulations){
    s_df <- df_i %>%
      slice_sample(n = N, replace=F)
    NNRIs <- c(NNRIs, NNRI(s_df))
    VDRIs <- c(VDRIs, VDRI(s_df))
  }
  output <- list("NNRI" = NNRIs, "VDRI" = VDRIs)
  return(output)
}
RI_simulator <- function(observed_sample, simulations = 100, subtype = '01_W3D1.1', diameter=200) {
  
  s_df <- observed_sample %>%
    filter(Prediction == subtype)
  
  x_0 <- mean(observed_sample$x) 
  y_0 <- mean(observed_sample$y) 
  N  <- nrow(s_df)
  NNRIs <- c()
  VDRIs <- c()
  for (i in 1:simulations){
    csr <- runif_circle(n=N, center_x=x_0, center_y=y_0, diameter=diameter)
    NNRIs <- c(NNRIs, NNRI(csr))
    VDRIs <- c(VDRIs, VDRI(csr))
  }
  output <- list("NNRI" = NNRIs, "VDRI" = VDRIs)
  return(output)
}



sample_regions <- find_sample_regions(updated_rgc_df, 
                                      max_diameter = 1000, 
                                      n_cells = 10, 
                                      n_regions = 30)
output <- list()
j = 1
N_s = 0
for (i in 1:length(sample_regions)){
  samp_i <- sample_regions[[i]] 
  df_i <- samp_i[['df_i']]
  N <- samp_i[['N']]
  d <- 2*samp_i[['r']] 
  s <- samp_i[['subtype']] 

  NNRI_obs <- NNRI(df_i)
  VDRI_obs <- VDRI(df_i)
  RIs_csr <- RI_simulator(df_i, subtype=s, diameter=d)
  RIs_boot <- RI_bootstrapper(df_i, N)
  output[[j]] <- list()
  output[[j]][['subtype']] <- s
  output[[j]][['N']] <- N
  output[[j]][['NNRI_obs']] <- NNRI_obs
  output[[j]][['VDRI_obs']] <- VDRI_obs
  output[[j]][['NNRIs_csr']] <- RIs_csr[['NNRI']]
  output[[j]][['NNRI_csr']] <- mean(RIs_csr[['NNRI']])
  output[[j]][['NNRIs_boot']] <- RIs_boot[['NNRI']]
  output[[j]][['NNRI_boot']] <- mean(RIs_boot[['NNRI']])
  
  output[[j]][['VDRIs_csr']] <- RIs_csr[['VDRI']]
  output[[j]][['VDRI_csr']] <- mean(RIs_csr[['VDRI']])
  output[[j]][['VDRIs_boot']] <- RIs_boot[['VDRI']]
  output[[j]][['VDRI_boot']] <- mean(RIs_boot[['VDRI']])
  j = j+1
  N_s = N_s+1
  print(paste("For sample", i , N, "cells observed", s, "NNRI =", round(NNRI_obs, 2), "VDRI =", round(VDRI_obs, 2)))
}




df_obs <- data.frame(
  subtype = sapply(output, function(x) x$subtype),
  N = sapply(output, function(x) x$N),
  NNRI = sapply(output, function(x) x$NNRI_obs),
  NNRI_sim = sapply(output, function(x) x$NNRI_csr),
  NNRI_boot = sapply(output, function(x) x$NNRI_boot),
  VDRI = sapply(output, function(x) x$VDRI_obs),
  VDRI_sim = sapply(output, function(x) x$VDRI_csr),
  VDRI_boot = sapply(output, function(x) x$VDRI_boot)
) 
df_obs_summary <- df_obs %>%
  group_by(subtype) %>%
  summarise(NNRI_sd = sd(NNRI),
            NNRI_sim_sd = sd(NNRI_sim),
            NNRI_boot_sd = sd(NNRI_boot),
            NNRI = mean(NNRI),
            NNRI_sim = mean(NNRI_sim),
            NNRI_boot = mean(NNRI_boot),
            VDRI_sd = sd(VDRI),
            VDRI_sim_sd = sd(VDRI_sim),
            VDRI_boot_sd = sd(VDRI_boot),
            VDRI = mean(VDRI),
            VDRI_sim = mean(VDRI_sim),
            VDRI_boot = mean(VDRI_boot)) %>%
  ungroup() %>%
  select(subtype, NNRI, NNRI_sd, NNRI_sim, NNRI_sim_sd, NNRI_boot, NNRI_boot_sd,
         VDRI, VDRI_sd, VDRI_sim, VDRI_sim_sd, VDRI_boot, VDRI_boot_sd)

(NNRI_plot <- df_obs %>%
    select(subtype, NNRI, NNRI_sim, NNRI_boot) %>%
  gather(key = 'Data', val = 'NNRI', -subtype) %>%
  mutate(Data = case_when(Data == "NNRI" ~ 'Observed',
                          Data == "NNRI_boot" ~ 'Bootstrapped',
                          Data == "NNRI_sim" ~ 'Simulated CSR')) %>%
  ggplot(aes(y=subtype, x = NNRI, color = Data)) +
  geom_boxplot(outliers=F) +
  ggtitle(paste('RGC Nearest Neighbor Regularity Index in', N_s_thresh, 'random', d, 'um diameter regions')) +
  theme_tufte())

(VDRI_plot <- df_obs %>%
    select(subtype, VDRI, VDRI_sim, VDRI_boot) %>%
    gather(key = 'Data', val = 'VDRI', -subtype) %>%
    mutate(Data = case_when(Data == "VDRI" ~ 'Observed',
                            Data == "VDRI_boot" ~ 'Bootstrapped',
                            Data == "VDRI_sim" ~ 'Simulated CSR')) %>%
    ggplot(aes(y=subtype, x = VDRI, color = Data)) +
    geom_boxplot(outliers=F) +
    ggtitle(paste('RGC Voronoi Domain Regularity Index in', N_s_thresh, 'random', d, 'um diameter regions')) +
    theme_tufte())


(NNRI_plot+ theme(legend.position="none") )+VDRI_plot

print(length(unique(df_obs$subtype)))

################################################################################




###################################################################################

updated_rgc_df %>%
  group_by(Prediction) %>%
  summarise(n=n()) %>%
  ggplot(aes(y=Prediction, x = n)) +
  geom_col()+
  theme_tufte()


df_i <- updated_rgc_df %>%
  filter(retina==1 & slice==10 & between(y, 6075, 6275) & between(x, 1200, 1400) ) 

# write.csv(df_i,'/home/sam/Downloads/R2.12_y5100TO5300_x7550TO7750.csv', row.names = FALSE)

df_i %>%
  # filter(Prediction == '01_W3D1.1') %>%
  ggplot(aes(x=x,y=y, color = Prediction)) +
  geom_point()

d = 200
sample_thresh = 5
output <- list()
j = 1
for (s in unique(rgc_df$Prediction)) {
  N <- nrow(filter(df_i, Prediction == s))
  if (N > sample_thresh){
    NNRI_obs <- NNRI(filter(df_i, Prediction==s))
    VDRI_obs <- VDRI(filter(df_i, Prediction==s))
    RIs_csr <- RI_simulator(df_i, subtype=s, diameter=d)
    RIs_boot <- RI_bootstrapper(df_i, N)
    output[[j]] <- list()
    output[[j]][['subtype']] <- s
    output[[j]][['N']] <- N
    output[[j]][['NNRI_obs']] <- NNRI_obs
    output[[j]][['VDRI_obs']] <- VDRI_obs
    output[[j]][['NNRIs_csr']] <- RIs_csr[['NNRI']]
    output[[j]][['NNRI_csr']] <- mean(RIs_csr[['NNRI']])
    output[[j]][['NNRIs_boot']] <- RIs_boot[['NNRI']]
    output[[j]][['NNRI_boot']] <- mean(RIs_boot[['NNRI']])
    
    output[[j]][['VDRIs_csr']] <- RIs_csr[['VDRI']]
    output[[j]][['VDRI_csr']] <- mean(RIs_csr[['VDRI']])
    output[[j]][['VDRIs_boot']] <- RIs_boot[['VDRI']]
    output[[j]][['VDRI_boot']] <- mean(RIs_boot[['VDRI']])
    j = j+1
    print(paste("For sample", N, "cells observed", s, "NNRI =", round(NNRI_obs, 2), "VDRI =", round(VDRI_obs, 2)))
  }
}

df_obs <- data.frame(
  subtype = sapply(output, function(x) x$subtype),
  N = sapply(output, function(x) x$N),
  NNRI = sapply(output, function(x) x$NNRI_obs),
  VDRI = sapply(output, function(x) x$VDRI_obs)
) 

# Add the simulation vectors as list columns
df_obs$NNRI_sim <- lapply(output, function(x) x$NNRIs_csr)
df_obs$VDRI_sim <- lapply(output, function(x) x$VDRIs_csr)
df_obs$NNRI_boot <- lapply(output, function(x) x$NNRIs_boot)
df_obs$VDRI_boot <- lapply(output, function(x) x$VDRIs_boot)

df_obs <- df_obs %>%
  unnest(cols = c(NNRI_sim, VDRI_sim,NNRI_boot, VDRI_boot))


NNRI_plot <- df_obs %>%
  select(subtype, NNRI, NNRI_sim, NNRI_boot) %>%
  gather(key = 'Data', val = 'NNRI', -subtype) %>%
  mutate(Data = case_when(Data == "NNRI" ~ 'Observed',
                          Data == "NNRI_boot" ~ 'Bootstrapped',
                          Data == "NNRI_sim" ~ 'Simulated CSR')) %>%
  ggplot(aes(y=subtype, x = NNRI, color = Data)) +
  geom_boxplot(outliers=F) +
  theme_tufte()

VDRI_plot <- df_obs %>%
    select(subtype, VDRI, VDRI_sim, VDRI_boot) %>%
    gather(key = 'Data', val = 'VDRI', -subtype) %>%
    mutate(Data = case_when(Data == "VDRI" ~ 'Observed',
                            Data == "VDRI_boot" ~ 'Bootstrapped',
                            Data == "VDRI_sim" ~ 'Simulated CSR')) %>%
    ggplot(aes(y=subtype, x = VDRI, color = Data)) +
    geom_boxplot(outliers=F) +
    theme_tufte()
  
  
(NNRI_plot+ theme(legend.position="none") )+VDRI_plot


updated_rgc_df %>%
  ggplot(aes(y = Prediction, x = volume)) +
  geom_violin() +
  theme_tufte()


updated_rgc_df %>%
  ggplot(aes(y = Prediction, x = nn_dist)) +
  geom_violin() +
  scale_x_log10() +
  theme_tufte()
