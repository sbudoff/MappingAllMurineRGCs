library(tidyverse)
library(dplyr)

root <- '/media/sam/Data2/baysor_rbpms_consolidated'

# Load the gene order for cuttlenet inference
gene_order <- read.csv('/media/sam/Data2/baysor_analysis/CuttleNet/input_order.csv', header = F) %>% pull()

# Find all processed segmentations
slides <- list.files(root)
all_paths <- c()
paths <- c()
for (slide in slides) {
  path <- paste0(root,"/",slide)
  for (p in list.files(path)){
    p_i <- paste0(root,"/",slide,"/",p)
    all_paths <- c(all_paths, p_i)
    for (file in list.files(p_i)){
      if (file == 'segmentation.csv') {
        
        #sample_root <- substr(path,1, nchar(path)-16)
        #sample_id <- substr(path,50, nchar(path)-24)
        #out_path <- paste0(sample_root, sample_id, 'expression_matrix.csv')
        
        # Construct the expected output path
        sample_root <- substr(p_i,1, nchar(p_i))
        sample_id <- substr(p_i,50, nchar(p_i)-7)
        expected_output <- paste0(sample_root, '/', sample_id, 'expression_matrix.csv')
        # Only add to paths if output doesn't exist
        if (!file.exists(expected_output)) {
          paths <- c(paths, paste0(p_i, '/segmentation.csv'))
          print(paste("segmentation found for", p, "- needs processing"))
        } else {
          print(paste("skipping", p, "- already processed"))
        }
      }
    }
  }
}


add_missing_columns <- function(df, required_cols) {
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    df[missing_cols] <- 0
  }
  return(df)
}

for (path in paths) {
  sample_root <- substr(path,1, nchar(path)-16)
  sample_id <- substr(path,50, nchar(path)-24)
  out_path <- paste0(sample_root, sample_id, 'expression_matrix.csv')
  
  print(paste("Working on", sample_id))
  
  df <- read.csv(path)
  
  df <- df %>%
    filter(is_noise == "false") %>% # Remove all genes predicted by Baysor to be noise
    filter(qv >= 20) %>% # Remove all poor quality transcripts
    filter(gene %in% gene_order) %>% # Remove all QC and unassigned probes
    dplyr::select(cell, x, y, z, gene) %>% # subset the df for calculations
    group_by(cell) %>%
    filter(n() >= 100) %>% # remove all "cells" with fewer than 100 transcripts
    ungroup()
  
  df_geom <- df %>%
    group_by(cell) %>%
    summarize( # Calculate geometric properties
      volume = (max(x) - min(x)) * (max(y) - min(y)) * (max(z) - min(z)),
      pca = list(prcomp(cbind(x, y, z))$sdev),
      elongation = map_dbl(pca, ~.x[1]/.x[2]),  # ratio of first two principal components
      flatness = map_dbl(pca, ~.x[2]/.x[3]),    # ratio of second two principal components
      compactness = map_dbl(pca, ~.x[3]/.x[1]),  # ratio of largest to smallest axis
      sphericity = map_dbl(pca, ~min(.x)/max(.x)),
      x = mean(x),
      y = mean(y),
      z = mean(z)
    ) %>%
    select(-pca)  %>%
    ungroup()
  
  df_genes <- df %>%
    group_by(cell) %>%
    count(cell, gene) %>% # Cont expression matrix
    pivot_wider(
      names_from = gene,
      values_from = n,
      values_fill = 0
    ) %>%
    ungroup() %>%
    add_missing_columns(gene_order)
  
  
  df_clean <- left_join(df_geom, df_genes,
      by = "cell"
    ) %>%
    dplyr::select(cell, x, y, z, volume, elongation, flatness, sphericity, compactness,all_of(gene_order)) # fix order of columns
  write.csv(df_clean,out_path, row.names = FALSE)
  print(paste(" Expression matrix computed and stored at", out_path))
}

print('All available expression matrices produced')
