library(tidyverse)
library(patchwork)

DATA_DIR <- '/home/sam/MappingAllMurineRGCs/Data/'
ROOT_DIR <- paste0(DATA_DIR, 'Figure6_Outputs/')

clusters <- read.csv(paste0(ROOT_DIR,'hc_consensus.csv')) %>%
  dplyr::select(x,y,consensus) %>%
  rename(cluster = consensus) %>%
  mutate(x=x, y =y)

volcano_path <- paste0(DATA_DIR, "Figure5_Outputs/volcan_expression_matrix.csv")
mouse <-read.csv(volcano_path) %>%
  dplyr::select(-X)


# We can use the closest grid point to assign clusters
# Add areacentralis column based on nearest neighbor in clusters grid
# Scale mouse coordinates to match clusters grid
mouse <- mouse %>%
  mutate(
    # Scale x from [-1,1] to [1,12]
    x_scaled = (x - min(x)) / (max(x) - min(x)) * (max(clusters$x) - min(clusters$x)) + min(clusters$x),
    # Scale y from [-1,1] to [1,11]
    y_scaled = (y - min(y)) / (max(y) - min(y)) * (max(clusters$y) - min(clusters$y)) + min(clusters$y)
  ) %>%
  rowwise() %>%
  mutate(areacentralis = {
    # Find closest grid point using scaled coordinates
    dists <- sqrt((clusters$x - x_scaled)^2 + (clusters$y - y_scaled)^2)
    closest_cluster <- clusters$cluster[which.min(dists)]
    # Convert to boolean 0/1
    as.numeric(closest_cluster == 1)
  }) %>%
  dplyr::select(-x_scaled, -y_scaled)  # Remove temporary scaling columns


# Define the columns to plot
bool_cols <- c("ipsi", "binocular", "peripheral", 
               "visual_sky", "visual_ground",  "visual_floor", 
               "binocular_sky", "binocular_ground", "peripheral_sky",
               "peripheral_ground", "peripheral_floor", "binocular_floor" ,
               "areacentralis", "w3"
               )

# Set the letter font size as a parameter for easy adjustment
letter_size <- 14  # Adjust this value to change letter sizes

# Create individual plots
plot_list <- map(seq_along(bool_cols), function(i) {
  col <- bool_cols[i]
  letter <- letters[i]
  
  ggplot(mouse, aes(x = x, y = y)) +
    geom_point(aes(fill = as.factor(!!sym(col))), 
               shape = 21,
               color = "black",
               size = 0.5,
               alpha = 0.4) +
    scale_fill_manual(values = c("white", "black")) +
    theme_void() +
    theme(legend.position = "none",
          plot.title = element_text(size = letter_size, hjust = 0, vjust = -1)) +
    labs(title = letter)
})

# Combine plots in a grid
# Adjust the number of columns as needed (e.g., 4 for a 4x4 grid)
(combined_plot <- wrap_plots(plot_list, ncol = 3))

ggsave(paste0(ROOT_DIR, "S12_visualSceneMasks.png"), 
       combined_plot,
       width = 8,
       height = 12, 
       dpi = 300) 
