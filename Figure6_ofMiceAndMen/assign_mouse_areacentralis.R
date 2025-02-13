library(tidyverse)

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
                            


ggplot(mouse, aes(x,y, color = areacentralis)) +
  geom_point(alpha=0.5)

write.csv(mouse, file = volcano_path)
