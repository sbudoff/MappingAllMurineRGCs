library(sf)
library(dplyr)


root <- '/media/sam/New Volume/Xenium_Data/IHC_experiments/Alignments'

# Search through root for subdirectory names
slides <- list.dirs(root, full.names = FALSE, recursive = FALSE)

# Initialize an empty dataframe to store all the data
polygons_df <- data.frame()

group <- 0
# Search through each slide for csv files
for (slide in slides) {
  path <- file.path(root, slide)
  # Find csv files in the subdirectory
  csvs <- list.files(path = path, pattern = "\\.csv$", full.names = TRUE)
  
  for (csv in csvs) {
    # Read the current CSV
    temp_df <- read.csv(csv, 
                        skip = 2,
                        col.names = c("Selection", "X", "Y"))
    
    # Add the slide ID column
    temp_df$slide_id <- slide
    temp_df$group <- group
    group = group + 1
    # Combine with the main dataframe
    polygons_df <- bind_rows(polygons_df, temp_df)
  }
}

polygons_df <- polygons_df %>%
  mutate(rectangle = paste(slide_id, group, Selection),
         rectangle = factor(rectangle),
         slide = factor(slide_id)) %>%
  select(-group, -Selection, -slide_id) %>%
  unique()


# Convert rectangles to sf polygons
rect_to_polygon <- function(rect_points) {
  # Add first point to close the polygon
  closed_points <- rbind(rect_points, rect_points[1,])
  # Create polygon from points
  st_polygon(list(as.matrix(closed_points[, c("X", "Y")])))
}

# Process each slide
process_slide <- function(slide_data, buffer_dist = 5) {
  # Convert rectangles to sf polygons
  polygons_list <- slide_data %>%
    group_by(rectangle) %>%
    group_split() %>%
    lapply(function(rect) rect_to_polygon(rect))
  
  # Convert to sf object
  polygons_sf <- st_sfc(polygons_list) %>%
    st_set_crs(NA)  # No CRS needed for this data
  
  # Buffer out
  buffered <- st_buffer(polygons_sf, dist = buffer_dist)
  
  # Union the buffered polygons
  united <- st_union(buffered)
  
  # Buffer back in
  final <- st_buffer(united, dist = -buffer_dist)
  
  return(final)
}

# Apply to all slides
final_polygons <- polygons_df %>%
  group_by(slide) %>%
  group_map(~process_slide(.x, buffer_dist = 5)) %>%
  setNames(unique(polygons_df$slide))



library(sf)
library(dplyr)
library(ggplot2)

# First create the original rectangles plot function
plot_original_rectangles <- function(slide_data, slide_name) {
  ggplot(slide_data) +
    geom_point(aes(x = X, y = Y, group = rectangle), color = "blue", size = 0.5) +
    geom_path(aes(x = X, y = Y, group = rectangle), color = "blue", alpha = 0.5) +
    coord_equal() +
    theme_minimal() +
    ggtitle(paste("Original Rectangles -", slide_name)) +
    theme(plot.title = element_text(hjust = 0.5))
}

# Then plot the merged polygons
plot_merged_polygons <- function(merged_polygon, slide_name) {
  # Convert to sf data frame for ggplot
  polygon_df <- st_as_sf(data.frame(geometry = merged_polygon))
  
  ggplot() +
    geom_sf(data = polygon_df, fill = "red", alpha = 0.3, color = "red") +
    theme_minimal() +
    ggtitle(paste("Merged Polygons -", slide_name)) +
    theme(plot.title = element_text(hjust = 0.5))
}

# Create comparison plots for each slide
plot_comparison <- function(slide_data, merged_result, slide_name) {
  original <- plot_original_rectangles(slide_data, slide_name)
  merged <- plot_merged_polygons(merged_result, slide_name)
  
  # Arrange plots side by side
  gridExtra::grid.arrange(original, merged, ncol = 2)
}

# Generate plots for each slide
for (current_slide in unique(polygons_df$slide)) {
  slide_data <- polygons_df %>% filter(slide == current_slide)
  plot_comparison(slide_data, 
                  final_polygons[[current_slide]], 
                  current_slide)
}

# Function to extract coordinates from a polygon and create a dataframe
polygon_to_df <- function(polygon) {
  # Extract coordinates from the geometry
  coords <- st_coordinates(polygon)
  
  # Create dataframe with X, Y coordinates
  data.frame(
    X = coords[, "X"],
    Y = coords[, "Y"],
    # L1 column from st_coordinates represents different polygons within the multipolygon
    L1 = coords[, "L1"],
    L2 = coords[, "L2"],
    L3 = coords[, "L3"]
  ) %>%
    mutate(polygon_id = paste(L1, L2, L3)) %>%
    select(-L1, -L2, -L3)
}

# Export each slide's polygons to a separate CSV
for (slide in names(final_polygons)) {
  # Convert the multipolygon to dataframe
  polygon_df <- polygon_to_df(final_polygons[[slide]])
  
  # Create filename
  filename <- file.path(root, paste0("merged_RBPMS_rectangles_", slide, ".csv"))
  
  # Write to CSV
  write.csv(polygon_df, filename, row.names = FALSE)
}
