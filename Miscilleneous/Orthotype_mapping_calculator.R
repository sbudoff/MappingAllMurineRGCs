ROOT_DIR <- '/home/sam/Primate_Comparison/'
orthotype_counts <- read.csv( paste0(ROOT_DIR, 'Mouse_Monkey_Orthotype.csv')) %>%
  mutate(Species = factor(Species)) %>%
  group_by(Species, Location, Datasource) %>%
  mutate(Subtype = factor(Subtype)) %>%
  ungroup()

library(tidyverse)
library(reshape2)

# Set working directory
ROOT_DIR <- '/home/sam/Primate_Comparison/'

# Custom function to extract numeric value from subtype names
extract_numeric <- function(x) {
  # Handle special cases first
  if (grepl("^(MG|PG)_(ON|OFF)", x)) {
    # Order: MG_ON, MG_OFF, PG_ON, PG_OFF
    ordering <- c("MG_OFF" = 1, "MG_ON" = 2, "PG_OFF" = 3, "PG_ON" = 4,
                  "PG_OFF_a" = 3.1, "PG_OFF_b" = 3.2)
    base_value <- ordering[sub("_[ab]$", "", x)]
    if (grepl("_[ab]$", x)) {
      return(base_value + ifelse(grepl("_a$", x), 0.1, 0.2))
    }
    return(base_value)
  }
  
  # Handle merged clusters (e.g., CRGC12.13)
  if (grepl("\\d+\\.\\d+$", x)) {
    nums <- as.numeric(strsplit(gsub("[^0-9.]", "", x), "\\.")[[1]])
    return(min(nums))
  }
  
  # Extract numeric part for regular cases
  num <- as.numeric(gsub("[^0-9]", "", x))
  if (is.na(num)) return(999)
  
  # Handle a/b variants
  if (grepl("[ab]$", x)) {
    num <- num + ifelse(grepl("a$", x), 0.1, 0.2)
  }
  
  return(num)
}

# Load and preprocess data
orthotype_counts <- read.csv(paste0(ROOT_DIR, 'Mouse_Monkey_Orthotype.csv')) %>%
  mutate(Species = factor(Species),
         Location = if_else(is.na(Location) | Location == "", "Unknown", Location)) %>%
  group_by(Species, Location) %>%
  mutate(Subtype = factor(Subtype)) %>%
  ungroup()

# Add numeric ordering value for subtypes
orthotype_counts <- orthotype_counts %>%
  mutate(subtype_order = sapply(as.character(Subtype), extract_numeric))

# Identify orthotype columns
orthotype_cols <- grep("^O[0-9]+$", names(orthotype_counts), value = TRUE)

# Convert to long format first, then calculate percentages
melted_data <- orthotype_counts %>%
  select(Species, Location, Subtype, Datasource, subtype_order, all_of(orthotype_cols)) %>%
  pivot_longer(cols = all_of(orthotype_cols),
               names_to = "Orthotype",
               values_to = "Count") %>%
  group_by(Species, Location, Subtype, Datasource) %>%
  mutate(Percentage = (Count / sum(Count)) * 100) %>%
  ungroup() %>%
  mutate(
    # Create reversed orthotype factor for proper top-to-bottom ordering
    Orthotype = factor(Orthotype, 
                       levels = rev(paste0("O", 1:21)),
                       labels = rev(paste0("oRGC", 1:21))),
    # Create combined subtype_datasource label
    SubtypeLabel = paste(Subtype, Datasource, sep="_")
  )

# Verify that each group sums to 100%
sum_check <- melted_data %>%
  group_by(Species, Location, Subtype, Datasource) %>%
  summarise(Total = sum(Percentage), .groups = 'drop') %>%
  filter(abs(Total - 100) > 0.01)  # Check for any that don't sum to 100 (allowing for small floating point differences)

if (nrow(sum_check) > 0) {
  warning("Some groups do not sum to 100%")
  print(sum_check)
}

# Create ordered factor for Species-Location combinations
location_order <- c("Fovea", "Periphery", "Unknown")
species_order <- c("Human", "Macaque", "Marmoset", "Mouse")

melted_data <- melted_data %>%
  filter(Species %in% species_order) %>%
  mutate(
    Location = factor(Location, levels = location_order),
    Species = factor(Species, levels = species_order),
    Panel = factor(paste(Species, Location),
                   levels = c(
                     paste("Human", location_order),
                     paste("Macaque", location_order),
                     paste("Marmoset", location_order),
                     "Mouse Unknown"
                   ))
  )

# Create the plot with ordered subtypes within each panel
p <- ggplot(melted_data, aes(x = reorder(SubtypeLabel, subtype_order), 
                             y = Orthotype, 
                             fill = Percentage)) +
  geom_tile() +
  scale_fill_gradient2(low = "white",
                       mid = "#FB8072",
                       high = "#67001F",
                       midpoint = 50,
                       limits = c(0, 100),
                       name = "Mapping (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 10),
        axis.title.x = element_blank()) +
  labs(y = "RGC orthotype") +
  facet_grid(. ~ Panel, scales = "free_x", space = "free")

# Display plot
print(p)

# Print some example percentages to verify normalization
example_check <- melted_data %>%
  group_by(Species, Location, Subtype, Datasource) %>%
  summarise(
    Total = sum(Percentage),
    N_Orthotypes = n(),
    .groups = 'drop'
  ) %>%
  head(10)

print("Example group sums:")
print(example_check)




#############################################################################################
library(tidyverse)
library(reshape2)

# Set working directory
ROOT_DIR <- '/home/sam/Primate_Comparison/'

# Custom function to extract numeric value from subtype names
extract_numeric <- function(x) {
  # Handle special cases first
  if (grepl("^(MG|PG)_(ON|OFF)", x)) {
    ordering <- c("MG_ON" = 1, "MG_OFF" = 2, "PG_ON" = 3, "PG_OFF" = 4,
                  "PG_OFF_a" = 4.1, "PG_OFF_b" = 4.2)
    base_value <- ordering[sub("_[ab]$", "", x)]
    if (grepl("_[ab]$", x)) {
      return(base_value + ifelse(grepl("_a$", x), 0.1, 0.2))
    }
    return(base_value)
  }
  
  # Handle merged clusters (e.g., CRGC12.13)
  if (grepl("\\d+\\.\\d+$", x)) {
    nums <- as.numeric(strsplit(gsub("[^0-9.]", "", x), "\\.")[[1]])
    return(min(nums))
  }
  
  # Extract numeric part for regular cases
  num <- as.numeric(gsub("[^0-9]", "", x))
  if (is.na(num)) return(999)
  
  # Handle a/b variants
  if (grepl("[ab]$", x)) {
    num <- num + ifelse(grepl("a$", x), 0.1, 0.2)
  }
  
  return(num)
}

# Load and preprocess data
orthotype_counts <- read.csv(paste0(ROOT_DIR, 'Mouse_Monkey_Orthotype.csv')) %>%
  mutate(Species = factor(Species),
         Location = if_else(is.na(Location) | Location == "", "Unknown", Location)) %>%
  group_by(Species, Location) %>%
  mutate(Subtype = factor(Subtype)) %>%
  ungroup()

# Add numeric ordering value for subtypes
orthotype_counts <- orthotype_counts %>%
  mutate(subtype_order = sapply(as.character(Subtype), extract_numeric))

# Identify orthotype columns
orthotype_cols <- grep("^O[0-9]+$", names(orthotype_counts), value = TRUE)

# Convert to long format and calculate percentages
melted_data <- orthotype_counts %>%
  select(Species, Location, Subtype, Datasource, subtype_order, all_of(orthotype_cols)) %>%
  pivot_longer(cols = all_of(orthotype_cols),
               names_to = "Orthotype",
               values_to = "Count") %>%
  group_by(Species, Location, Subtype, Datasource) %>%
  mutate(Percentage = (Count / sum(Count)) * 100) %>%
  ungroup()

# Find the most likely orthotype for each subtype within each Species-Location
max_orthotypes <- melted_data %>%
  group_by(Species, Location, Subtype, Datasource) %>%
  slice_max(order_by = Percentage, n = 1) %>%
  mutate(
    orthotype_num = as.numeric(gsub("O", "", Orthotype)),
    # Create composite sort value: orthotype_num * 1000 + subtype_order
    # This ensures primary sort by orthotype, secondary by subtype number
    composite_order = orthotype_num * 1000 + subtype_order
  ) %>%
  ungroup()

# Join back to get the composite order for all rows
melted_data <- melted_data %>%
  left_join(
    max_orthotypes %>% select(Species, Location, Subtype, Datasource, composite_order),
    by = c("Species", "Location", "Subtype", "Datasource")
  ) %>%
  mutate(
    # Create reversed orthotype factor for proper top-to-bottom ordering
    Orthotype = factor(Orthotype, 
                       levels = rev(paste0("O", 1:21)),
                       labels = rev(paste0("oRGC", 1:21))),
    # Create combined subtype_datasource label
    SubtypeLabel = paste(Subtype, Datasource, sep="_")
  )

# Create ordered factor for Species-Location combinations
location_order <- c("Fovea", "Periphery", "Unknown")
species_order <- c("Human", "Macaque", "Marmoset", "Mouse")

melted_data <- melted_data %>%
  filter(Species %in% species_order) %>%
  mutate(
    Location = factor(Location, levels = location_order),
    Species = factor(Species, levels = species_order),
    Panel = factor(paste(Species, Location),
                   levels = c(
                     paste("Human", location_order),
                     paste("Macaque", location_order),
                     paste("Marmoset", location_order),
                     "Mouse Unknown"
                   ))
  )

# Create the plot with new ordering
p <- ggplot(melted_data, aes(x = reorder(SubtypeLabel, composite_order), 
                             y = Orthotype, 
                             fill = Percentage)) +
  geom_tile() +
  scale_fill_gradient2(low = "white",
                       mid = "#FB8072",
                       high = "#67001F",
                       midpoint = 50,
                       limits = c(0, 100),
                       name = "Mapping (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 10),
        axis.title.x = element_blank()) +
  labs(y = "RGC orthotype") +
  facet_grid(. ~ Panel, scales = "free_x", space = "free")

# Display plot
print(p)

# Print some example mappings to verify the ordering
example_mappings <- max_orthotypes %>%
  arrange(Species, Location, composite_order) %>%
  select(Species, Location, Subtype, Datasource, 
         Most_Likely_Orthotype = Orthotype, 
         Max_Percentage = Percentage) %>%
  head(20)

print("Example mappings and their ordering:")
print(example_mappings)


##################################################################################3
library(tidyverse)
library(reshape2)

# Set working directory
ROOT_DIR <- '/home/sam/Primate_Comparison/'

# Custom function to extract numeric value from subtype names
extract_numeric <- function(x) {
  # Handle special cases first
  if (grepl("^(MG|PG)_(ON|OFF)", x)) {
    ordering <- c("MG_ON" = 1, "MG_OFF" = 2, "PG_OFF" = 3, "PG_OFN" = 4,
                  "PG_OFF_a" = 3.1, "PG_OFF_b" = 3.2)
    base_value <- ordering[sub("_[ab]$", "", x)]
    if (grepl("_[ab]$", x)) {
      return(base_value + ifelse(grepl("_a$", x), 0.1, 0.2))
    }
    return(base_value)
  }
  
  # Handle merged clusters (e.g., CRGC12.13)
  if (grepl("\\d+\\.\\d+$", x)) {
    nums <- as.numeric(strsplit(gsub("[^0-9.]", "", x), "\\.")[[1]])
    return(min(nums))
  }
  
  # Extract numeric part for regular cases
  num <- as.numeric(gsub("[^0-9]", "", x))
  if (is.na(num)) return(999)
  
  # Handle a/b variants
  if (grepl("[ab]$", x)) {
    num <- num + ifelse(grepl("a$", x), 0.1, 0.2)
  }
  
  return(num)
}

# Load and preprocess data
orthotype_counts <- read.csv(paste0(ROOT_DIR, 'Mouse_Monkey_Orthotype.csv')) %>%
  mutate(Species = factor(Species),
         Location = if_else(is.na(Location) | Location == "", "Unknown", Location)) %>%
  group_by(Species, Location) %>%
  mutate(Subtype = factor(Subtype)) %>%
  ungroup()

# Add numeric ordering value for subtypes
orthotype_counts <- orthotype_counts %>%
  mutate(subtype_order = sapply(as.character(Subtype), extract_numeric))

# Identify orthotype columns
orthotype_cols <- grep("^O[0-9]+$", names(orthotype_counts), value = TRUE)

# Sum across datasets and convert to percentages
merged_data <- orthotype_counts %>%
  group_by(Species, Location, Subtype) %>%
  summarise(across(all_of(orthotype_cols), sum),
            subtype_order = first(subtype_order),
            .groups = 'drop') %>%
  pivot_longer(cols = all_of(orthotype_cols),
               names_to = "Orthotype",
               values_to = "Count") %>%
  group_by(Species, Location, Subtype) %>%
  mutate(Percentage = (Count / sum(Count)) * 100) %>%
  ungroup()

# Find the most likely orthotype for each subtype within each Species-Location
max_orthotypes <- merged_data %>%
  group_by(Species, Location, Subtype) %>%
  slice_max(order_by = Percentage, n = 1) %>%
  mutate(
    orthotype_num = as.numeric(gsub("O", "", Orthotype)),
    # Create composite sort value: orthotype_num * 1000 + subtype_order
    composite_order = orthotype_num * 1000 + subtype_order
  ) %>%
  ungroup()

# Join back to get the composite order for all rows
merged_data <- merged_data %>%
  left_join(
    max_orthotypes %>% select(Species, Location, Subtype, composite_order),
    by = c("Species", "Location", "Subtype")
  ) %>%
  mutate(
    # Create reversed orthotype factor for proper top-to-bottom ordering
    Orthotype = factor(Orthotype, 
                       levels = rev(paste0("O", 1:21)),
                       labels = rev(paste0("oRGC", 1:21)))
  )

# Create ordered factor for Species-Location combinations
location_order <- c("Fovea", "Periphery", "Unknown")
species_order <- c("Human", "Macaque", "Marmoset", "Mouse")

merged_data <- merged_data %>%
  filter(Species %in% species_order) %>%
  mutate(
    Location = factor(Location, levels = location_order),
    Species = factor(Species, levels = species_order),
    Panel = factor(paste(Species, Location),
                   levels = c(
                     paste("Human", location_order),
                     paste("Macaque", location_order),
                     paste("Marmoset", location_order),
                     "Mouse Unknown"
                   ))
  )

# Create the plot with new ordering
p <- ggplot(merged_data, aes(x = reorder(Subtype, composite_order), 
                             y = Orthotype, 
                             fill = Percentage)) +
  geom_tile() +
  scale_fill_gradient2(low = "white",
                       mid = "#FB8072",
                       high = "#67001F",
                       midpoint = 50,
                       limits = c(0, 100),
                       name = "Mapping (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 10),
        axis.title.x = element_blank()) +
  labs(y = "RGC orthotype") +
  facet_grid(. ~ Panel, scales = "free_x", space = "free")

# Display plot
print(p)

# Print some example mappings to verify the ordering
example_mappings <- max_orthotypes %>%
  arrange(Species, Location, composite_order) %>%
  select(Species, Location, Subtype, 
         Most_Likely_Orthotype = Orthotype, 
         Max_Percentage = Percentage) %>%
  head(20)

print("Example mappings and their ordering:")
print(example_mappings)

write_csv(max_orthotypes, file = paste0(ROOT_DIR, 'Orthotype_mapping.csv'))

mouse_orthos <- max_orthotypes %>%
  filter(Species=='Mouse') %>%
  select(orthotype_num,Subtype)
