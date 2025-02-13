library(tidyverse)
library(patchwork)
library(ggpubr)
library(ggrepel)
library(dplyr)


DATA_DIR <- '/home/sam/MappingAllMurineRGCs/Data/'
ROOT_DIR <- paste0(DATA_DIR, 'Figure6_Outputs_supW3/')
orthotype_key <- read.csv(paste0(ROOT_DIR, 'Orthotype_mapping.csv')) %>%
  dplyr::select(Species, Location, Subtype, Orthotype)

mouse_Locations <- read.csv(paste0(DATA_DIR,'Figure4_Outputs/cluster_assignments.csv') )%>%
  mutate(Species = "Mouse") %>%
  rename(Subtype = Label,
         Location = Cluster) %>%
  dplyr::select(-Method)
mouse_ortho2Location <- orthotype_key %>%
  filter(Species == 'Mouse') %>%
  dplyr::select(-Location) %>%
  left_join(mouse_Locations, by=c("Species", "Subtype"))
orthotype_key <- orthotype_key %>%
  filter(Species != 'Mouse') %>%
  rbind(mouse_ortho2Location)

write.csv(orthotype_key, file = paste0(ROOT_DIR, 'mouse_primate_orthotypes_by_Location.csv'))


ortholog_key <- read.csv( paste0(DATA_DIR, 'mouse2primate_xeniumOrthologs_verified.csv'))
  

human <- read.csv(paste0(ROOT_DIR, 'human_spatial_scRNAseq_xeniumSubset.csv'))
huamn_ortholog_mapping <- ortholog_key %>%
  dplyr::select(xenium_target, human_id) %>%
  drop_na() %>%
  dplyr::distinct()

# Identify which xenium_targets will have duplicates
duplicate_counts <- huamn_ortholog_mapping %>%
  group_by(xenium_target) %>%
  mutate(count = n()) %>%
  arrange(xenium_target) %>%
  mutate(suffix = if_else(count > 1, 
                          paste0(".", row_number()), 
                          ""))

# Get current column names
old_colnames <- colnames(human)

# Create new column names using a for loop
new_colnames <- old_colnames
for(i in seq_along(old_colnames)) {
  # Check if this column name exists in human_id
  mapping <- duplicate_counts %>%
    filter(human_id == old_colnames[i])
  
  # If we found a mapping, replace the name
  if(nrow(mapping) == 1) {
    # Combine xenium_target with suffix (if any)
    new_colnames[i] <- paste0(mapping$xenium_target, mapping$suffix)
  }
}

# Apply the new column names
colnames(human) <- new_colnames

# Print summary of changes and show any duplicates that were handled
duplicates_handled <- duplicate_counts %>%
  filter(count > 1) %>%
  dplyr::select(human_id, xenium_target, suffix) %>%
  arrange(xenium_target)

print("The following duplicates were handled with suffixes:")
print(duplicates_handled)

# Print total number of changes
n_changed <- sum(new_colnames != old_colnames)
print(sprintf("Changed %d column names", n_changed))

# Verify no duplicate names remain
if(any(duplicated(colnames(human)))) {
  warning("There are still duplicate column names!")
} else {
  print("No duplicate column names remain in the dataset")
}

human_ortho <- orthotype_key %>%
  filter(Species=="Human") %>%
  mutate(Subtype = case_when(Subtype == "MG_ON" ~ "HRGC1",
                             Subtype == "MG_OFF" ~ "HRGC2",
                             Subtype == "PG_OFF" ~ "HRGC3",
                             Subtype == "PG_ON" ~ "HRGC4",
                             T~Subtype),
         Subtype = if_else(Subtype == "HRGC4" & Orthotype == "O2",NA, Subtype)) %>% # this is a duplication mistake so i am making a judgment call from the plot
  filter(Location != "Unknown") %>%
  drop_na()%>%
  # Get all combinations of existing Subtypes with both Locations
  group_by(Subtype) %>%
  complete(Location = c("Fovea", "Periphery")) %>%
  # Fill in Orthotype within each Subtype group
  group_by(Subtype) %>%
  fill(Orthotype, .direction = "updown") %>%
  ungroup()

print(unique(human$organ))

human_FoveaMaculaPool <- human %>%
  rename(Subtype=rgc_cluster) %>% 
  mutate(Location = case_when(organ == "peripheral region of retina" ~ "Periphery",
                              TRUE ~ "Fovea")) %>%
  left_join(human_ortho, by = c('Subtype', 'Location'))

print(unique(human_FoveaMaculaPool$organ))

human_FoveaVSPeripheryMaculaPool <- human %>%
  rename(Subtype=rgc_cluster) %>% 
  mutate(Location = case_when(organ == "fovea centralis" ~ "Fovea",
                              TRUE ~ "Periphery")) %>%
  left_join(human_ortho, by = c('Subtype', 'Location'))

print(unique(human_FoveaVSPeripheryMaculaPool$organ))

human_MaculaPoolvsPeriphery <- human %>%
  filter(organ != "fovea centralis") %>%
  rename(Subtype=rgc_cluster) %>% 
  mutate(Location = case_when(organ == "peripheral region of retina" ~ "Periphery",
                              TRUE ~ "Fovea")) %>%
  left_join(human_ortho, by = c('Subtype', 'Location'))

print(unique(human_MaculaPoolvsPeriphery$organ))

human_MaculaPropervsPeriphery <- human %>%
  filter(organ != "fovea centralis", organ != "macula lutea") %>%
  rename(Subtype=rgc_cluster) %>% 
  mutate(Location = case_when(organ == "peripheral region of retina" ~ "Periphery",
                              TRUE ~ "Fovea")) %>%
  left_join(human_ortho, by = c('Subtype', 'Location'))

print(unique(human_MaculaPropervsPeriphery$organ))


human_FoveavsPeriphery <- human %>%
  filter(organ != "macula lutea proper", organ != "macula lutea") %>%
  rename(Subtype=rgc_cluster) %>% 
  mutate(Location = case_when(organ == "peripheral region of retina" ~ "Periphery",
                              TRUE ~ "Fovea")) %>%
  left_join(human_ortho, by = c('Subtype', 'Location'))


print(unique(human_FoveavsPeriphery$organ))

  





macaque_ortho <- orthotype_key %>%
  filter(Species=="Macaque") %>%
  mutate(Subtype = case_when(Subtype == "MG_ON" ~ "MRGC1",
                             Subtype == "MG_OFF" ~ "MRGC2",
                             Subtype == "PG_OFF" ~ "MRGC3",
                             Subtype == "PG_ON" ~ "MRGC4",
                             T~Subtype)         ) %>% 
  filter(Location != "Unknown") %>%
  drop_na()%>%
  # Get all combinations of existing Subtypes with both Locations
  group_by(Subtype) %>%
  complete(Location = c("Fovea", "Periphery")) %>%
  # Fill in Orthotype within each Subtype group
  group_by(Subtype) %>%
  fill(Orthotype, .direction = "updown") %>%
  ungroup()


macaque <- read.csv(paste0(ROOT_DIR, 'macaque_spatial_scRNAseq_xeniumSubset.csv'))%>%
  mutate(
    Subtype = case_when(Subtype == "MG_ON" ~ "MRGC1",
                        Subtype == "MG_OFF" ~ "MRGC2",
                        Subtype == "PG_OFF" ~ "MRGC3",
                        Subtype == "PG_ON" ~ "MRGC4",
                        startsWith(Subtype, "pRGC") ~ paste0("M",substr(Subtype, 2, nchar(Subtype))),
                        startsWith(Subtype, "fRGC") ~ paste0("M",substr(Subtype, 2, nchar(Subtype))),
                        TRUE ~ Subtype),
    Subtype = factor(Subtype, 
                     levels = c(paste0("MRGC", sort(as.numeric(gsub("MRGC", "", 
                                                                    unique(Subtype[grepl("MRGC", Subtype)])))))),
                     ordered = TRUE),
    Location = case_when(Location == "peripheral" ~ "Periphery",
                         TRUE ~ "Fovea"),
    Location = factor(Location, levels = c("Fovea", "Periphery"))
  ) %>% 
  left_join(macaque_ortho, by = c('Subtype', 'Location')) %>%
  mutate(
    Orthotype = factor(Orthotype, 
                       levels = c(paste0("O", sort(as.numeric(gsub("O", "", 
                                                                   unique(Orthotype[grepl("O", Orthotype)])))))),
                       ordered = TRUE),
    Subtype = factor(Subtype, 
                     levels = c(paste0("MRGC", sort(as.numeric(gsub("MRGC", "", 
                                                                    unique(Subtype[grepl("MRGC", Subtype)])))))),
                     ordered = TRUE))

mouse_ortho <- orthotype_key %>%
  filter(Species=="Mouse") %>%
  filter(Location != "Unknown") %>%
  drop_na() %>%
  dplyr::select(-Location)
# mutate(Location = case_when(Location==4~ 'Fovea',
#                             T ~ 'Periphery'),
#        Location = factor(Location, levels = c("Fovea", "Periphery"))
#        )

mouse <-read.csv(paste0(DATA_DIR, "Figure5_Outputs/volcan_expression_matrix.csv"))%>%
  dplyr::select(-c(X, x, y, peripheral, visual_sky, visual_ground, visual_floor,     
                   binocular, binocular_sky,binocular_ground, peripheral_sky, peripheral_ground, 
                   peripheral_floor, binocular_floor, art, art_thin,horizon   )  ) %>%
  rename(Location = w3,
         Subtype = Prediction) %>%
  mutate(Location = case_when(Location == 1 ~ "Fovea",
                              Location == 0 ~ "Periphery"),
         Location = factor(Location, levels = c("Fovea", "Periphery")),
         Subtype = paste0("T", as.numeric(sub("_.*", "", Subtype))),
  ) %>%
  left_join(mouse_ortho, by = c('Subtype')) %>%
  mutate(
    Orthotype = factor(Orthotype, 
                       levels = c(paste0("O", sort(as.numeric(gsub("O", "", 
                                                                   unique(Orthotype[grepl("O", Orthotype)])))))),
                       ordered = TRUE),
    Subtype = factor(Subtype, 
                     levels = c(paste0("T", sort(as.numeric(gsub("T", "", 
                                                                 unique(Subtype[grepl("T", Subtype)])))))),
                     ordered = TRUE))

plot_gene_distribution_with_stats <- function(data_list, species_names, gene_name) {
  # Input validation
  if (length(data_list) != length(species_names)) {
    stop("Number of datasets must match number of species names")
  }
  
  # Helper function to check if gene exists in dataset
  check_gene_presence <- function(df, gene) {
    gene_cols <- grep(paste0("^", gene, "(\\.\\d+)?$"), names(df), value = TRUE)
    return(length(gene_cols) > 0)
  }
  
  # Check gene presence in each dataset and filter accordingly
  valid_indices <- which(sapply(data_list, check_gene_presence, gene = gene_name))
  
  if (length(valid_indices) == 0) {
    stop(paste("Gene", gene_name, "not found in any dataset"))
  }
  
  # Filter data_list and species_names to only include valid datasets
  data_list <- data_list[valid_indices]
  species_names <- species_names[valid_indices]
  
  # Notify about skipped species
  skipped_species <- setdiff(seq_along(species_names), valid_indices)
  if (length(skipped_species) > 0) {
    warning(paste("Gene", gene_name, "not found in dataset(s) for species:",
                  paste(species_names[skipped_species], collapse = ", "),
                  "\nProceeding with available data."))
  }
  
  # Helper function to handle duplicates and create simple dataset
  merge_duplicates <- function(df, gene, group_col = NULL) {
    gene_cols <- grep(paste0("^", gene, "(\\.\\d+)?$"), names(df), value = TRUE)
    
    if (length(gene_cols) > 1) {
      gene_data <- df %>%
        dplyr::select(Location, if(!is.null(group_col)) all_of(group_col), all_of(gene_cols)) %>%
        pivot_longer(
          cols = all_of(gene_cols),
          names_to = "gene_version",
          values_to = gene
        )
      return(gene_data %>% 
               dplyr::select(-gene_version) )#%>%
      # filter(.data[[gene]] != 0)) # This drops 0 valued observations #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    } else if (length(gene_cols) == 1) {
      return(df %>% 
               dplyr::select(Location, if(!is.null(group_col)) all_of(group_col), all_of(gene_cols)) %>%
               rename_with(~gene, all_of(gene_cols))  )# %>%
      # filter(.data[[gene]] != 0))  # This drops 0 valued observations #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }
  }
  
  # Process all datasets
  processed_data <- list()
  for (i in seq_along(data_list)) {
    species <- species_names[i]
    df <- data_list[[i]]
    
    processed_data[[species]] <- list(
      ortho = merge_duplicates(df, gene_name, "Orthotype"),
      sub = merge_duplicates(df, gene_name, "Subtype"),
      simple = merge_duplicates(df, gene_name)
    )
  }
  
  # Function to perform statistical analysis for each dataset
  perform_stats <- function(data, gene_name, group_col = NULL) {
    if (!is.null(group_col)) {
      groups <- unique(data[[group_col]])
      n_groups <- length(groups)
    } else {
      groups <- "Overall"
      n_groups <- 1
    }
    
    Locations <- unique(data$Location)
    results <- data.frame()
    
    if (is.null(group_col)) {
      # Simple t-test for overall comparison
      loc1_data <- data[[gene_name]][data$Location == Locations[1]]
      loc2_data <- data[[gene_name]][data$Location == Locations[2]]
      
      if(length(loc1_data) >= 2 && length(loc2_data) >= 2) {
        t_test_result <- t.test(loc1_data, loc2_data)
        
        results <- data.frame(
          Group = "Overall",
          group1 = Locations[1],
          group2 = Locations[2],
          statistic = t_test_result$statistic,
          p.value = t_test_result$p.value,
          p.adj = t_test_result$p.value,  # No correction needed for single test
          GroupType = "None"
        )
      }
    } else {
      # Group-wise analysis
      for(grp in groups) {
        subset_data <- data %>% filter(.data[[group_col]] == grp)
        
        if(nrow(subset_data) < 3) next
        
        loc1_data <- subset_data[[gene_name]][subset_data$Location == Locations[1]]
        loc2_data <- subset_data[[gene_name]][subset_data$Location == Locations[2]]
        
        if(length(loc1_data) >= 2 && length(loc2_data) >= 2) {
          t_test_result <- t.test(loc1_data, loc2_data)
          p_adj <- t_test_result$p.value * n_groups
          p_adj[is.na(p_adj)] <- 1
          
          results <- rbind(results, data.frame(
            Group = grp,
            group1 = Locations[1],
            group2 = Locations[2],
            statistic = t_test_result$statistic,
            p.value = t_test_result$p.value,
            p.adj = p_adj,
            GroupType = group_col
          ))
        }
      }
    }
    
    return(results)
  }
  
  # Perform statistical analysis for all datasets
  stats_data <- list()
  for (species in species_names) {
    stats_data[[species]] <- list(
      ortho = perform_stats(processed_data[[species]]$ortho, gene_name, "Orthotype"),
      sub = perform_stats(processed_data[[species]]$sub, gene_name, "Subtype"),
      simple = perform_stats(processed_data[[species]]$simple, gene_name)
    )
  }
  
  # Function to create plot with statistical annotations
  create_plot_with_stats <- function(data, stats, title, group_col = NULL) {
    y_max <- max(data[[gene_name]], na.rm = TRUE)
    y_range <- diff(range(data[[gene_name]], na.rm = TRUE))
    y_pos <- y_max + (y_range * 0.1)
    
    if (!is.null(group_col)) {
      p <- data %>%
        mutate(
          Group = factor(.data[[group_col]]),
          Location = factor(Location)
        ) %>%
        ggplot(aes(x = Group, y = .data[[gene_name]], color = Location)) +
        geom_boxplot(outlier.alpha = 0.1)
    } else {
      p <- data %>%
        mutate(Location = factor(Location)) %>%
        ggplot(aes(x = Location, y = .data[[gene_name]], color = Location)) +
        geom_boxplot(outlier.alpha = 0.1)
    }
    
    p <- p +
      ggtitle(title) +
      theme_minimal() +
      ylab(gene_name) +
      xlab(if(!is.null(group_col)) group_col else "Location")
    
    # Always remove legend - color scheme is shown in overall plots
    p <- p + theme(legend.position = "none")
    
    # Add slanted x-axis labels for Subtype plots
    if (!is.null(group_col) && group_col == "Subtype") {
      p <- p + theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
      )
    }
    
    if(nrow(stats) > 0) {
      for(i in 1:nrow(stats)) {
        if(stats$p.adj[i] <= 0.05) {
          sig_symbol <- ifelse(stats$p.adj[i] <= 0.001, "***",
                               ifelse(stats$p.adj[i] <= 0.01, "**", "*"))
          
          if (!is.null(group_col)) {
            x_pos <- which(levels(factor(data[[group_col]])) == stats$Group[i])
          } else {
            x_pos <- 1.5  # Center position between the two boxes
          }
          
          p <- p + annotate("text", 
                            x = x_pos,
                            y = y_pos,
                            label = sig_symbol,
                            color = "black")
        }
      }
    }
    
    p <- p + coord_cartesian(ylim = c(NA, y_max + y_range * 0.2))
    
    return(p)
  }
  
  # Create plots for each species
  plot_list <- list()
  for (species in species_names) {
    plot_list[[species]] <- list(
      overall = create_plot_with_stats(
        processed_data[[species]]$simple, 
        stats_data[[species]]$simple, 
        paste(species, "- Overall")
      ),
      subtype = create_plot_with_stats(
        processed_data[[species]]$sub, 
        stats_data[[species]]$sub, 
        paste(species, "- Subtype"),
        "Subtype"
      ),
      orthotype = create_plot_with_stats(
        processed_data[[species]]$ortho, 
        stats_data[[species]]$ortho, 
        paste(species, "- Orthotype"),
        "Orthotype"
      )
    )
  }
  
  # Create the combined plot based on number of species
  plot_rows <- list()
  for (species in species_names) {
    row <- plot_list[[species]]$overall + 
      plot_list[[species]]$subtype + 
      plot_list[[species]]$orthotype
    plot_rows[[species]] <- row
  }
  
  # Combine all rows
  combined_plot <- Reduce(`/`, plot_rows) +
    plot_annotation(
      title = paste("Gene Distribution:", gene_name),
      theme = theme_minimal()
    )
  
  # Create summary dataframe with all statistical results
  summary_df <- bind_rows(
    lapply(species_names, function(species) {
      bind_rows(
        stats_data[[species]]$simple %>% mutate(Species = species),
        stats_data[[species]]$sub %>% mutate(Species = species),
        stats_data[[species]]$ortho %>% mutate(Species = species)
      )
    })
  ) %>%
    mutate(gene = gene_name)
  
  # Print plot
  # print(combined_plot)
  
  # Return summary
  return(summary_df)
}






volcanoes <- read_csv(paste0(DATA_DIR, 'Figure5_Outputs/volcanoes/all_volcano_results.csv') )

genes <- unique(volcanoes$gene) 

gois <- volcanoes %>%
  filter(abs(log2FC) >=0.26 & adj_p_value <=0.00001 & 
           cell_type !='44_Novel' & mask == 'binocular_ground') %>%
  pull(gene) %>%
  unique() %>%
  sort()







  
species_data <- list(human_FoveaMaculaPool, 
                       human_FoveaVSPeripheryMaculaPool, 
                       human_MaculaPoolvsPeriphery, 
                       human_MaculaPropervsPeriphery, 
                       human_FoveavsPeriphery,
                     macaque,mouse)  
species_names <- c("human_FoveaMaculaPool", 
                   "human_FoveaVSPeripheryMaculaPool", 
                   "human_MaculaPoolvsPeriphery", 
                   "human_MaculaPropervsPeriphery", 
                   "human_FoveavsPeriphery",
                   "Macaque", "Mouse")
# Create directory for results
# Create directory for results
dir.create(file.path(ROOT_DIR, "HumanRegionsPooled_ttest"), showWarnings = FALSE)

# Initialize progress bar
pb <- txtProgressBar(min = 0, max = length(genes), style = 3)

# Initialize empty list to store all results
all_results <- list()

for (i in seq_along(genes)) {
  gene <- genes[i]
  setTxtProgressBar(pb, i)
  message("\nProcessing gene: ", gene)
  
  # Try-catch block for the main analysis
  tryCatch({
    # Run t-test analysis first
    t_test_results <- plot_gene_distribution_with_stats(
      data_list = species_data,
      species_names = species_names,
      gene_name = gene
    )
    
    # If we get valid t_test_results, proceed with mean calculations
    if (!is.null(t_test_results) && nrow(t_test_results) > 0) {
      # Calculate means for each species
      means_list <- list()
      
      for (j in seq_along(species_data)) {
        species <- species_data[[j]]
        species_name <- species_names[j]
        
        # Only calculate means for groups present in t_test_results for this species
        species_results <- t_test_results %>% 
          filter(Species == species_name)
        
        # Process each group type present in the results
        for (group_type in unique(species_results$GroupType)) {
          for (group_name in unique(species_results$Group[species_results$GroupType == group_type])) {
            if (group_type == "None") {
              # Calculate overall means (AreaCentralis)
              mu <- species %>%
                dplyr::select(all_of(gene), Location) %>%
                group_by(Location) %>%
                summarise(Mean = mean(get(gene), na.rm = TRUE)) %>%
                ungroup() %>%
                t() %>%
                data.frame() %>%
                rename(Fovea_mean = X1,
                       Periphery_mean = X2) %>%
                remove_rownames() %>%
                tail(1) %>%
                mutate(Species = species_name,
                       Group = 'Overall',
                       GroupType = 'None',
                       gene = gene)
            } else {
              # Calculate means for specific group (Subtype or Orthotype)
              mu <- species %>%
                filter(get(group_type) == group_name) %>%
                dplyr::select(all_of(gene), Location) %>%
                group_by(Location) %>%
                summarise(Mean = mean(get(gene), na.rm = TRUE)) %>%
                ungroup() %>%
                t() %>%
                data.frame() %>%
                rename(Fovea_mean = X1,
                       Periphery_mean = X2) %>%
                remove_rownames() %>%
                tail(1) %>%
                mutate(Species = species_name,
                       Group = group_name,
                       GroupType = group_type,
                       gene = gene)
            }
            means_list[[length(means_list) + 1]] <- mu
          }
        }
      }
      
      # Combine all means
      combined_means <- bind_rows(means_list)
      
      # Join with t-test results
      final_results <- t_test_results %>%
        left_join(combined_means, 
                  by = c("Species", "Group", "GroupType", "gene"))
      
      # Add to comprehensive results list
      all_results[[gene]] <- final_results
    }
  }, error = function(e) {
    message("Error processing gene ", gene, ": ", conditionMessage(e))
    return(NULL)
  })
}

# Close progress bar
close(pb)

# Combine all results and remove NULL entries
final_comprehensive_results <- bind_rows(all_results)  
  

# Save results
dir.create(file.path(ROOT_DIR,  "AreaCentralis_ttest"), showWarnings = FALSE)

write.csv(final_comprehensive_results, 
          file = file.path(ROOT_DIR, "AreaCentralis_ttest", "comprehensive_results.csv"),
          row.names = FALSE)

message("\nAnalysis complete. Results saved to: ", 
        file.path(ROOT_DIR, "AreaCentralis_ttest", "comprehensive_results.csv"))


# Combine all results and remove NULL entries
ac_comprehensive_results <- final_comprehensive_results %>%
  filter(Group == "Overall") %>%
  dplyr::select(-Group, -GroupType)

# Save results
write.csv(ac_comprehensive_results, 
          file = file.path(ROOT_DIR, "AreaCentralis_ttest", "ACvPeriphery_results.csv"),
          row.names = FALSE)

message("\nAnalysis complete. Results saved to: ", 
        file.path(ROOT_DIR, "AreaCentralis_ttest", "ACvPeriphery_results.csv"))




analyze_species_correlations <- function(final_comprehensive_results, 
                                         x_axis_species = "Mouse",
                                         target_genes = NULL,
                                         group_by = "Overall",
                                         count_thresh = 5) {
  
  # Initial filtering
  if(group_by == "Overall") {
    filtered_results <- final_comprehensive_results %>%
      filter(Group == "Overall") %>%
      mutate(Fovea_mean = as.numeric(Fovea_mean),
             Periphery_mean = as.numeric(Periphery_mean)) %>%
      # Calculate species-specific means for threshold scaling
      group_by(Species) %>%
      mutate(species_mean = mean(Fovea_mean + Periphery_mean, na.rm = TRUE),
             relative_sum = (Fovea_mean + Periphery_mean) / species_mean) %>%
      ungroup() %>%
      # Filter genes where ALL species are below threshold
      group_by(gene) %>%
      filter(
        gene %in% target_genes,  # Separate the conditions
        !all(relative_sum < count_thresh)  # Keep if NOT all species are below threshold
      ) %>%
      ungroup() %>%
      # Clean up temporary columns
      dplyr::select(-species_mean, -relative_sum)
      
  } else if(group_by == "Orthotype") {
    filtered_results <- final_comprehensive_results %>%
      filter(str_detect(Group, "^O") & Group != "Overall") %>%
      mutate(Fovea_mean = as.numeric(Fovea_mean),
             Periphery_mean = as.numeric(Periphery_mean) ) 
  } else {
    stop("group_by must be either 'Overall' or 'Orthotype'")
  }
  
  # Calculate log fold changes with proper numeric conversion
  results_with_lfc <- filtered_results %>%
    mutate(log2FC = log2(Fovea_mean / Periphery_mean)) %>%
    filter(log2FC != abs(Inf))
  
  # Handle grouping logic
  if(group_by == "Overall") {
    results_processed <- results_with_lfc %>%
      dplyr::select(Species, gene, log2FC)
    
    if(!is.null(target_genes)) {
      results_processed <- results_processed %>%
        filter(gene %in% target_genes)
    }
    id_col <- "gene"
    
  } else {
    results_processed <- results_with_lfc %>%
      dplyr::select(Species, gene, log2FC, Group)
    
    if(!is.null(target_genes)) {
      results_processed <- results_processed %>%
        filter(gene %in% target_genes) %>%
        group_by(Group, Species) %>%
        summarize(log2FC = mean(log2FC)) %>%
        ungroup()
    }
    id_col <- "Group"
  }

  
  # Pivot and clean data
  wide_data <- results_processed %>%
    pivot_wider(
      id_cols = all_of(id_col),
      names_from = Species,
      values_from = log2FC
    ) %>%
    drop_na()
  
  y_axis_species <- setdiff(unique(filtered_results$Species), x_axis_species)
  regression_list <- list()
  plots_list <- list()
  
  for(species in y_axis_species) {
    analysis_data <- wide_data %>%
      dplyr::select(all_of(c(id_col, x_axis_species, species))) %>%
      drop_na() %>%
      filter(is.finite(get(species)) & is.finite(get(x_axis_species)))
    
    if(nrow(analysis_data) > 0) {
      tryCatch({
        model <- lm(get(species) ~ get(x_axis_species), data = analysis_data)
        summary_stats <- summary(model)
        
        regression_list[[species]] <- data.frame(
          y_axis_species = species,
          x_axis_species = x_axis_species,
          slope = coef(model)[2],
          intercept = coef(model)[1],
          r_squared = summary_stats$r.squared,
          p_value = summary_stats$coefficients[2,4],
          n_points = nrow(analysis_data)
        )
        
        p <- ggplot(analysis_data, aes_string(x = x_axis_species, y = species)) +
          geom_point() +
          geom_smooth(method = "lm", formula = y ~ x, color = "red") +
          theme_minimal() +
          labs(
            title = paste(species, "vs", x_axis_species),
            x = paste(x_axis_species, "log2FC"),
            y = paste(species, "log2FC")
          ) +
          annotate(
            "text",
            x = min(analysis_data[[x_axis_species]]),
            y = max(analysis_data[[species]]),
            label = sprintf(
              "R² = %.3f\np = %.2e\nn = %d",
              summary_stats$r.squared,
              summary_stats$coefficients[2,4],
              nrow(analysis_data)
            ),
            hjust = 0,
            vjust = 1
          )
        
        if(group_by == "Orthotype") {
          p <- p + geom_text_repel(aes_string(label = id_col), size = 3)
        } else {
          p <- p + geom_text_repel(aes_string(label = id_col), size = 3)
        }
        
        plots_list[[species]] <- p
      }, error = function(e) {
        warning(sprintf("Error processing species %s: %s", species, e$message))
        NULL
      })
    }
  }
  
  regression_df <- bind_rows(regression_list)
  combined_plot <- wrap_plots(plots_list, ncol = 2)
  
  return(list(
    regression_results = regression_df,
    plots = combined_plot
  ))
}


calcium_channels <- c("CAchannels", "Cacna1a", "Cacna1b", "Cacna1e", "Cacna1g", "Cacna1i", "Cacna2d1", "Cacna2d2", "Cacna2d3", "Cacna2d4", "Cacnb3", "Cacnb4", "Cacng2", "Cacng3", "Cacng4", "Cacng5", "Cacng7")

potassium_channels <- c("Kchannels", "Hcn1", "Hcn2", "Kcna1", "Kcna2", "Kcna4", "Kcna5", "Kcna6", "Kcnab1", "Kcnab2", "Kcnab3", "Kcnb1", "Kcnb2", "Kcnc1", "Kcnc2", "Kcnc3", "Kcnc4", "Kcnd2", "Kcnd3", "Kcnh1", "Kcnh2", "Kcnip1", "Kcnip2", "Kcnip3", "Kcnip4", "Kcnj12", "Kcnj3", "Kcnj9", "Kcnk1", "Kcnn1", "Kcnn2", "Kcnq1ot1", "Kcnq2", "Kcnq3")

sodium_channels <- c("NAchannels", "Scn1a", "Scn1b", "Scn2b", "Scn3a", "Scn3b", "Scn4b", "Scn7a", "Scn8a", "Scn9a")

glutamate_receptors <- c("Glutamate", "Gria1", "Gria2", "Gria3", "Gria4", "Grik1", "Grik5", "Grin1", "Grin2b", "Grin3a", "Grm3", "Grm4", "Grm5", "Grm6", "Grm8")

gaba_glycine_receptors <- c("GABA_Gly", "Gabra1", "Gabra2", "Gabra3", "Gabra4", "Gabrb1", "Gabrb3", "Gabrg1", "Gabrg2", "Gabrg3", "Gabrr1", "Gabrr2", "Gabrr3", "Glra1", "Glra2", "Glrb")
gaba_receptors <- c("GABA", "Gabra1", "Gabra2", "Gabra3", "Gabra4", "Gabrb1", "Gabrb3", "Gabrg1", "Gabrg2", "Gabrg3", "Gabrr1", "Gabrr2", "Gabrr3")# "Glra1", "Glra2", "Glrb")


for (gene_set in list(sodium_channels, 
                      gaba_receptors,
                      gaba_glycine_receptors, 
                      glutamate_receptors,
                      potassium_channels, 
                      calcium_channels)) {
  
  # Example usage:
  mouse_reg_result <- analyze_species_correlations(
    final_comprehensive_results,
    x_axis_species = "Mouse",
    target_genes = gene_set,
    count_thresh = 0.5
  )
    write.csv(mouse_reg_result$regression_results, 
              file = file.path(ROOT_DIR, "AreaCentralis_ttest", paste0(gene_set[1], "_Mouse_regression_results.csv")),
              row.names = FALSE)
    
    print(mouse_reg_result$plots)
  
  # (macaque_reg_result <- analyze_species_correlations(
  #   final_comprehensive_results,
  #   x_axis_species = "Macaque",
  #   target_genes = gaba_glycine_receptors,
  #   count_thresh = 0.5
  # ))
  # 
  #   write.csv(macaque_reg_result$regression_results, 
  #             file = file.path(ROOT_DIR, "AreaCentralis_ttest", paste0(gene_set[1], "_Maqaque_regression_results.csv")),
  #             row.names = FALSE)
  #   print(macaque_reg_result$plots)

}


functional_measures <- c("Functional_Group","RF_size","sup_PeakFiringRate",
                         "sup_ResponseDuration","sup_BLFiringRate","sup_SuppressionIndex",     
                         "sup_OnOffIndex", "sup_PeakResponseLatency")
unif_df <- read.csv(paste0(DATA_DIR, 'RGC Subtype Unification.csv')) %>%
  filter(Goetz2022_PatchSeq == 1) %>%
  mutate(Subtype = paste0("T", Tran2019_Clusters),
         Functional_Group = factor(Functional_Group),
         RF_size = as.numeric( Spot_Size_.um.)) %>%
  dplyr::select(Subtype, all_of(functional_measures)) %>%
  dplyr::rename(PeakFiringRate = sup_PeakFiringRate,
         ResponseDuration = sup_ResponseDuration,
         BLFiringRate = sup_BLFiringRate,
         SuppressionIndex = sup_SuppressionIndex,
         OnOffIndex = sup_OnOffIndex,
         PeakResponseLatency = sup_PeakResponseLatency) %>%
  left_join(mouse %>%
              dplyr::select(Subtype, Location),
            by = "Subtype") %>%
  mutate(Subtype = factor(Subtype),
         Location = factor(Location)
          ) 


# Apply Bonferroni correction
n_tests <- length(unique(unif_df$Functional_Group))
alpha <- 0.05
bonferroni_threshold <- alpha / n_tests

# Chi-square test for Functional_Group
func_table <- table(unif_df$Location, unif_df$Functional_Group)
# Convert raw counts to proportions within each region
fovea_total <- sum(mouse$Location=="Fovea")
periphery_total <- sum(mouse$Location=="Fovea")
fovea_props <- func_table["Fovea",] / fovea_total
periphery_props <- func_table["Periphery",] / periphery_total



# First perform overall chi-square test
chisq_result <- chisq.test(func_table)

# Create function for post-hoc analysis of individual functional groups
analyze_functional_group <- function(group, data, bonferroni_threshold) {
  # Subset data for specific functional group
  group_data <- data %>%
    mutate(is_group = Functional_Group == group) %>%
    with(table(Location, is_group))
  
  # Perform chi-square test for this group
  group_test <- chisq.test(group_data)
  
  # Calculate proportions in each location
  props <- data %>%
    group_by(Location) %>%
    summarize(prop = mean(Functional_Group == group))
  
  return(list(
    group = group,
    p_value = group_test$p.value,
    significant = group_test$p.value < bonferroni_threshold,
    fovea_prop = props$prop[props$Location == "Fovea"],
    periphery_prop = props$prop[props$Location == "Periphery"]
  ))
}

# Apply analysis to each functional group
functional_groups <- levels(unif_df$Functional_Group)
results <- lapply(functional_groups, function(group) {
  analyze_functional_group(group, unif_df, bonferroni_threshold)
})

# Create summary dataframe
results_df <- do.call(rbind, lapply(results, function(x) {
  data.frame(
    Functional_Group = x$group,
    P_Value = x$p_value,
    Significant = x$significant,
    Fovea_Proportion = x$fovea_prop,
    Periphery_Proportion = x$periphery_prop,
    Fold_Change = log2( x$fovea_prop/ x$periphery_prop)
  )
}))



# Create a custom color scale based on p-values
# Convert p-values to a gray scale where:
# p > alpha -> white with black outline
# alpha > p > alpha/1000 -> graduating gray scale
# p < alpha/1000 -> black
results_df <- results_df %>%
  mutate(
    color_value = case_when(
      P_Value > bonferroni_threshold ~ "white",
      P_Value < bonferroni_threshold/1000 ~ "black",
      TRUE  ~ "grey"
    )
  )


ggplot(results_df, aes(x = reorder(Functional_Group, -Fold_Change))) +
  geom_point(aes(y = Fold_Change, fill = color_value), 
             shape = 21, # Use filled circle with outline
             size = 3,
             color = "black") + # Always black outline
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_identity() + # Use the actual colors we specified
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Functional Group",
    y = "Log2 Fold Change (ART/Periphery)",
    title = "Distribution of Functional Groups Between Locations",
    subtitle = paste("Bonferroni-corrected α =", bonferroni_threshold)
  )

write_csv(results_df, file = file.path(ROOT_DIR, "AreaCentralis_ttest", "Functional_chiSquare.csv"))




# Get the receptor columns (excluding Subtype and Location)
receptor_cols <- setdiff(names(dplyr::select(mouse, any_of(gaba_glycine_receptors))), c("Subtype", "Location"))

# Calculate both raw means and indices
GABA_analysis_mouse <- mouse %>%
  dplyr::select(Subtype, Location, any_of(gaba_glycine_receptors)) %>%
  # Create index columns: (x-Gabra1)/(x+Gabra1)
  mutate(across(all_of(receptor_cols), 
                list(index = ~(.-Gabra1)/(. + Gabra1)),
                .names = "{.col}_index")) %>%
  # Group and calculate means for both original values and indices
  group_by(Subtype, Location) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

## Create long format data for plotting
GABA_long <- mouse %>%
  dplyr::select(Subtype, Location, any_of(gaba_glycine_receptors)) %>%
  # Calculate indices for each receptor
  mutate(across(all_of(receptor_cols), 
                list(index = ~(.-Gabra1)/(. + Gabra1)),
                .names = "{.col}_index")) %>%
  # Keep only index columns, Subtype, and Location
  dplyr::select(Subtype, Location, ends_with("_index")) %>%
  # Reshape to long format
  pivot_longer(cols = ends_with("_index"),
               names_to = "Receptor",
               values_to = "Index") %>%
  # Clean up receptor names by removing "_index" suffix
  mutate(Receptor = str_remove(Receptor, "_index"))

# Create faceted boxplots
ggplot(GABA_long, aes(x = Index, y = Subtype, fill = Location)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  # Facet by receptor
  facet_wrap(~Receptor, scales = "free_x") +
  # Add theme elements
  theme_bw() +
  theme(axis.text.y = element_text(size = 7),
        strip.text = element_text(size = 10),
        legend.position = "top") +
  # Add labels
  labs(x = "Index relative to Gabra1: (x-Gabra1)/(x+Gabra1)",
       y = "Subtype",
       title = "Distribution of GABA/Glycine Receptor Expression Indices") +
  # Use colorblind-friendly colors
  scale_fill_brewer(palette = "Set2")



gaba_receptors <- c("Gabra1", "Gabra2", "Gabra3", "Gabra4", "Gabrb1", "Gabrb3", "Gabrg1", "Gabrg2", "Gabrg3", "Gabrr1", "Gabrr2", "Gabrr3")
gaba_a_receptors <- c("Gabra1", "Gabra2", "Gabra3", "Gabra4", "Gabrb1", "Gabrb3", "Gabrg1", "Gabrg2", "Gabrg3")
gaba_r_receptors <- c("Gabrr1", "Gabrr2", "Gabrr3")
glycine_receptors <- c("Glra1", "Glra2", "Glrb")

GABA_a_mouse <- mouse %>%
  dplyr::select(Subtype, Location, any_of(gaba_a_receptors)) %>%
  # Calculate indices for each receptor
  mutate(alpha = Gabra1+ Gabra2 + Gabra3 + Gabra4,
         beta = Gabrb1+ Gabrb3,
         gamma = Gabrg1+ Gabrg2 + Gabrg3) %>%
  group_by(Subtype, Location) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  mutate(across(all_of(c( "alpha", "beta", "gamma", gaba_a_receptors) ), 
                list(prop = ~./alpha),
                .names = "{.col}_prop")) 


GABA_r_mouse <- mouse %>%
  dplyr::select(Subtype, Location, any_of(gaba_r_receptors)) %>%
  # Calculate indices for each receptor
  mutate(rho = Gabrr1+ Gabrr2 + Gabrr3) %>%
  group_by(Subtype, Location) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  mutate(across(all_of(c( "rho", gaba_r_receptors) ), 
                list(prop = ~./rho),
                .names = "{.col}_prop")) 

GABA_mouse <- mouse %>%
  dplyr::select(Subtype, Location, any_of(gaba_receptors)) %>%
  # Calculate indices for each receptor
  mutate(alpha = Gabra1+ Gabra2 + Gabra3 + Gabra4,
         beta = Gabrb1+ Gabrb3,
         gamma = Gabrg1+ Gabrg2 + Gabrg3,
         rho = Gabrr1+ Gabrr2 + Gabrr3) %>%
  group_by(Subtype, Location) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  mutate(across(all_of(c( "alpha", "rho", "beta", "gamma", gaba_receptors) ), 
                list(prop = ~./alpha),
                .names = "{.col}_prop")) 

GABA_mouse_stoichiometries <- GABA_mouse %>%
  dplyr::select(Subtype, Location, alpha, beta, gamma, rho) %>%
  group_by(Location, Subtype) %>%
  mutate(alpha =round(alpha,2),
         beta=round(beta),
         gamma=round(gamma,2),
         rho=round(rho,2),
         min_val = min(alpha, beta, gamma, rho),
         factor = 1/min_val,
         alpha =round(factor*alpha,2),
         beta=round(factor*beta),
         gamma=round(factor*gamma,2),
         rho=round(factor*rho,2),
         max_val = max(alpha, beta, gamma, rho)         
         ) %>%
  dplyr::select(-factor) %>%


unif_GABA <- read.csv(paste0(DATA_DIR, 'RGC Subtype Unification.csv')) %>%
  filter(Goetz2022_PatchSeq == 1) %>%
  mutate(Subtype = paste0("T", Tran2019_Clusters),
         Functional_Group = factor(Functional_Group) ) %>%
  dplyr::select(Subtype, Functional_Group) %>%
  left_join(GABA_mouse ,
            by = "Subtype") %>%
  mutate(Subtype = factor(Subtype),
         Location = factor(Location)
  ) %>%
  dplyr::select(-Subtype) %>%
  group_by(Functional_Group, Location) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

  
gaba_receptors_h <- c( "Gabra1",        "Gabra2",        "Gabra3",        "Gabra4",        
                       "Gabrb1",        "Gabrb3",      
                       "Gabrg1",        "Gabrg2",        "Gabrg3" ,      
                       "Gabrr1",        "Gabrr2",        "Gabrr3" )

GABA_human <- human %>%
  dplyr::select(organ, any_of(gaba_receptors_h)) %>%
  # Calculate indices for each receptor
  mutate(alpha = Gabra1+ Gabra2 + Gabra3 + Gabra4,
         beta = Gabrb1+ Gabrb3,
         gamma = Gabrg1+ Gabrg2 + Gabrg3,
         rho = Gabrr1+ Gabrr2 + Gabrr3) %>%
  group_by(organ) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  mutate(across(all_of(c( "alpha", "rho", "beta", "gamma", gaba_receptors_h) ), 
                list(prop = ~./alpha),
                .names = "{.col}_prop")) 

GABA_mouse_sum <- mouse %>%
  dplyr::select(Location, any_of(gaba_receptors)) %>%
  # Calculate indices for each receptor
  mutate(alpha = Gabra1+ Gabra2 + Gabra3 + Gabra4,
         beta = Gabrb1+ Gabrb3,
         gamma = Gabrg1+ Gabrg2 + Gabrg3,
         rho = Gabrr1+ Gabrr2 + Gabrr3) %>%
  group_by(Location) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  mutate(across(all_of(c( "alpha", "rho", "beta", "gamma", gaba_receptors) ), 
                list(prop = ~./alpha),
                .names = "{.col}_prop")) 



####################################################################################
# CSV table output of GABA stoichiometries
##################################################################################3
# Function to calculate receptor stoichiometries
calculate_stoichiometries <- function(data, species_name, group_cols = NULL) {
  # Set default group columns based on species
  if (is.null(group_cols)) {
    if (species_name == "Human") {
      group_cols <- c("organ", "rgc_cluster")
    } else {
      group_cols <- c("Location", "Subtype")
    }
  }
  
  # Define receptor groups
  gaba_receptors <- c("Gabra1", "Gabra2", "Gabra3", "Gabra4", 
                      "Gabrb1", "Gabrb3",
                      "Gabrg1", "Gabrg2", "Gabrg3",
                      "Gabrr1", "Gabrr2", "Gabrr3")
  
  glycine_receptors <- c("Glra1", "Glra2", "Glrb")
  
  # Add missing columns with NA values
  missing_cols <- setdiff(c(gaba_receptors, glycine_receptors), names(data))
  if (length(missing_cols) > 0) {
    for (col in missing_cols) {
      data[[col]] <- NA_real_
    }
  }
  
  # Calculate combined metrics
  results <- data %>%
    # Select relevant columns
    dplyr::select(all_of(group_cols), any_of(c(gaba_receptors, glycine_receptors))) %>%
    
    # Calculate combined metrics
    mutate(
      # GABA receptor subunits
      alpha = rowSums(across(c("Gabra1", "Gabra2", "Gabra3", "Gabra4"), ~replace_na(., 0))),
      beta = rowSums(across(c("Gabrb1", "Gabrb3"), ~replace_na(., 0))),
      gamma = rowSums(across(c("Gabrg1", "Gabrg2", "Gabrg3"), ~replace_na(., 0))),
      rho = rowSums(across(c("Gabrr1", "Gabrr2", "Gabrr3"), ~replace_na(., 0))),
      
      # Glycine receptor subunits
      Glyr_alpha = rowSums(across(c("Glra1", "Glra2"), ~replace_na(., 0))),
      Glyr_beta = replace_na(Glrb, 0),
      
      # Combined metrics
      GABA_r_all = alpha + beta + gamma + rho,
      Glyr_all = Glyr_alpha + Glyr_beta
    ) %>%
    
    # Group by location and subtype
    group_by(across(all_of(group_cols))) %>%
    summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
    
    # Calculate GABA stoichiometries
    mutate(
      # Round initial values
      across(c(alpha, beta, gamma, rho), ~round(., 2)),
      
      # Calculate GABA stoichiometry factors
      min_val_gaba = pmin(alpha, beta, gamma, rho, na.rm = TRUE),
      gaba_factor = ifelse(min_val_gaba > 0, 1/min_val_gaba, NA_real_),
      
      # Calculate normalized GABA stoichiometries as integers
      alpha_stoich = round(alpha * gaba_factor),
      beta_stoich = round(beta * gaba_factor),
      gamma_stoich = round(gamma * gaba_factor),
      rho_stoich = round(rho * gaba_factor),
      
      # Calculate Glycine stoichiometries
      min_val_gly = pmin(Glyr_alpha, Glyr_beta, na.rm = TRUE),
      gly_factor = ifelse(min_val_gly > 0, 1/min_val_gly, NA_real_),
      Glyr_alpha_stoich = round(Glyr_alpha * gly_factor),
      Glyr_beta_stoich = round(Glyr_beta * gly_factor),
      
      # Calculate GABA vs Glycine stoichiometry
      min_val_total = pmin(GABA_r_all, Glyr_all, na.rm = TRUE),
      total_factor = ifelse(min_val_total > 0, 1/min_val_total, NA_real_),
      GABA_total_stoich = round(GABA_r_all * total_factor),
      Glyr_total_stoich = round(Glyr_all * total_factor)
    ) %>%
    
    # Add species identifier
    mutate(Species = species_name) %>%
    
    # Select and arrange columns
    dplyr::select(
      Species,
      all_of(group_cols),
      # Raw subunit values
      all_of(gaba_receptors),
      all_of(glycine_receptors),
      # Combined metrics
      alpha, beta, gamma, rho,
      Glyr_alpha, Glyr_beta,
      GABA_r_all, Glyr_all,
      # Stoichiometries
      ends_with("stoich")
    )
  
  # Calculate aggregate (ALL) values
  all_results <- results %>%
    ungroup() %>%
    dplyr::select(-all_of(group_cols[2])) %>%  # Remove Subtype/rgc_cluster
    group_by(Species, across(all_of(group_cols[1]))) %>%  # Group by Species and Location/organ
    summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
    mutate(!!sym(group_cols[2]) := "ALL") %>%  # Add back Subtype/rgc_cluster column
    dplyr::select(names(results))  # Reorder columns to match
  
  # Combine regular and ALL results
  bind_rows(results, all_results)
}

# Process each species dataset
mouse_results <- calculate_stoichiometries(mouse, "Mouse")
human_results <- calculate_stoichiometries(human, "Human") %>%
  rename(Location = organ, Subtype = rgc_cluster)
macaque_results <- calculate_stoichiometries(macaque, "Macaque")

# Combine all results
all_results <- bind_rows(mouse_results, human_results, macaque_results) %>%
  dplyr::select(Species, Location, Subtype, 
         GABA_total_stoich,  Glyr_total_stoich,
         alpha_stoich, beta_stoich, gamma_stoich, rho_stoich,
         Glyr_alpha_stoich, Glyr_beta_stoich,
         alpha, Gabra1, Gabra2, Gabra3, Gabra4,
         beta, Gabrb1,Gabrb3,           
         gamma, Gabrg1, Gabrg2, Gabrg3, 
         rho, Gabrr1, Gabrr2, Gabrr3,
         Glyr_alpha, Glra1, Glra2, 
         Glyr_beta, Glrb,             
         GABA_r_all, Glyr_all, 
         )

# Write output
write.csv(all_results, file.path(ROOT_DIR, "GABA_stoichiometries.csv"), row.names = FALSE)




# Function to normalize stoichiometry values
normalize_stoich <- function(df) {
  df %>%
    group_by(Species, Location, Subtype) %>%
    mutate(total = sum(stoichiometry),
           stoichiometry = stoichiometry / total) %>%
    ungroup()
}

# Function to calculate max stoichiometric range for GABA receptors
calc_gaba_range <- function(data) {
  data %>%
    dplyr::mutate(
      max_ratio = pmax(
        alpha_stoich/beta_stoich,
        alpha_stoich/gamma_stoich,
        beta_stoich/gamma_stoich,
        na.rm = TRUE
      ),
      min_ratio = pmin(
        alpha_stoich/beta_stoich,
        alpha_stoich/gamma_stoich,
        beta_stoich/gamma_stoich,
        na.rm = TRUE
      ),
      stoich_range = max_ratio/min_ratio
    )
}

# Prepare base dataset with locations renamed
base_data <- all_results %>%
  filter(Species %in% c("Human", "Mouse")) %>%
  filter(
    (Species == "Human" & Location %in% c("peripheral region of retina", "fovea centralis")) |
      (Species == "Mouse" & Location %in% c("Fovea", "Periphery"))) %>%
  dplyr::mutate(
    Location = case_when(
      Location == "peripheral region of retina" ~ "Periphery",
      Location == "fovea centralis" ~ "Fovea",
      Location == "Fovea" ~ "ART",
      TRUE ~ Location
    )
  )

# Find extreme GABA stoichiometry subtypes
extreme_subtypes <- base_data %>%
  group_by(Species, Subtype) %>%
  calc_gaba_range() %>%
  dplyr::summarise(max_range = max(stoich_range, na.rm = TRUE), .groups = "drop") %>%
  group_by(Species) %>%
  dplyr::slice(which.max(max_range), which.min(max_range)) %>%
  dplyr::pull(Subtype)

# Prepare GABA visualization data
gaba_viz_data <- base_data %>%
  filter(Subtype %in% c("ALL", extreme_subtypes)) %>%
  dplyr::select(Species, Location, Subtype, 
                alpha_stoich, beta_stoich, gamma_stoich, rho_stoich) %>%
  # Add control data
  bind_rows(data.frame(
    Species = "Control",
    Location = "Expected",
    Subtype = "Expected",
    alpha_stoich = 20,
    beta_stoich = 20,
    gamma_stoich = 10,
    rho_stoich = 1
  )) %>%
  # Create ordered display groups
  dplyr::mutate(
    display_group = case_when(
      Subtype == "Expected" ~ paste0("A_", Subtype),
      Subtype == "ALL" & Species == "Mouse" & Location == "ART" ~ "B_Mouse_ALL_ART",
      Subtype == "ALL" & Species == "Mouse" & Location == "Periphery" ~ "C_Mouse_ALL_Periphery",
      Subtype == "ALL" & Species == "Human" & Location == "Fovea" ~ "D_Human_ALL_Fovea",
      Subtype == "ALL" & Species == "Human" & Location == "Periphery" ~ "E_Human_ALL_Periphery",
      TRUE ~ paste0(
        LETTERS[6 + as.numeric(factor(Subtype))],
        "_",
        Species,
        "_",
        Subtype,
        "_",
        Location
      )
    ),
    display_group = factor(display_group, levels = sort(unique(display_group)))
  ) %>%
  # Reshape for plotting
  pivot_longer(cols = ends_with("stoich"),
               names_to = "subunit",
               values_to = "stoichiometry") %>%
  dplyr::mutate(subunit = factor(subunit, 
                                 levels = c("alpha_stoich", "beta_stoich", 
                                            "gamma_stoich", "rho_stoich"))) %>%
  # Normalize values
  normalize_stoich()

# Find extreme Glycine stoichiometry subtypes
gly_extreme_subtypes <- base_data %>%
  group_by(Species, Subtype) %>%
  dplyr::mutate(gly_ratio = Glyr_alpha_stoich/Glyr_beta_stoich) %>%
  dplyr::summarise(ratio_range = max(gly_ratio, na.rm = TRUE)/min(gly_ratio, na.rm = TRUE), 
                   .groups = "drop") %>%
  group_by(Species) %>%
  dplyr::slice(which.max(ratio_range), which.min(ratio_range)) %>%
  dplyr::pull(Subtype)

# Prepare Glycine visualization data
gly_viz_data <- base_data %>%
  filter(Subtype %in% c("ALL", gly_extreme_subtypes)) %>%
  dplyr::select(Species, Location, Subtype, 
                Glyr_alpha_stoich, Glyr_beta_stoich) %>%
  # Add control data
  bind_rows(data.frame(
    Species = "Control",
    Location = "Expected",
    Subtype = "Expected",
    Glyr_alpha_stoich = 3,
    Glyr_beta_stoich = 2
  )) %>%
  # Create ordered display groups (matching GABA order)
  dplyr::mutate(
    display_group = case_when(
      Subtype == "Expected" ~ paste0("A_", Subtype),
      Subtype == "ALL" & Species == "Mouse" & Location == "ART" ~ "B_Mouse_ALL_ART",
      Subtype == "ALL" & Species == "Mouse" & Location == "Periphery" ~ "C_Mouse_ALL_Periphery",
      Subtype == "ALL" & Species == "Human" & Location == "Fovea" ~ "D_Human_ALL_Fovea",
      Subtype == "ALL" & Species == "Human" & Location == "Periphery" ~ "E_Human_ALL_Periphery",
      TRUE ~ paste0(
        LETTERS[6 + as.numeric(factor(Subtype))],
        "_",
        Species,
        "_",
        Subtype,
        "_",
        Location
      )
    ),
    display_group = factor(display_group, levels = sort(unique(display_group)))
  ) %>%
  pivot_longer(cols = ends_with("stoich"),
               names_to = "subunit",
               values_to = "stoichiometry") %>%
  dplyr::mutate(subunit = factor(subunit, 
                                 levels = c("Glyr_alpha_stoich", "Glyr_beta_stoich"))) %>%
  normalize_stoich()

# Prepare ratio visualization data using the same subtypes as GABA
ratio_viz_data <- base_data %>%
  filter(Subtype %in% c("ALL", extreme_subtypes)) %>%
  dplyr::select(Species, Location, Subtype, 
                GABA_total_stoich, Glyr_total_stoich) %>%
  # Add control data
  bind_rows(data.frame(
    Species = "Control",
    Location = "Expected",
    Subtype = "Expected",
    GABA_total_stoich = 55,
    Glyr_total_stoich = 45
  )) %>%
  # Create ordered display groups (matching GABA order)
  dplyr::mutate(
    display_group = case_when(
      Subtype == "Expected" ~ paste0("A_", Subtype),
      Subtype == "ALL" & Species == "Mouse" & Location == "ART" ~ "B_Mouse_ALL_ART",
      Subtype == "ALL" & Species == "Mouse" & Location == "Periphery" ~ "C_Mouse_ALL_Periphery",
      Subtype == "ALL" & Species == "Human" & Location == "Fovea" ~ "D_Human_ALL_Fovea",
      Subtype == "ALL" & Species == "Human" & Location == "Periphery" ~ "E_Human_ALL_Periphery",
      TRUE ~ paste0(
        LETTERS[6 + as.numeric(factor(Subtype))],
        "_",
        Species,
        "_",
        Subtype,
        "_",
        Location
      )
    ),
    display_group = factor(display_group, levels = sort(unique(display_group))),
    total = GABA_total_stoich + Glyr_total_stoich,
    GABA_prop = GABA_total_stoich/total,
    Gly_prop = Glyr_total_stoich/total
  ) %>%
  pivot_longer(
    cols = c(GABA_prop, Gly_prop),
    names_to = "receptor_type",
    values_to = "proportion"
  ) %>%
  dplyr::mutate(
    receptor_type = factor(receptor_type, 
                           levels = c("GABA_prop", "Gly_prop"))
  )



# Function to create nice labels
create_labels <- function(display_group) {
  gsub("_", " ", substr(display_group, 3, nchar(display_group)))
}

# GABA Plot
p1 <- ggplot(gaba_viz_data, 
             aes(x = "1", y = stoichiometry, fill = subunit)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~display_group, labeller = labeller(display_group = create_labels)) +
  scale_fill_brewer(palette = "Set1", 
                    labels = c("α", "β", "γ", "ρ")) +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  labs(title = "GABA Receptor Stoichiometry",
       y = "Relative Proportion",
       x = NULL,
       fill = "Subunit") +
  theme(axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.text.x = element_text(angle = 45))

# Glycine Plot
p2 <- ggplot(gly_viz_data, 
             aes(x = "1", y = stoichiometry, fill = subunit)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~display_group, labeller = labeller(display_group = create_labels)) +
  scale_fill_brewer(palette = "Set2", 
                    labels = c("α", "β")) +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  labs(title = "Glycine Receptor Stoichiometry",
       y = "Relative Proportion",
       x = NULL,
       fill = "Subunit") +
  theme(axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.text.x = element_text(angle = 45))

# Ratio Plot
p3 <- ggplot(ratio_viz_data, 
             aes(x = "1", y = proportion, fill = receptor_type)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~display_group, labeller = labeller(display_group = create_labels)) +
  scale_fill_brewer(palette = "Set3", 
                    labels = c("GABA", "Glycine")) +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  labs(title = "GABA/Glycine Receptor Ratio",
       y = "Relative Proportion",
       x = NULL,
       fill = "Receptor") +
  theme(axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.text.x = element_text(angle = 45))


# Combine plots
combined_plot <- p1 / p2 / p3 +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(
    title = "Receptor Stoichiometry Comparison",
    theme = theme_minimal()
  )

# Save the plot
ggsave(file.path(ROOT_DIR, "receptor_stoichiometry_comparison.pdf"), 
       combined_plot, 
       width = 20, 
       height = 12)

# Prepare GABA data for export
gaba_export <- gaba_viz_data %>%
  dplyr::select(Species, Location, Subtype, display_group, subunit, stoichiometry) %>%
  tidyr::pivot_wider(
    names_from = subunit,
    values_from = stoichiometry
  ) %>%
  dplyr::arrange(display_group)

# Prepare Glycine data for export
glycine_export <- gly_viz_data %>%
  dplyr::select(Species, Location, Subtype, display_group, subunit, stoichiometry) %>%
  tidyr::pivot_wider(
    names_from = subunit,
    values_from = stoichiometry
  ) %>%
  dplyr::arrange(display_group)

# Prepare ratio data for export
ratio_export <- ratio_viz_data %>%
  dplyr::select(Species, Location, Subtype, display_group, receptor_type, proportion) %>%
  tidyr::pivot_wider(
    names_from = receptor_type,
    values_from = proportion
  ) %>%
  dplyr::arrange(display_group)

# Save the files with clear names
write.csv(gaba_export, 
          file.path(ROOT_DIR, "gaba_receptor_stoichiometry_proportions.csv"), 
          row.names = FALSE)

write.csv(glycine_export, 
          file.path(ROOT_DIR, "glycine_receptor_stoichiometry_proportions.csv"), 
          row.names = FALSE)

write.csv(ratio_export, 
          file.path(ROOT_DIR, "gaba_glycine_receptor_proportions.csv"), 
          row.names = FALSE)
