library(tidyverse)
library(data.table)

root <- '/home/sam/Primate_Comparison/'

macaque_targets <- read.csv(paste0(root,'mouse2primate_xeniumOrthologs_verified.csv'))%>%
  select(macaque_id, macaque_name,xenium_target) %>%
  distinct() 

setwd(paste0(root,'Macaque_Spatial_Gene_Analysis'))
# Function to standardize cell IDs (replace dots with dashes)
standardize_cell_ids <- function(ids) {
  gsub("\\.", "-", ids)
}

# Load expression matrices
per_expression <- fread("Macaque_per_RGC_expression.txt")
fov_expression <- fread("Macaque_fov_RGC_expression.txt")

# Load coordinate files (skip first row for peripheral)
per_coords <- fread("Macaque_per_RGC_coordinates.txt", skip = 1)
fov_coords <- fread("Macaque_fov_RGC_coordinates.txt")

# Load metadata
metadata <- fread("Macaque_NN_RGC_AC_BC_HC_PR_metadata_3.txt")

# Standardize cell IDs in expression matrices
colnames(per_expression) <- c("GENE", standardize_cell_ids(colnames(per_expression)[-1]))
colnames(fov_expression) <- c("GENE", standardize_cell_ids(colnames(fov_expression)[-1]))

# Create a list containing all data
scRNA_data <- list(
  peripheral = list(
    expression = per_expression,
    coordinates = per_coords
  ),
  foveal = list(
    expression = fov_expression,
    coordinates = fov_coords
  ),
  metadata = metadata
)

scRNAseq_genes <- scRNA_data$peripheral$expression$GENE[scRNA_data$peripheral$expression$GENE %in% 
      scRNA_data$foveal$expression$GENE]

sum(macaque_targets$macaque_name %in% scRNAseq_genes)


# Let's look at the format of genes in both datasets
head(scRNAseq_genes, 10)  # Show first 10 genes from scRNA-seq
head(macaque_targets, 10)  # Show your target gene information

# Also let's see what format types we have
print("scRNA-seq gene name pattern:")
print(head(grep("^ENS", scRNAseq_genes, value=TRUE), 5))  # Any Ensembl IDs?
print(head(grep("^ENS", scRNAseq_genes, invert=TRUE, value=TRUE), 5))  # Non-Ensembl

print("\nTarget gene name pattern:")
print(head(grep("^ENS", macaque_targets$macaque_id, value=TRUE), 5))  # Ensembl IDs
print(head(macaque_targets$macaque_name, 5))  # Gene symbols

# 
# get_symbols_from_ensembl <- function(ensembl_ids, batch_size=50) {
#   library(httr)
#   library(jsonlite)
#   
#   base_url <- "https://rest.ensembl.org/lookup/id"
#   all_symbols <- data.frame(ensembl_id=character(), symbol=character(), stringsAsFactors=FALSE)
#   
#   # Process in batches
#   for(i in seq(1, length(ensembl_ids), by=batch_size)) {
#     end_idx <- min(i+batch_size-1, length(ensembl_ids))
#     batch <- ensembl_ids[i:end_idx]
#     
#     # Prepare POST request
#     response <- POST(
#       url = base_url,
#       body = list(ids = batch),
#       encode = "json",
#       add_headers(
#         "Content-Type" = "application/json",
#         "Accept" = "application/json"
#       )
#     )
#     
#     if(status_code(response) == 200) {
#       result <- fromJSON(rawToChar(response$content))
#       
#       # Extract IDs and symbols, handling missing display_names
#       batch_ids <- names(result)
#       batch_symbols <- sapply(result, function(x) {
#         if(!is.null(x$display_name)) x$display_name else NA
#       })
#       
#       # Combine into dataframe
#       if(length(batch_ids) > 0) {
#         batch_df <- data.frame(
#           ensembl_id = batch_ids,
#           symbol = batch_symbols,
#           stringsAsFactors = FALSE
#         )
#         all_symbols <- rbind(all_symbols, batch_df)
#       }
#     }
#     
#     # Print progress and response status
#     cat(sprintf("Batch %d-%d: Status %d, Found %d IDs\n", 
#                 i, end_idx, status_code(response), 
#                 length(names(result))))
#     
#     # Be nice to the API
#     Sys.sleep(1)
#   }
#   
#   return(all_symbols)
# }
# 
# # Clean macaque_targets
# clean_targets <- macaque_targets %>%
#   filter(macaque_id != "") %>%  # Remove empty IDs
#   filter(!is.na(macaque_id)) %>%  # Remove NA IDs
#   distinct(macaque_id, .keep_all = TRUE)  # Remove duplicates
# 
# # Check original vs cleaned counts
# print("Original targets:")
# print(nrow(macaque_targets))
# print("After cleaning:")
# print(nrow(clean_targets))
# 
# # Now try API with clean IDs
# id_symbol_map <- get_symbols_from_ensembl(clean_targets$macaque_id)
# 
# # Analyze response coverage
# print("Response analysis:")
# print(sprintf("Total mappings found: %d out of %d", nrow(id_symbol_map), nrow(clean_targets)))
# print("Proportion of NAs in symbols:")
# print(sum(is.na(id_symbol_map$symbol)) / nrow(id_symbol_map))

# Get overlapping genes between foveal and peripheral data
scRNAseq_genes <- per_expression$GENE[per_expression$GENE %in% fov_expression$GENE]
# testscRNAseq_genes <- fov_expression$GENE[fov_expression$GENE %in% per_expression$GENE]

# Check overlap between scRNA-seq genes and our mapped symbols
# mapped_overlap <- id_symbol_map$symbol %in% scRNAseq_genes

# Print summary statistics
# print("Overlap analysis:")
# print(sprintf("Total mapped symbols: %d", nrow(id_symbol_map)))
# print(sprintf("Symbols found in scRNA-seq: %d", sum(mapped_overlap, na.rm=TRUE)))
# 
# # Look at examples of matches and non-matches
# print("\nSample of matching genes:")
# head(id_symbol_map$symbol[mapped_overlap])
# 
# print("\nSample of non-matching genes:")
# head(id_symbol_map$symbol[!mapped_overlap])
# 
# # Check for case sensitivity issues
# case_sensitive_check <- id_symbol_map$symbol %in% toupper(scRNAseq_genes)
# if(sum(case_sensitive_check) > sum(mapped_overlap)) {
#   print("\nWarning: Case sensitivity might be affecting matching")
#   print(sprintf("Additional matches when ignoring case: %d", 
#                 sum(case_sensitive_check) - sum(mapped_overlap)))
# }
# 




# 1. Update target CSV with matched symbols
original_targets <- read.csv(paste0(root,'mouse2primate_xeniumOrthologs_verified.csv'))
# Merge with our id_symbol_map to update macaque_name
# updated_targets <- original_targets %>%
#   left_join(id_symbol_map, by=c("macaque_id" = "ensembl_id")) %>%
#   mutate(macaque_name = symbol) %>%
#   select(-symbol)  # Remove the temporary symbol column
# # Save updated targets
# write.csv(updated_targets, paste0(root,'mouse2primate_xeniumOrthologs.csv'), row.names=FALSE)

# 2. Subset expression matrices with target symbols and update gene names
# Create mapping dictionary
gene_name_map <- setNames(original_targets$xenium_target, original_targets$macaque_name)


# Function to handle duplicate names
make_unique_names <- function(names) {
  counts <- table(names)
  duplicates <- names[duplicated(names)]
  
  if(length(duplicates) > 0) {
    for(dup in unique(duplicates)) {
      # Find all instances of this name
      idx <- which(names == dup)
      # Add numbered suffix starting from 1
      names[idx] <- paste0(dup, ".", seq_along(idx))
    }
  }
  return(names)
}

# 2. Subset expression matrices with target symbols and update gene names
# Create mapping dictionary
gene_name_map <- setNames(original_targets$xenium_target, original_targets$macaque_name)


# Subset and prepare foveal
fov_subset <- fov_expression %>%
  select(GENE, all_of(intersect(colnames(.), fov_coords$NAME))) %>%
  filter(GENE %in% original_targets$macaque_name)

# Create new gene name column first
fov_subset$new_gene <- gene_name_map[fov_subset$GENE]
# Handle duplicates
fov_subset$new_gene <- make_unique_names(fov_subset$new_gene)
# Now transpose with gene names as columns
fov_t <- data.frame(t(fov_subset %>% select(-GENE, -new_gene)))
colnames(fov_t) <- fov_subset$new_gene

# Repeat for peripheral
per_subset <- per_expression %>%
  select(GENE, all_of(intersect(colnames(.), per_coords$TYPE))) %>%
  filter(GENE %in% original_targets$macaque_name)

per_subset$new_gene <- gene_name_map[per_subset$GENE]
per_subset$new_gene <- make_unique_names(per_subset$new_gene)
per_t <- data.frame(t(per_subset %>% select(-GENE, -new_gene)))
colnames(per_t) <- per_subset$new_gene

# Add cell IDs as explicit column before join
fov_t$cell_id <- rownames(fov_t)
per_t$cell_id <- rownames(per_t)

# Add Location and Subtype for foveal data
fov_final <- fov_t %>%
  left_join(fov_coords %>% select(NAME, Cluster), by = c("cell_id" = "NAME")) %>%
  mutate(Location = "fovea") %>%
  rename(Subtype = Cluster)

# Add Location and Subtype for peripheral data
per_final <- per_t %>%
  left_join(per_coords %>% select(TYPE, group), by = c("cell_id" = "TYPE")) %>%
  mutate(Location = "peripheral") %>%
  rename(Subtype = group)

# Check results
print("Foveal data:")
print(head(select(fov_final, cell_id, Location, Subtype)))
print("Peripheral data:")
print(head(select(per_final, cell_id, Location, Subtype)))



# 6. Combine and save
combined_data <- bind_rows(
  fov_final,
  per_final
) %>%
  column_to_rownames("cell_id") 

#drop all na columns
combined_data[is.na(combined_data)] <- 0



# Save final dataset
write.csv(combined_data, 
          paste0(root, 'macaque_spatial_scRNAseq_xeniumSubset.csv'), 
          row.names=TRUE)  # Use row.names=TRUE to keep cell IDs as row names

# Print summary stats
print("Dataset summary:")
print(sprintf("Total cells: %d", nrow(combined_data)))
print(sprintf("Number of genes: %d", ncol(combined_data) - 2))  # -2 for Location and Subtype
print("\nCells per location:")
print(table(combined_data$Location))
print("\nSubtypes per location:")
print(table(combined_data$Location, combined_data$Subtype))

