# Load necessary libraries
library(Seurat)
library(hdf5r)
library(Matrix)
library(tidyverse)

# Read the h5ad file directly with hdf5r
root <- '/home/sam/Primate_Comparison/'
h5_file <- H5File$new(paste0(root,"Human_Spatial_Gene_Analysis/HRCA_snRNA_RGC.h5ad"), "r")

# Get the sparse matrix components
x_data <- h5_file[["X/data"]][]
x_indices <- h5_file[["X/indices"]][]
x_indptr <- h5_file[["X/indptr"]][]

# Create the sparse matrix
counts_matrix <- sparseMatrix(
  i = x_indices + 1,  # Convert to 1-based indexing for R
  p = x_indptr,
  x = x_data,
  dims = c(length(h5_file[["var/_index"]][]),  # number of genes
           length(h5_file[["obs/_index"]][]))    # number of cells
)


# Get the gene names from the h5ad file
gene_names <- h5_file[["var/_index"]][]

# Convert gene names to a more manageable format if they're in a special HDF5 format
if (is.array(gene_names)) {
  gene_names <- apply(gene_names, 1, function(x) paste(rawToChar(x[x != 0]), collapse = ""))
}

# Create a proper matrix with gene names as rownames
rownames(counts_matrix) <- gene_names
colnames(counts_matrix) <- h5_file[["obs/_index"]][]

human_targets <- read.csv(paste0(root,'mouse2primate_xeniumOrthologs.csv'))%>%
  select(human_id, human_name,xenium_target) %>%
  distinct() %>%
  drop_na()
  

# Find which genes in our counts matrix match our targets
matching_genes <- rownames(counts_matrix) %in% human_targets$human_id

# Subset the counts matrix to only include these genes
counts_matrix_subset <- counts_matrix[matching_genes, ]

# Let's verify the dimensions
dim(counts_matrix_subset)

# We can also check how many of our target genes we found
sum(matching_genes)
n_missing <- sum(!human_targets$human_id %in% rownames(counts_matrix))
if(n_missing > 0) {
  cat("Missing genes:", n_missing, "\n")
  # To see which genes are missing:
  missing_genes <- human_targets$human_id[!human_targets$human_id %in% rownames(counts_matrix)]
  missing_names <- human_targets$human_name[!human_targets$human_id %in% rownames(counts_matrix)]
  missing_df <- data.frame(
    ensembl_id = missing_genes,
    gene_name = missing_names
  )
  print(missing_df)
}









# First, let's check the structure of the RGC_cluster group
h5_file[["obs/RGC_cluster"]]$ls()

# And check the organ ontology label group structure
h5_file[["obs/organ__ontology_label"]]$ls()

# Get the categories and codes for RGC clusters
rgc_categories <- h5_file[["obs/RGC_cluster/categories"]][]
rgc_codes <- h5_file[["obs/RGC_cluster/codes"]][]

# Get the categories and codes for organ labels
organ_categories <- h5_file[["obs/organ__ontology_label/categories"]][]
organ_codes <- h5_file[["obs/organ__ontology_label/codes"]][]

# Create a data frame with the mapped values
metadata_df <- data.frame(
  cell_id = h5_file[["obs/_index"]][],
  rgc_cluster = rgc_categories[rgc_codes + 1],  # Adding 1 because R is 1-based indexing
  organ = organ_categories[organ_codes + 1]
)

# Look at the unique values
cat("Unique RGC clusters:\n")
print(table(metadata_df$rgc_cluster))
cat("\nUnique organ labels:\n")
print(table(metadata_df$organ))

# Function to safely get categorical data
get_categorical_data <- function(h5_file, field_path) {
  tryCatch({
    categories <- h5_file[[paste0("obs/", field_path, "/categories")]][]
    codes <- h5_file[[paste0("obs/", field_path, "/codes")]][]
    return(categories[codes + 1])
  }, error = function(e) {
    return(NULL)
  })
}

# Get multiple metadata fields
additional_metadata <- data.frame(
  disease = get_categorical_data(h5_file, "disease"),
  disease_ontology = get_categorical_data(h5_file, "disease__ontology_label"),
  development_stage = get_categorical_data(h5_file, "developmentStage"),
  donor_age = get_categorical_data(h5_file, "donor_age"),
  sex = get_categorical_data(h5_file, "sex"),
  tissue_type = get_categorical_data(h5_file, "tissue_type")
)

# Print unique values for each field
for(col in names(additional_metadata)) {
  cat("\nUnique values in", col, ":\n")
  print(table(additional_metadata[[col]], useNA = "ifany"))
}

# Add these to your existing metadata
metadata_df <- cbind(metadata_df, additional_metadata)


# Create a logical vector of cells to keep
cells_to_keep <- colnames(counts_matrix_subset) %in% metadata_df$cell_id

# Now subset the counts matrix
counts_matrix_final <- counts_matrix_subset[, cells_to_keep]

# Let's verify the dimensions and ordering
dim(counts_matrix_final)

# Make sure the metadata is ordered the same way as our counts matrix
metadata_df_ordered <- metadata_df[match(colnames(counts_matrix_final), metadata_df$cell_id), ]

# Verify alignment
all(colnames(counts_matrix_final) == metadata_df_ordered$cell_id)

# Print summary of what we kept
cat("Final dimensions of count matrix:", dim(counts_matrix_final), "\n")
cat("\nBreakdown of cells by RGC cluster:\n")
print(table(metadata_df_ordered$rgc_cluster))
cat("\nBreakdown by organ:\n")
print(table(metadata_df_ordered$organ))
cat("\nBreakdown by sex:\n")
print(table(metadata_df_ordered$sex))
cat("\nBreakdown by development stage:\n")
print(table(metadata_df_ordered$development_stage))



# First, let's convert the sparse matrix to a regular matrix and transpose it
# so cells are rows and genes are columns
counts_df <- as.data.frame(t(as.matrix(counts_matrix_final)))

# Add the metadata columns
final_df <- cbind(
  metadata_df_ordered[, c("cell_id", "rgc_cluster", "organ", "sex", "development_stage")],
  counts_df
)

# Verify the dimensions
dim(final_df)

# Look at the first few rows and columns to verify structure
head(final_df[, 1:10])  # First 10 columns only for readability

# Verify column names
cat("\nMetadata columns:\n")
head(colnames(final_df)[1:5])
cat("\nFirst few gene ID columns:\n")
head(colnames(final_df)[6:10])


# Extract age and create adult boolean
final_df <- within(final_df, {
  # Extract age from development_stage
  age <- as.numeric(gsub("([0-9]+).*", "\\1", 
                         ifelse(development_stage == "human adult stage" | 
                                  development_stage == "80 year-old and over human stage",
                                NA, 
                                development_stage)))
  
  # Create adult boolean
  adult <- ifelse(age >= 18 | development_stage == "human adult stage" | 
                    development_stage == "80 year-old and over human stage", 1, 0)
  
  # Remove development_stage
  development_stage <- NULL
})

final_df <- final_df %>%
  mutate(age = ifelse(is.na(age), -1, age))

# Verify the change
summary(final_df$age)

# Save the dataframe to CSV
write.csv(final_df, 
          file = paste0(root, "human_spatial_scRNAseq_xeniumSubset.csv"), 
          row.names = FALSE, 
          quote = TRUE)

# Verify the file was created and check its size
file.info(paste0(root, "human_spatial_scRNAseq_xeniumSubset.csv"))$size / 1e9  # Size in GB
