library(tidyverse)


rgcs <-read.csv("/home/sam/FinalRGC_xenium/volcan_expression_matrix.csv") %>%
  select(-all_of(c("x", "y" , "ipsi" , "binocular","peripheral","visual_sky",
             "visual_ground","visual_floor","binocular_sky","binocular_ground", 
              "peripheral_sky","peripheral_ground", "peripheral_floor",
             "binocular_floor","areacentralis", "X"       
             ))
         ) %>%
  mutate(Class = 'RGC')

all_RGC_layer <- read.csv('/media/sam/Data2/baysor_rbpms_consolidated/OldModel/all_retinas_prediction_expmat.csv')  %>%
  filter(Class != "RGC") %>% # remove the unfiltered and unmerged RGCs
  select(all_of(names(rgcs))) %>%
  rbind(rgcs) # Add back in the cleaned up RGC dataset

print(all_RGC_layer %>%
  group_by(Class) %>%
  summarise(n()))


mean_matrix <- all_RGC_layer %>%
  group_by(Class, Prediction) %>%
  summarize_all(~mean(., na.rm = T)) %>%
  ungroup() %>%
  mutate(DataSource = 'Xenium')

write.csv(mean_matrix, file="/home/sam/MappingAllMurineRGCs/Data/Supplements_Outputs/average_xenium_expression_all_cells.csv",row.names = FALSE)

scRNAseq <- read.csv('/home/sam/scRNAseq/ComprehensiveAtlas/retina_average_expression.csv') %>%
  select(-c(ct_idx, atlas, X)) %>%
  rename(Prediction = cell_type) 

name_mapping <- c(
  # RGC types - add leading 0
  "1_W3D1.1" = "01_W3D1.1",
  "2_W3D1.2" = "02_W3D1.2",
  "3_FminiON" = "03_FminiON",
  "4_FminiOFF" = "04_FminiOFF",
  "5_J-RGC" = "05_J-RGC",
  "6_W3B" = "06_W3B",
  "7_Novel" = "07_Novel",
  "8_Novel" = "08_Novel",
  "9_Tbr1_Novel" = "09_Tbr1_Novel",
  
  # Cell type prefix changes
  "BC1A" = "0BC1A",
  "BC1B" = "0BC1B",
  "BC2" = "0BC2",
  "BC3A" = "0BC3A",
  "BC3B" = "0BC3B",
  "BC4" = "0BC4",
  "BC5A" = "0BC5A (Cone Bipolar cell 5A)",
  "BC5B" = "0BC5B",
  "BC5C" = "0BC5C",
  "BC5D" = "0BC5D",
  "BC6" = "0BC6",
  "BC7" = "0BC7 (Cone Bipolar cell 7)",
  "BC8" = "0BC8/9 (mixture of BC8 and BC9)",
  "BC9" = "0BC8/9 (mixture of BC8 and BC9)",
  
  # Full name changes
  "Cone" = "0Cone Photoreceptors",
  "Rod" = "0Rod Photoreceptors",
  "Endotheliocyte" = "0Endothelial",
  "HC" = "0Horizontal Cell",
  "MG" = "0MG (Mueller Glia)",
  "Microglia" = "0Microglia",
  "Pericyte" = "0Pericyte",
  "RBC" = "0RBC (Rod Bipolar cell)"
)

# Apply the mapping
scRNAseq$Prediction_aligned <- as.character(scRNAseq$Prediction)
for(old_name in names(name_mapping)) {
  scRNAseq$Prediction_aligned[scRNAseq$Prediction == old_name] <- name_mapping[old_name]
}


scRNAseq <- scRNAseq %>%
  mutate(Prediction = Prediction_aligned,
         DataSource = 'scRNAseq') %>%
  select(-Prediction_aligned) %>%
  left_join(select(mean_matrix, Prediction, Class), by = "Prediction") %>%
  drop_na()

final_df <- scRNAseq %>%
  rbind(mean_matrix)
write.csv(final_df, file="/home/sam/MappingAllMurineRGCs/Data/Supplements_Outputs/average_xeniumANDscrnaseq_expression_all_cells.csv",row.names = FALSE)


# 1. Calculate z-scores for gene columns
gene_columns <- final_df %>% 
  select(-Class, -Prediction, -DataSource) %>% 
  names()

z_score_matrix <- final_df %>%
  group_by(Class, DataSource) %>%
  mutate(across(all_of(gene_columns), scale)) %>%
  ungroup()

# 2. Reshape to long format
long_data <- z_score_matrix %>%
  pivot_longer(
    cols = all_of(gene_columns),
    names_to = "Gene",
    values_to = "Z_score"
  )

# 3. Calculate mean z-scores for RGC class to order genes - MODIFIED THIS PART
gene_order <- long_data %>%
  filter(Class == "RGC", DataSource == 'scRNAseq') %>%
  group_by(Gene) %>%
  summarize(z = mean(Z_score)) %>%
  arrange(desc(z)) %>%  # This ensures highest to lowest ordering
  pull(Gene)

# Make Gene a factor with levels in the correct order
long_data$Gene <- factor(long_data$Gene, levels = gene_order)

# 4. Create factor for ordering classes
class_order <- c("RGC", "AC", "BC", "Ph", "Other")
long_data$Class <- factor(long_data$Class, levels = class_order)

# 5. Create the heatmap with subtypes
ggplot(long_data, aes(x = Prediction, y = Gene, fill = Z_score)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "purple",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-2.5, 7.5)
  ) +
  facet_grid(. ~ Class*DataSource, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    strip.text = element_text(size = 10),
    panel.spacing = unit(0.1, "lines")
  ) +
  labs(
    title = "Gene Expression Z-scores by Cell Class",
    x = "Cell Subtype",
    y = "Genes",
    fill = "Z-score"
  )
