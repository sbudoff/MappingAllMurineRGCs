library(tidyverse)
library(umap)
library(patchwork)

ROOT_DIR <- '/home/sam/Primate_Comparison/'

calcium_channels <- c("CAchannels", "Cacna1a", "Cacna1b", "Cacna1e", "Cacna1g", "Cacna1i", "Cacna2d1", "Cacna2d2", "Cacna2d3", "Cacna2d4", "Cacnb3", "Cacnb4", "Cacng2", "Cacng3", "Cacng4", "Cacng5", "Cacng7")

potassium_channels <- c("Kchannels", "Hcn1", "Hcn2", "Kcna1", "Kcna2", "Kcna4", "Kcna5", "Kcna6", "Kcnab1", "Kcnab2", "Kcnab3", "Kcnb1", "Kcnb2", "Kcnc1", "Kcnc2", "Kcnc3", "Kcnc4", "Kcnd2", "Kcnd3", "Kcnh1", "Kcnh2", "Kcnip1", "Kcnip2", "Kcnip3", "Kcnip4", "Kcnj12", "Kcnj3", "Kcnj9", "Kcnk1", "Kcnn1", "Kcnn2", "Kcnq1ot1", "Kcnq2", "Kcnq3")

sodium_channels <- c("NAchannels", "Scn1a", "Scn1b", "Scn2b", "Scn3a", "Scn3b", "Scn4b", "Scn7a", "Scn8a", "Scn9a")

glutamate_receptors <- c("Glutamate", "Gria1", "Gria2", "Gria3", "Gria4", "Grik1", "Grik5", "Grin1", "Grin2b", "Grin3a", "Grm3", "Grm4", "Grm5", "Grm6", "Grm8")

gaba_glycine_receptors <- c("GABA", "Gabra1", "Gabra2", "Gabra3", "Gabra4", "Gabrb1", "Gabrb3", "Gabrg1", "Gabrg2", "Gabrg3", "Gabrr1", "Gabrr2", "Gabrr3")# "Glra1", "Glra2", "Glrb")

mouse <-read.csv("/home/sam/FinalRGC_xenium/volcan_expression_matrix.csv")%>%
  dplyr::select(-c(X, x, y, peripheral, visual_sky, visual_ground, visual_floor,     
                   binocular, binocular_sky,binocular_ground, peripheral_sky, peripheral_ground, 
                   peripheral_floor, binocular_floor)  ) %>%
  rename(Location = areacentralis,
         Subtype = Prediction) %>%
  mutate(Location = case_when(Location == 1 ~ "Fovea",
                              Location == 0 ~ "Periphery"),
         organ = factor(Location, levels = c("Fovea", "Periphery")),
         Subtype = paste0("T", as.numeric(sub("_.*", "", Subtype))),
  ) %>%
  dplyr::select(organ, 
                any_of(calcium_channels),
                any_of(potassium_channels),
                any_of(sodium_channels),
                # any_of(glutamate_receptors),
                # any_of(gaba_glycine_receptors),
                "Rbpms"
  ) %>%
  mutate(across(!c(organ, Rbpms), ~./Rbpms)) %>%
  dplyr::select(-c(Rbpms))#, Kcnb2, Kcnj12, Kcnq1ot1))

human <- read.csv( paste0(ROOT_DIR, "human_gois.csv")) %>%
  mutate(across(!c(organ, Rbpms), ~./Rbpms)) %>%
  dplyr::select(-c(X, Rbpms))





N = 100000




human_gois <- human 
# %>%
#   slice_sample(n=N)
# Create UMAP projection
# First, remove the organ column for the UMAP calculation
# If there are issues, clean the data:
umap_input  <- human_gois %>%
  dplyr::select(-organ) %>%
  mutate(across(everything(), ~replace(., is.infinite(.), NA))) %>%
  mutate(across(everything(), ~replace(., is.na(.), mean(., na.rm = TRUE)))) %>%
  as.matrix()

# Run UMAP
set.seed(18)  # for reproducibility
umap_result <- umap(umap_input)

# Create a dataframe with UMAP coordinates and organ information
umap_df <- data.frame(
  UMAP1 = umap_result$layout[,1],
  UMAP2 = umap_result$layout[,2],
  organ = human_gois$organ
)

# Create the plot
p1 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = organ)) +
  geom_point(size = 0.1, alpha = 0.2) +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Human UMAP Projection of Gene Sets",
    x = "UMAP1",
    y = "UMAP2"
  ) +
  scale_color_brewer(palette = "Set1") 


# Create UMAP projection
mouse_gois <- mouse %>%
  slice_sample(n=10000)
# First, remove the organ column for the UMAP calculation
umap_input_m <- mouse_gois %>%
  dplyr::select(-organ) %>%
  mutate(across(everything(), ~replace(., is.infinite(.), NA))) %>%
  mutate(across(everything(), ~replace(., is.na(.), mean(., na.rm = TRUE)))) %>%
  as.matrix() 

# Run UMAP
umap_result_m <- umap(umap_input_m)


umap_df_m <- data.frame(
  UMAP1 = umap_result_m$layout[,1],
  UMAP2 = umap_result_m$layout[,2],
  organ = mouse_gois$organ
)

# Create the plot
ggplot(umap_df_m, aes(x = UMAP1, y = UMAP2, color = organ)) +
  geom_point(size = 1, alpha = 0.2) +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Human UMAP Projection of Gene Sets",
    x = "UMAP1",
    y = "UMAP2"
  ) +
  scale_color_brewer(palette = "Set1") +
  facet_grid(~organ)




mouse_full <- read.csv("/home/sam/FinalRGC_xenium/volcan_expression_matrix.csv") %>%
  rename(Location = areacentralis,
       Subtype = Prediction) %>%
  mutate(Location = case_when(Location == 1 ~ "Fovea",
                              Location == 0 ~ "Periphery"),
         Location = factor(Location, levels = c("Fovea", "Periphery")),
         Subtype = paste0("T", as.numeric(sub("_.*", "", Subtype))))

# First create a clean unification key with unique entries
unification_key <- read.csv('/home/sam/FinalRGC_xenium/RGC Subtype Unification.csv') %>%
  filter(Goetz2022_PatchSeq == 1) %>%
  mutate(Subtype = paste0("T", Tran2019_Clusters)) %>%
  distinct(Subtype, Functional_Group)

# Perform the join
mouse_full <- mouse_full %>%
  left_join(unification_key,
            by = "Subtype",
            relationship = "many-to-one")

# Create a dataframe with UMAP coordinates and organ information
# Create the UMAP dataframe with separate columns
umap_df_m <- data.frame(
  UMAP1 = umap_result_m$layout[,1],
  UMAP2 = umap_result_m$layout[,2],
  Functional_Group = mouse_full$Functional_Group,
  Location = mouse_full$Location
) %>%
  drop_na() %>%
  mutate(Funcional_Group = factor(Functional_Group))

# Create the plot
(p2 <- ggplot(umap_df_m, aes(x = UMAP1, y = UMAP2, 
                             color = Functional_Group)) +
    geom_point(size = 0.5, alpha = 0.5) +
    theme_bw() +
    theme(
      legend.position = "right",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "Mouse UMAP Projection of Gene Sets",
      x = "UMAP1",
      y = "UMAP2"
    ) +
    scale_color_brewer(palette = "Set1") +
    facet_grid(~Location)
)

p1+p2

p1+
  facet_grid(~organ) +
p2 + 
facet_grid(~organ) 


# Save UMAP coordinates for human data
human_umap_coords <- data.frame(
  UMAP1 = umap_result$layout[,1],
  UMAP2 = umap_result$layout[,2]
)

# Save UMAP coordinates for mouse data
mouse_umap_coords <- data.frame(
  UMAP1 = umap_result_m$layout[,1],
  UMAP2 = umap_result_m$layout[,2]
)

# Save to CSV files
write.csv(human_umap_coords, file = paste0(ROOT_DIR, "human_functionalGOI_umap_coordinates.csv"), row.names = FALSE)
write.csv(mouse_umap_coords, file = paste0(ROOT_DIR, "mouse_functionalGOI_umap_coordinates.csv"), row.names = FALSE)

