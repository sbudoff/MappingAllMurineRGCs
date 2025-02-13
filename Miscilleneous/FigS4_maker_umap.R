library(tidyverse)
library(patchwork)
library(umap)

load('/home/sam/scRNAseq/Xenium/Retina_expMatrix_clean.RData')


topNgenes_old <- read.csv('/home/sam/scRNAseq/Xenium/top_225_genes.csv', head=FALSE) %>%
  pull(V1)


Retina_interesting_genes <- read.csv('/home/sam/scRNAseq/Xenium/Xenium Gene List_final.csv') %>%
  mutate_all(~str_replace_all(., " ", "")) %>%  # Remove spaces in all columns
  mutate(gene = str_replace_all(gene, "-", "."))
topNgenes <- Retina_interesting_genes$gene

Retina_300 <- Retina_expMatrix_candidateGenes %>%
  dplyr::select(Cluster, any_of(topNgenes))

# round_df <- function(df, digits) {
#   nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
#   
#   df[,nums] <- round(df[,nums], digits = digits)
#   
#   (df)
# }

# Retina_300 <- round_df(Retina_300, digits=2)


retina_umap <- Retina_300 %>%
  dplyr::select(-Cluster) %>%
  umap()

retina_umap_coordinates <- retina_umap[['layout']] %>%
  data.frame() %>%
  mutate(Cluster = Retina_300$Cluster,
         X = X1,
         Y = X2,
         Class = case_when(substr(Cluster, 1, 2) == "AC" ~ "AC",
                           substr(Cluster, 1, 2) == "BC" ~ "BC",
                           substr(Cluster, 1, 2) == "RB" ~ "BC",
                           substr(Cluster, 1, 2) == "Ro" ~ "Ph",
                           substr(Cluster, 1, 2) == "Co" ~ "Ph",
                           substr(Cluster, 1, 2) == "MG" ~ "Glia",
                           T~ "RGC")) %>%
  dplyr::select(-X1,-X2) 

distinct_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#393b79", "#31a354", "#756bb1", "#252525", "#fed9a6", 
  "#fdd0a2", "#525252", "#084081", "#b30000", "#7f0000", 
  "#bdbdbd", "#969696", '#000075', "#f768a1", "#599861", 
  "#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", 
  "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", 
  "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", 
  "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", 
  "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", 
  '#3cb44b', '#ffe119', '#911eb4', '#46f0f0', '#f032e6', 
  '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', 
  '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', 
  "#808080"
)
alpha_bkgrnd = 0.05
alpha_frgrnd = 0.4


retina_umap_colors <- retina_umap_coordinates %>%
  group_by(Class) %>%
  dplyr::select(Cluster) %>%
  distinct() %>%
  mutate(n = 1:n(),
         Color = distinct_colors[n]) %>%
  dplyr::select(-n) %>%
  ungroup()

retina_umap_coordinates <- left_join(retina_umap_coordinates,
                                     retina_umap_colors, by=c("Cluster", "Class"))

RGC_umap <- Retina_300 %>%
  mutate(Class = case_when(substr(Cluster, 1, 2) == "AC" ~ "AC",
                           substr(Cluster, 1, 2) == "BC" ~ "BC",
                           substr(Cluster, 1, 2) == "RB" ~ "BC",
                           substr(Cluster, 1, 2) == "Ro" ~ "Ph",
                           substr(Cluster, 1, 2) == "Co" ~ "Ph",
                           substr(Cluster, 1, 2) == "MG" ~ "Glia",
                           T~ "RGC")) %>%
  filter(Class == "RGC") 

RGC_umap_labels <- RGC_umap$Cluster

RGC_umap <- RGC_umap %>%
  dplyr::select(-Cluster, -Class) %>%
  umap()

RGC_umap_coordinates <- RGC_umap[['layout']] %>%
  data.frame() %>%
  mutate(Cluster = RGC_umap_labels,
         X = X1,
         Y = X2) %>%
  dplyr::select(-X1,-X2) 




###################################
# Create a mapping for RGC labels from C[d]_[string] to T[d]
rgc_label_mapping <- retina_umap_coordinates %>%
  filter(Class == "RGC") %>%
  distinct(Cluster) %>%
  arrange(Cluster) %>%
  mutate(new_label = paste0("T", row_number()))

# Update the coordinates with new labels for RGCs
retina_umap_coordinates <- retina_umap_coordinates %>%
  left_join(rgc_label_mapping, by = "Cluster") %>%
  mutate(Cluster = ifelse(Class == "RGC", new_label, Cluster)) %>%
  select(-new_label)

# Extract the RGC-specific color mapping
rgc_colors <- retina_umap_coordinates %>%
  filter(Class == "RGC") %>%
  distinct(Cluster, Color)

# Update plotting data
RGC_umap <- retina_umap_coordinates %>%
  filter(Class == "RGC")
BC_umap <- retina_umap_coordinates %>%
  filter(Class == "BC")
AC_umap <- retina_umap_coordinates %>%
  filter(Class == "AC")

# Create plots with consistent theme
plot_theme <- theme_void() +
  theme(
    legend.position = "none",
    plot.tag = element_text(size = 12, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

BC_plot <- retina_umap_coordinates %>% 
  filter(Class != "Glia", Class != "BC") %>%
  ggplot(aes(x=X, y=Y)) +
  geom_point(alpha = alpha_bkgrnd) +
  geom_point(data = BC_umap, aes(x=X, y=Y, color = Cluster), alpha = alpha_frgrnd) +
  xlim(-15,15)+
  plot_theme +
  labs(tag = "a")

AC_plot <- retina_umap_coordinates %>% 
  filter(Class != "Glia", Class != "AC") %>%
  ggplot(aes(x=X, y=Y)) +
  geom_point(alpha = alpha_bkgrnd) +
  geom_point(data = AC_umap, aes(x=X, y=Y, color = Cluster), alpha = alpha_frgrnd) +
  xlim(-15,15)+
  plot_theme +
  labs(tag = "b")

RGC_plot <- retina_umap_coordinates %>% 
  filter(Class != "Glia", Class != "RGC") %>%
  ggplot(aes(x=X, y=Y)) +
  geom_point(alpha = alpha_bkgrnd) +
  geom_point(data = RGC_umap, aes(x=X, y=Y, color = Cluster), alpha = alpha_frgrnd) +
  xlim(-15,15)+
  plot_theme +
  labs(tag = "c")

# Update RGC-only UMAP with consistent labels
# Update RGC labels in RGC_umap_coordinates
RGC_umap_coordinates <- RGC_umap_coordinates %>%
  left_join(rgc_label_mapping, by = c("Cluster" = "Cluster")) %>%
  mutate(Cluster = ifelse(!is.na(new_label), new_label, Cluster)) %>%
  left_join(rgc_colors, by = "Cluster")

# Use single color assignment
rgc_colors <- retina_umap_colors %>% 
  filter(Class == "RGC")

RGC_only_plot <- RGC_umap_coordinates %>%
  filter(Class == "RGC") %>%
  ggplot(aes(x=X, y=Y, color = Cluster)) +
  geom_point(alpha = alpha_frgrnd) +
  scale_color_manual(values = setNames(rgc_colors$Color, rgc_colors$Cluster)) +
  ylim(-12,10) +
  plot_theme +
  theme(legend.position = "right",
        legend.title = element_blank()) +
  labs(tag = "d")

# Combine plots
combined_plot <- (BC_plot + AC_plot) / (RGC_plot + RGC_only_plot)
print(combined_plot)

# Save the plot
ggsave('/home/sam/MappingAllMurineRGCs/Data/Supplements_Outputs/S4_UMAP_retina.png',
       combined_plot, height = 10, width = 15, dpi = 300)
