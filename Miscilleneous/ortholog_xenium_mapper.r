library(dplyr)

root <- '/home/sam/Primate_Comparison/'
target_genes <- read.csv('/media/sam/Data2/CuttleNet_atlases/OriginalRetina/input_order.csv', header =  F, col.names='MouseName')

# Download human mapping
download.file("https://ftp.ensembl.org/pub/release-113/tsv/homo_sapiens/Homo_sapiens.GRCh38.113.entrez.tsv.gz",
              "human_genes.tsv.gz")

# Download orthology relationships
download.file("https://ftp.ensembl.org/pub/release-110/tsv/mus_musculus/Mus_musculus.GRCm39.110.homologies.tsv.gz",
              "mouse_human_orthologs.tsv.gz")

# Read and process files
mouse_genes <- read.delim("mouse_genes.tsv.gz")
human_genes <- read.delim("human_genes.tsv.gz")
orthologs <- read.delim(paste0(root,'Human_Spatial_Gene_Analysis/mouse2_human_mart_export.txt'))

# Clean up the orthology data
clean_orthologs <- orthologs %>%
  # Remove empty mappings and duplicates
  filter(Human.gene.stable.ID != "") %>%
  # Select just the essential columns
  select(mouse_id = Gene.stable.ID,
         human_id = Human.gene.stable.ID,
         human_name = Human.gene.name) %>%
  distinct()


mouse_symbols <- read.delim(paste0(root,'Human_Spatial_Gene_Analysis/mouse_sym_mart_export.txt'))  %>%
  rename(mouse_id = Gene.stable.ID,
         MouseName = Gene.name,
         MouseSynonym = Gene.Synonym) 

# Now filter for matches in either MouseName or MouseSynonym
target_mouse_symbols <- mouse_symbols %>%
  filter(MouseName %in% target_genes$MouseName | 
           MouseSynonym %in% target_genes$MouseName) %>%
  mutate(xenium_target = case_when(
    MouseName %in% target_genes$MouseName ~ MouseName,
    MouseSynonym %in% target_genes$MouseName ~ MouseSynonym,
    TRUE ~ NA_character_
  )) %>%
  select(xenium_target, mouse_id) %>%
  distinct()

target_mouse2human <- target_mouse_symbols %>%
  left_join(clean_orthologs, by='mouse_id')

monkey_orthologs <- read.delim(paste0(root,'Macaque_Spatial_Gene_Analysis/mouse2_macaque_mart_export.txt')) %>%
  rename(mouse_id = Gene.stable.ID,
         macaque_id = Macaque.gene.stable.ID, 
         macaque_name = Macaque.gene.name
  )

target_mouse2monkeys <- target_mouse2human %>%
  left_join(monkey_orthologs, by='mouse_id')

write.csv(target_mouse2monkeys, file = paste0(root,'mouse2primate_xeniumOrthologs.csv'))
