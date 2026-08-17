library(dplyr)
library(tidyr)
library(ComplexUpset)
library(readxl)
library(ggplot2)
set.seed(123)

# Read DEG files
all_de <- read.csv("pseudobulk_shrinkage_0.3res_Salmonella_vs_Control.csv")
wb <- read_excel("Results_DGE_Crl02vSal02.xlsx", sheet = 1, range = NULL)
[, c("GeneName", "baseMean", "log2FoldChange", "padj")]

mono_clusters <- c("0", "4", "5", "9")
mono <- all_de[all_de$cluster %in% mono_clusters, ]
rm(all_de)

## Cluster names
mono_cluster_names <- c(
  "0" = "Mono_1",
  "4" = "Mono_2",
  "5" = "Mono_3",
  "9" = "Mono_4"
)

# #####################################
# Use just Mono
monos <- mono %>%
  filter(padj < 0.05) %>%
  mutate(cluster = unname(mono_cluster_names[as.character(cluster)]))
 
# Use just Mono UP
mono_up <- mono %>%
  filter(padj < 0.05, log2FoldChange >= 1) %>%
  mutate(cluster = unname(mono_cluster_names[as.character(cluster)]))

# Prep data
mono_list_up <- mono_up %>%
  dplyr::select(gene, cluster) %>%
  distinct()

deg_binary_mono_up <- mono_list_up %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = cluster,
    values_from = value,
    values_fill = 0
  )

contrast_cols_mono_up <- setdiff(colnames(deg_binary_mono_up), "gene")
deg_binary_mono_up[contrast_cols_mono_up] <- lapply(deg_binary_mono_up[contrast_cols_mono_up], as.logical)

colSums(as.data.frame(deg_binary_mono_up[contrast_cols_mono_up]))

# Plot
upset(
  deg_binary_mono_up,
  intersect = contrast_cols_mono_up,
  # sort_sets = FALSE,
  name = "Monocyte UP (log2FC>1) DEGs comparison",
  base_annotations = list(
    "Intersection size" = intersection_size(text = list(size = 2.8))
  ),
  set_sizes = upset_set_size(),
  width_ratio = 0.2
)
#####################################
# Use just Mono  DOWN
mono_down <- mono %>%
  filter(padj < 0.05, log2FoldChange =< -1) %>%
  mutate(cluster = unname(mono_cluster_names[as.character(cluster)]))

# Prep data
mono_list_down <- mono_down %>%
  dplyr::select(gene, cluster) %>%
  distinct()

deg_binary_mono_down <- mono_list_down %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = cluster,
    values_from = value,
    values_fill = 0
  )

contrast_cols_mono_down <- setdiff(colnames(deg_binary_mono_down), "gene")
deg_binary_mono_down[contrast_cols_mono_down] <- lapply(deg_binary_mono_down[contrast_cols_mono_down], as.logical)

colSums(as.data.frame(deg_binary_mono_down[contrast_cols_mono_down]))

# Plot
upset(
  deg_binary_mono_down,
  intersect = contrast_cols_mono_down,
  name = "Monocyte DOWN (log2FC<-1) DEGs comparison",
  base_annotations = list(
    "Intersection size" = intersection_size(text = list(size = 2.8))
  ),
  set_sizes = upset_set_size(),
  width_ratio = 0.2
)


#####################################
# Use WB list and plot both
wb_filtered <- wb %>%
  rename(gene = GeneName) %>%
  filter(padj < 0.05, log2FoldChange >= 1) %>%
  mutate(cluster = "WB")

# prep
wb_list <- bind_rows(monos, wb_filtered) %>%
  dplyr::select(gene, cluster) %>%
  distinct()

deg_binary_wb <- wb_list %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = cluster,
    values_from = value,
    values_fill = 0
  )

# Keep only WB genes
deg_binary_wb <- deg_binary_wb %>%
  filter(gene %in% wb_filtered$gene)

contrast_cols_wb <- setdiff(colnames(deg_binary_wb), "gene")
deg_binary_wb[contrast_cols_wb] <- lapply(deg_binary_wb[contrast_cols_wb], as.logical)

colSums(as.data.frame(deg_binary_wb[contrast_cols_wb]))

# Plot
upset(
  deg_binary_wb,
  intersect = contrast_cols_wb,
  name = "WB genes, Mono DEG comparison",
  base_annotations = list(
    "Intersection size" = intersection_size(text = list(size = 2.8))
  ),
  set_sizes = upset_set_size(),
  width_ratio = 0.2
)
