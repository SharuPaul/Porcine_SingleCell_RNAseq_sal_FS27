library(ggplot2)
library(Seurat)
library(SeuratDisk)
library(dplyr)
set.seed(123)


DefaultAssay(seu) <- 'SCT'
table(seu$seurat_clusters)

seu <- PrepSCTFindMarkers(seu, assay = "SCT", verbose = TRUE)

## Find All Markers for 1 vs all cluster comparisons
all_markers <- FindAllMarkers(seu, pseudocount.use=1e-9,
                              assay = "SCT",
                              group.by = "seurat_clusters",
                              min.pct = 0.25, logfc.threshold = 0.25)

significant_markers <- all_markers %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))

## Get top markers per cluster based on positive avg_log2FC
top_markers_per_cluster <- significant_markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) %>%
  ungroup()

#####################################################################
#####################################################################
## Find Markers for treatment comparison within cluster

## For wilcoxon test on SCT assay
DefaultAssay(seu) <- 'SCT'
seu <- PrepSCTFindMarkers(seu, assay = "SCT", verbose = TRUE)

## For wilcoxon test on RNA assay or
## For negbinom test on RNA assay
DefaultAssay(seu) <- 'RNA'

table(seu$orig.ident, seu$treatment)

clusters <- levels(seu$seurat_clusters)
de_cluster <- setNames(vector("list", length(clusters)), clusters)
Idents(seu) <- seu$seurat_clusters

seu$clus_trt <- paste(seu$seurat_clusters, seu$treatment, sep = '_')  
Idents(seu) <- seu$clus_trt


# Initialize a list to store DE results
de_list <- list()

# Loop through each cluster and perform differential expression analysis
for (i in 0:18) {
  de <- FindMarkers(seu, pseudocount.use=5e-4, 
                    ident.1 = paste0(i, "_Salmonella"), ident.2 = paste0(i, "_Control"),
                    test.use = "wilcox")

  de$cluster <- as.character(i)
  de$gene <- rownames(de)
  de_list[[i + 1]] <- de
}


de <- do.call(rbind, de_list)
colnames(de) <- c('p_val', 'avg_log2FC', 'pct_Salmonella', 
                  'pct_Control', 'p_val_adj', 'cluster', 'gene')

head(de)


# Function to plot DE results
plot_DE_comparison <- function(de_results, comparison) {
  ggplot(de_results, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = p_val_adj < 0.05 & abs(avg_log2FC) > 1), alpha = 0.5) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey50"),
                       name = "Significant DE") +
    labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value", 
         title = paste("Differentially Expressed Genes (0.3res, Wilcox, pseudo 5e_4, SCT) -", comparison),
         subtitle = "Red dots: p_val_adj < 0.05 & abs(avg_log2FC) > 1") +
    theme_minimal() +
    coord_cartesian(xlim = c(-2, 2)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black")
}

for (cl in 0:18) {
  p <- plot_DE_comparison(
    subset(de, cluster == as.character(cl)),
    paste("Cluster", cl)
  )
  print(p)

}

de$cluster <- factor(de$cluster, levels = as.character(0:18))
ggplot(de, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = p_val_adj < 0.05 & abs(avg_log2FC) > 1), alpha = 0.5) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey50"),
                     name = "Significant DE") +
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value",
       title = "Differentially Expressed Genes by Cluster (0.3res, Wilcox, pseudo 5e_4, SCT)",
       subtitle = "Red = p_adj < 0.05 & |log2FC| > 1") +
  theme_minimal() +
  coord_cartesian(xlim = c(-2, 2)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~ cluster, ncol = 5)
