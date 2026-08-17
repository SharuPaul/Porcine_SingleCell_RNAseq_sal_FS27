library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
set.seed(123)


table(seu$treatment, seu$seurat_clusters)
table(seu$orig.ident, seu$seurat_clusters)

# Define colorblind-friendly colors
control_color <- "#56B4E9"  # Sky blue
salmonella_color <- "#E69F00"  # Orange

# Create a table of treatments by clusters
treatment_by_cluster <- table(seu$treatment, seu$seurat_clusters)

# Generate the stacked bar plot for proportions
ggplot(as.data.frame(as.table(prop.table(treatment_by_cluster, margin = 2))), 
       aes(x = Var2, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c(control_color, salmonella_color)) +
  labs(x = "Cluster", y = "Proportion", fill = "Treatment",
       title = "Proportion of cells across Clusters (0.3 Res)") +
  theme_minimal()

# Generate the stacked bar plot for counts
ggplot(as.data.frame(as.table(treatment_by_cluster)), aes(x = Var2, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), size = 3) + 
  scale_fill_manual(values = c(control_color, salmonella_color)) +
  labs(x = "Cluster", y = "Number", fill = "Treatment",
       title = "Cell counts across Clusters (0.3 Res)") +
  theme_minimal()


# Create a table of treatments by clusters
sample_by_cluster <- table(seu$orig.ident, seu$seurat_clusters)

# Generate the stacked bar plot for proportions
ggplot(as.data.frame(as.table(prop.table(sample_by_cluster, margin = 2))), 
       aes(x = Var2, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c('orange', 'red3', 'dodgerblue4', 'indianred1', 
                             'mediumorchid', 'cornflowerblue', 'chartreuse4', 'hotpink')) +
  labs(x = "Cluster", y = "Proportion", fill = "Treatment",
       title = "Proportion of cells per Sample across Clusters (0.3 Res)") +
  theme_minimal()

# Generate the stacked bar plot for counts
ggplot(as.data.frame(as.table(sample_by_cluster)), aes(x = Var2, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), size = 3) + 
  scale_fill_manual(values=c('orange', 'red3', 'dodgerblue4', 'indianred1', 
                             'mediumorchid', 'cornflowerblue', 'chartreuse4', 'hotpink')) +
  labs(x = "Cluster", y = "Number", fill = "Treatment",
       title = "Cell counts per Sample across Clusters (0.3 Res)") +
  theme_minimal()

