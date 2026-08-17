library(Seurat)
library(SeuratDisk)
library(miloR)
library(SingleCellExperiment)
library(data.table)
library(ggplot2)
library(readxl)
library(dplyr)
library(writexl)
set.seed(123)


## Order clusters
Idents(seu) <- seu$seurat_clusters
cluster_order <- c("0", "4", "5", "9", "12", "16", "3", "10", "17", "1", 
                   "7", "8", "6", "14", "2", "13", "11", "18", "15")
Idents(seu) <- factor(Idents(seu), levels = cluster_order)

## Cell types named
CellTypes <- c('Mono_1', 'Mono_2', 'Mono_3', 'Mono_4', 
               'pDC', 'cDC', 'Bcell_1', 'Bcell_2', 'ASC', 
               'abT_CD4_1', 'abT_CD4_2', 'abT_NK_1', 'abT_NK_2',
               'abT_prolif', 'gdT_CD2n', 'gdT_CD2p', 
               'Mixed_Tcells_1', 'Mixed_Tcells_2', 'Mixed_T_Mono')

length(CellTypes)
names(CellTypes) <- cluster_order 
seu <- RenameIdents(seu, CellTypes)
seu$celltypes <- Idents(seu)

CellType_plot <- DimPlot(seu, label = TRUE) +
  ggtitle("Cell Types at 0.3 Res")
CellType_plot


## Milo object
milo <- as.SingleCellExperiment(seu, assay = 'SCT')
milo <- Milo(milo)


# PCs # Calculate PCdims: 1:pcs
k_param = 40

#find nearest neighbors
milo <- buildGraph(milo, k = k_param, d = pcs, reduced.dim = "PCA")

### Create & visualize cell neighborhoods
#generate the neighborhoods
milo <- makeNhoods(milo,
                   prop = 0.2, 
                   k = k_param, 
                   d = pcs, 
                   refined = TRUE, 
                   refinement_scheme="graph")

plotNhoodSizeHist(milo) 
milo <- buildNhoodGraph(milo)
plotNhoodGraph(milo, layout = 'UMAP')

### Count cells in each neighborhood
milo <- countCells(milo, meta.data = data.frame(colData(milo)), sample="orig.ident")
head(nhoodCounts(milo))

### Create experimental design
milo_design <- data.frame(colData(milo))[,c("orig.ident", "treatment")]
milo_design <- distinct(milo_design)
rownames(milo_design) <- milo_design$orig.ident
milo_design

### Perform DA testing
da_results <- testNhoods(milo,
                         design = ~ treatment,
                         design.df = milo_design,
                         fdr.weighting = 'graph-overlap', 
                         reduced.dim = 'PCA')
head(da_results)

#Make a histogram of p-values found across cell neighborhoods:
ggplot(da_results, aes(PValue)) + geom_histogram(bins=100)

#Make a volcano plot of DA. Each dot is one cell neighborhood:
ggplot(da_results, aes(logFC, -log10(FDR))) +
  geom_point() +
  geom_hline(yintercept = 1)

#Overlay logFC scores onto cell neighborhood central coordinates on t-SNE & UMAP plots:
plotNhoodGraphDA(milo, da_results, layout="UMAP",alpha=0.1)

#And we can also look at all cell neighborhoods on a bee swarm plot:
plotDAbeeswarm(da_results, alpha = 0.1)


### Annotate cell neighborhoods
da_results <- annotateNhoods(milo, da_results, coldata_col = "celltypes")
head(da_results)
#Create a histogram to look at the largest percentages for a single cell type within each cell neighborhood:

ggplot(da_results, aes(celltypes_fraction)) + geom_histogram(bins=100)

da_results$celltypes <- ifelse(da_results$celltypes_fraction < 0.7, "Mixed", da_results$celltypes)

### Plot DA across annotated cell neighborhoods:
table(da_results$celltypes)
#Make a bee swarm plot:
da_results$celltypes <- factor(
  da_results$celltypes,
  levels = rev(c(CellTypes, "Mixed")))
plotDAbeeswarm(da_results, group.by = "celltypes", alpha = 0.1)
