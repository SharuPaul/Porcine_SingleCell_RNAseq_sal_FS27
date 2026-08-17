library(Matrix)
library(SeuratDisk)
library(Seurat)
library(ggplot2)
library(clustree)
set.seed(5)

PCdims <- 1:16 # from previous step
DefaultAssay(seu) <- 'integrated'

## Set treatments
seu$treatment <- seu$orig.ident
Idents(seu) <- seu$treatment
levels(seu)
trt <- c('Control', 'Salmonella', 'Control', 'Salmonella',
         'Control', 'Salmonella', 'Control', 'Salmonella')
names(trt) <- levels(seu)
seu <- RenameIdents(seu, trt)
seu$treatment <- Idents(seu)
table(seu$orig.ident, seu$treatment)

## Cluster cells
seu <- FindNeighbors(seu, dims = PCdims, reduction = 'pca', assay = 'integrated')
resolution <- 0.3
seu <- FindClusters(seu, 
                    resolution = resolution)

DimPlot(seu, group.by = "seurat_clusters", label = TRUE) +
  ggtitle(paste("Clusters at", resolution, "Res"))

DimPlot(seu, group.by = "treatment", label = TRUE) +
  ggtitle(paste("Treatments"))

DimPlot(seu, group.by = "seurat_clusters", split.by = "treatment", label = TRUE) +
  ggtitle(paste("Clusters at", resolution, "Res"))

levels(seu$treatment)
table(seu$treatment)

# Identify cell types in each cluster
DefaultAssay(seu) <- 'SCT'
genes_1 <- c('CD19', 'CD79A', 'CD79B', 'PAX5', 'JCHAIN', 'PRDM1', 
             'IRF4', 'AIF1', 'CSF1R', 'FLT3', 'FCGR3A', 'CD14', 'CST3', 
             'CD163', 'LYZ', 'CD3E', 'CD3G', 'CD247', 'ZAP70', 'CD4', 
             'CD8A', 'CD8B', 'TRDC', 'CD2', 'NCR1', 'AHSP', 'HBM',
             'PCLAF', 'MKI67', 'PCNA', 'UBE2C')
DotPlot(seu,
        features = genes_1,
                        cols = c('gold', 'purple3')) +
  ggtitle(paste("Feature plot at", resolution, "res (SCT)")) & RotatedAxis()

## new gene list
genes_2 <- c('CD19', 'CD79A', 'PAX5', 'JCHAIN', 'IRF4', 'AIF1', 
             'SIRPA', 'FLT3', 'FCGR3A', 'CD14', 'CST3', 'CD163', 
             'LYZ', 'CD3E', 'CD3G', 'CD2', 'TRDC', 'CD5', 'CD4', 
             'CD8A', 'CD8B', 'NCR1', 'KLRK1', 'KLRB1', 'GZMA', 
             'STMN1', 'FCER1G', 'NK67', 'PCLAF')
DotPlot(seu,
        features = genes_2,
        cols = c('gold', 'purple3')) +
  ggtitle(paste("Feature plot at", resolution, "res (SCT)")) & RotatedAxis()


##################
## plot by treatment groups

Idents(seu) <- "treatment"

DotPlot(seu,
        features = genes_1,
        group.by = "seurat_clusters", idents = "Control",
        cols = c('gold', 'purple3')) +
  ggtitle(paste("Control Feature plot at", resolution, "res (SCT)")) & RotatedAxis()
DotPlot(seu,
        features = genes_1, 
        group.by = "seurat_clusters", idents = "Salmonella",
        cols = c('gold', 'purple3')) +
  ggtitle(paste("Salmonella Feature plot at", resolution, "res (SCT)")) & RotatedAxis()
DotPlot(seu,
        features = genes_2,
        group.by = "seurat_clusters", idents = "Control",
        cols = c('gold', 'purple3')) +
  ggtitle(paste("Control Feature plot at", resolution, "res (SCT)")) & RotatedAxis()
DotPlot(seu,
        features = genes_2, 
        group.by = "seurat_clusters", idents = "Salmonella",
        cols = c('gold', 'purple3')) +
  ggtitle(paste("Salmonella Feature plot at", resolution, "res (SCT)")) & RotatedAxis()

Idents(seu) <- "seurat_clusters"
########################################################

### Clustree

DefaultAssay(seu) <- 'integrated'

# Find clusters at resolution 0.1
seu <- FindClusters(seu, resolution = 0.1) 
# Add to metadata
seu[['res.0.1']] <- Idents(seu)

# Find clusters at resolution 0.2
seu <- FindClusters(seu, resolution = 0.2)
seu[['res.0.2']] <- Idents(seu)

# Find clusters at resolution 0.3
seu <- FindClusters(seu, resolution = 0.3)
seu[['res.0.3']] <- Idents(seu)

# Find clusters at resolution 0.4
seu <- FindClusters(seu, resolution = 0.4) 
seu[['res.0.4']] <- Idents(seu)

# Find clusters at resolution 0.8
seu <- FindClusters(seu, resolution = 0.8)
seu[['res.0.8']] <- Idents(seu)

# Find clusters at resolution 1.0
seu <- FindClusters(seu, resolution = 1.0) 
seu[['res.1.0']] <- Idents(seu)

# Generate clustree
tree <- clustree(seu@meta.data, prefix = "res.") # prop_filter = 0.0
tree
