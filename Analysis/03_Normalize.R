library(Matrix)
library(SeuratDisk)
library(Seurat)
library(ggplot2)
set.seed(5)


# Split Seurat object into list by sample ID
seu.list <- SplitObject(seu, split.by = "orig.ident")
seu.list

## Normalize data using SCTransform method
#Perform SCTransform on each sample from list
for (i in 1:length(seu.list)) { 
  seu.list[[i]] <- SCTransform(seu.list[[i]], 
                               return.only.var.genes = FALSE, 
                               verbose = TRUE) 
}


# Find variable features and integrate data
features <- SelectIntegrationFeatures(seu.list, verbose = TRUE) 
seu.list <- PrepSCTIntegration(seu.list,
                               anchor.features = features,
                               verbose = TRUE)
# Identify anchors for integration
anchors <- FindIntegrationAnchors(seu.list,
                                  normalization.method = "SCT", 
                                  anchor.features = features, 
                                  dims = 1:30)
# Integrate 
seu <- IntegrateData(anchors, 
                                normalization.method = "SCT", 
                                dims = 1:30)
rm(seu.list, features, anchors)


#### Calculate principle components
seu <- RunPCA(seu, npcs = 50, verbose = TRUE) 

#Quantitiatively calculate PC cutoff, which we will use to set our data dimensions in most of our subsequent analyses:
pct <- seu[["pca"]]@stdev / sum(seu[["pca"]]@stdev) * 100
cumu <- cumsum(pct) 
co1 <- which(cumu > 90 & pct < 5)[1] 
co1 
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 
co2 
pcs <- min(co1, co2) 
pcs 
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
PCdims <- 1:pcs 
rm(pct, cumu, co1, co2, pcs)

# Create UMAP
seu <- RunUMAP(seu, dims = PCdims,
               reduction = "pca", 
               min.dist = 0.5,
               spread = 0.2,
               assay = "SCT")

# Visualize UMAP
DimPlot(seu,
        reduction = 'umap',
        group.by = 'orig.ident')+ 
  ggtitle("UMAP Integrated")


# Create tSNE
seu <- RunTSNE(seu, dims = PCdims, 
               reduction = "pca", assay = "SCT") 

# Visualize tSNE
DimPlot(seu,
        reduction = 'tsne',
        group.by = 'orig.ident')+ 
  ggtitle("tSNE Integrated")

#### Add normalized/scaled data to RNA assay
dim(seu[["RNA"]]@scale.data) 
seu <- NormalizeData(seu, 
                     normalization.method = "LogNormalize", 
                     scale.factor = 10000, 
                     assay = "RNA")
seu <- ScaleData(seu, assay = "RNA")
