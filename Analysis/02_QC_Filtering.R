library(scDblFinder)
library(SingleCellExperiment)
library(DropletUtils) 
library(Matrix)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(readxl)
library(ggplot2)
library(scales)
set.seed(5)


## Perform doublet calculation within each sample
# Recommendation: best to perform this on data that has had empty drops removed but has not been filtered any further: https://github.com/plger/scDblFinder

# Create Seurat object
scRNA_data <- Read10X(data.dir = datadir)
seu = CreateSeuratObject(counts = scRNA_data)

# Plot genes/cell and UMIs/cell detected in each sample
meta <- seu@meta.data
ggplot(meta, aes(x=nCount_RNA,y=..density..)) + 
  geom_histogram(fill="white",color="black",bins=500) + 
  scale_x_continuous(breaks = seq(0, 10000, 2000), lim = c(0, 10000)) + 
  facet_wrap(~orig.ident) +
  geom_vline(aes(xintercept=500),color="red",lty="longdash") + 
  RotatedAxis() + 
  ggtitle('nCount_RNA')
## From the above plot, we determine most of our cells have sufficiently high total transcript counts to perform doublet calculations, but we do have a handful of low-read cells.
## We opt to get rid of cells with <500 total UMIs before performing doublet removal

# Filter out exceptionally low-transcript cells
select <- WhichCells(seu, expression = nCount_RNA > 500)
seu <- subset(seu, cells = select)

# Convert filtered Seurat object to SingleCellExperiment object
sce <- as.SingleCellExperiment(seu)

# Calculate doublets:
sce <- scDblFinder(sce, 
                   samples="orig.ident", #run each sample individually
                   clusters=TRUE) # use clustering approach; 
#clusters will be generated here since they have not been pre-calculated; 
#clustering approach suggested for well-segmented data
#non-clustering recommended for data with poor segregation, like a developmental trajectory: https://github.com/plger/scDblFinder
table(sce$scDblFinder.class, sce$orig.ident) # see table of predicted doublets

# Convert SingleCellExperiment back to a Seurat object that includes doublet information
seu <- as.Seurat(sce, counts = "counts", data = "logcounts")
seu <- AddMetaData(seu, metadata = sce$scDblFinder.class, col.name = "scDblFinder.class")

## At this point, doublets have NOT been removed, only calculated

#We need to calculate the percentage of mitochondrial genes expressed within each cell.
# extract mitochondrial gene Ensembl IDs from annotation file
mitoGenes <- annotKey[annotKey$ENSID %in% mitoGenes,] 
mitoGenes <- mitoGenes$FinalList
length(mitoGenes) # make sure the length is 37, the same number of mitochondrial genes found in pigs as in humans

# Now we need to calculate the percentage of mitochondrial read counts in each cell:
counts <- GetAssayData(object = seu, slot = "counts")
mitoCounts <- counts[rownames(counts) %in% mitoGenes,] 
pctMito <- ((colSums(mitoCounts))/(colSums(counts)))*100
seu <- AddMetaData(seu, pctMito, col.name = "percent_mito")
rm(counts, mitoCounts,pctMito,mitoGenes)

### Plot QC metrics
#Violin plots:
VlnPlot(seu, 
        features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), 
        split.by = 'orig.ident', 
        pt.size = 0.00, 
        ncol = 3)

#Number of genes vs. percent mitochondrial reads
meta <- seu@meta.data
ggplot(meta, aes(x=percent_mito, y=nFeature_RNA, color = orig.ident))+
  geom_point(size = 0.75) + 
  facet_wrap(~orig.ident, nrow =1)+
  theme_get()

#Number of genes vs. percent mitochondrial reads (color = # UMIs)
ggplot(meta, aes(x=percent_mito, y=nFeature_RNA, color = nCount_RNA))+
  geom_point(size = 0.75) + 
  facet_wrap(~orig.ident, nrow =1)+
  theme_get() +
  scale_colour_gradient(low = "gold", high = "red", limits=c(0, 4000), oob=squish)

#Number of genes vs. percent mitochondrial reads (color = doublet classification)
ggplot(meta, aes(x=percent_mito, y=nFeature_RNA, color = scDblFinder.class))+
  geom_point(size = 0.75) + 
  facet_wrap(~orig.ident, nrow =1)+
  theme_get()

# Number of genes vs. number reads (color = mitochondrial reads)
ggplot(meta, aes(x=nCount_RNA, y=nFeature_RNA, color = percent_mito))+
  geom_point(size = 0.75) + 
  facet_wrap(~orig.ident, nrow =1)+
  theme_get() + 
  xlim(0,20000) +
  ylim(0,4000) +
  scale_colour_gradient(low = "gold", high = "red", limits=c(0, 20), oob=squish)

# Number of genes vs. number reads (color = doublet classification)
ggplot(meta, aes(x=nCount_RNA, y=nFeature_RNA, color = scDblFinder.class))+
  geom_point(size = 0.75) + 
  facet_wrap(~orig.ident, nrow =1)+
  theme_get() + 
  xlim(0,20000) +
  ylim(0,4000) 

#Histograms with our desired thresholds shown from looking at all plots
ggplot(meta, aes(x=percent_mito,y=..density..)) + 
  geom_histogram(fill="white",color="black",bins=500) + 
  scale_y_continuous(breaks = seq(0, 1, 0.2), lim = c(0, 1)) + 
  scale_x_continuous(breaks = seq(0, 50, 5), lim = c(0, 50)) + 
  facet_wrap(~orig.ident) +
  geom_vline(aes(xintercept=25),color="red",lty="longdash") + 
  RotatedAxis() + 
  ggtitle('percent_mito')

ggplot(meta, aes(x=nFeature_RNA,y=..density..)) + 
  geom_histogram(fill="white",color="black",bins=500) + 
  scale_x_continuous(breaks = seq(0, 4000, 250), lim = c(0, 4000)) + 
  facet_wrap(~orig.ident) +
  geom_vline(aes(xintercept=400),color="red",lty="longdash") + 
  RotatedAxis() + 
  ggtitle('nFeature_RNA')

ggplot(meta, aes(x=nCount_RNA,y=..density..)) + 
  geom_histogram(fill="white",color="black",bins=500) + 
  scale_x_continuous(breaks = seq(0, 20000, 2000), lim = c(0, 20000)) + 
  facet_wrap(~orig.ident) +
  geom_vline(aes(xintercept=500),color="red",lty="longdash") + 
  RotatedAxis() + 
  ggtitle('nCount_RNA')

# Remake some xy plots with cutoffs
#Number of genes vs. percent mitochondrial reads (color = # UMIs):
ggplot(meta, aes(x=percent_mito, y=nFeature_RNA, color = nCount_RNA))+
  geom_point(size = 0.75) + 
  facet_wrap(~orig.ident, nrow =1)+
  theme_get() +
  scale_colour_gradient(low = "gold", high = "red", limits=c(0, 4000), oob=squish) + 
  geom_hline(aes(yintercept=400),color="blue",lty="longdash") + 
  geom_vline(aes(xintercept=25),color="blue",lty="longdash")

#Number of genes vs. percent mitochondrial reads (color = doublet classification):
ggplot(meta, aes(x=percent_mito, y=nFeature_RNA, color = scDblFinder.class))+
  geom_point(size = 0.75) + 
  facet_wrap(~orig.ident, nrow =1)+
  theme_get() + 
  geom_hline(aes(yintercept=400),color="blue",lty="longdash") + 
  geom_vline(aes(xintercept=25),color="blue",lty="longdash")

# Number of genes vs. number reads (color = mitochondrial reads):
ggplot(meta, aes(x=nCount_RNA, y=nFeature_RNA, color = percent_mito))+
  geom_point(size = 0.75) + 
  facet_wrap(~orig.ident, nrow =1)+
  theme_get() + 
  xlim(0,20000) +
  ylim(0,4000) +
  scale_colour_gradient(low = "gold", high = "red", limits=c(0, 20), oob=squish) + 
  geom_hline(aes(yintercept=400),color="blue",lty="longdash") + 
  geom_vline(aes(xintercept=500),color="blue",lty="longdash")

# Number of genes vs. number reads (color = doublet classification):
ggplot(meta, aes(x=nCount_RNA, y=nFeature_RNA, color = scDblFinder.class))+
  geom_point(size = 0.75) + 
  facet_wrap(~orig.ident, nrow =1)+
  theme_get() + 
  xlim(0,20000) +
  ylim(0,4000) + 
  geom_hline(aes(yintercept=400),color="blue",lty="longdash") + 
  geom_vline(aes(xintercept=500),color="blue",lty="longdash")



### After testing, filters used were:  Keep cells with <25% mitochondrial reads, >400 genes detected, >500 total reads, not called as doublets
##Filter out poor quality cells
#Identify cells passing each/every QC filter:
keepMito <- WhichCells(seu, expression = percent_mito < 25) 
keepGenes <- WhichCells(seu, expression = nFeature_RNA > 400) 
keepUMI <- WhichCells(seu, expression = nCount_RNA > 500)
keepDub <- WhichCells(seu, expression = scDblFinder.class == 'singlet') 
keep <- Reduce(intersect, list(keepMito, keepGenes, keepUMI, keepDub)) 
rm(keepMito, keepGenes, keepUMI, keepDub) 

# Create new Seurat object with only cells passing all QC filters:
seuKeep <- subset(seu, cells = keep)


### Look over QC plots to see which cells did or did not pass QC
# Number of genes vs. percent mitochondrial reads
metaKeep <- seuKeep@meta.data
ggplot(meta, aes(x=percent_mito, y=nFeature_RNA))+
  geom_point(color = 'black', size = 0.75) + 
  facet_wrap(~orig.ident, nrow =1)+
  theme_get() +
  geom_point(data=metaKeep, aes(x=percent_mito, y=nFeature_RNA), color = 'red', size = 0.75)

# Number of genes vs. number reads:
ggplot(meta, aes(x=nCount_RNA, y=nFeature_RNA))+
  geom_point(color = 'black', size = 0.75) + 
  facet_wrap(~orig.ident, nrow =1)+
  theme_get() + 
  xlim(0,20000) +
  ylim(0,4000) +
  geom_point(data=metaKeep, aes(x=nCount_RNA, y=nFeature_RNA), color = 'red', size = 0.75)

## Save counts of cells passing QC from each sample
counts <- seuKeep@assays$RNA@counts
seu <- CreateSeuratObject(counts = counts, 
                         min.cells = 1) 
