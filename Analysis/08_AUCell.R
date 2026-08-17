library(Matrix)
library(Seurat)
library(SeuratDisk)
library(readxl)
library(AUCell)
library(ggplot2)
library(BiocParallel)
library(dplyr)
library(scales)
set.seed(123)

## AUCell analysis
## Done with complete dataset or subset to just control or sal objects

## Exclude cluster 18 from seu object
seu <- subset(seu, subset=seurat_clusters != "18")

## DGE Crl02 v Sal02  
##  https://github.com/jwiarda/FS27_WholeBlood_BulkRNAseq/blob/main/Results_DGE_DESeq2/Results_DGE_Crl02vSal02.xlsx
tbl <- read_excel("Results_DGE_Crl02vSal02.xlsx",
                  sheet = 1, range = NULL)[, c("GeneName","log2FoldChange", "padj")]

## Subset
ord <- order(tbl$log2FoldChange, decreasing = TRUE, na.last = NA)
keep <- tbl$log2FoldChange[ord] > 1 & tbl$padj[ord] < 0.05 & !is.na(tbl$padj[ord])
sig_genes <- tbl$GeneName[ ord[keep] ]

## Get data from seu object
mat <- GetAssayData(seu, assay = "RNA", layer = "counts")
gene_set <- intersect(sig_genes, rownames(mat))

## build ranking in AUCell
rankings  <- AUCell_buildRankings(mat, 
                                  splitByBlocks = TRUE,          # TRUE needed to use BPPARAM
                                  BPPARAM = BiocParallel::MulticoreParam(16))
aucMaxRank <- ceiling(0.05 * nrow(rankings))  # 0.05 means top 5%

auc <- AUCell_calcAUC(list(genes = gene_set), 
                      rankings, 
                      aucMaxRank = aucMaxRank)

## AUC score to seurat object
seu$top_AUC <- as.numeric(getAUC(auc)["genes", ])

## Plot
FeaturePlot(seu, features = "top_AUC", reduction = "umap")+
  ggtitle("AUC scores (DEGs Day2 Con vs Sal)") 

# hist(seu$top_AUC, breaks=50) 

## Pick threshold
thr_info <- AUCell_exploreThresholds(auc["genes", ], plotHist = TRUE)
thr <- thr_info[[1]]$aucThr$selected


seu$AUC_thres <- seu$top_AUC > thr
DimPlot(seu, group.by="AUC_thres") +
  ggtitle(paste0("Cells above AUC threshold ", round(thr, 3), " (DEGs Day2 Con vs Sal)"))

#################
#### Make UMAP with colors
# get UMAP coordinates
umap_df <- as.data.frame(Embeddings(seu, "umap"))
colnames(umap_df) <- c("UMAP_1", "UMAP_2")

# add AUC values
umap_df$AUC <- seu$top_AUC


min_auc <- min(umap_df$AUC, na.rm = TRUE)
max_auc <- max(umap_df$AUC, na.rm = TRUE)


umap_df$AUC_scaled <- ifelse(
  umap_df$AUC <= thr,
  rescale(umap_df$AUC, to = c(-1, 0), from = c(min_auc, thr)),
  rescale(umap_df$AUC, to = c(0, 1), from = c(thr, max_auc))
)

## Plot
ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = AUC_scaled)) +
  geom_point(size = 0.8) +
  scale_color_gradient2(
    name = "AUC",
    low = "darkblue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1),
    breaks = c(-1, 0, 1),
    labels = c("low", "threshold", "high")
  ) +
  theme_classic() +
  labs(x = NULL, y = NULL)
#######

seu$AUC_z <- scale(seu$AUC_thres)[,1]
cluster_means <- tapply(seu$AUC_z, Idents(seu), mean, na.rm = TRUE)
cluster_means <- sort(cluster_means, decreasing = TRUE)
cluster_means


# active cells
pct_active <- mean(seu$top_AUC > thr) * 100
pct_active
##################################
df <- data.frame(
  cluster = as.character(Idents(seu)),
  top_AUC = seu$top_AUC,
  AUC_thres = seu$AUC_thres
) %>%
  group_by(cluster) %>%
  summarise(
    pct_pos = mean(AUC_thres, na.rm = TRUE) * 100,
    mean_auc = mean(top_AUC, na.rm = TRUE),
    .groups = "drop"
  )
df$cluster <- factor(df$cluster, levels = df$cluster[order(df$mean_auc, decreasing = TRUE)])
ggplot(df, aes(x = cluster, y = mean_auc, size = pct_pos)) +
  geom_point() +
  labs(x = "Cluster", y = "Mean AUC", size = "% AUC+ cells") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  ggtitle(paste0("Mean AUC per cluster is All thres ", round(thr, 3), " (DEGs Day2 Con vs Sal)"))

##################################
## Second threshold
thr2 <- 0.16 #0.06
seu$AUC_thres2 <- seu$top_AUC > thr2
DimPlot(seu, group.by="AUC_thres2") +
  ggtitle(paste0("Cells above AUC threshold ", round(thr2, 3), " (DEGs Day2 Con vs Sal)"))

#################
#### Make UMAP with colors

## Scaled plot
umap_df$AUC_scaled2 <- ifelse(
  umap_df$AUC <= thr2,
  rescale(umap_df$AUC, to = c(-1, 0), from = c(min_auc, thr2)),
  rescale(umap_df$AUC, to = c(0, 1), from = c(thr2, max_auc))
)

ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = AUC_scaled2)) +
  geom_point(size = 0.8) +
  scale_color_gradient2(
    name = "AUC",
    low = "darkblue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1),
    breaks = c(-1, 0, 1),
    labels = c("low", "threshold", "high")
  ) +
  theme_classic() +
  labs(x = NULL, y = NULL)
#######


seu$AUC_z <- scale(seu$AUC_thres2)[,1]
cluster_means <- tapply(seu$AUC_z, Idents(seu), mean, na.rm = TRUE)
cluster_means <- sort(cluster_means, decreasing = TRUE)
cluster_means


# active cells
pct_active <- mean(seu$top_AUC > thr2) * 100
pct_active
##################################
df <- data.frame(
  cluster = as.character(Idents(seu)),
  top_AUC = seu$top_AUC,
  AUC_thres = seu$AUC_thres2
) %>%
  group_by(cluster) %>%
  summarise(
    pct_pos = mean(AUC_thres, na.rm = TRUE) * 100,
    mean_auc = mean(top_AUC, na.rm = TRUE),
    .groups = "drop"
  )
df$cluster <- factor(df$cluster, levels = df$cluster[order(df$mean_auc, decreasing = TRUE)])
ggplot(df, aes(x = cluster, y = mean_auc, size = pct_pos)) +
  geom_point() +
  labs(x = "Cluster", y = "Mean AUC", size = "% AUC+ cells") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  ggtitle(paste0("Mean AUC per cluster is All thres ", round(thr2, 3), " (DEGs Day2 Con vs Sal)"))


########## Histogram with both thresholds
# build histogram
h <- hist(seu$top_AUC, breaks = 50, plot = FALSE)

# color bars depending on whether the bin is above thr1
bar_cols <- ifelse(h$mids > thr, "#2166AC", "#C6DBEF")

# plot histogram
plot(h,
  col = bar_cols,
  border = "white",
  xlab = "AUC score",
  ylab = "Number of cells"
)

# add threshold lines
abline(v = thr, col = "blue", lwd = 3, lty = 1)
abline(v = thr2, col = "red",  lwd = 2, lty = 2)

# add threshold labels
usr <- par("usr")

text(
  x = thr,
  y = usr[4] * 0.60,
  labels = paste0("AUC > ", round(thr, 3)),
  col = "blue",
  pos = 4,
  cex = 1.1,
  xpd = TRUE
)

text(
  x = thr2,
  y = usr[4] * 0.60,
  labels = paste0("AUC > ", round(thr2, 3)),
  col = "red",
  pos = 4,
  cex = 1,
  xpd = TRUE
)

