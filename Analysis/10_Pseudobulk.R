library(DESeq2)
library(edgeR)
library(ggplot2)
library(Seurat)
library(SeuratDisk)
library(apeglm)
library(dplyr)
library(tibble)
library(scales)
set.seed(123)

DefaultAssay(seu) <- "RNA" 

seu$orig.ident <- factor(seu$orig.ident)
seu$treatment <- factor(seu$treatment)

Idents(seu) <- "seurat_clusters" 
levels(Idents(seu)) 

## Pseudobulk
pseudo <- AggregateExpression(
  object = seu, 
  assays = "RNA",
  slot = "counts", 
  group.by = c("treatment", "orig.ident", "seurat_clusters"),
  return.seurat = FALSE)

##cpm filter
cnt <- pseudo$RNA
nrow(cnt)
keep_genes <- rowSums(cpm(cnt) > 0.25) >= 3
sum(keep_genes)
cnt <- cnt[keep_genes, , drop = FALSE]

## build metadata
cn    <- colnames(cnt)
parts <- do.call(rbind, strsplit(cn, "_"))
meta_from_names <- as.data.frame(parts, stringsAsFactors = FALSE)
colnames(meta_from_names) <- c("treatment", "orig.ident", "seurat_clusters")
rownames(meta_from_names) <- cn

meta_from_names$treatment       <- factor(meta_from_names$treatment)
meta_from_names$orig.ident      <- factor(meta_from_names$orig.ident)
meta_from_names$seurat_clusters <- factor(meta_from_names$seurat_clusters)

coldata <- meta_from_names


# DESeq2 DE for each cluster
all_clusters <- sort(levels(coldata$seurat_clusters))
results_list <- list()

cond_ref  <- "Control"
cond_test <- "Salmonella"

for (cl in all_clusters) {
  message("Cluster: ", cl)
  
  idx_cl  <- coldata$seurat_clusters == cl
  cnt_cl  <- cnt[, idx_cl, drop = FALSE]
  meta_cl <- coldata[idx_cl, , drop = FALSE]
  meta_cl$treatment <- relevel(meta_cl$treatment, ref = cond_ref)
  
  # need both conditions present
  if (!all(c(cond_ref, cond_test) %in% meta_cl$treatment)) {
    message("Skipping: one condition missing")
    next
  }
  
  # require at least 3 pseudobulk samples per group
  tab_treat <- table(meta_cl$treatment)
  n_ref  <- tab_treat[cond_ref]
  n_test <- tab_treat[cond_test]
  if (n_ref < 3 || n_test < 3) {
    message("Skipping: not enough replicates (", n_ref, " vs ", n_test, ")")
    next
  }
  
  # DESeq2
  dds <- DESeqDataSetFromMatrix(
    countData = round(cnt_cl),
    colData   = meta_cl,
    design    = ~ treatment
  )
  
  # QC check: sample-level separation before DE
  rld <- rlog(dds, blind = TRUE)
  
  dds <- DESeq(dds)
  
  # QC check: dispersion trend fit
  plotDispEsts(dds, main = paste0("Dispersion Fit - Cluster ", cl))
  
  coef_name <- paste0("treatment_", cond_test, "_vs_", cond_ref)
  if (!coef_name %in% resultsNames(dds)) {
    stop("Coefficient not found: ", coef_name,
         " | available: ", paste(resultsNames(dds), collapse = ", "))
  }
  
  res_raw <- results(dds, name = coef_name)
  DESeq2::plotMA(res_raw, alpha = 0.05, main = paste0("MA (Raw LFC) - Cluster ", cl))
  res_shr <- lfcShrink(dds, coef = coef_name, type = "apeglm")
  DESeq2::plotMA(res_shr, alpha = 0.05, main = paste0("MA (Shrunken LFC) - Cluster ", cl))
  
  # keep p-values from Wald test, replace only LFC/SE with shrunken estimates
  res_cl <- res_raw
  res_cl$log2FoldChange <- res_shr$log2FoldChange
  res_cl$lfcSE <- res_shr$lfcSE
  
  
  res_cl <- as.data.frame(res_cl) %>%
    rownames_to_column("gene") %>%
    mutate(
      cluster  = cl,
      contrast = paste0(cond_test, "_vs_", cond_ref)
    )
  
  results_list[[cl]] <- res_cl
}

# Combine all clusters
all_de <- bind_rows(results_list)

all_de_sig <- all_de %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)


# Save results
write.csv(all_de,
  file = paste0("pseudobulk_shrinkage_0.3res_Salmonella_vs_Control.csv"),
  row.names = FALSE
)

## plot
all_de$cluster <- factor(all_de$cluster, levels = as.character(0:17))
ggplot(all_de, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1), alpha = 0.5) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey50"),
                     name = "Significant DE") +
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value",
       title = "Differentially Expressed Genes by Cluster (0.3res, DESeq2, filtered)",
       subtitle = "Red = p_adj < 0.05 & |log2FC| > 1") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 30), oob = squish) +
  scale_x_continuous(limits = c(-4, 4), oob = squish) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~ cluster, ncol = 6)
