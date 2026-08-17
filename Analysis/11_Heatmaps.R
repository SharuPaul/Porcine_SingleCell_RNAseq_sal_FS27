library(ComplexHeatmap)
library(circlize)
library(biomaRt)
library(readxl)
set.seed(123)

# shrunk LFC
all_de <- read.csv("pseudobulk_shrinkage_0.3res_Salmonella_vs_Control.csv")
wb <- read_excel("Results_DGE_Crl02vSal02.xlsx", sheet = 1,  range = NULL)
[, c("GeneName", "baseMean", "log2FoldChange", "padj")]

# gene list
genes <- c(
  "ACOD1", "ACSL1", "ADGRE1", "ADGRG3", "AIF1", "ALOX5AP", "ANKRD22",
  "AOAH", "B3GALNT1", "B3GNT5", "BATF2", "BCL6", "C15orf48", "C1QC",
  "C1R", "CALHM6", "CCR1", "CCR5", "CD180", "CD68", "CHIT1", "CMKLR1",
  "CSTA", "DOK3", "ENSSSCG00000005209", "ENSSSCG00000006919",
  "ENSSSCG00000014618", "ENSSSCG00000014997", "ENSSSCG00000028802",
  "ENSSSCG00000030882", "ENSSSCG00000031087", "ENSSSCG00000032135",
  "ENSSSCG00000039758", "FAM111A", "GBP1", "GBP2", "HAVCR2", "HK3", "HLX",
  "HPS5", "HSD3B7", "IDH2", "IFI35", "IL21R", "IRF1", "KMO", "KPNA2",
  "LAMTOR2", "LAP3", "LCN2", "LPAR6", "LPCAT2", "MARVELD2", "MCEMP1",
  "MICALL2", "MRPL24", "MS4A7", "MS4A8", "NSMF", "OSCAR", "P2RY13", "P2RY14",
  "PADI4", "PDCD1LG2", "PEAK3", "PSMB10", "PSMB9", "PSTPIP2", "RETN", "RIPK3",
  "S100A11", "S100A12", "S100A8", "S100A9", "SERPING1", "SLC2A5", "SLC39A8",
  "SOD2", "SOWAHD", "TCEA3", "TCN1", "TEC", "TESC", "TFEC", "TGM1", "TICAM2",
  "TMEM92", "TNF", "TSPO", "UPP1", "WARS", "ZCWPW1"
)

all_de$gene <- as.character(all_de$gene)
all_de$cluster <- as.character(all_de$cluster)
wb$GeneName <- as.character(wb$GeneName)

cluster_order <- c("0", "4", "5", "9", "12", "16", "3", "10", "17", "1",
                   "7", "8", "6", "14", "2", "13", "11", "18", "15")
cluster_labels <- c("Mono_1", "Mono_2", "Mono_3", "Mono_4",
                    "pDC", "cDC", "Bcell_1", "Bcell_2", "ASC",
                    "abT_CD4_1", "abT_CD4_2", "abT_NK_1", "abT_NK_2",
                    "abT_prolif", "gdT_CD2n", "gdT_CD2p",
                    "Mixed_Tcells_1", "Mixed_Tcells_2", "Mixed_T_Mono")
cluster_map <- setNames(cluster_labels, cluster_order)

cluster_levels <- unique(all_de$cluster)
mono_clusters <- c("0", "4", "5", "9")
mono <- all_de[all_de$cluster %in% mono_clusters, ]
genes_use <- unique(genes)

makeStars <- function(p) {
  s <- rep("", length(p))
  s[!is.na(p) & p < 0.05] <- "*"
  s[!is.na(p) & p < 0.01] <- "**"
  s[!is.na(p) & p < 0.001] <- "***"
  s[!is.na(p) & p < 1e-4] <- "****"
  s
}

# Plot settings
gene_font <- 8
title_font <- 10
star_font <- 7

width_wb <- unit(8, "mm")
width_mono <- unit(28, "mm")
width_all <- unit(max(28, 7 * length(cluster_levels)), "mm")
col_dend_height <- unit(4, "mm")

make_mono_heatmap <- function(mat, stars = NULL, col_fun, title,
                              cluster_rows = FALSE,
                              cluster_columns = TRUE,
                              width_mono = unit(28, "mm"),
                              col_dend_height = unit(4, "mm"),
                              gene_font = 8, title_font = 10, star_font = 7) {
  Heatmap(
    mat,
    name = "log2FC",
    col = col_fun,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    column_dend_height = col_dend_height,
    width = width_mono,
    row_names_gp = gpar(fontsize = gene_font),
    column_title = title,
    column_title_gp = gpar(fontsize = title_font),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = title_font),
      labels_gp = gpar(fontsize = gene_font)
    ),
    na_col = "grey90",
    show_row_names = TRUE,
    cell_fun = if (is.null(stars)) NULL else function(j, i, x, y, w, h, fill) {
      s <- stars[i, j]
      if (nzchar(s)) grid.text(s, x, y, gp = gpar(fontsize = star_font))
    }
  )
}

make_wb_heatmap <- function(mat, stars = NULL, col_fun, title = "WB",
                            cluster_rows = FALSE,
                            cluster_columns = FALSE,
                            width_wb = unit(8, "mm"),
                            gene_font = 8, title_font = 10, star_font = 7) {
  Heatmap(
    mat,
    name = "log2FC",
    col = col_fun,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    width = width_wb,
    show_heatmap_legend = FALSE,
    column_title = title,
    row_names_gp = gpar(fontsize = gene_font),
    column_title_gp = gpar(fontsize = title_font),
    na_col = "grey90",
    show_row_names = TRUE,
    cell_fun = if (is.null(stars)) NULL else function(j, i, x, y, w, h, fill) {
      s <- stars[i, j]
      if (nzchar(s)) grid.text(s, x, y, gp = gpar(fontsize = star_font))
    }
  )
}

make_cluster_mats <- function(de_df, clusters, genes_use) {
  de_ord <- de_df[order(abs(de_df$log2FoldChange), decreasing = TRUE), ]
  de_ord <- de_ord[!duplicated(paste(de_ord$gene, de_ord$cluster)), ]
  
  mat <- matrix(
    NA_real_,
    nrow = length(genes_use),
    ncol = length(clusters),
    dimnames = list(genes_use, clusters)
  )
  
  stars <- matrix(
    "",
    nrow = length(genes_use),
    ncol = length(clusters),
    dimnames = list(genes_use, clusters)
  )
  padj_mat <- matrix(
    NA_real_,
    nrow = length(genes_use),
    ncol = length(clusters),
    dimnames = list(genes_use, clusters)
  )
  for (i in seq_len(nrow(de_ord))) {
    g <- de_ord$gene[i]
    cl <- de_ord$cluster[i]
    
    if (!(g %in% genes_use) || !(cl %in% clusters)) next
    
    mat[g, cl] <- de_ord$log2FoldChange[i]
    stars[g, cl] <- makeStars(de_ord$padj[i])
    padj_mat[g, cl] <- de_ord$padj[i]
  }
  
  list(mat = mat, stars = stars, padj = padj_mat)
}

make_wb_mats <- function(wb, genes_use) {
  wb_ord <- wb[order(abs(wb$log2FoldChange), decreasing = TRUE), ]
  wb_ord <- wb_ord[!duplicated(wb_ord$GeneName), ]
  
  mat <- matrix(
    NA_real_,
    nrow = length(genes_use),
    ncol = 1,
    dimnames = list(genes_use, "WB")
  )
  
  stars <- matrix(
    "",
    nrow = length(genes_use),
    ncol = 1,
    dimnames = list(genes_use, "WB")
  )
  padj_mat <- matrix(
    NA_real_,
    nrow = length(genes_use),
    ncol = 1,
    dimnames = list(genes_use, "WB")
  )
  
  for (i in seq_len(nrow(wb_ord))) {
    g <- wb_ord$GeneName[i]
    
    if (!(g %in% genes_use)) next
    
    mat[g, "WB"] <- wb_ord$log2FoldChange[i]
    stars[g, "WB"] <- makeStars(wb_ord$padj[i])
    padj_mat[g, "WB"] <- wb_ord$padj[i]
  }
  
  list(mat = mat, stars = stars, padj = padj_mat)
}

rename_heatmap_rows <- function(mat, name_map) {
  mat2 <- mat
  rn <- rownames(mat2)
  
  rownames(mat2) <- ifelse(
    rn %in% names(name_map) & !is.na(name_map[rn]) & name_map[rn] != "",
    name_map[rn],
    rn
  )
  
  mat2
}

draw_wb_plus_clusters <- function(wb_mat, cluster_mat,
                                  wb_stars = NULL, cluster_stars = NULL,
                                  cluster_title = "Monocytes",
                                  cluster_width = width_mono,
                                  col_fun = col_fun) {
  clust_mat <- cluster_mat
  clust_mat[is.na(clust_mat)] <- 0
  row_dend <- as.dendrogram(hclust(dist(clust_mat)))
  
  draw(
    make_wb_heatmap(
      wb_mat,
      stars = wb_stars,
      col_fun = col_fun,
      title = "WB",
      cluster_rows = row_dend, ##  row_dend, or FALSE
      cluster_columns = FALSE,
      width_wb = width_wb,
      gene_font = gene_font,
      title_font = title_font,
      star_font = star_font
    ) +
      make_mono_heatmap(
        cluster_mat,
        stars = cluster_stars,
        col_fun = col_fun,
        title = cluster_title,
        cluster_rows = row_dend,
        cluster_columns = TRUE,
        width_mono = cluster_width,
        col_dend_height = col_dend_height,
        gene_font = gene_font,
        title_font = title_font,
        star_font = star_font
      ),
    heatmap_legend_side = "right"
  )
}

##########################################
##########################################

# Build matrices for ALL-CLUSTER heatmap
all_mats <- make_cluster_mats(all_de, cluster_levels, genes_use_all)
wb_mats_all <- make_wb_mats(wb, genes_use_all)

# Build matrices for MONOCYTE heatmap
mono_mats <- make_cluster_mats(mono, mono_clusters, genes_use_mono)
wb_mats_mono <- make_wb_mats(wb, genes_use_mono)

all_mat <- all_mats$mat

mono_mat <- mono_mats$mat
mono_stars <- mono_mats$stars

wb_mat_all <- wb_mats_all$mat
# wb_stars_all <- wb_mats_all$stars

wb_mat_mono <- wb_mats_mono$mat
wb_stars_mono <- wb_mats_mono$stars

all_padj <- all_mats$padj
mono_padj <- mono_mats$padj

wb_padj_all <- wb_mats_all$padj
wb_padj_mono <- wb_mats_mono$padj

# Rename cluster columns for heatmap display only
colnames(all_mat) <- unname(cluster_map[colnames(all_mat)])
colnames(mono_mat) <- unname(cluster_map[colnames(mono_mat)])
colnames(mono_stars) <- colnames(mono_mat)

# Common color scale across all plots
all_vals <- c(as.vector(wb_mat_all), as.vector(wb_mat_mono),
  as.vector(all_mat),  as.vector(mono_mat))
# all_vals <- c(as.vector(wb_mat), as.vector(all_mat), as.vector(mono_mat))
lim <- quantile(abs(all_vals), 0.95, na.rm = TRUE)
col_fun <- colorRamp2(c(-lim, 0, lim), c("blue", "white", "red"))

# # Positive-only log2FC color cap for UP-gene heatmap
col_fun <- colorRamp2(c(0, 2, 4), c("white", "pink", "red"))


## Plot Heatmaps

# Heatmap 1 (WB + All Clusters)
wb_mat = wb_mat_all
# wb_stars = wb_stars_all
cluster_mat = all_mat

draw_wb_plus_clusters(
  wb_mat = wb_mat,
  cluster_mat = all_mat,
  wb_stars = NULL,
  cluster_stars = NULL,
  cluster_title = "All clusters",
  cluster_width = width_all,
  col_fun = col_fun
)

# Heatmap 2/3 (WB + Monocytes)
wb_mat = wb_mat_mono
wb_stars = wb_stars_mono
cluster_mat = mono_mat

draw_wb_plus_clusters(
  wb_mat = wb_mat,
  cluster_mat = mono_mat,
  wb_stars = NULL,
  cluster_stars = NULL,
  cluster_title = "Monocytes",
  cluster_width = width_mono,
  col_fun = col_fun
)

# with asterisks
draw_wb_plus_clusters(
  wb_mat = wb_mat,
  cluster_mat = mono_mat,
  wb_stars = wb_stars,
  cluster_stars = mono_stars,
  cluster_title = "Monocytes",
  cluster_width = width_mono,
  col_fun = col_fun
)
