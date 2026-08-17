library(ggplot2)
library(clusterProfiler)
library(biomaRt)
library(readr)
library(dplyr)
library(GO.db)
library(Matrix)
library(Seurat)
library(SeuratDisk)
library(readxl)
set.seed(123)


DefaultAssay(seu) <- "RNA"
bg <- rownames(seu)
rm(seu)
## load annotations
#annot <- read_excel()

# Convert
sym2ens <- setNames(annot$ENSID, annot$FinalList)

bg_ens <- sym2ens[bg]          # lookup
bg_ens <- bg_ens[!is.na(bg_ens)]  # drop NAs
bg_ens <- unique(bg_ens)   # remove duplicates

## mart
mart <- useEnsembl("ensembl", dataset = "sscrofa_gene_ensembl")

# Build BP map
go_map <- getBM(
  attributes = c("ensembl_gene_id", "go_id", "name_1006", "namespace_1003"),
  filters    = "ensembl_gene_id",
  values     = bg_ens,
  mart       = mart
)

go_map_bp <- go_map[go_map$namespace_1003 == "biological_process" &
                      !is.na(go_map$go_id) & nzchar(go_map$go_id), ]

# TERM2GENE 
TERM2GENE_BP <- unique(data.frame(
  GO      = go_map_bp$go_id,
  ENSEMBL = go_map_bp$ensembl_gene_id,
  stringsAsFactors = FALSE
))

# GO term dictionary 
GO_TERMS_BP <- unique(data.frame(
  ID          = go_map_bp$go_id,
  Description = go_map_bp$name_1006,
  stringsAsFactors = FALSE
))

bg_ens <- intersect(bg_ens, unique(TERM2GENE_BP$ENSEMBL))


# Per-cluster GO BP enrichment
go_results_by_cluster <- list()
clusters <- unique(deg_data$cluster)

############################################
####### do up and down reg degs separately
for (cl in clusters) {
  ## subset once per cluster
  de_cl <- deg_data[deg_data$cluster == cl, ]
  
  ## split into up- and down-regulated
  genes_cl_up   <- unique(de_cl$gene[de_cl$log2FoldChange > 0])
  genes_cl_down <- unique(de_cl$gene[de_cl$log2FoldChange < 0])
  
  ## --- UP ---
  genes_cl_up_ens <- sym2ens[genes_cl_up]
  genes_cl_up_ens <- unique(genes_cl_up_ens[!is.na(genes_cl_up_ens)])
  genes_cl_up_ens <- intersect(genes_cl_up_ens, bg_ens)
  
  res_up <- NULL
  if (length(genes_cl_up_ens) > 0) {
    res_up <- enricher(
      gene          = genes_cl_up_ens,
      universe      = bg_ens,
      TERM2GENE     = TERM2GENE_BP,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05
    )
  }
  
  ## --- DOWN ---
  genes_cl_down_ens <- sym2ens[genes_cl_down]
  genes_cl_down_ens <- unique(genes_cl_down_ens[!is.na(genes_cl_down_ens)])
  genes_cl_down_ens <- intersect(genes_cl_down_ens, bg_ens)
  
  res_down <- NULL
  if (length(genes_cl_down_ens) > 0) {
    res_down <- enricher(
      gene          = genes_cl_down_ens,
      universe      = bg_ens,
      TERM2GENE     = TERM2GENE_BP,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05
    )
  }
  
  ## store both for this cluster
  go_results_by_cluster[[as.character(cl)]] <- list(
    up   = res_up,
    down = res_down
  )
}


# Combine all clusters (up + down) into one data.frame and add BP term names
all_results <- do.call(rbind, lapply(as.character(clusters), function(cl) {
  res_cl <- go_results_by_cluster[[cl]]
  if (is.null(res_cl)) return(NULL)
  
  out_list <- list()
  
  ## ---- UP ----
  if (!is.null(res_cl$up)) {
    df_up <- as.data.frame(res_cl$up)
    if (nrow(df_up)) {
      df_up <- merge(df_up, GO_TERMS_BP, by = "ID", all.x = TRUE, sort = FALSE)
      df_up$cluster   <- cl
      df_up$direction <- "up"
      out_list[["up"]] <- df_up
    }
  }
  
  ## ---- DOWN ----
  if (!is.null(res_cl$down)) {
    df_down <- as.data.frame(res_cl$down)
    if (nrow(df_down)) {
      df_down <- merge(df_down, GO_TERMS_BP, by = "ID", all.x = TRUE, sort = FALSE)
      df_down$cluster   <- cl
      df_down$direction <- "down"
      out_list[["down"]] <- df_down
    }
  }
  
  if (!length(out_list)) return(NULL)
  do.call(rbind, out_list)
}))

## Calculate Fold Enrichment from GeneRatio and BgRatio
if (!is.null(all_results) && nrow(all_results)) {
  # GeneRatio: "a/b" -> numeric a/b
  gr_parts <- strsplit(all_results$GeneRatio, "/")
  br_parts <- strsplit(all_results$BgRatio, "/")
  
  all_results$GeneRatio_num <- sapply(gr_parts, function(x) as.numeric(x[1]) / as.numeric(x[2]))
  all_results$BgRatio_num   <- sapply(br_parts, function(x) as.numeric(x[1]) / as.numeric(x[2]))
  
  all_results$FoldEnrichment <- all_results$GeneRatio_num / all_results$BgRatio_num
}

## Reverse lookup: Ensembl to Symbol 
ens2sym <- setNames(annot$FinalList, annot$ENSID)

all_results$geneID_symbol <- vapply(strsplit(all_results$geneID, "/"), function(ids) {
  syms <- ens2sym[ids]
  syms[is.na(syms)] <- ids[is.na(syms)]
  paste(unique(syms), collapse = "/")
}, FUN.VALUE = character(1), USE.NAMES = FALSE)


####################################################################
## Plots for all clusters, up and down separately
for (cl in clusters) {
  cl_chr <- as.character(cl)
  res_cl <- go_results_by_cluster[[cl_chr]]
  
  if (is.null(res_cl)) next
  
  term_map <- setNames(GO_TERMS_BP$Description, GO_TERMS_BP$ID)
  
  for (direction in c("up", "down")) {
    res <- res_cl[[direction]]
    if (is.null(res)) next
    
    df_res <- as.data.frame(res)
    if (!nrow(df_res)) next
    
    # Attach BP term names
    res@result$Description <- term_map[res@result$ID]
    
    # Plot
    p1 <- dotplot(res, showCategory = 20) +
      ggtitle(paste("GO BP 0.25 cpm, 0.3 res cluster", cl_chr, "-", direction))
    print(p1)
  }
}

