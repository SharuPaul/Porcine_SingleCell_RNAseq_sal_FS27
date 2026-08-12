# Methods

## Data Analysis

10x Genomics Cell Ranger (v. 7.2.0) was used for aligning the reads to a reference generating feature matrix and bam files for each sample. Feature matrix was used for further analysis.

R libraries from the tidyverse (Wickham et al., 2023) and related packages such as dplyr (version 1.1.4), stringr (version 1.5.1), readxl (version 1.4.3), ggplot2 (version 3.5.1), patchwork (version 1.2.0), and scales (version 1.3.0) were used for data wrangling and plotting.

The data analysis pipeline was run on R (version 4.3.2; R Core Team, 2021) and Seurat (version 4.3.0.1; Hao et al., 2021). SeuratObject (version 5.0.1) and SeuratDisk (version 0.0.0.9021) were used for loading and saving Seurat data.

Ambient RNA estimation and removal was performed using SoupX (v. 1.6.2; Young and Behjati, 2020). SingleCellExperiment (version 1.24.0; Amezquita et al., 2020) was used to convert Seurat object into a single cell experiment object for finding doublets using scDblFinder (version 1.16.0; Germain et al., 2022). QC plots were made using ggplot2 and cells with <25% mitochondrial reads, >400 genes detected, >500 total reads, not called as doublets were kept.

Data were normalized using SCTransform and integrated for visualization on a UMAP. PCA was done using 1:16 pcs. Clustree was made using clustree (version 0.5.1; Zappia et al., 2018) and ggraph (version 2.1.0; Pedersen T, 2022) to choose clustering resolution. After exploration, a final clustering resolution of 0.3 was selected.  

Marker genes were found using FindAllMarkers() function. DotPlot() was used to plot the marker gene expression per cluster. Clusters were annotated after manual curation of marker gene lists.

GO enrichment analysis was performed using clusterProfiler (version 4.10.1; Yu et al., 2012) and biomaRt (version 2.58.2; Durinck et al., 2009) to query the *S. scrofa* ensembl dataset. GO.db (version 3.18.0) was used to map GO IDs to the human readable descriptions.

Gene-set activity scores were computed per cell using AUCell (v. 1.24.0; Aibar et al., 2017). Differential abundance analysis was done using miloR (v. 1.10.0; Dann et al., 2022). Differential Expression analysis was performed using DESeq2 (v. 1.42.1; Love et al., 2014). DE data were visualized as volcano plots made using ggplot and as heatmaps made using ComplexHeatmap (v. 2.18.0; Gu et al., 2016) packages. Upsetplots were made using ComplexUpset (v. 1.3.3; Krassowski, 2020).

### References
* Aibar S, Bravo Gonzalez-Blas C, Moerman T, Huynh-Thu V, Imrichova H, Hulselmans G, Rambow F, Marine J, Geurts P, Aerts J, van den Oord J, Kalender Atak Z, Wouters J, Aerts S (2017). “SCENIC: Single-Cell Regulatory Network Inference And Clustering.” Nature Methods, 14, 1083-1086. doi:10.1038/nmeth.4463.
* Amezquita R, Lun A, Becht E, Carey V, Carpp L, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pages H, Smith M, Huber W, Morgan M, Gottardo R, Hicks S (2020). “Orchestrating single-cell analysis with Bioconductor.” Nature Methods, 17, 137–145. doi: 10.1038/s41592-019-0654-x.
* Dann, E., Henderson, N. C., Teichmann, S. A., Morgan, M. D., & Marioni, J. C. (2022). Differential abundance testing on single-cell data using k-nearest neighbor graphs. Nature biotechnology, 40(2), 245-253.
* Durinck S, Spellman PT, Birney E, Huber W. Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt. Nat Protoc. 2009;4(8):1184-91. doi: 10.1038/nprot.2009.97.
* Germain P, Lun A, Garcia Meixide C, Macnair W, Robinson M (2022). “Doublet identification in single-cell sequencing data using scDblFinder.” f1000research. doi: 10.12688/f1000research.73600.2.
* Gu Z, Eils R, Schlesner M (2016). “Complex heatmaps reveal patterns and correlations in multidimensional genomic data.” Bioinformatics. doi:10.1093/bioinformatics/btw313.
* Hao, Y., Hao, S., Andersen-Nissen, E., Mauck, W.M., Zheng, S., Butler, A., et al. (2021). Integrated analysis of multimodal single-cell data. Cell 184(13), 3573-3587.e3529. doi: 10.1016/j.cell.2021.04.048.
* Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi:10.1186/s13059-014-0550-8.
* Michał Krassowski. (2020). krassowski/complex-upset. Zenodo. http://doi.org/10.5281/zenodo.3700590
* Pedersen T (2022). ggraph: An Implementation of Grammar of Graphics for Graphs and Networks. R package version 2.1.0, https://github.com/thomasp85/ggraph.
* R Core Team (2021) R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna. https://www.R-project.org.
* Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to the tidyverse.” Journal of Open Source Software, 4(43), 1686. doi: 10.21105/joss.01686.
* Matthew D Young, Sam Behjati, SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data, GigaScience, Volume 9, Issue 12, December 2020, giaa151, https://doi.org/10.1093/gigascience/giaa151
* Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS. 2012 May;16(5):284-7. doi: 10.1089/omi.2011.0118
* Zappia, L., and Oshlack, A. (2018). Clustering trees: a visualization for evaluating clusterings at multiple resolutions. GigaScience 7(7). doi: 10.1093/gigascience/giy083.
