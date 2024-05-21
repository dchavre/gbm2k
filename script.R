# Dataset: https://www.10xgenomics.com/datasets/2-k-sorted-cells-from-human-glioblastoma-multiforme-3-v-3-1-3-1-standard-6-0-0

# Import Statements
library(dplyr)
library(Seurat)
library(patchwork)

# Setting up the Seurat Object

# Loading Dataset
gbm.data <- Read10X(data.dir = 'C:/Users/dchav/Documents/10X_Genomics/Brain_Tumor_3p/filtered_feature_bc_matrix/') # Read10X() is able to read the output of the cellranger pipeline

# Initialize Seurat Object with Raw (Non-Normalized) Data
gbm <- CreateSeuratObject(counts=gbm.data, project = 'gbm2k', min.cells = 3, min.features = 200)

gbm

# Standard Pre-Processing Workflow

# Quality Control & Finding Cells for Further Analysis

# [[ can add columns to object metadata (to store QC stats)
# PercentageFeatureSet() calculates percentage
# of counts originating from a set of features
gbm[['percent.mt']] <- PercentageFeatureSet(gbm, pattern = '^MT-') # uses genes that start with MT- (mitochondrial genes)

# Visualizing Quality Control via Violin Plot
VlnPlot(gbm, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)

# FeatureScatter() allows you to view feature-feature relationships,
# but can be used for anything calculated by the object,
# like columns, in object metadata, PC scores, etc.
plot1 <- FeatureScatter(gbm, feature1 = 'nCount_RNA', feature2 = 'percent.mt')
plot2 <- FeatureScatter(gbm, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')

plot1 + plot2

gbm <- subset(gbm, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25)

# Data Normalization

# Using LogNormalize, which normalizes feature expression for each cell
# via the total expression, multiples by a scale factor (10k default),
# and then log-transforms the result

gbm <- NormalizeData(gbm) # LogNormalize via scale factor of 10k by default

# Identifying highly variable features (Feature Selection)

# Finding a subset of features that show
# high cell-to-cell variation in the dataset
# (highly expressed in some cells, lowly expressed in others)
# this helps to find biological signals in scRNA-seq datasets

# This is done via modeling mean-variance relationship in the data
# with FindVariableFeatures(), returning 2k feaures by default 
# (per dataset), and these can be used in later analyses, like PCA

gbm <- FindVariableFeatures(gbm, selection.method = 'vst', nfeatures = 2000)

# Identifying top 10 most variable genes
top_10 <- head(VariableFeatures(gbm), 10)

top_10

# Plotting variable feature with and without labels
plot1 <- VariableFeaturePlot(gbm)
plot2 <- LabelPoints(plot = plot1, points = top_10, repel = TRUE)
plot1 + plot2

# Scaling the Data

# Applying a linear transformation, which is standard pre-processing prior to
# PCA or other dimensionality reduction techniques.
# ScaleData() shifts the expression of each gene, so that mean expression is 0
# Scales expression of each gene so that variance is 1, for equal weighting in
# downstream analysis, only variable features are scaled by default, or specify.

all.genes <- rownames(gbm)
gbm <- ScaleData(gbm, features = all.genes)
 
# Performing Linear Dimensional Reduction

# Performing PCA on scaled data; variable features from before are the input,
# but it is possible to define specific features as well.

# For the first principal components, Seurat outputs a list of genes with the
# most positive and negative loadings, which represent modules of genes that 
# exhibit correlation/anti-correlation across single-cells in the dataset.

gbm <- RunPCA(gbm, features = VariableFeatures(object = gbm))

# Can examine cells and features that define PCA via VizDimReduction(),
# DimPlot(), and DimHeatMap()

print(gbm[['pca']], dim = 1:5, nfeatures = 5)

VizDimLoadings(gbm, dims = 1:2, reduction = 'pca')

DimPlot(gbm, reduction = 'pca') + NoLegend()

# DimHeatmap() allows for examination of primary sources of heterogeneity,
# and is helpful when trying to determine what PCs to analyze further.
# Both cells and features are ordered according to their PCA scores.
# Setting cells to a number plots the extreme cells on both ends of the spectrum,
# while makes plotting faster for larger datasets. 

DimHeatmap(gbm)

# Determining the Dimensionality of the Dataset

# To overcome extensive tehcnical noise in any feature for scRNA-seq data,
# Seurat clusters cells based on their PCA scores, and each PC represents a 
# metafeature that combines information across a correlated feature set

# An elbow plot ranks PCs based on the percentage of variance explained by each.
# An 'elbow' can be seen at around PC7, meaning that the majority of true signal
# lies withn the first 7 PCs

ElbowPlot(gbm)

# The true dimensionality of a dataset can also be found via exploring PCs to 
# find sources of heterogeneity, with the use of GSEA, and the use of a
# heuristic, which is commonly used, and calculated instantly (also elbow plot).

# Cell Clustering

# Seurat uses a graph-based clustering approach, with KNN, and other techniques.

# FindNeighbors() constructs a KNN graph based on euclidian distance in 
# PCA space, and then refines edge weights between any two cells based on
# the shared overlap in their local neighborhoods (takes in first 10 PCs).

gbm <- FindNeighbors(gbm, dims = 1:10)

# FindClusters() uses modularity optimization technqiues, such as
# the Louvain algorithm (default), or SLM, which iteratively groups cells together,
# with the goal of optimizing the standard modularity function.

# The resolution parameter sets the granularity of the clustering, where with
# a greater resolution, more clusters are found. Best to keep parameter 0.4-1.2,
# with datasets that are 3K in size (higher resolution = better for large data)

gbm <- FindClusters(gbm, resolution = 0.5)

# Clusters can then be found via the Idents() function
head(Idents(gbm), 5)

# Running Non-Linear Dimensional Reduction (UMAP/tSNE)

# tSNE and UMAP can be used to visualize and explore datasets. These algorithms
# assist with trying to find an underlying structure within the dataset, 
# in order to place similar cells together, in a lo-dimensional space.
# Cells that are grouped together within graph-based clusters from previous
# analyses should co-localize on these dimension reduction plots.

gbm <- RunUMAP(gbm, dims = 1:10)

DimPlot(gbm, reduction = 'umap', label = TRUE) # you can set label = TRUE or LabelClusters() to label individual clusters

# It is recommended to save the Seurat object to avoid redundancy
saveRDS(gbm, file = 'C:/Users/dchav/Documents/10X_Genomics/Brain_Tumor_3p/output/gbm_Seurat_object.rds')

# Finding Differentially Expressed Features (Cluster Biomarkers)

# Seurat can find markers that define clusters via differential expression.
# By default, it identifies positive and negative markers of a single cluster
# (specified in ident.1), compared to all other cells.
# FindAllMarkers() automates tis process for all clusters.


# Finding all markers of cluster 0

cluster0.markers <- FindMarkers(gbm, ident.1 = 0)
head(cluster0.markers, n = 5)

# Find all markers distinguishing cluster 10 from cluster 6 and 7

cluster10.markers <- FindMarkers(gbm, ident.1 = 10, ident.2 = c(6, 7))
head(cluster10.markers, n = 5)

# Find markers for every cluster compared to all remaining cells

gbm.markers <- FindAllMarkers(gbm, only.pos = TRUE) # only positive markers
gbm.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# Seurat has multiple tests for DE, which can be set via the test.use parameter.
# Ex. the ROC test returns the classification power for any marker, that ranges
# from 0 (random) to 1 (perfect)

cluster1.markers <- FindMarkers(gbm, ident.1 = 1, logfc.threshold = 0.25, test.use = 'roc', only.pos = TRUE)

# There are numerous tools fr visualizing marker expression. VlnPlot() shows
# expression probability distributions across clusters. FeaturePlot() visualizes 
# feature expression on a tSNE or PCA plot. RidgePlot(), CellScatter(),
# and DotPlot() are other ways to view a dataset.

VlnPlot(gbm, features = c('FCER1G'))
VlnPlot(gbm, features = c('FCER1G'), slot = 'counts', log = TRUE) # to plot raw counts

FeaturePlot(gbm, features = c('CHL1', 'SLC4A4', 'PTPRZ1', 'F3', 'FABP7', 'LUZP2', 'FAM181B', 'POU3F2', 'SOX9', 'NPAS3'))

# DoHeatmap() creates an expression heatmap for given cells and features.

gbm.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top_10
DoHeatmap(gbm, features = top_10$gene) + NoLegend()

# Assigning Cell Type Identity to Clusters (needs canonical marker genes)
# to match the unbiased clustering to known cell types.


