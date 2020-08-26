install.packages("ggfortify")
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
set.seed(1234)
library(ggfortify)

counts <- Read10X_h5("atac_v1_pbmc_5k_filtered_peak_bc_matrix.h5")
pbmc <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'pbmc5k',
  min.cells = 1,
)

counts
pbmc

# Log-normalization
# By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.
pbmcCounts <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
# Calculating variances
# vst: First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess). Then standardizes the feature values using the observed mean and expected variance (given by the fitted line). Feature variance is then calculated on the standardized values after clipping to a maximum (see clip.max parameter).
pbmcCounts <- FindVariableFeatures(pbmcCounts, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmcCounts), 10)

# plot variable features with and without labels
#plot1 <- VariableFeaturePlot(pbmcCounts)
#plot2 <- LabelPoints(plot = plot1, points = top10)
#plot1 + plot2
all.genes <- rownames(pbmcCounts)
pbmcCounts <- ScaleData(pbmcCounts, features = all.genes)
pbmcCounts <- RunPCA(pbmcCounts, features = VariableFeatures(object = pbmcCounts))

# Examine and visualize PCA results a few different ways
print(pbmcCounts[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmcCounts, dims = 1:2, reduction = "pca")
#DimPlot(pbmcCounts,  reduction = "pca")

# Exploration of the primary sources of heterogeneity in a dataset, and can be usefu, when trying to decide which PCs to include for further downstream analyses
#DimHeatmap(pbmcCounts, dims = 1, cells = 500, balanced = TRUE)
# PC1-15
#DimHeatmap(pbmcCounts, dims = 1:15, cells = 500, balanced = TRUE)


# Determine the 'dimensionality' of the dataset
# Seurat clusters cells based on their PCA scores, with each PC essentially representing a ‘metafeature’ that combines information across a correlated feature set. The top principal components therefore represent a robust compression of the dataset.
# We randomly permute a subset of the data (1% by default) and rerun PCA, constructing a ‘null distribution’ of feature scores, and repeat this procedure. We identify ‘significant’ PCs as those who have a strong enrichment of low p-value features.
pbmcCounts <- JackStraw(pbmcCounts, num.replicate = 100)
pbmcCounts <- ScoreJackStraw(pbmcCounts, dims = 1:20)
# Visualisation
# ‘Significant’ PCs will show a strong enrichment of features with low p-values (solid curve above the dashed line).
JackStrawPlot(pbmcCounts, dims = 1:15)
ElbowPlot(pbmcCounts)

# Cluster the cells
# We chose 15 here
pbmcCounts <- FindNeighbors(pbmcCounts, dims = 1:15)
pbmcCounts <- FindClusters(pbmcCounts, resolution = 0.5)
#head(Idents(pbmcCounts), 5)


# UMAP
# to learn the underlying manifold of the data in order to place similar cells together in low-dimensional space
pbmcCounts <- RunUMAP(pbmcCounts, dims = 1:15)
DimPlot(pbmcCounts, label = TRUE, reduction = "umap")


#saveRDS(pbmcCounts, "pbmcCounts.rds")
#write.csv(pbmcCounts@active.ident, "pbmcCounts_ident.csv")











