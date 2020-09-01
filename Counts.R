install.packages("ggfortify")
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(ggfortify)
library(patchwork)
set.seed(1234)


counts <- Read10X_h5(filename = "atac_v1_pbmc_5k_filtered_peak_bc_matrix.h5")

metadata <- read.csv(
  file = "atac_v1_pbmc_5k_singlecell.csv",
  header = TRUE,
  row.names = 1
)
chromatinassay <- CreateChromatinAssay(counts = counts, genome = "hg19")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = 'atac_v1_pbmc_5k_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

#pbmc
#pbmc[['peaks']]
#granges(pbmc)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg19"

# add the gene information to the object
Annotation(pbmc) <- annotations


# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()

pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')

VlnPlot(
  object = pbmc,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

pbmc <- subset(
  x = pbmc,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

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


saveRDS(pbmcCounts, "pbmcCounts.rds")
write.csv(pbmcCounts@active.ident, "pbmcCounts_ident.csv")











