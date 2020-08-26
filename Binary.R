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

# Convert to binary data
counts[counts > 0] = 1

pbmc <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'pbmc5k',
  min.cells = 1,
)

# counts
# pbmc


# Normalization: Signac performs term frequency-inverse document frequency (TF-IDF) normalization.
# This is a two-step normalization procedure, that both normalizes across cells to correct for differences in cellular sequencing depth, and across peaks to give higher values to more rare peaks.
pbmcBinary <- RunTFIDF(pbmc)

# Feature selection: use only the top n% of features (peaks) for dimensional reduction, or remove features present in less that n cells with the FindTopFeatures function. 
# min.cutoff = 'q0' -> all features    min.cutoff to ‘q75’ -> the top 25% all peaks
pbmcBinary <- FindTopFeatures(pbmcBinary, min.cutoff = 'q0')

# Dimensional reduction: We next run a singular value decomposition (SVD) on the TD-IDF normalized matrix, using the features (peaks) selected above. This returns a low-dimensional representation of the object
pbmcBinary <- RunSVD(
  object = pbmcBinary,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)

# The first LSI component often captures sequencing depth (technical variation) rather than biological variation. If this is the case, the component should be removed from downstream analysis. 
# We can assess the correlation between each LSI component and sequencing depth using the DepthCor function:
# DepthCor(pbmcBinary)

# Here we see there is a very strong correlation between the first LSI component and the total number of counts for the cell, so we will perform downstream steps without this component.

pbmcBinary <- RunUMAP(object = pbmcBinary, reduction = 'lsi', dims = 2:30)
pbmcBinary <- FindNeighbors(object = pbmcBinary, reduction = 'lsi', dims = 2:30)
pbmcBinary <- FindClusters(object = pbmcBinary, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmcBinary, label = TRUE) 


# saveRDS(pbmcBinary, "pbmcBinary.rds")
# write.csv(pbmcBinary$seurat_clusters, "pbmcBinary_ident.csv")



# pbmcBinary@active.ident <- pbmcCounts@active.ident
# DimPlot(pbmcBinary, label = TRUE)





