library(Seurat)
library(Signac)
library(ggplot2)
library(patchwork)

#pbmcBinary <- readRDS("pbmcBinary.rds")

# Gene activity quantification
peaks <- Read10X_h5("atac_v1_pbmc_5k_filtered_peak_bc_matrix.h5")

# convert to binary
# peaks[peaks > 0] = 1


# Object setup
pbmc.atac <- CreateSeuratObject(counts = peaks, assay = "ATAC", project = "10x_ATAC", min.cells = 1)

pbmc.atac <- RunTFIDF(pbmc.atac)
pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = 'q0')
pbmc.atac <- RunSVD(
  object = pbmc.atac,
  assay = 'ATAC',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)
pbmc.atac <- RunUMAP(pbmc.atac, reduction = "lsi", dims = 2:30)
p1 <- DimPlot(pbmc.atac, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq")
p1
# create a gene activity matrix from the peak matrix and GTF, using chromosomes 1:22, X, and Y.
# Peaks that fall within gene bodies, or 2kb upstream of a gene, are considered
activity.matrix <- CreateGeneActivityMatrix(peak.matrix = peaks, annotation.file = "Homo_sapiens.GRCh38.100.gtf", 
                                            seq.levels = c(1:22, "X", "Y"), upstream = 2000, verbose = TRUE)
pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
meta <- read.table("atac_v1_pbmc_5k_singlecell.csv", sep = ",", header = TRUE, row.names = 1, 
                   stringsAsFactors = FALSE)
meta <- meta[colnames(pbmc.atac), ]
pbmc.atac <- AddMetaData(pbmc.atac, metadata = meta)
# pbmc.atac <- subset(pbmc.atac, subset = nCount_ATAC > 5000)
pbmc.atac$tech <- "atac"

# Data processing
DefaultAssay(pbmc.atac) <- "ACTIVITY"
pbmc.atac <- FindVariableFeatures(pbmc.atac)
pbmc.atac <- NormalizeData(pbmc.atac)
pbmc.atac <- ScaleData(pbmc.atac)

DefaultAssay(pbmc.atac) <- "ATAC"
VariableFeatures(pbmc.atac) <- names(which(Matrix::rowSums(pbmc.atac) > 100))

pbmc.rna <- readRDS("pbmc_10k_v3.rds")
pbmc.rna$tech <- "rna"


p2 <- DimPlot(pbmc.rna, group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scRNA-seq")
p1 + p2


transfer.anchors <- FindTransferAnchors(reference = pbmc.rna, query = pbmc.atac, features = VariableFeatures(object = pbmc.rna), 
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = pbmc.rna$celltype, 
                                     weight.reduction = pbmc.atac[["lsi"]])
pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions)


hist(pbmc.atac$prediction.score.max)
abline(v = 0.3, col = "red")

table(pbmc.atac$prediction.score.max > 0.3)


pbmc.atac.filtered <- subset(pbmc.atac, subset = prediction.score.max > 0.3)
pbmc.atac.filtered$predicted.id <- factor(pbmc.atac.filtered$predicted.id, levels = levels(pbmc.rna))  # to make the colors match
p1 <- DimPlot(pbmc.atac.filtered, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + 
  NoLegend() + scale_colour_hue(drop = FALSE)
p2 <- DimPlot(pbmc.rna, group.by = "celltype", label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") + 
  NoLegend()
p1 + p2



# Co-embedding
# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(pbmc.rna)
refdata <- GetAssayData(pbmc.rna, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = pbmc.atac[["lsi"]])

# this line adds the imputed data matrix to the pbmc.atac object
pbmc.atac[["RNA"]] <- imputation
coembed <- merge(x = pbmc.rna, y = pbmc.atac)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)
coembed$celltype <- ifelse(!is.na(coembed$celltype), coembed$celltype, coembed$predicted.id)

p1 <- DimPlot(coembed, group.by = "tech")
p2 <- DimPlot(coembed, group.by = "celltype", label = TRUE, repel = TRUE)
p1 + p2

DimPlot(coembed, split.by = "tech", group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend()

write.csv(coembed$celltype, "coembed_ct.csv")
write.csv(pbmcBinary$seurat_clusters, "pbmcBinary_clu.csv")

table(coembed$celltype)
table(coembed$tech)
table(coembed$celltype, coembed$tech)

#DimPlot(pbmcBinary, label = TRUE)


