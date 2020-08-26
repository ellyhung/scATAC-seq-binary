library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
set.seed(1234)


counts <- Read10X_h5(filename = "atac_v1_pbmc_5k_filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "atac_v1_pbmc_5k_singlecell.csv",
  header = TRUE,
  row.names = 1
)

# Convert to binary data
##counts[counts > 0] = 1



pbmc <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'ATAC',
  min.cells = 1,
  meta.data = metadata
)

fragment.path <- 'atac_v1_pbmc_5k_fragments.tsv.gz'

pbmcDAR <- SetFragments(
  object = pbmc,
  file = fragment.path
)


## Computing QC Metrics
pbmcDAR <- NucleosomeSignal(object = pbmcDAR)

pbmcDAR$pct_reads_in_peaks <- pbmcDAR$peak_region_fragments / pbmcDAR$passed_filters * 100
pbmcDAR$blacklist_ratio <- pbmcDAR$blacklist_region_fragments / pbmcDAR$peak_region_fragments

p1 <- VlnPlot(pbmcDAR, c('pct_reads_in_peaks', 'peak_region_fragments'), pt.size = 0.1)
p2 <- VlnPlot(pbmcDAR, c('blacklist_ratio', 'nucleosome_signal'), pt.size = 0.1) & scale_y_log10()
p1 | p2

# group by cells with high or low nucleosomal signal strength
# You can see that cells which are outliers for the mononucleosomal/ nucleosome-free ratio (based off the plots above) have different banding patterns. 
# The remaining cells exhibit a pattern that is typical for a successful ATAC-seq experiment.
pbmcDAR$nucleosome_group <- ifelse(pbmcDAR$nucleosome_signal > 10, 'NS > 10', 'NS < 10')
FragmentHistogram(object = pbmcDAR, group.by = 'nucleosome_group')

# calculate the TSS enrichment score for each cell
# The enrichment of Tn5 integration events at transcriptional start sites (TSSs) can also be an important quality control metric to assess the targeting of Tn5 in ATAC-seq experiments.
# create granges object with TSS positions
gene.ranges <- genes(EnsDb.Hsapiens.v75)
seqlevelsStyle(gene.ranges) <- 'UCSC'
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')

tss.ranges <- resize(gene.ranges, width = 1, fix = "start")
seqlevelsStyle(tss.ranges) <- 'UCSC'
tss.ranges <- keepStandardChromosomes(tss.ranges, pruning.mode = 'coarse')

# to save time use the first 2000 TSSs
pbmcDAR <- TSSEnrichment(object = pbmcDAR, tss.positions = tss.ranges[1:2000])

pbmcDAR$high.tss <- ifelse(pbmcDAR$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(pbmcDAR, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()

# remove cells that are outliers for these QC metric
pbmcDAR <- subset(
  x = pbmcDAR,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 10 &
    TSS.enrichment > 2
)


## Normalization and linear dimensional reduction
pbmcDAR <- RunTFIDF(pbmcDAR)
pbmcDAR <- FindTopFeatures(pbmcDAR, min.cutoff = 'q0')
pbmcDAR <- RunSVD(
  object = pbmcDAR,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)
DepthCor(pbmcDAR)


## Non-linear dimension reduction and clustering
pbmcDAR <- RunUMAP(object = pbmcDAR, reduction = 'lsi', dims = 2:30)
pbmcDAR <- FindNeighbors(object = pbmcDAR, reduction = 'lsi', dims = 2:30)
pbmcDAR <- FindClusters(object = pbmcDAR, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmcDAR, label = TRUE) + NoLegend()


## Create a gene activity matrix
# Extend coordinates upstream to include the promoter
genebodyandpromoter.coords <- Extend(x = gene.ranges, upstream = 2000, downstream = 0)

# create a gene by cell matrix
gene.activities <- FeatureMatrix(
  fragments = fragment.path,
  features = genebodyandpromoter.coords,
  cells = colnames(pbmcDAR),
  chunk = 20
)

# convert rownames from chromsomal coordinates into gene names
gene.key <- genebodyandpromoter.coords$gene_name
names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
rownames(gene.activities) <- gene.key[rownames(gene.activities)]

# add the gene activity matrix to the Seurat object as a new assay, and normalize it
pbmcDAR[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmcDAR <- NormalizeData(
  object = pbmcDAR,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmcDAR$nCount_RNA)
)

DefaultAssay(pbmcDAR) <- 'RNA'

FeaturePlot(
  object = pbmcDAR,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)


## Integrating with scRNA-seq data
# Load the pre-processed scRNA-seq data for PBMCs
pbmc_rna <- readRDS("pbmc_10k_v3.rds")

transfer.anchors <- FindTransferAnchors(
  reference = pbmc_rna,
  query = pbmcDAR,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc_rna$celltype,
  weight.reduction = pbmcDAR[['lsi']],
  dims = 2:30
)

pbmcDAR <- AddMetaData(object = pbmcDAR, metadata = predicted.labels)

plot1 <- DimPlot(
  object = pbmc_rna,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = pbmcDAR,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot1 + plot2

pbmcDAR <- subset(pbmcDAR,idents = 14, invert = TRUE)
pbmcDAR <- RenameIdents(
  object = pbmcDAR,
  '0' = 'CD14 Mono',
  '1' = 'CD4 Memory',
  '2' = 'CD4 Naive',
  '3' = 'CD14 Mono',
  '4' = 'CD8 Effector',
  '5' = 'CD14 Mono',
  '6' = 'DN T',
  '7' = 'CD8 Naive',
  '8' = 'NK',
  '9' = 'pre-B',
  '10' = 'CD16 Mono',
  '11' = 'pro-B',
  '12' = 'Dendritic',
  '13' = 'pDC'
)


## Find differentially accessible peaks between clusters
# switch back to working with peaks instead of gene activities
DefaultAssay(pbmcDAR) <- 'peaks'

da_peaks <- FindMarkers(
  object = pbmcDAR,
  ident.1 = "CD14 Mono",
  ident.2 = "CD4 Naive", #for comparing ident1 with another cell type, if not all remaining cells
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

head(da_peaks)

plot1 <- VlnPlot(
  object = pbmcDAR,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = c("CD4 Naive","CD14 Mono")
)
plot2 <- FeaturePlot(
  object = pbmcDAR,
  features = rownames(da_peaks)[1],
  pt.size = 0.1
)

plot1 | plot2

open_cd4naive <- rownames(da_peaks[da_peaks$avg_logFC > 0.5, ])
open_cd14mono <- rownames(da_peaks[da_peaks$avg_logFC < -0.5, ])

closest_genes_cd4memory<- ClosestFeature(
  regions = open_cd4naive,
  annotation = gene.ranges,
  sep = c(':', '-')
)
closest_genes_cd14mono <- ClosestFeature(
  regions = open_cd14mono,
  annotation = gene.ranges,
  sep = c(':', '-')
)

head(closest_genes_cd4memory)
head(closest_genes_cd14mono)

# set plotting order
levels(pbmcDAR) <- c("CD4 Naive","CD4 Memory","CD8 Naive","CD8 Effector","DN T","NK","pre-B",
                     'pro-B',"pDC","Dendritic","CD14 Mono",'CD16 Mono')

CoveragePlot(
  object = pbmcDAR,
  region = rownames(da_peaks)[c(1,5)],
  sep = c(":", "-"),
  peaks = StringToGRanges(rownames(pbmcDAR), sep = c(":", "-")),
  annotation = gene.ranges,
  extend.upstream = 20000,
  extend.downstream = 20000,
  ncol = 1
)








