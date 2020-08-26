pbmc_Counts <- readRDS("pbmcCounts.rds")
pbmc_Binary <- readRDS("pbmcBinary.rds")

pbmc_Counts@active.ident <- pbmc_Binary@active.ident
DimPlot(pbmc_Counts, label = TRUE)

pbmc_Binary@active.ident <- pbmc_Counts@active.ident
DimPlot(pbmc_Binary, label = TRUE)
