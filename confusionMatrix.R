install.packages("caret")
library("caret")
library("dplyr")


pbmc_Counts <- readRDS("pbmcCounts.rds")
pbmc_Binary <- readRDS("pbmcBinary.rds")
cM <- confusionMatrix(pbmc_Counts$seurat_clusters,pbmc_Binary$seurat_clusters)
cM


# extract the confusion matrix values as data.frame
values <- cM$table
sel_values <- values[1:9,1:15]
cM_d <- as.data.frame(sel_values)
# confusion matrix statistics as data.frame
cM_st <-data.frame(cM$overall)
# round the values
cM_st$cM.overall <- round(cM_st$cM.overall,2)

# here we also have the rounded percentage values
cM_p <- as.data.frame(prop.table(cM$table))
cM_d$Perc <- round(cM_p$Freq*100,2)


library(ggplot2)     # to plot
library(gridExtra)   # to put more
library(grid)        # plot together

# plotting the matrix
cM_d_p <-  ggplot(data = cM_d, aes(x = Prediction , y =  Reference, fill = Freq))+
  geom_tile() +
  geom_text(aes(label = paste("",Freq)), color = 'white', size = 4) +
  theme_light() +
  guides(fill=FALSE) 

# plotting the stats
cM_st_p <-  tableGrob(cM_st)

# all together
grid.arrange(cM_d_p, cM_st_p,nrow = 1, ncol = 2, 
             top=textGrob("Confusion Matrix and Statistics",gp=gpar(fontsize=25,font=1)))
