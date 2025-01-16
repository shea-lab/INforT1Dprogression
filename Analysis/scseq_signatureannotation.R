#29 October 2024
# Jessica King
# scRNAseq data annotation w/ signature genes

library(pacman)
  pacman::p_load(tidyverse, plyr, magrittr, stats, dplyr, limma, RColorBrewer, gplots, glmnet, 
                 biomaRt, colorspace, ggplot2, fmsb, car, mixOmics, DESeq2, apeglm, rgl, qpcR,
                 boot, caret, ggvenn, grid, devtools, reshape2, gridExtra, factoextra, edgeR, 
                 cowplot, pheatmap, coefplot, randomForest, ROCR, genefilter, Hmisc, rdist, 
                 factoextra, ggforce, NormqPCR, ggpubr, matrixStats, GSEAmining, ggrepel,
		     Seurat,SeuratData, SeuratDisk, convert)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")

BiocManager::install("convert")
library(BiocManager)
pacman::p_load(SeuratData, SeuratDisk)

library(SeuratDisk)

pacman::p_load(Seurat,SeuratData, SeuratDisk, Convert)

message("Reading counts...")
x <- read.csv(".\\Data\\scraw-counts.csv",header=TRUE)
rownames(x) <- x[,1]
x[,1] <- NULL
print(dim(x))
print(x[1:5,1:5])

message("Reading metadata...")
m <- read.csv("scmetadata.csv",header=TRUE)
rownames(m) <- m[,1]
colnames(m)[1] <- "sample"
print(dim(m))
print(head(m))

message("Writing seurat object...")
saveRDS(
  CreateSeuratObject(counts=t(x),meta.data=m,project="seurat",min.cells=0,min.features=0),
  "seurat.Rds"
)
seurat_seq <- readRDS('seurat.Rds')

p <- VlnPlot(object = seurat_seq, features =c(wk6genes))
p$data %>% group_by(ident) %>% plyr::summarize(counts = sum(wk6genes, na.rm = TRUE))


## Percent of all cells, and total number of cells expressing the week 6 signature genes
prcts_celltype = PrctCellExpringGene(seurat_seq, wk6genes, group.by = "celltype")
prcts_time = PrctCellExpringGene(seurat_seq, wk6genes, group.by = "Time")
total_celltype = PrctCellExpringGene(seurat_seq, wk6genes, group.by = "celltype")
total_time = PrctCellExpringGene(seurat_seq, wk6genes, group.by = c("Time")

write.csv(total_celltype, "total_celltypes.csv")
write.csv(total_time, "total_times.csv")

VlnPlot(seurat_seq, features = wk6genes)
DotPlot(seurat_seq, features = wk6genes, split.by = "celltype", cols="RdBu") + RotatedAxis()
DoHeatmap(subset(seurat_seq, downsample = 100), features = wk6genes, size = 3)

sum(GetAssayData(object = seurat_seq, slot = "data")[wk6genes[2],]>0)

#Progressor signature
progenes = g_EN.85

## Percent of all cells, and total number of cells expressing the progressor signature genes
prcts_celltype2 = PrctCellExpringGene(seurat_seq, progenes, group.by = "celltype")
prcts_time2 = PrctCellExpringGene(seurat_seq, progenes, group.by = "Time")
total_celltype2 = PrctCellExpringGene(seurat_seq, progenes, group.by = "celltype")
total_time2 = PrctCellExpringGene(seurat_seq, progenes, group.by = "Time")

write.csv(total_celltype2, "total_celltypes_prog.csv")
write.csv(total_time2, "total_times_prog.csv")

VlnPlot(seurat_seq, features = progenes)
DotPlot(seurat_seq, features = progenes, split.by = "celltype", cols="RdBu") + RotatedAxis()
DoHeatmap(subset(seurat_seq, downsample = 100), features = wk6genes, size = 3)
