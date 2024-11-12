# Jessica King
# 18 July 2024
# Batch correction and signature analysis of original and validation NOD cohorts

setwd("C:/Users/jlkin/OneDrive/Documents/Research/T1D sensor/6872-JK")

library("pacman")
  pacman::p_load(tidyverse, plyr, magrittr, stats, dplyr, limma, RColorBrewer, gplots, glmnet, 
                 biomaRt, colorspace, ggplot2, fmsb, car, mixOmics, DESeq2, apeglm, rgl, qpcR,
                 boot, caret, ggvenn, grid, devtools, reshape2, gridExtra, factoextra, edgeR, 
                 cowplot, pheatmap, coefplot, randomForest, ROCR, genefilter, Hmisc, rdist, 
                 factoextra, ggforce, NormqPCR, ggpubr, matrixStats, GSEAmining, ggrepel)

#Load data
  #Original cohort
  og_file = read.csv("gene_expected_count.annot.csv", header=T)
  dim(og_file) #55492 genes
  r_file <- og_file[ grep("Gm", og_file[,3], invert = T) , ] # Remove Gm genes
  dim(r_file) #27821 genes
  r_file <- r_file[ grep("Rik", r_file[,3], invert = T) , ] # Remove Riken genes  
  dim(r_file) #24843 genes
  r_file <- r_file[ grep("pseudogene", vr_file[,4], invert = T) , ] # Remove pseudogenes  
  dim(r_file) #23096 genes

#Vector of gene names
genesog = r_file[,3]
genesog = make.names(genesog, unique = T)
rownames(r_file) = genesog

sample_og = c("NODt0m1","NODt0m2","NODt0m3","NODt0m4","NODt0m5","NODt1m1",
	"NODt1m2","NODt1m3","NODt1m4","NODt1m5","NODt2m1","NODt2m2","NODt2m3",
	"NODt2m4","NODt2m5","NODt3m1","NODt3m2","NODt3m3","NODt3m4","NODt3m5",
	"NODt4m1","NODt4m2","NODt4m3","NODt4m4","NODt4m5","NORt0m1","NORt0m2",
	"NORt0m4","NORt0m5","NORt1m1","NORt1m2","NORt1m3","NORt1m4","NORt1m5",
	"NORt2m1","NORt2m2","NORt2m3","NORt2m4","NORt2m5","NORt3m1","NORt3m3",
	"NORt3m4","NORt3m5","NORt4m1","NORt4m2","NORt4m3","NORt4m4","NORt4m5",
	"NSt0m1","NSt0m2","NSt0m3","NSt0m4","NSt1m1","NSt1m2","NSt1m3",
	"NSt1m4","NSt2m1","NSt2m2","NSt2m3","NSt2m4","NSt3m1","NSt3m2","NSt3m3",
	"NSt3m4","NSt4m1","NSt4m2","NSt4m3","NSt4m4","AHMt1m1","AHMt1m2",
	"AHMt2m1","AHMt2m2","AHMt3m1","AHMt3m2")
og_order = c(1,12,22,33,72,53,64,73,74,2,3,4,5,6,7,8,9,10,11,13,14,15,16,17,
	18,19,20,21,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38,39,40,41,42,43,
	44,45,46,47,48,49,50,51,52,54,55,56,57,58,59,60,61,62,63,65,66,67,68,69,
	70,71)

og_counts = r_file[,-(1:4)]
og_counts = og_counts[,og_order]
colnames(og_counts) = sample_og
sample_numog = length(og_counts[1,])

  #Validation cohort
  setwd("C:/Users/jlkin/OneDrive/Documents/Research/T1D sensor/11067-JK")
  vc_file = read.csv("gene_expected_count.annot.csv", header=T)
  dim(vc_file) #55491 genes
  vr_file <- vc_file[ grep("Gm", vc_file[,3], invert = T) , ] # Remove Gm genes
  dim(vr_file) #27820 genes
  vr_file <- vr_file[ grep("Rik", vr_file[,3], invert = T) , ] # Remove Riken genes  
  dim(vr_file) #24842 genes
  vr_file <- vr_file[ grep("pseudogene", vr_file[,4], invert = T) , ] # Remove pseudogenes  
  dim(vr_file) #23096 genes

#Vector of gene names
genesvc = vr_file[,3]
genesvc = make.names(genesvc, unique = T)
rownames(vr_file) = genesvc

sample_vc = c("NODt7bm1","NODt7bm2","NODt7bm3","NODt7bm4",
	"NSt0m5","NSt0m6","NSt0m7","NSt0m8","NSt0m9","NSt0m10","NSt0m11",
	"NSt2m5","NSt2m6","NSt2m7","NSt2m8","NSt2m9","NSt2m10","NODt0m6",
	"NODt0m7","NODt0m8","NODt0m9","NODt0m10","NODt7bm6","NODt7bm7",
	"NODt7bm9","NODt7bm10","NODt1m6","NODt1m7","NODt1m8","NODt1m9",
	"NODt1m10","NODt2m6","NODt2m7","NODt2m8","NODt2m9","NODt2m10","NODt3m6",
	"NODt3m7","NODt3m8","NODt3m9","NODt3m10","NODt4m6","NODt4m7","NODt4m8",
	"NODt4m9","NODt4m10")
vc_order = c(1,12,23,34,42,43,44,45,46,2,3,4,5,6,7,8,9,10,11,13,14,15,16,17,
	18,19,20,21,22,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41)

vc_counts = vr_file[,-(1:4)]
vc_counts = vc_counts[,vc_order]
colnames(vc_counts) = sample_vc
sample_numvc = length(vc_counts[1,])

#Combine datasets
 genes <- Reduce(intersect, list(genesog,genesvc)) # exclude genes not in both
  sample_num <- sample_numog + sample_numvc # Number of total samples
  r_counts <- cbind(og_counts[genes,],vc_counts[genes,])
  r_counts <- r_counts[rowSums(r_counts)>0,] # Remove 0 expression genes
  r_counts <- r_counts[rowSums(r_counts == 0) <= sample_num*(3/4),] # Remove mostly 0 genes
  r_counts <- r_counts[rowSums(r_counts) > sample_num, ] # Remove low expression genes
  genes <- rownames(r_counts)
  samples <- colnames(r_counts)

#Define factors
  batch <- c(rep("Original cohort", 74), rep("Validation Cohort", 46))
  grouped <- c(rep("NOD pre-symptomatic", 20), rep("NOD symptomatic", 5),
	rep("NOR healthy", 23), rep("NOD never-sick", 26), 
	rep("NOD pre-symptomatic", 4), rep("NOD never-sick", 13), 
	rep("NOD pre-symptomatic", 24), rep("NOD symptomatic", 5))
  time <- c(rep("t0", 5), rep("t1", 5), rep("t2", 5), rep("t3", 5), 
	rep("t4", 5), rep("t0", 4), rep("t1", 5), rep("t2", 5), rep("t3", 5),
	rep("t4", 4), rep("t0", 4), rep("t1", 4), rep("t2", 4), rep("t3", 4),
	rep("t4", 4), rep("AMt1", 2), rep("AMt2", 2), rep("AMt3", 2),
	rep("t0.5", 4), rep("t0", 7),rep("t2", 6), rep("t0", 5), rep("t0.5", 4),
	rep("t1", 5), rep("t2", 5), rep("t3", 5), rep("t4", 5))
  mouse <- c(rep(c("NOD1","NOD2","NOD3","NOD4","NOD5"), 5), "NOR1","NOR2",
	"NOR4","NOR5",rep(c("NOR1","NOR2","NOR3","NOR4","NOR5"), 2),"NOR1",
	"NOR3","NOR4","NOR5","NOR1","NOR2","NOR3","NOR4","NOR5",rep(c("NS1",
	"NS2","NS3","NS4"), 5), rep(c("AM1","AM2"),3), "NOD1","NOD2","NOD3",
	"NOD4","NS5","NS6","NS7","NS8","NS9","NS10","NS11","NS5","NS6","NS7",
	"NS8","NS9","NS10","NOD6","NOD7","NOD8","NOD9","NOD10","NOD6","NOD7",
	"NOD9","NOD10",rep(c("NOD6","NOD7","NOD8","NOD9","NOD10"), 4))
  simple <- c(rep("Euglycemic", 20), rep("Hyperglycemic", 5),
	 rep("Euglycemic", 90), rep("Hyperglycemic", 5))
  strain <- c(rep("NOD", 25), rep("NOR", 23), rep("NOD", 72))
  strain_subsets <- c(rep("NOD progressors", 25), rep("NOR", 23), 
	rep("NOD non-progressors", 26), rep("NOD progressors", 4), 
	rep("NOD non-progressors", 13), rep("NOD progressors", 29))

clin_info <- cbind.data.frame( batch, time, mouse, simple, grouped, strain,
	strain_subsets)
rownames(clin_info) <- samples


#Batch correction
 sum(is.na(r_counts)) #0
fac_DEseq <- cbind(grouped, batch) # Combine groups and batch 
    coldata <- data.frame(cbind(fac_DEseq))
    row.names(coldata) <- samples # DESeq2 uses RAW! counts; rows:genes & cols:samples
    dds <- DESeqDataSetFromMatrix(countData=r_counts, colData=coldata, design = ~ batch + grouped) # Used to specify which variables to use in analysis
    paste(nrow(dds), " genes input into DESeq2 for batch correction", sep="")
   dds <- DESeq(dds)
 # deseq_counts <- get_counts(dds)
  vsd <- vst(dds, blind=FALSE)
 # mat <- assay(vsd)
 # mat <- limma::removeBatchEffect(mat, vsd$batch, design = model.matrix( ~grouped))
 # assay(vsd) <- mat
  bc_counts <- assay(vsd)
  bc_counts <- as.data.frame(bc_counts)
  sample.diff(bc_counts, samples) # Likeness of samples, function in 30AugRRU

#Add filters
ps_counts <- log(bc_counts + 1, base = 2) # Transform to pseudo-normalized
filt_low <- 2 # These three filters are per sample per gene
filt_high <- 16
keep1 <- rowSums(ps_counts) > (filt_low*sample_num) # Low count filter, required
ps_filt <- ps_counts[keep1,] # Apply filter to pseudonorm counts
dim(ps_filt) #2:18269, 5:17179
keep2 <- rowSums(ps_filt) < (filt_high*sample_num) # High count filter, optional
ps_filt <- ps_filt[keep2,] # Apply filter to pseudonorm counts
dim(ps_filt) #16:2773, 30:3082, n=25 11525, n=36 high filts 14985/low filts 7104

#16499 bc2
bc_filt <- bc_counts[rownames(ps_filt),]
 counts_filt = bc_filt
 # Filter by variance across all samples to exclude genes that don't change much
    mn_cut <- 0.2 # Mean cutoff
    var_cut <- 0.01 # Variance cutoff
    df <- cbind.data.frame("x" = rowMeans(counts_filt), # Calculate Mean & Var/Mean
                           "y" = rowVars(counts_filt) / rowMeans(counts_filt)) 
    g_var_all <- rownames(subset(df, x > mn_cut & y > var_cut)) # Genes to keep
    # Plot variance and mean
    par(mfrow = c(1, 1)) # Set the layout to a single plot in a 1x1 grid
    with(df, plot(x, y, main = paste("All Var vs Mean, ", length(g_var_all), 
                                     " genes w/ Mean > ", mn_cut, " & Var > ", var_cut, sep = ""), 
                  pch = 20, cex = 1, xlab = "Mean across all samples", 
                  ylab = "Var/Mean across All samples"))
    with(subset(df, x > mn_cut & y > var_cut), 
         points(x, y, pch = 20, col = "red3", cex = 1)) # Genes above cutoff are red
    abline(v = mn_cut, col = "red3", lty = 2, lwd = 1.5) # Line for mean cutoff
    abline(h = var_cut, col = "red3", lty = 2, lwd = 1.5) # Line for variance cutoff
    removed <- nrow(counts_filt) - length(g_var_all) # Number of filtered genes
    print(paste(removed, "genes removed due to expression below the variance filter"))
    counts_filt <- counts_filt[g_var_all,]  # Update counts_filt 
#    if (length(unique(groups)) == 2) { # Additional variance by group 
#      num_groups <- length(unique(groups)) # Determine the number of groups 
      # ++++ This needs to be finished still
      
#      removed <- nrow(counts_filt) - length(g_var_group) # Number of filtered genes 
#      print(paste(removed, " genes removed by the pairwise group variance filter")) 
#      counts_filt <- counts_filt[g_var_group,]  # Update counts_filt 

bc_filt = counts_filt
#Elastic net

#Figure 2: NOR
#Factors
groups = factor(strain, labels = c(0,1)) #0: NOD, 1: NOR

#Colors for grouped
coul <- colorRampPalette(brewer.pal(11, "PuOr"))(25)
ecolor <- c("skyblue", "seagreen3")
colSide <- clin_info$strain_subsets
for (x in 1:sample_num) {
  if (strain[x] == "NOR") {colSide[x] <- "seagreen3"  
  } else if (strain[x] == "NOD") {colSide[x] <- "skyblue"
  }}
colind = colSide
names(colind) = c(rep("NOD", 25), rep("NOR", 23), rep("NOD", 72))




#FlexiDEG
NORsamples = colnames(bc_filt)
metadata = cbind(NORsamples, groups)
colnames(metadata) = c("Samples", "Group")
metadata <- as.data.frame(metadata)
pdf(file= "Flexi, NOD vs NOR, bc2, var filt.pdf" )
NORflexi <- flexiDEG.function4(bc_filt, metadata, validation_option = 2)

EN1 <- na.omit(bc_filt[unique(rownames(bc_filt)[as_vector(NORflexi[[1]])]), ]) 
EN2 <- na.omit(bc_filt[unique(rownames(bc_filt)[as_vector(NORflexi[[2]])]), ]) 
EN3 <- na.omit(bc_filt[unique(rownames(bc_filt)[as_vector(NORflexi[[3]])]), ]) 

heatmap.2(as.matrix(EN1), scale="row", col=coul, key= T, xlab="", ylab="", 
          margins=c(7,7), ColSideColors=colind, trace="none", key.title=NA, 
          key.ylab=NA, keysize=0.8, dendrogram="both")
ggbiplot(prcomp(t(EN1), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
         var.scale=1, circle=T) + 
  theme_classic() + scale_color_manual(name="Group", values=ecolor) 
heatmap.2(as.matrix(EN2), scale="row", col=coul, key= T, xlab="", ylab="", 
          margins=c(7,7), ColSideColors=colind, trace="none", key.title=NA, 
          key.ylab=NA, keysize=0.8, dendrogram="both")
ggbiplot(prcomp(t(EN2), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
         var.scale=1, circle=T) + 
  theme_classic() + scale_color_manual(name="Group", values=ecolor) 

heatmap.2(as.matrix(EN3), scale="row", col=coul, key= T, xlab="", ylab="", 
          margins=c(7,7), ColSideColors=colind, trace="none", key.title=NA, 
          key.ylab=NA, keysize=0.8, dendrogram="both")
ggbiplot(prcomp(t(EN3), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
         var.scale=1, circle=T) + 
  theme_classic() + scale_color_manual(name="Group", values=ecolor) 
dev.off()

#Grouped EN, NOR vs NOD and NS
  xfactors <- cbind.data.frame(groups) # SELECT THE FACTOR YOU WANT!!!!!
samplesbc = colnames(bc_filt)
  rownames(xfactors) <- samplesbc
  dataset <- t(bc_filt) # Load the data. Should it be raw or normalized!?!?!?!
  dataset <- as.data.frame(cbind(xfactors, dataset)) # Add factor as the first column
  dataset <- dataset[complete.cases(dataset), ] # For cases without missed items
  set.seed(123) # Split the data into training and test set
  training.samples <- dataset$groups %>% createDataPartition(p = 1.0, list = F)
  train.data  <- dataset[training.samples, ]
  test.data <- dataset[-training.samples, ]
  x <- model.matrix(groups~., train.data)[,-1] # Predictor variables; need reduced gene set to run, runs w/ 9814
  y <- train.data$groups # Outcome variable
  genes_EN <- binom_EN( x, y, 100) # number following multinom is number of groups
NODNORbc2var = genes_EN
batchcorr_counts = bc_filt
  g_EN1 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`1`])
  g_EN.95 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.95`])
  g_EN.9 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.9`])
  g_EN.85 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.85`])
  g_EN.8 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.8`])
  g_EN.75 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.75`])
  g_EN.7 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.7`])
  g_EN.65 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.65`])
  g_EN.6 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.6`])
  g_EN.55 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.55`])
  g_EN.5 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.5`])
  g_EN.45 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.45`])
  g_EN.4 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.4`])
  g_EN.35 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.35`])
  g_EN.3 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.3`])
  g_EN.25 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.25`])
  g_EN.2 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.2`])
  g_EN.15 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.15`])
  g_EN.1 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.1`])
  g_EN.05 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.0499999999999999`])
  g_EN0 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0`])
  n_EN1 <- na.omit(batchcorr_counts[g_EN1, ])
  n_EN.95 <- na.omit(batchcorr_counts[g_EN.95, ])
  n_EN.9 <- na.omit(batchcorr_counts[g_EN.9, ])
  n_EN.85 <- na.omit(batchcorr_counts[g_EN.85, ])
  n_EN.8 <- na.omit(batchcorr_counts[g_EN.8, ])
  n_EN.75 <- na.omit(batchcorr_counts[g_EN.75, ])
  n_EN.7 <- na.omit(batchcorr_counts[g_EN.7, ])
  n_EN.65 <- na.omit(batchcorr_counts[g_EN.65, ])
  n_EN.6 <- na.omit(batchcorr_counts[g_EN.6, ])
  n_EN.55 <- na.omit(batchcorr_counts[g_EN.55, ])
  n_EN.5 <- na.omit(batchcorr_counts[g_EN.5, ])
  n_EN.45 <- na.omit(batchcorr_counts[g_EN.45, ])
  n_EN.4 <- na.omit(batchcorr_counts[g_EN.4, ])
  n_EN.35 <- na.omit(batchcorr_counts[g_EN.35, ])
  n_EN.3 <- na.omit(batchcorr_counts[g_EN.3, ])
  n_EN.25 <- na.omit(batchcorr_counts[g_EN.25, ])
  n_EN.2 <- na.omit(batchcorr_counts[g_EN.2, ])
  n_EN.15 <- na.omit(batchcorr_counts[g_EN.15, ])
  n_EN.1 <- na.omit(batchcorr_counts[g_EN.1, ])
  n_EN.05 <- na.omit(batchcorr_counts[g_EN.05, ])
  n_EN0 <- na.omit(batchcorr_counts[g_EN0, ])
pdf(file= "Grouped EN, NOD, NOR combined, bc2, var filt.pdf" )
  heatmap.2(as.matrix(n_EN1), scale = "row", col = coul, key = T, 
            xlab = "Sample", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN1), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle("STx EN alpha=1 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.95), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.95), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle("STx EN alpha=.95 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.9), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.9), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle("STx EN alpha=.9 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.85), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.85), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle("STx EN alpha=.85 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.8), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.8), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle("STx EN alpha=.8 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.75), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.75), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle("STx EN alpha=.75 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.7), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.7), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle("STx EN alpha=.7 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.65), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.65), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle("STx EN alpha=.65 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.6), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.6), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle("STx EN alpha=.6 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.55), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.55), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle("STx EN alpha=.55 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.5), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.5), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle("STx EN alpha=.5 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.45), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.45), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle("STx EN alpha=.45 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.4), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.4), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle("STx EN alpha=.4 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.35), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.35), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle("STx EN alpha=.35 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.3), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.3), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle("STx EN alpha=.3 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.25), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.25), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle("STx EN alpha=.25 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.2), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.2), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle("STx EN alpha=.2 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.15), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.15), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle("STx EN alpha=.15 Geneset") + 
    geom_point(size=3, color=colind) + 
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.1), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.1), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle("STx EN alpha=.1 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.05), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.05), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle("STx EN alpha=.05 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  g_Groupedfull <- g_EN.05
  g_Groupedlasso <- g_EN1
dev.off() #stop adding to PDF file

#Factors
NOR = bc_filt[ ,strain_subsets=="NOR"]
NS = bc_filt[ ,strain_subsets=="NOD non-progressors"]
NOR_NS = cbind.data.frame(NOR, NS)
  clin_info_nopd <- clin_info[ grep("NOD progressors", clin_info[,"strain_subsets"], invert = T) , ] #Remove progressors
groupedfNORNS = factor(clin_info_nopd[,"strain_subsets"], labels = c(0,1)) #0: NOD, 1:NS, 2:NOR
sample_num_nopd = length(colnames(NOR_NS))


#Colors for grouped
coul <- colorRampPalette(brewer.pal(11, "PuOr"))(25)
ecolor <- c("gold", "seagreen3")
colSide <- clin_info_nopd$strain_subsets
for (x in 1:sample_num_nopd) {
  if (colSide[x] == "NOR") {colSide[x] <- "seagreen3"  
  } else if (colSide[x] == "NOD non-progressors") {colSide[x] <- "gold"
  }}
colind = colSide
names(colind) = c(rep("NOR", 23), rep("NOD non-progressors", 39))

#Colors for test group
nopd_sub = clin_info_nopd[NORNS.test, ]
colSide <- nopd_sub$strain_subsets
for (x in 1:length(colSide)) {
  if (colSide[x] == "NOR") {colSide[x] <- "seagreen3"  
  } else if (colSide[x] == "NOD non-progressors") {colSide[x] <- "gold"
  }}
colind = colSide
names(colind) = c(rep("NOR", 11), rep("NOD non-progressors", 19))

#Grouped EN, NOR vs NS

NORNSsamples = colnames(NOR_NS)
groups = groupedfNORNS
metadata = cbind(NORNSsamples, groups)
colnames(metadata) = c("Samples", "Group")
metadata <- as.data.frame(metadata)
NORNSflexi <- flexiDEG.function4(NOR_NS, metadata, validation_option = 2)
EN4 <- na.omit(NOR_NS[unique(rownames(NOR_NS)[as_vector(NORNSflexi[[1]])]), ]) 
EN5 <- na.omit(NOR_NS[unique(rownames(NOR_NS)[as_vector(NORNSflexi[[2]])]), ]) 
EN6 <- na.omit(NOR_NS[unique(rownames(NOR_NS)[as_vector(NORNSflexi[[3]])]), ]) 

heatmap.2(as.matrix(EN6), scale="row", col=coul, key= T, xlab="", ylab="", 
          margins=c(7,7), ColSideColors=colind, trace="none", key.title=NA, 
          key.ylab=NA, keysize=0.8, dendrogram="both")
ggbiplot(prcomp(t(EN6), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
         var.scale=1, circle=T) + 
  theme_classic() + scale_color_manual(name="Group", values=ecolor)

  xfactors <- cbind.data.frame(groupedfNORNS) 
samplesbc = colnames(NOR_NS)
  rownames(xfactors) <- samplesbc
  dataset <- t(NOR_NS) # Load the data. Should it be raw or normalized!?!?!?!
  dataset <- as.data.frame(cbind(xfactors, dataset)) # Add factor as the first column
  dataset <- dataset[complete.cases(dataset), ] # For cases without missed items
  set.seed(123) # Split the data into training and test set
  training.samples <- dataset$groupedfNORNS %>% createDataPartition(p = 0.5, list = F)
  train.data  <- dataset[training.samples, ]
  test.data <- dataset[-training.samples, ]
  x <- model.matrix(groupedfNORNS~., train.data)[,-1] # Predictor variables; need reduced gene set to run, runs w/ 9814
  y <- train.data$groupedfNORNS # Outcome variable
  genes_EN <- binom_EN( x, y, 100) # number following multinom is number of groups
NORvNSENboot = genes_EN
NORNS.test = rownames(test.data)
test = NOR_NS[ , NORNS.test]
batchcorr_counts = test
  g_EN1 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`1`])
  g_EN.95 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.95`])
  g_EN.9 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.9`])
  g_EN.85 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.85`])
  g_EN.8 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.8`])
  g_EN.75 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.75`])
  g_EN.7 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.7`])
  g_EN.65 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.65`])
  g_EN.6 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.6`])
  g_EN.55 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.55`])
  g_EN.5 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.5`])
  g_EN.45 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.45`])
  g_EN.4 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.4`])
  g_EN.35 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.35`])
  g_EN.3 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.3`])
  g_EN.25 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.25`])
  g_EN.2 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.2`])
  g_EN.15 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.15`])
  g_EN.1 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.1`])
  g_EN.05 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.0499999999999999`])
  g_EN0 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0`])
  n_EN1 <- na.omit(batchcorr_counts[g_EN1, ])
  n_EN.95 <- na.omit(batchcorr_counts[g_EN.95, ])
  n_EN.9 <- na.omit(batchcorr_counts[g_EN.9, ])
  n_EN.85 <- na.omit(batchcorr_counts[g_EN.85, ])
  n_EN.8 <- na.omit(batchcorr_counts[g_EN.8, ])
  n_EN.75 <- na.omit(batchcorr_counts[g_EN.75, ])
  n_EN.7 <- na.omit(batchcorr_counts[g_EN.7, ])
  n_EN.65 <- na.omit(batchcorr_counts[g_EN.65, ])
  n_EN.6 <- na.omit(batchcorr_counts[g_EN.6, ])
  n_EN.55 <- na.omit(batchcorr_counts[g_EN.55, ])
  n_EN.5 <- na.omit(batchcorr_counts[g_EN.5, ])
  n_EN.45 <- na.omit(batchcorr_counts[g_EN.45, ])
  n_EN.4 <- na.omit(batchcorr_counts[g_EN.4, ])
  n_EN.35 <- na.omit(batchcorr_counts[g_EN.35, ])
  n_EN.3 <- na.omit(batchcorr_counts[g_EN.3, ])
  n_EN.25 <- na.omit(batchcorr_counts[g_EN.25, ])
  n_EN.2 <- na.omit(batchcorr_counts[g_EN.2, ])
  n_EN.15 <- na.omit(batchcorr_counts[g_EN.15, ])
  n_EN.1 <- na.omit(batchcorr_counts[g_EN.1, ])
  n_EN.05 <- na.omit(batchcorr_counts[g_EN.05, ])
  n_EN0 <- na.omit(batchcorr_counts[g_EN0, ])
pdf(file= "Grouped NS and NOR all times, .5 data partition.pdf" )
  heatmap.2(as.matrix(n_EN1), scale = "row", col = coul, key = T, 
            xlab = "Sample", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN1), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=1 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.95), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.95), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.95 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.9), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.9), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.9 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.85), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.85), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.85 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.8), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.8), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.8 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.75), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.75), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.75 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.7), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.7), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.7 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.65), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.65), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.65 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.6), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.6), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.6 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.55), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.55), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.55 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.5), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.5), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.5 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.45), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.45), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.45 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.4), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.4), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.4 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.35), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.35), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.35 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.3), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.3), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.3 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.25), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.25), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.25 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.2), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.2), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.2 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.15), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.15), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.15 Geneset") + 
    geom_point(size=3, color=colind) + 
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.1), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.1), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.1 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.05), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.05), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.05 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  g_Groupedfull <- g_EN.05
  g_Groupedlasso <- g_EN1
dev.off() #stop adding to PDF file

#NOR timecourse
sample_num_NOR = length(colnames(NOR))
  clin_info_NOR <- clin_info[ grep("NOD", clin_info[,"strain_subsets"], invert = T) , ] #Remove NOD
groupedfNOR = factor(clin_info_NOR[,"time"], labels = c(0,1,2,3,4)) 

#Colors for grouped
coul <- colorRampPalette(brewer.pal(11, "PuOr"))(25)
ecolor <- c("palegreen","lightseagreen","mediumspringgreen","seagreen3","seagreen4")
colSide <- clin_info_NOR$time
for (x in 1:sample_num_NOR) {
  if (colSide[x] == "t0") {colSide[x] <- "palegreen"  
  } else if (colSide[x] == "t1") {colSide[x] <- "lightseagreen"
  } else if (colSide[x] == "t2") {colSide[x] <- "mediumspringgreen"
  } else if (colSide[x] == "t3") {colSide[x] <- "seagreen3"
  } else if (colSide[x] == "t4") {colSide[x] <- "seagreen4"
  }}
colind = colSide
names(colind) = c(rep("t0", 4), rep("t1", 5), rep("t2", 5), rep("t3", 5), rep("t4", 4))


#Grouped EN, NOR timecourse
  xfactors <- cbind.data.frame(groupedfNOR) 
samplesbc = colnames(NOR)
  rownames(xfactors) <- samplesbc
  dataset <- t(NOR) # Load the data. Should it be raw or normalized!?!?!?!
  dataset <- as.data.frame(cbind(xfactors, dataset)) # Add factor as the first column
  dataset <- dataset[complete.cases(dataset), ] # For cases without missed items
  set.seed(123) # Split the data into training and test set
  training.samples <- dataset$groupedfNOR %>% createDataPartition(p = 1.0, list = F)
  train.data  <- dataset[training.samples, ]
  test.data <- dataset[-training.samples, ]
  x <- model.matrix(groupedfNOR~., train.data)[,-1] # Predictor variables; need reduced gene set to run, runs w/ 9814
  y <- train.data$groupedfNOR # Outcome variable
  genes_EN <- multinom5_EN( x, y, 100) # number following multinom is number of groups
NORtime = genes_EN
batchcorr_counts = NOR
  g_EN1 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`1`])
  g_EN.95 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.95`])
  g_EN.9 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.9`])
  g_EN.85 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.85`])
  g_EN.8 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.8`])
  g_EN.75 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.75`])
  g_EN.7 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.7`])
  g_EN.65 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.65`])
  g_EN.6 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.6`])
  g_EN.55 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.55`])
  g_EN.5 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.5`])
  g_EN.45 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.45`])
  g_EN.4 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.4`])
  g_EN.35 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.35`])
  g_EN.3 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.3`])
  g_EN.25 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.25`])
  g_EN.2 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.2`])
  g_EN.15 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.15`])
  g_EN.1 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.1`])
  g_EN.05 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.0499999999999999`])
  g_EN0 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0`])
  n_EN1 <- na.omit(batchcorr_counts[g_EN1, ])
  n_EN.95 <- na.omit(batchcorr_counts[g_EN.95, ])
  n_EN.9 <- na.omit(batchcorr_counts[g_EN.9, ])
  n_EN.85 <- na.omit(batchcorr_counts[g_EN.85, ])
  n_EN.8 <- na.omit(batchcorr_counts[g_EN.8, ])
  n_EN.75 <- na.omit(batchcorr_counts[g_EN.75, ])
  n_EN.7 <- na.omit(batchcorr_counts[g_EN.7, ])
  n_EN.65 <- na.omit(batchcorr_counts[g_EN.65, ])
  n_EN.6 <- na.omit(batchcorr_counts[g_EN.6, ])
  n_EN.55 <- na.omit(batchcorr_counts[g_EN.55, ])
  n_EN.5 <- na.omit(batchcorr_counts[g_EN.5, ])
  n_EN.45 <- na.omit(batchcorr_counts[g_EN.45, ])
  n_EN.4 <- na.omit(batchcorr_counts[g_EN.4, ])
  n_EN.35 <- na.omit(batchcorr_counts[g_EN.35, ])
  n_EN.3 <- na.omit(batchcorr_counts[g_EN.3, ])
  n_EN.25 <- na.omit(batchcorr_counts[g_EN.25, ])
  n_EN.2 <- na.omit(batchcorr_counts[g_EN.2, ])
  n_EN.15 <- na.omit(batchcorr_counts[g_EN.15, ])
  n_EN.1 <- na.omit(batchcorr_counts[g_EN.1, ])
  n_EN.05 <- na.omit(batchcorr_counts[g_EN.05, ])
  n_EN0 <- na.omit(batchcorr_counts[g_EN0, ])
pdf(file= "Grouped NOR time course.pdf" )
  heatmap.2(as.matrix(n_EN1), scale = "row", col = coul, key = T, 
            xlab = "Sample", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN1), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=1 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.95), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.95), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.95 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.9), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.9), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.9 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.85), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.85), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.85 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.8), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.8), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.8 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.75), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.75), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.75 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.7), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.7), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.7 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.65), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.65), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.65 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.6), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.6), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.6 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.55), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.55), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.55 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.5), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.5), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.5 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.45), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.45), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.45 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.4), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.4), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.4 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.35), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.35), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.35 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.3), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.3), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.3 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.25), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.25), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.25 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.2), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.2), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.2 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.15), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.15), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.15 Geneset") + 
    geom_point(size=3, color=colind) + 
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.1), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.1), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.1 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.05), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.05), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.05 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  g_Groupedfull <- g_EN.05
  g_Groupedlasso <- g_EN1
dev.off() #stop adding to PDF file

#Filter NOR time course genes out of NOD data frame
NORtc_genes = g_EN.9 #39 genes (.85 in bc1)
NOD <- bc_filt[ ,strain=="NOD"]
NOD_genes = rownames(NOD)
NOD_filt = NOD[setdiff(NOD_genes, NORtc_genes),]
clin_info_NOD = clin_info[ grep("NOR", clin_info[,"strain_subsets"], invert = T) , ] #Remove NOR

#NOD week 6 analysis
tNOD = clin_info_NOD$time
wk6 = NOD_filt[ ,tNOD=="t0"] 
sample_num_wk6 = length(colnames(wk6))
clin_info_wk6 = clin_info_NOD[colnames(wk6) ,]

#Colors for all samples
coul <- colorRampPalette(brewer.pal(11, "PuOr"))(25)
ecolor <- c("gold", "skyblue")
colSide <- clin_info_wk6$strain_subsets
for (x in 1:sample_num_wk6) {
  if (colSide[x] == "NOD progressors") {colSide[x] <- "skyblue" 
  } else if (colSide[x] == "NOD non-progressors") {colSide[x] <- "gold"
  }}
colind = colSide
names(colind) = c(rep("NOD progressors", 5), rep("NOD non-progressors", 11), 
	rep("NOD progressors", 5))

#Factors
groups = factor(clin_info_wk6$strain_subsets, labels = c(0,1)) #1: NOD, 0:NS
#groups= simplef
NODwk6samples = colnames(wk6)
metadata = cbind(NODwk6samples, groups)
colnames(metadata) = c("Samples", "Group")
metadata <- as.data.frame(metadata)
pdf(file= "NOD week 6, bc2, 50x.pdf" )
NODwk6flexi <- flexiDEG.function4(wk6, metadata, validation_option = 2)
EN7 <- na.omit(wk6[unique(rownames(wk6)[as_vector(NODwk6flexi[[1]])]), ]) 
EN8 <- na.omit(wk6[unique(rownames(wk6)[as_vector(NODwk6flexi[[2]])]), ]) 
EN9 <- na.omit(wk6[unique(rownames(wk6)[as_vector(NODwk6flexi[[3]])]), ]) 

#50x
EN7.2 <- na.omit(wk6[unique(rownames(wk6)[as_vector(NODwk6flexi[[1]])]), ]) 
EN8.2 <- na.omit(wk6[unique(rownames(wk6)[as_vector(NODwk6flexi[[2]])]), ]) 
EN9.2 <- na.omit(wk6[unique(rownames(wk6)[as_vector(NODwk6flexi[[3]])]), ]) 

#5x
EN10.3 <- na.omit(wk6[unique(rownames(wk6)[as_vector(NODwk6flexi[[1]])]), ]) 
EN11.3 <- na.omit(wk6[unique(rownames(wk6)[as_vector(NODwk6flexi[[2]])]), ]) 
EN12.3 <- na.omit(wk6[unique(rownames(wk6)[as_vector(NODwk6flexi[[3]])]), ]) 


heatmap.2(as.matrix(EN9.3), scale="row", col=coul, key= T, xlab="", ylab="", 
          margins=c(7,7), ColSideColors=colind, trace="none", key.title=NA, 
          key.ylab=NA, keysize=0.8, dendrogram="both")
ggbiplot(prcomp(t(EN9.3), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
         var.scale=1, circle=T) + 
  theme_classic() + scale_color_manual(name="Group", values=ecolor)

heatmap.2(as.matrix(EN8.2), scale="row", col=coul, key= T, xlab="", ylab="", 
          margins=c(7,7), ColSideColors=colind, trace="none", key.title=NA, 
          key.ylab=NA, keysize=0.8, dendrogram="both")
ggbiplot(prcomp(t(EN8.2), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
         var.scale=1, circle=T) + 
  theme_classic() + scale_color_manual(name="Group", values=ecolor)

heatmap.2(as.matrix(EN9.2), scale="row", col=coul, key= T, xlab="", ylab="", 
          margins=c(7,7), ColSideColors=colind, trace="none", key.title=NA, 
          key.ylab=NA, keysize=0.8, dendrogram="both")
ggbiplot(prcomp(t(EN9.2), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
         var.scale=1, circle=T) + 
  theme_classic() + scale_color_manual(name="Group", values=ecolor)

dev.off()


#Colors for test group
wk6.test = rownames(test.data)
wk6_sub = clin_info_wk6[wk6.test, ]
colSide <- wk6_sub$strain_subsets
for (x in 1:length(colSide)) {
  if (colSide[x] == "NOD progressors") {colSide[x] <- "skyblue"  
  } else if (colSide[x] == "NOD non-progressors") {colSide[x] <- "gold"
  }}
colind = colSide
names(colind) = c(rep("NOD progressors", 2), rep("NOD non-progressors", 4),
	rep("NOD progressors", 2))

#Simple EN, NOD vs NS week 6 all run in signature
  xfactors <- cbind.data.frame(groups) 
samplesbc = colnames(wk6)
  rownames(xfactors) <- samplesbc
  dataset <- t(wk6) # Load the data. Should it be raw or normalized!?!?!?!
  dataset <- as.data.frame(cbind(xfactors, dataset)) # Add factor as the first column
  dataset <- dataset[complete.cases(dataset), ] # For cases without missed items
  set.seed(123) # Split the data into training and test set
  training.samples <- dataset$groups %>% createDataPartition(p = 1, list = F)
  train.data  <- dataset[training.samples, ]
  test.data <- dataset[-training.samples, ]
  x <- model.matrix(groups~., train.data)[,-1] # Predictor variables; need reduced gene set to run, runs w/ 9814
  y <- train.data$groups # Outcome variable
  genes_EN <- binom_EN( x, y, 100) # number following multinom is number of groups
NODvNSwk6bc2 = genes_EN
 wk6.test = rownames(test.data)
 test = wk6[ , wk6.test]
 batchcorr_counts = test
batchcorr_counts = wk6
  g_EN1 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`1`])
  g_EN.95 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.95`])
  g_EN.9 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.9`])
  g_EN.85 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.85`])
  g_EN.8 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.8`])
  g_EN.75 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.75`])
  g_EN.7 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.7`])
  g_EN.65 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.65`])
  g_EN.6 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.6`])
  g_EN.55 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.55`])
  g_EN.5 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.5`])
  g_EN.45 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.45`])
  g_EN.4 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.4`])
  g_EN.35 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.35`])
  g_EN.3 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.3`])
  g_EN.25 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.25`])
  g_EN.2 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.2`])
  g_EN.15 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.15`])
  g_EN.1 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.1`])
  g_EN.05 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.0499999999999999`])
  g_EN0 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0`])
  n_EN1 <- na.omit(batchcorr_counts[g_EN1, ])
  n_EN.95 <- na.omit(batchcorr_counts[g_EN.95, ])
  n_EN.9 <- na.omit(batchcorr_counts[g_EN.9, ])
  n_EN.85 <- na.omit(batchcorr_counts[g_EN.85, ])
  n_EN.8 <- na.omit(batchcorr_counts[g_EN.8, ])
  n_EN.75 <- na.omit(batchcorr_counts[g_EN.75, ])
  n_EN.7 <- na.omit(batchcorr_counts[g_EN.7, ])
  n_EN.65 <- na.omit(batchcorr_counts[g_EN.65, ])
  n_EN.6 <- na.omit(batchcorr_counts[g_EN.6, ])
  n_EN.55 <- na.omit(batchcorr_counts[g_EN.55, ])
  n_EN.5 <- na.omit(batchcorr_counts[g_EN.5, ])
  n_EN.45 <- na.omit(batchcorr_counts[g_EN.45, ])
  n_EN.4 <- na.omit(batchcorr_counts[g_EN.4, ])
  n_EN.35 <- na.omit(batchcorr_counts[g_EN.35, ])
  n_EN.3 <- na.omit(batchcorr_counts[g_EN.3, ])
  n_EN.25 <- na.omit(batchcorr_counts[g_EN.25, ])
  n_EN.2 <- na.omit(batchcorr_counts[g_EN.2, ])
  n_EN.15 <- na.omit(batchcorr_counts[g_EN.15, ])
  n_EN.1 <- na.omit(batchcorr_counts[g_EN.1, ])
  n_EN.05 <- na.omit(batchcorr_counts[g_EN.05, ])
  n_EN0 <- na.omit(batchcorr_counts[g_EN0, ])
pdf(file= "NOD week 6, bc2.pdf" )
  heatmap.2(as.matrix(n_EN1), scale = "row", col = coul, key = T, 
            xlab = "Sample", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN1), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=1 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.95), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.95), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.95 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.9), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.9), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.9 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.85), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.85), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.85 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.8), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.8), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.8 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.75), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.75), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.75 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.7), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.7), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.7 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.65), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.65), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.65 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.6), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.6), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.6 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.55), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.55), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic(base_size=15) + ggtitle(" EN alpha=.55 Geneset") + 
    theme(legend.position= c(0.75, 0.85)) + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.5), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.5), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.5 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.45), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.45), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.45 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.4), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.4), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.4 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.35), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.35), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.35 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.3), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.3), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.3 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.25), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.25), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.25 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.2), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.2), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.2 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.15), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.15), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.15 Geneset") + 
    geom_point(size=3, color=colind) + 
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.1), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.1), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.1 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.05), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.05), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.05 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  g_Groupedfull <- g_EN.05
  g_Groupedlasso <- g_EN1
dev.off() #stop adding to PDF file

#NOD timecourse
strainNOD = clin_info_NOD$strain_subsets
NODp = NOD_filt[,strainNOD=="NOD progressors"]
  clin_info_NODp <- clin_info[ grep("NOD progressors", clin_info[,"strain_subsets"], invert = F) , ] #Remove NOD
times = clin_info_NODp$time
mID = clin_info_NODp$mouse 
m1 = NODp[, mID=="NOD1"]
dm1 = m1-m1[,1]
m2 = NODp[, mID=="NOD2"]
dm2 = m2-m2[,1]
m3 = NODp[, mID=="NOD3"]
dm3 = m3-m3[,1]
m4 = NODp[, mID=="NOD4"]
dm4 = m4-m4[,1]
m5 = NODp[, mID=="NOD5"]
dm5 = m5-m5[,1]
m6 = NODp[, mID=="NOD6"]
dm6 = m6-m6[,1]
m7 = NODp[, mID=="NOD7"]
dm7 = m7-m7[,1]
m8 = NODp[, mID=="NOD8"]
dm8 = m8-m8[,1]
m9 = NODp[, mID=="NOD9"]
dm9 = m9-m9[,1]
m10 = NODp[, mID=="NOD10"]
dm10 = m10-m10[,1]

dmall = cbind.data.frame(dm1,dm2,dm3,dm4,dm5,dm6,dm7,dm8,dm9,dm10)
mtIDs = colnames(dmall)

NODnp = NOD_filt[,strainNOD=="NOD non-progressors"]
  clin_info_NODnp <- clin_info[ grep("NOD non-progressors", clin_info[,"strain_subsets"], invert = F) , ] #Remove NOD progressors
times = clin_info_NODnp$time
mID = clin_info_NODnp$mouse 
m1 = NODnp[, mID=="NS1"]
ndm1 = m1-m1[,1]
m2 = NODnp[, mID=="NS2"]
ndm2 = m2-m2[,1]
m3 = NODnp[, mID=="NS3"]
ndm3 = m3-m3[,1]
m4 = NODnp[, mID=="NS4"]
ndm4 = m4-m4[,1]
m5 = NODnp[, mID=="NS5"]
ndm5 = m5-m5[,1]
m6 = NODnp[, mID=="NS6"]
ndm6 = m6-m6[,1]
m7 = NODnp[, mID=="NS7"]
ndm7 = m7-m7[,1]
m8 = NODnp[, mID=="NS8"]
ndm8 = m8-m8[,1]
m9 = NODnp[, mID=="NS9"]
ndm9 = m9-m9[,1]
m10 = NODnp[, mID=="NS10"]
ndm10 = m10-m10[,1]

ndmall = cbind.data.frame(ndm1,ndm2,ndm3,ndm4,ndm5,ndm6,ndm7,ndm8,ndm9,ndm10)
ndmtIDs = colnames(ndmall)

#Delta relative to t0
dD.5 = dmall[,grep("t7b", mtIDs, invert = F)]
dD1 = dmall[,grep("t1", mtIDs, invert = F)]
dD2 = dmall[,grep("t2", mtIDs, invert = F)]
dD3 = dmall[,grep("t3", mtIDs, invert = F)]
dD4 = dmall[,grep("t4", mtIDs, invert = F)]
deltas = cbind.data.frame(dD.5,dD1,dD2,dD3,dD4)
deltimes= c(rep("t0.5",length(dD.5[1,])),rep("t1",length(dD1[1,])),
	rep("t2",length(dD2[1,])),rep("t3",length(dD3[1,])),rep("t4",length(dD4[1,])))
timegroups = c(rep("early",length(dD.5[1,])),rep("early",length(dD1[1,])),
	rep("late",length(dD2[1,])),rep("late",length(dD3[1,])),rep("diabetic",length(dD4[1,])))
groupedfNODdel = factor(deltimes, labels = c(0,1,2,3,4)) 

dND1 = ndmall[,grep("t1", ndmtIDs, invert = F)]
dND2 = ndmall[,grep("t2", ndmtIDs, invert = F)]
dND3 = ndmall[,grep("t3", ndmtIDs, invert = F)]
dND4 = ndmall[,grep("t4", ndmtIDs, invert = F)]
nsdeltas = cbind.data.frame(dND1,dND2,dND3,dND4)
nsdeltimes= c(rep("t1",length(dND1[1,])),rep("t2",length(dND2[1,])),
	rep("t3",length(dND3[1,])),rep("t4",length(dND4[1,])))
nsgroups = factor(nsdeltimes, labels = c(0,1,2,3))

nsearlyd = cbind.data.frame(dND1,dND2,dND3,dND4,dD.5,dD1,dD4)
nsvsevsdgroups = c(rep("Non-progressor",length(nsdeltas[1,])),
	rep("early",length(dD.5[1,])), rep("early",length(dD1[1,])), 
	rep("diabetic",length(dD4[1,])))
groupednsvsevsd = factor(nsvsevsdgroups, labels = c(2,1,0))

nswdeltas = cbind.data.frame(nsdeltas, deltas)
nswdtimes = c(rep("Non-progressor",length(nsdeltas[1,])),
	rep("t0.5",length(dD.5[1,])),rep("t1",length(dD1[1,])),
	rep("t2",length(dD2[1,])),rep("t3",length(dD3[1,])),
	rep("t4",length(dD4[1,])))
nswdgroups = factor(nswdtimes, labels = c(0,1,2,3,4,5))

earlyvsd = cbind.data.frame(dD.5,dD1,dD4)
evsdtimes= c(rep("t0.5",length(dD.5[1,])),rep("t1",length(dD1[1,])),
	rep("t4",length(dD4[1,])))
evsdgroups = c(rep("early",length(dD.5[1,])),rep("early",length(dD1[1,])),
	rep("diabetic",length(dD4[1,])))
groupedevsd = factor(evsdgroups, labels = c(0,1))

evsl = cbind.data.frame(dD.5,dD1,dD3,dD4)
evsltimes = c(rep("t0.5",length(dD.5[1,])),rep("t1",length(dD1[1,])),
	rep("t3",length(dD3[1,])),rep("t4",length(dD4[1,])))
evslgroups = c(rep("early",length(dD.5[1,])),rep("early",length(dD1[1,])),
	rep("late",length(dD3[1,])),rep("late",length(dD4[1,])))
groupedevsl = factor(evslgroups, labels = c(0,1))
grouped75late = factor(deltimes, labels = c(0,1,2,2,2))

t.5t1 = cbind.data.frame(dD.5,dD1)
t.5t1times= c(rep("t0.5",length(dD.5[1,])),rep("t1",length(dD1[1,])))
t.5t1groups = c(rep("t0.5",length(dD.5[1,])),rep("t1",length(dD1[1,])))
groupedt.5t1 = factor(t.5t1groups, labels = c(0,1))

#Colors for grouped
coul <- colorRampPalette(brewer.pal(11, "PuOr"))(25)
ecolor <- c("dodgerblue","darkorchid","mediumvioletred","orangered","red4")
colSide <- deltimes
for (x in 1:length(deltimes)) {
  if (colSide[x] == "t0.5") {colSide[x] <- "dodgerblue"
  } else if (colSide[x] == "t1") {colSide[x] <- "darkorchid"
  } else if (colSide[x] == "t2") {colSide[x] <- "mediumvioletred"
  } else if (colSide[x] == "t3") {colSide[x] <- "orangered"
  } else if (colSide[x] == "t4") {colSide[x] <- "red4"
  }}
colind = colSide
names(colind) = deltimes

#Colors for grouped w/ NS
coul <- colorRampPalette(brewer.pal(11, "PuOr"))(25)
ecolor <- c("gold","dodgerblue","darkorchid","mediumvioletred","orangered","red4")
colSide <- nswdtimes
for (x in 1:length(nswdtimes)) {
  if (colSide[x] == "Non-progressor") {colSide[x] <- "gold"
  } else if (colSide[x] == "t0.5") {colSide[x] <- "dodgerblue"
  } else if (colSide[x] == "t1") {colSide[x] <- "darkorchid"
  } else if (colSide[x] == "t2") {colSide[x] <- "mediumvioletred"
  } else if (colSide[x] == "t3") {colSide[x] <- "orangered"
  } else if (colSide[x] == "t4") {colSide[x] <- "red4"
  }}
colind = colSide
names(colind) = nswdtimes

#Colors for NS grouped
coul <- colorRampPalette(brewer.pal(11, "PuOr"))(25)
ecolor <- c("gold","goldenrod4","darksalmon","darkorange")
colSide <- nsdeltimes
for (x in 1:length(nsdeltimes)) {
  if (colSide[x] == "t1") {colSide[x] <- "gold"
  } else if (colSide[x] == "t2") {colSide[x] <- "goldenrod4"
  } else if (colSide[x] == "t3") {colSide[x] <- "darksalmon"
  } else if (colSide[x] == "t4") {colSide[x] <- "darkorange"
  }}
colind = colSide
names(colind) = nsdeltimes

#Colors for t0.5 vs t1 grouped
coul <- colorRampPalette(brewer.pal(11, "PuOr"))(25)
ecolor <- c("dodgerblue","darkorchid")
colSide <- t.5t1groups
for (x in 1:length(t.5t1groups)) {
  if (colSide[x] == "t0.5") {colSide[x] <- "dodgerblue"
  } else if (colSide[x] == "t1") {colSide[x] <- "darkorchid"
  }}
colind = colSide
names(colind) = t.5t1groups

#Colors for early vs diabetic grouped
coul <- colorRampPalette(brewer.pal(11, "PuOr"))(25)
ecolor <- c("dodgerblue","darkorchid","red4")
colSide <- evsdtimes
for (x in 1:length(evsdtimes)) {
  if (colSide[x] == "t0.5") {colSide[x] <- "dodgerblue"
  } else if (colSide[x] == "t1") {colSide[x] <- "darkorchid"
  } else if (colSide[x] == "t4") {colSide[x] <- "red4"
  }}
colind = colSide
names(colind) = evsdtimes

#Colors for early (t0.5/t1) vs late (t3/t4) grouped
coul <- colorRampPalette(brewer.pal(11, "PuOr"))(25)
ecolor <- c("dodgerblue","darkorchid","orangered","red4")
colSide <- evsltimes
for (x in 1:length(evsltimes)) {
  if (colSide[x] == "t0.5") {colSide[x] <- "dodgerblue"
  } else if (colSide[x] == "t1") {colSide[x] <- "darkorchid"
  } else if (colSide[x] == "t3") {colSide[x] <- "orangered"
  } else if (colSide[x] == "t4") {colSide[x] <- "red4"
  }}
colind = colSide
names(colind) = evsltimes

groups = nsvsevsdgroups
ecolor <- c("red4","dodgerblue","gold")
colSide <- groups
for (x in 1:length(groups)) {
  if (groups[x] == "Non-progressor") {colSide[x] <- "gold" 
  } else if (groups[x] == "early") {colSide[x] <- "skyblue"  
  } else if (groups[x] == "diabetic") {colSide[x] <- "red4"  
  }}
colind = colSide
names(colind) = nsvsevsdgroups


#FlexiDEG
groups = groupedevsd
NODevsdsamples = colnames(earlyvsd)
metadata = cbind(NODevsdsamples, groups)
colnames(metadata) = c("Samples", "Group")
metadata <- as.data.frame(metadata)
NODevsdflexi <- flexiDEG.function4(earlyvsd, metadata, validation_option = 2)
EN13 <- na.omit(earlyvsd[unique(rownames(earlyvsd)[as_vector(NODevsdflexi[[1]])]), ]) 
EN14 <- na.omit(earlyvsd[unique(rownames(earlyvsd)[as_vector(NODevsdflexi[[2]])]), ]) 
EN15 <- na.omit(earlyvsd[unique(rownames(earlyvsd)[as_vector(NODevsdflexi[[3]])]), ]) 
heatmap.2(as.matrix(EN14), scale="row", col=coul, key= T, xlab="", ylab="", 
          margins=c(7,7), ColSideColors=colind, trace="none", key.title=NA, 
          key.ylab=NA, keysize=0.8, dendrogram="both")
ggbiplot(prcomp(t(EN14), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
         var.scale=1, circle=T) + 
  theme_classic() + scale_color_manual(name="Group", values=ecolor)


#Grouped EN, NOD timecourse
  xfactors <- cbind.data.frame(groupedfNODdel) 
samplesbc = colnames(deltas)
  rownames(xfactors) <- samplesbc
  dataset <- t(deltas) # Load the data. Should it be raw or normalized!?!?!?!
  dataset <- as.data.frame(cbind(xfactors, dataset)) # Add factor as the first column
  dataset <- dataset[complete.cases(dataset), ] # For cases without missed items
  set.seed(123) # Split the data into training and test set
  training.samples <- dataset$groupedfNODdel %>% createDataPartition(p = 1.0, list = F)
  train.data  <- dataset[training.samples, ]
  test.data <- dataset[-training.samples, ]
  x <- model.matrix(groupedfNODdel~., train.data)[,-1] # Predictor variables; need reduced gene set to run, runs w/ 9814
  y <- train.data$groupedfNODdel # Outcome variable
  genes_EN <- multinom5_EN( x, y, 100) # number following multinom is number of groups
#NODdelevsd = genes_EN
NODdelevsd = genes_EN
batchcorr_counts = deltas

#Grouped EN, NOD timecourse w/ NS
  xfactors <- cbind.data.frame(nswdgroups) 
samplesbc = colnames(nswdeltas)
  rownames(xfactors) <- samplesbc
  dataset <- t(nswdeltas) # Load the data. Should it be raw or normalized!?!?!?!
  dataset <- as.data.frame(cbind(xfactors, dataset)) # Add factor as the first column
  dataset <- dataset[complete.cases(dataset), ] # For cases without missed items
  set.seed(123) # Split the data into training and test set
  training.samples <- dataset$nswdgroups %>% createDataPartition(p = 1.0, list = F)
  train.data  <- dataset[training.samples, ]
  test.data <- dataset[-training.samples, ]
  x <- model.matrix(nswdgroups~., train.data)[,-1] # Predictor variables; need reduced gene set to run, runs w/ 9814
  y <- train.data$nswdgroups # Outcome variable
  genes_EN <- multinom6_EN( x, y, 100) # number following multinom is number of groups
NODdelwNS = genes_EN
batchcorr_counts = nswdeltas

#Grouped EN, NS timecourse 
  xfactors <- cbind.data.frame(nsgroups) 
samplesbc = colnames(nsdeltas)
  rownames(xfactors) <- samplesbc
  dataset <- t(nsdeltas) # Load the data. Should it be raw or normalized!?!?!?!
  dataset <- as.data.frame(cbind(xfactors, dataset)) # Add factor as the first column
  dataset <- dataset[complete.cases(dataset), ] # For cases without missed items
  set.seed(123) # Split the data into training and test set
  training.samples <- dataset$nsgroups %>% createDataPartition(p = 1.0, list = F)
  train.data  <- dataset[training.samples, ]
  test.data <- dataset[-training.samples, ]
  x <- model.matrix(nsgroups~., train.data)[,-1] # Predictor variables; need reduced gene set to run, runs w/ 9814
  y <- train.data$nsgroups # Outcome variable
  genes_EN <- multinom4_EN( x, y, 100) # number following multinom is number of groups
NSdel = genes_EN
batchcorr_counts = nsdeltas

#Grouped EN, NOD t0.5, t1, late (t2-4) timecourse 
  xfactors <- cbind.data.frame(grouped75late) 
samplesbc = colnames(deltas)
  rownames(xfactors) <- samplesbc
  dataset <- t(deltas) # Load the data. Should it be raw or normalized!?!?!?!
  dataset <- as.data.frame(cbind(xfactors, dataset)) # Add factor as the first column
  dataset <- dataset[complete.cases(dataset), ] # For cases without missed items
  set.seed(123) # Split the data into training and test set
  training.samples <- dataset$grouped75late %>% createDataPartition(p = 1.0, list = F)
  train.data  <- dataset[training.samples, ]
  test.data <- dataset[-training.samples, ]
  x <- model.matrix(grouped75late~., train.data)[,-1] # Predictor variables; need reduced gene set to run, runs w/ 9814
  y <- train.data$grouped75late # Outcome variable
  genes_EN <- multinom3_EN( x, y, 100) # number following multinom is number of groups
NODdel75late = genes_EN
batchcorr_counts = deltas

#Grouped EN, NOD t0.5, t1 
  xfactors <- cbind.data.frame(groupedt.5t1) 
samplesbc = colnames(t.5t1)
  rownames(xfactors) <- samplesbc
  dataset <- t(t.5t1) # Load the data. Should it be raw or normalized!?!?!?!
  dataset <- as.data.frame(cbind(xfactors, dataset)) # Add factor as the first column
  dataset <- dataset[complete.cases(dataset), ] # For cases without missed items
  set.seed(123) # Split the data into training and test set
  training.samples <- dataset$groupedt.5t1 %>% createDataPartition(p = 1.0, list = F)
  train.data  <- dataset[training.samples, ]
  test.data <- dataset[-training.samples, ]
  x <- model.matrix(groupedt.5t1~., train.data)[,-1] # Predictor variables; need reduced gene set to run, runs w/ 9814
  y <- train.data$groupedt.5t1 # Outcome variable
  genes_EN <- binom_EN( x, y, 100) # number following multinom is number of groups
NODt.5t1sig = genes_EN
batchcorr_counts = t.5t1


  g_EN1 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`1`])
  g_EN.95 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.95`])
  g_EN.9 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.9`])
  g_EN.85 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.85`])
  g_EN.8 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.8`])
  g_EN.75 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.75`])
  g_EN.7 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.7`])
  g_EN.65 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.65`])
  g_EN.6 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.6`])
  g_EN.55 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.55`])
  g_EN.5 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.5`])
  g_EN.45 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.45`])
  g_EN.4 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.4`])
  g_EN.35 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.35`])
  g_EN.3 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.3`])
  g_EN.25 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.25`])
  g_EN.2 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.2`])
  g_EN.15 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.15`])
  g_EN.1 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.1`])
  g_EN.05 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.0499999999999999`])
  g_EN0 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0`])
  n_EN1 <- na.omit(batchcorr_counts[g_EN1, ])
  n_EN.95 <- na.omit(batchcorr_counts[g_EN.95, ])
  n_EN.9 <- na.omit(batchcorr_counts[g_EN.9, ])
  n_EN.85 <- na.omit(batchcorr_counts[g_EN.85, ])
  n_EN.8 <- na.omit(batchcorr_counts[g_EN.8, ])
  n_EN.75 <- na.omit(batchcorr_counts[g_EN.75, ])
  n_EN.7 <- na.omit(batchcorr_counts[g_EN.7, ])
  n_EN.65 <- na.omit(batchcorr_counts[g_EN.65, ])
  n_EN.6 <- na.omit(batchcorr_counts[g_EN.6, ])
  n_EN.55 <- na.omit(batchcorr_counts[g_EN.55, ])
  n_EN.5 <- na.omit(batchcorr_counts[g_EN.5, ])
  n_EN.45 <- na.omit(batchcorr_counts[g_EN.45, ])
  n_EN.4 <- na.omit(batchcorr_counts[g_EN.4, ])
  n_EN.35 <- na.omit(batchcorr_counts[g_EN.35, ])
  n_EN.3 <- na.omit(batchcorr_counts[g_EN.3, ])
  n_EN.25 <- na.omit(batchcorr_counts[g_EN.25, ])
  n_EN.2 <- na.omit(batchcorr_counts[g_EN.2, ])
  n_EN.15 <- na.omit(batchcorr_counts[g_EN.15, ])
  n_EN.1 <- na.omit(batchcorr_counts[g_EN.1, ])
  n_EN.05 <- na.omit(batchcorr_counts[g_EN.05, ])
  n_EN0 <- na.omit(batchcorr_counts[g_EN0, ])
pdf(file= "NOD t.5, t1, bc2.pdf" )
  heatmap.2(as.matrix(n_EN1), scale = "row", col = coul, key = T, 
            xlab = "Sample", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN1), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=1 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.95), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.95), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.95 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.9), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.9), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.9 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.85), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.85), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.85 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.8), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.8), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.8 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.75), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.75), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.75 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.7), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.7), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.7 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.65), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.65), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.65 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.6), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.6), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.6 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.55), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.55), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.55 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.5), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.5), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.5 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.45), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.45), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.45 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.4), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.4), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.4 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.35), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.35), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.35 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.3), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.3), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.3 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.25), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.25), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.25 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.2), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.2), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.2 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.15), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.15), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.15 Geneset") + 
    geom_point(size=3, color=colind) + 
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.1), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.1), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.1 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.05), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.05), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.05 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  g_Groupedfull <- g_EN.05
  g_Groupedlasso <- g_EN1
dev.off() #stop adding to PDF file

#Grouped EN, NOD early vs late with NS
  xfactors <- cbind.data.frame(groupednsvsevsd) 
samplesbc = colnames(nsearlyd)
  rownames(xfactors) <- samplesbc
  dataset <- t(nsearlyd) # Load the data. Should it be raw or normalized!?!?!?!
  dataset <- as.data.frame(cbind(xfactors, dataset)) # Add factor as the first column
  dataset <- dataset[complete.cases(dataset), ] # For cases without missed items
  set.seed(123) # Split the data into training and test set
  training.samples <- dataset$groupednsvsevsd %>% createDataPartition(p = 1.0, list = F)
  train.data  <- dataset[training.samples, ]
  test.data <- dataset[-training.samples, ]
  x <- model.matrix(groupednsvsevsd~., train.data)[,-1] # Predictor variables; need reduced gene set to run, runs w/ 9814
  y <- train.data$groupednsvsevsd # Outcome variable
  genes_EN <- multinom3_EN( x, y, 100) # number following multinom is number of groups
delNSearlylate = genes_EN
batchcorr_counts = nsearlyd
  g_EN1 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`1`])
  g_EN.95 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.95`])
  g_EN.9 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.9`])
  g_EN.85 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.85`])
  g_EN.8 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.8`])
  g_EN.75 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.75`])
  g_EN.7 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.7`])
  g_EN.65 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.65`])
  g_EN.6 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.6`])
  g_EN.55 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.55`])
  g_EN.5 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.5`])
  g_EN.45 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.45`])
  g_EN.4 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.4`])
  g_EN.35 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.35`])
  g_EN.3 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.3`])
  g_EN.25 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.25`])
  g_EN.2 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.2`])
  g_EN.15 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.15`])
  g_EN.1 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.1`])
  g_EN.05 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.0499999999999999`])
  g_EN0 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0`])
  n_EN1 <- na.omit(batchcorr_counts[g_EN1, ])
  n_EN.95 <- na.omit(batchcorr_counts[g_EN.95, ])
  n_EN.9 <- na.omit(batchcorr_counts[g_EN.9, ])
  n_EN.85 <- na.omit(batchcorr_counts[g_EN.85, ])
  n_EN.8 <- na.omit(batchcorr_counts[g_EN.8, ])
  n_EN.75 <- na.omit(batchcorr_counts[g_EN.75, ])
  n_EN.7 <- na.omit(batchcorr_counts[g_EN.7, ])
  n_EN.65 <- na.omit(batchcorr_counts[g_EN.65, ])
  n_EN.6 <- na.omit(batchcorr_counts[g_EN.6, ])
  n_EN.55 <- na.omit(batchcorr_counts[g_EN.55, ])
  n_EN.5 <- na.omit(batchcorr_counts[g_EN.5, ])
  n_EN.45 <- na.omit(batchcorr_counts[g_EN.45, ])
  n_EN.4 <- na.omit(batchcorr_counts[g_EN.4, ])
  n_EN.35 <- na.omit(batchcorr_counts[g_EN.35, ])
  n_EN.3 <- na.omit(batchcorr_counts[g_EN.3, ])
  n_EN.25 <- na.omit(batchcorr_counts[g_EN.25, ])
  n_EN.2 <- na.omit(batchcorr_counts[g_EN.2, ])
  n_EN.15 <- na.omit(batchcorr_counts[g_EN.15, ])
  n_EN.1 <- na.omit(batchcorr_counts[g_EN.1, ])
  n_EN.05 <- na.omit(batchcorr_counts[g_EN.05, ])
  n_EN0 <- na.omit(batchcorr_counts[g_EN0, ])
pdf(file= "Grouped NS and NOD early vs diabetic delta groups.pdf" )
  heatmap.2(as.matrix(n_EN1), scale = "row", col = coul, key = T, 
            xlab = "Sample", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN1), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=1 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.95), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.95), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.95 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.9), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.9), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.9 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.85), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.85), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.85 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.8), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.8), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.8 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.75), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.75), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.75 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.7), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.7), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.7 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.65), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.65), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.65 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.6), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.6), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.6 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.55), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.55), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.55 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.5), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.5), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.5 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.45), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.45), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.45 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.4), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.4), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.4 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.35), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.35), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.35 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.3), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.3), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.3 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.25), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.25), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.25 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.2), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.2), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.2 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.15), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.15), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.15 Geneset") + 
    geom_point(size=3, color=colind) + 
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.1), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.1), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.1 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  heatmap.2(as.matrix(n_EN.05), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
  ggbiplot(prcomp(t(n_EN.05), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.05 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)
  g_Groupedfull <- g_EN.05
  g_Groupedlasso <- g_EN1
dev.off() #stop adding to PDF file

