# Jessica King
# 18 July 2024
# Batch correction and signature analysis of original and validation NOD cohorts


library("pacman")
  pacman::p_load(tidyverse, plyr, magrittr, stats, dplyr, limma, RColorBrewer, gplots, glmnet, 
                 biomaRt, colorspace, ggplot2, fmsb, car, mixOmics, DESeq2, apeglm, rgl, qpcR,
                 boot, caret, ggvenn, grid, devtools, reshape2, gridExtra, factoextra, edgeR, 
                 cowplot, pheatmap, coefplot, randomForest, ROCR, genefilter, Hmisc, rdist, 
                 factoextra, ggforce, NormqPCR, ggpubr, matrixStats, GSEAmining, ggrepel)

#Load data
  #Original cohort
  og_file = read.csv(".\\Data\\og_gene_expected_count.annot.csv", header=T)
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
  vc_file = read.csv(".\\Data\\vc_gene_expected_count.annot.csv", header=T)
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
    row.names(coldata) <- 
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
filt_low <- 2 # Filters are per sample per gene
filt_high <- 16
keep1 <- rowSums(ps_counts) > (filt_low*sample_num) # Low count filter
ps_filt <- ps_counts[keep1,] 
dim(ps_filt) 
keep2 <- rowSums(ps_filt) < (filt_high*sample_num) # High count filter
ps_filt <- ps_filt[keep2,] 
dim(ps_filt) 

bc_filt <- bc_counts[rownames(ps_filt),]
 counts_filt = bc_filt
 # Filter by variance across all samples to exclude genes that don't change much
    mn_cut <- 0.2 # Mean cutoff
    var_cut <- 0.01 # Variance cutoff
    df <- cbind.data.frame("x" = rowMeans(counts_filt), # Calculate Mean & Var/Mean
                           "y" = rowVars(counts_filt) / rowMeans(counts_filt)) 
    g_var_all <- rownames(subset(df, x > mn_cut & y > var_cut)) # Genes to keep

    # Plot variance and mean Figure 1E
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

#Figure 2B, ROC for top 3 alphas
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

#Figure 2A
heatmap.2(as.matrix(EN3), scale="row", col=coul, key= T, xlab="", ylab="", 
          margins=c(7,7), ColSideColors=colind, trace="none", key.title=NA, 
          key.ylab=NA, keysize=0.8, dendrogram="both")
#Figure 2B
ggbiplot(prcomp(t(EN3), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
         var.scale=1, circle=T) + 
  theme_classic() + scale_color_manual(name="Group", values=ecolor) 
dev.off()


#NOR timecourse
NOR = bc_filt[ ,strain_subsets=="NOR"]
NS = bc_filt[ ,strain_subsets=="NOD non-progressors"]

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
  g_EN.9 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.9`])
  n_EN.9 <- na.omit(batchcorr_counts[g_EN.9, ])

#Figure 2D
  heatmap.2(as.matrix(n_EN.9), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "colu  heatmap.2(as.matrix(n_EN.85), scale = "row", col = coul, key = T, 
            xlab = "Scaffold", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
#Figure 2E
  ggbiplot(prcomp(t(n_EN.9), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.9 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)

#GSVA
library(GSVA)
library(org.Mm.eg.db)
library(qusage)
library(magrittr)

genes_EN = NORtime
  g_EN.9tc <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.9`])
vst_df <- NOR[g_EN.9tc,] %>%
  as.data.frame() %>% # Make into a data frame
  tibble::rownames_to_column("ensembl_id") # Make Gene IDs into their own column

cgp_gene_sets = msigdbr::msigdbr(
  species = "mouse", category = "C2", 
  subcategory = "CGP"
)
cgp_list <- split(
  cgp_gene_sets$entrez_gene, # The genes we want split into pathways
  cgp_gene_sets$gs_name # The pathways made as the higher levels of the list
)

keytypes(org.Mm.eg.db)

mapped_df = data.frame(
  "entrez_id" = mapIds(
    org.Mm.eg.db,
    keys = vst_df$ensembl_id,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )
) %>%
  # If an Ensembl gene identifier doesn't map to a Entrez gene identifier,
  # drop that from the data frame
  dplyr::filter(!is.na(entrez_id)) %>%
  # Make an `Ensembl` column to store the row names
  tibble::rownames_to_column("Symbol") %>%
  # Join the rest of the expression data
  dplyr::inner_join(vst_df, by = c("Symbol" = "ensembl_id"))

filtered_mapped_matrix <- mapped_df %>%
  # GSVA can't the Ensembl IDs so we should drop this column as well as the means
  dplyr::select(-Symbol) %>%
  # We need to store our gene identifiers as row names
  tibble::column_to_rownames("entrez_id") %>%
  # Now we can convert our object into a matrix
  as.matrix()

gsva_results <- gsva(
  filtered_mapped_matrix,
  cgp_list,
  method = "gsva",
  # Appropriate for our vst transformed data
  kcdf = "Gaussian",
  # Minimum gene set size
  min.sz = 5,
  # Maximum gene set size
  max.sz = 500,
  # Compute Gaussian-distributed scores
  mx.diff = TRUE,
  # Don't print out the progress bar
  verbose = T
)

gsva_results %>%
  as.data.frame() %>%
  tibble::rownames_to_column("pathway") %>%
  readr::write_tsv(
    "NORtc_gsva_results_ENfilt.tsv"
  )

(clin_info_NOR$time)

annot_df <- clin_info_NOR %>%
  # We need the sample IDs and the main column that contains the metadata info
  dplyr::select(
    time
  ) 

pathway_heatmap = pheatmap::pheatmap(gsva_results,
    annotation_col = annot_df,
    show_colnames = F,
    fontsize = 8
    )

heatmap(DEpwys_es, ColSideColors=sample.color.map, xlab="samples",
        ylab="Pathways", margins=c(2, 20),
        labRow=substr(gsub("_", " ", gsub("^KEGG_|^REACTOME_|^BIOCARTA_", "",
                                          rownames(DEpwys_es))), 1, 35),
        labCol="", scale="row", Colv=as.dendrogram(sampleClustering),
        Rowv=as.dendrogram(geneSetClustering))

#Figure 2F
png("gsva_heatmap.png", units = "in", width = 10, height = 8, res = 400)
pathway_heatmap
dev.off()


#Filter NOR time course genes out of NOD data frame
NORtc_genes = g_EN.9 #39 genes 
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

#Figure 3D, first ROC for 5 iterations
NODwk6flexi <- flexiDEG.function4(wk6, metadata, validation_option = 2) 

#Simple EN, NOD vs NS week 6 
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
batchcorr_counts = wk6
  g_EN.55 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.55`])
  n_EN.55 <- na.omit(batchcorr_counts[g_EN.55, ])

#Figure 3A
  heatmap.2(as.matrix(n_EN.55), scale = "row", col = coul, key = T, 
            xlab = "", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")
#Figure 3B
  ggbiplot(prcomp(t(n_EN.55), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic(base_size=15) + ggtitle(" EN alpha=.55 Geneset") + 
    theme(legend.position= c(0.75, 0.85)) + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)

#NOD week 6 analysis w/ all NS/NP time points
NODp = NOD_filt[,strainNOD=="NOD progressors"]
  clin_info_NODp <- clin_info[ grep("NOD progressors", clin_info[,"strain_subsets"], invert = F) , ] #Remove NOD
times = clin_info_NODp$time

NODnp = NOD_filt[,strainNOD=="NOD non-progressors"]
clin_info_NODnp <- clin_info[ grep("NOD non-progressors", clin_info[,"strain_subsets"], invert = F) , ] #Remove NOD progressors

wk6p = NODp[ , times=="t0"]
wk6pnpall = cbind.data.frame(NODnp, wk6p)

wk6allgroups = c(rep("Non-progressor",39),(rep("Progressor wk 6",10)))
groupswk6all = factor(wk6allgroups, labels = c(0,1))

#Colors for all samples
coul <- colorRampPalette(brewer.pal(11, "PuOr"))(25)
ecolor <- c("gold", "skyblue")
colSide <- wk6allgroups
for (x in 1:length(wk6allgroups)) {
  if (colSide[x] == "Non-progressor") {colSide[x] <- "gold" 
  } else if (colSide[x] == "Progressor wk 6") {colSide[x] <- "skyblue"
  }}
colind = colSide
names(colind) = wk6allgroups

#FlexiDEG
groups = groupswk6all #1: NOD, 0:NS
wk6samples = colnames(wk6pnpall)
metadata = cbind(wk6samples, groups)
colnames(metadata) = c("Samples", "Group")
metadata <- as.data.frame(metadata)

#Figure 4B, ROC of top 3 alphas for 50 iterations
NODwk6flexi <- flexiDEG.function4(wk6pnpall, metadata, validation_option = 2)

#Figure 4A
  heatmap.2(as.matrix(EN8.3), scale = "row", col = coul, key = T, 
            xlab = "", ylab="", labCol=FALSE, # density.info="none",
            margins = c(7, 9), ColSideColors = colSide, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")

#Figure 4C
  ggbiplot(prcomp(t(EN8.3), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic(base_size=15) +  
    theme(legend.position= c(0.75, 0.09)) + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)


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

#Colors for all times
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

#Grouped EN, NOD timecourse
  xfactors <- cbind.data.frame(groupedfNODdel) 
samplesbc = colnames(deltas)
  rownames(xfactors) <- samplesbc
  dataset <- t(deltas) 
  dataset <- as.data.frame(cbind(xfactors, dataset)) # Add factor as the first column
  dataset <- dataset[complete.cases(dataset), ] # For cases without missed items
  set.seed(123) # Split the data into training and test set
  training.samples <- dataset$groupedfNODdel %>% createDataPartition(p = 1.0, list = F)
  train.data  <- dataset[training.samples, ]
  test.data <- dataset[-training.samples, ]
  x <- model.matrix(groupedfNODdel~., train.data)[,-1] # Predictor variables
  y <- train.data$groupedfNODdel # Outcome variable
  genes_EN <- multinom5_EN( x, y, 100) # number following multinom is number of groups
#NODdelevsd = genes_EN
NODdelevsd = genes_EN
batchcorr_counts = deltas
  g_EN.9 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.9`])
  n_EN.9 <- na.omit(batchcorr_counts[g_EN.9, ])

#Figure 5A
  ggbiplot(prcomp(t(n_EN.9), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic() + ggtitle(" EN alpha=.9 Geneset") + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)


#-7 to -5 signature
t.5t1 = cbind.data.frame(dD.5,dD1)
t.5t1times= c(rep("t0.5",length(dD.5[1,])),rep("t1",length(dD1[1,])))
t.5t1groups = c(rep("t0.5",length(dD.5[1,])),rep("t1",length(dD1[1,])))
groupedt.5t1 = factor(t.5t1groups, labels = c(0,1))

#Colors for t0.5 vs t1 
coul <- colorRampPalette(brewer.pal(11, "PuOr"))(25)
ecolor <- c("dodgerblue","darkorchid")
colSide <- t.5t1groups
for (x in 1:length(t.5t1groups)) {
  if (colSide[x] == "t0.5") {colSide[x] <- "dodgerblue"
  } else if (colSide[x] == "t1") {colSide[x] <- "darkorchid"
  }}
colind = colSide
names(colind) = t.5t1groups

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
  g_EN.85 <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.85`])
  n_EN.85 <- na.omit(batchcorr_counts[g_EN.85, ])

#Figure 5B
  heatmap.2(as.matrix(n_EN.85), scale = "row", col = coul, key = T, 
            xlab = "", ylab="", # density.info="none",
            margins = c(7, 7), ColSideColors = colind, trace="none", # Colv = "Na",
            key.title=NA, key.ylab=NA, keysize = 0.8, dendrogram = "column")

#Figure 5C
  ggbiplot(prcomp(t(n_EN.85), scale.=T), ellipse=T, groups=names(colind), var.axes=F, 
           var.scale=1, circle=T) + theme_classic(base_size=15) + ggtitle(" EN alpha=.85 Geneset") +  
    theme(legend.position= c(0.75, 0.8)) + 
    geom_point(size=3, color=colind) +
    scale_color_manual(name="Group", values = ecolor)

 
