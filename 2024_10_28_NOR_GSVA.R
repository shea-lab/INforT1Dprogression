#GSEA and/or GSVA
# 28 October 2024

  # * 10.1 DE Analysis ----
  factor_DEseq <- c(strain)
  coldata <- data.frame(cbind(factor_DEseq))
  row.names(coldata) <- colnames(r_counts)
  dds <- DESeqDataSetFromMatrix(countData=round(r_counts), colData=coldata, 
                                design = ~ factor_DEseq)
  paste(nrow(dds), " genes input into DESeq2 from these filters", sep="")
  dds <- DESeq(dds)
  res <- results(dds) # Table of log2 fold changes, p-values, & p-adj values
#  dds <- estimateSizeFactors(dds)
# dds <- estimateDispersionsGeneEst(dds)
#  dispersions(dds) <- mcols(dds)$dispGeneEst



  # * 10.2 Organize Results ----
  resultsordered <- res[order(res$padj),]   # Order results by p-adj
  sum(res$padj < 0.1, na.rm=T)   # How many genes have p-adj < 0.1
  resSig <- subset(resultsordered, padj < 0.1) # Subset by p-adj < 0.1, 8398 genes
#  head(resSig[order(resSig$log2FoldChange),]) # Sig genes w/ strongest down-regulation
#  tail(resSig[order(resSig$log2FoldChange),]) # Sig genes w/ strongest up-regulation
  
  dds_genes <- rownames(resSig) # Differentialy-expressed genes
  r_DEG <- r_counts[dds_genes,] 

 # 10.3 GSEA ----
library(GOSemSim)
library(AnnotationHub)
library(org.Mm.eg.db)

#NOR vs NOD
ENgenes = EN3
egoNORNOD <- enrichGO(gene         = as.character(ENgenes),
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.1,
                readable  = TRUE)
goplot(egoNORNOD)
plot(barplot(egoNORNOD, showCategory = 10))



#NOR timecourse
genes_EN = NORtime
  g_EN.9tc <- unique(rownames(batchcorr_counts)[genes_EN$ENgenes$`0.9`])

egoNORtc <- enrichGO(gene         = g_EN.9tc,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                readable  = TRUE)
goplot(egoNORtc)
plot(barplot(egoNORtc, showCategory = 10))

#GSVA (from Laila)
library(GSVA)
library(org.Mm.eg.db)
library(qusage)
library(magrittr)
#library(DESeq2)

dds = DESeqDataSetFromMatrix(
  countData = count.data_filt,
  colData = sample_info,
  design = ~1
)

dds_norm = vst(dds)

# Retrieve the normalized data from the `DESeqDataSet`
#vst_df <- assay(dds_norm) %>%
vst_df <- NOR[g_EN.9tc,] %>%
  as.data.frame() %>% # Make into a data frame
  tibble::rownames_to_column("ensembl_id") # Make Gene IDs into their own column


hallmark_gene_sets = msigdbr::msigdbr(
  species = "Mus musculus",
  category = "H"
)

head(hallmark_gene_sets)

hallmarks_list <- split(
  hallmark_gene_sets$entrez_gene, # The genes we want split into pathways
  hallmark_gene_sets$gs_name # The pathways made as the higher levels of the list
)

head(hallmarks_list, n = 2)

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
  # Now let's join the rest of the expression data
  dplyr::inner_join(vst_df, by = c("Symbol" = "ensembl_id"))


#write.csv(mapped_df2,paste(path,"/EAE_BI_cohort1genesigngenes.csv", sep = ""), row.names = T)

head(mapped_df)
sum(duplicated(mapped_df$entrez_id))

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

head(gsva_results[, 1:10])

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


png("gsva_heatmap.png", units = "in", width = 10, height = 8, res = 400)
pathway_heatmap
dev.off()
 

dds_g = DESeqDataSetFromMatrix(countData = gsva_results,
                               colData = sample_info,
                               design = ~Condition_t2)




