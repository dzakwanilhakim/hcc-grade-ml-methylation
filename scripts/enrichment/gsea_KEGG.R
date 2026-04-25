# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
gsea_kegg <- function(df_sigs){
  original_gene_list <- df_sigs$logFC
  names(original_gene_list) <- df_sigs$gene
  gene_list = sort(original_gene_list, decreasing = TRUE)
  organism = "org.Hs.eg.db"
  

  ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
  dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
  df2 = df_sigs[df_sigs$gene %in% dedup_ids$SYMBOL,]
  
  # Create a new column in df2 with the corresponding ENTREZ IDs
  df2$Y = dedup_ids$ENTREZID
  
  # Create a vector of the gene unuiverse
  kegg_gene_list <- df2$logFC
  
  # Name vector with ENTREZ ids
  names(kegg_gene_list) <- df2$Y
  
  # omit any NA values 
  kegg_gene_list<-na.omit(kegg_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
  
  kegg_organism = "hsa"
  kk2 <- gseKEGG(geneList     = kegg_gene_list,
                 organism     = kegg_organism,
                 nPerm        = 10000000,
                 minGSSize    = 5,
                 maxGSSize    = 800,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = 'none',
                 keyType       = "ncbi-geneid")
  
  return(kk2)
  

}

gokegg_NG1 <- gsea_kegg(sigs_NG1)
dotplot(gokegg_NG1, showCategory = 10, title = "G1" , split=".sign") + facet_grid(.~.sign)

gokegg_NG2 <- gsea_kegg(sigs_NG2)
dotplot(gokegg_NG2, showCategory = 10, title = "G2" , split=".sign") + facet_grid(.~.sign)

gokegg_NG3 <- gsea_kegg(sigs_NG3)
dotplot(gokegg_NG3, showCategory = 10, title = "G3" , split=".sign") + facet_grid(.~.sign)

gokegg_G2G3 <- gsea_kegg(sigs_G2G3)
dotplot(, showCategory = 10, title = "G3, baseline G2" , split=".sign") + facet_grid(.~.sign)


original_gene_list <- sigs_NG2$logFC
names(original_gene_list) <- sigs_NG2$gene
gene_list = sort(original_gene_list, decreasing = TRUE)
organism = "org.Hs.eg.db"

gsea_kegg <- f


ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]


# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = sigs_NG1[sigs_NG1$gene %in% dedup_ids$SYMBOL,]

library(dplyr)
anyDuplicated(dedup_ids$SYMBOL) > 0
anyDuplicated(dedup_ids$ENTREZID) > 0
dedup_ids <- dedup_ids %>% distinct(SYMBOL, .keep_all = TRUE)
dedup_ids <- dedup_ids %>% distinct(ENTREZID, .keep_all = TRUE)
# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$logFC

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 1000000,
               minGSSize    = 5,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = 'none',
               keyType       = "ncbi-geneid")

dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
