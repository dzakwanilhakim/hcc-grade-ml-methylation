BiocManager::install("pathview")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
library(org.Dm.eg.db)

NG1 <- filter(ph_label, c(0,1), beta_avg_p56_01)
NG2 <- filter(ph_label, c(0,2), beta_avg_p56_01)
NG3 <- filter(ph_label, c(0,3), beta_avg_p56_01)
G2G3 <- filter(ph_label, c(2,3), beta_avg_p56_01)
#======================================================

gene2cpg <- function(bins, beta){
  first_cg_mapping <- sapply(bins, `[`, 1)
  # Replace row names in beta_NG1
  rownames(beta) <- first_cg_mapping[rownames(beta)]
  return(beta)
}

beta_NG1 <- NG1$beta_filter
ph_NG1 <- NG1$ph_filter

beta_NG2 <- NG2$beta_filter
ph_NG2 <- NG2$ph_filter

beta_NG3 <- NG3$beta_filter
ph_NG3 <- NG3$ph_filter

beta_G2G3 <- G2G3$beta_filter
ph_G2G3 <- G2G3$ph_filter


beta_NG1_cpg <- gene2cpg(bin_p_56_01, beta_NG1)
beta_NG2_cpg <- gene2cpg(bin_p_56_01, beta_NG2)
beta_NG3_cpg <- gene2cpg(bin_p_56_01, beta_NG3)
beta_G2G3_cpg <- gene2cpg(bin_p_56_01, beta_G2G3)

DMP_NG1 <- champ.DMP(beta = beta_NG1_cpg, pheno = ph_NG1$label, adjPVal = 0.01)
DMP_NG2 <- champ.DMP(beta = beta_NG2_cpg, pheno = ph_NG2$label, adjPVal = 0.01)
DMP_NG3 <- champ.DMP(beta = beta_NG3_cpg, pheno = ph_NG3$label, adjPVal = 0.01)
DMP_G2G3 <- champ.DMP(beta = beta_G2G3_cpg, pheno = ph_G2G3$label, adjPVal = 0.01)

sigs_NG1 <- DMP_NG1[[1]]
sigs_NG2 <- DMP_NG2[[1]]
sigs_NG3 <- DMP_NG3[[1]]
sigs_G2G3 <- DMP_G2G3[[1]]

prior_genelist <- colnames(beta_avg_p56_01)
length(prior_genelist)
all(sigs_NG1$gene %in% prior_genelist)
all(sigs_NG2$gene %in% prior_genelist)
all(sigs_NG3$gene %in% prior_genelist)
all(sigs_G2G3$gene %in% prior_genelist)



go_enrichment <- function(df_sigs){
  original_gene_list <- df_sigs$logFC
  names(original_gene_list) <- df_sigs$gene
  gene_list = sort(original_gene_list, decreasing = TRUE)
  
  organism = "org.Hs.eg.db"
  
  gse_fgsea <- gseGO(geneList=gene_list, 
                     ont ="ALL", 
                     keyType = "SYMBOL", 
                     nPermSimple = 1000000,
                     minGSSize = 5, 
                     maxGSSize = 800, 
                     pvalueCutoff = 0.05, 
                     verbose = TRUE, 
                     OrgDb = organism, 
                     pAdjustMethod = "BH",
                     by = 'fgsea')
  return(gse_fgsea)
}

gsea_NG1 <- go_enrichment(sigs_NG1)
gsea_NG2 <- go_enrichment(sigs_NG2)
gsea_NG3 <- go_enrichment(sigs_NG3)
gsea_G2G3 <- go_enrichment(sigs_G2G3)

dotplot(gsea_NG3, showCategory=10, split=".sign",
        font.size=8, title='G3',) + facet_grid(.~.sign)

dotplot(gsea_NG1_none, showCategory=10, split=".sign",
        font.size=8, title='G1',) + facet_grid(.~.sign)

dotplot(gsea_NG2_none, showCategory=10, split=".sign",
        font.size=8, title='G2',) + facet_grid(.~.sign)

dotplot(gsea_NG3_none, showCategory=10, split=".sign",
        font.size=8, title='G3',) + facet_grid(.~.sign)
emapplot(gsea_NG1)    

original_gene_list <- sigs_NG1$logFC
names(original_gene_list) <- sigs_NG1$gene
gene_list = sort(original_gene_list, decreasing = TRUE)

organism = "org.Hs.eg.db"
#======================================================



sigs_G2G3 <- DMP_G2G3[[1]]
G2G3_genelist <- sigs_G2G3$gene


write.table(G2G3_genelist, "G2G3_genelist.txt", sep = " ", row.names = FALSE)

#======================================================
original_gene_list <- sigs_NG1$logFC
names(original_gene_list) <- sigs_NG1$gene
gene_list = sort(original_gene_list, decreasing = TRUE)

organism = "org.Hs.eg.db"

gse_fgsea <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPermSimple = 1000000,
             minGSSize = 5, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "BH",
             by = 'fgsea')

require(DOSE)
require(fgsea)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
dotplot(gse_fgsea, showCategory=5, split=".sign",
        font.size=8, title='G1',) + facet_grid(.~.sign)
emapplot(gse)
# Return the modified DataFrame
beta_NG1


# we want the log2 fold change 
original_gene_list <- test$log2FoldChange

# name the vector
names(original_gene_list) <- rownames(test)

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

organism = "org.Dm.eg.db"

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
