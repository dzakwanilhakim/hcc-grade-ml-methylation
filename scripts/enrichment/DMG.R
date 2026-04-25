library(DESeq2)
library(tidyverse)
library(limma)
library(lumi)
library(org.Hs.eg.db)
library(methylGSA)
library(clusterProfiler)
filter <- function(ph, listpair, beta) {
  # Step 1: Filter ph to create ph_filter
  ph_filter <- ph[ph$label %in% listpair, ]
  
  # Step 2: Map ph_filter$Sample_Name to the rownames of beta and filter out unmapped rows
  sample_names <- ph_filter$Sample_Name
  beta_filter <- beta[colnames(beta) %in% sample_names, ]
  
  # Ensure the row names of beta are identical to the Sample_Name in ph_filter
  identical_check_before <- identical(colnames(beta_filter), sample_names)
  print(paste("Row names identical before transposition:", identical_check_before))
  
  # Transpose the filtered beta data frame
  #beta_transpose <- t(beta_filter)
  
  # Set the new column names to be the old row names
  #colnames(beta_transpose) <- sample_names
  
  # Check if the new column names of the transposed beta match the Sample_Name in ph_filter
  #identical_check_after <- identical(colnames(beta_transpose), sample_names)
  #print(paste("Column names identical after transposition:", identical_check_after))
  
  return(list(ph_filter = ph_filter, beta_filter = beta_transpose))
}
#================================================================================================
NG1 <- filter(ph_label, c(0,1), beta_avg_p56_01)
NG2 <- filter(ph_label, c(0,2), beta_avg_p56_01)
NG3 <- filter(ph_label, c(0,3), beta_avg_p56_01)
G2G3 <- filter(ph_label, c(2,3), beta_avg_p56_01)
#============================================================================
phNG1 <- NG1$ph_filter
beta_NG1 <- NG1$beta_filter
gene_NG1 <- rownames(beta_NG1)
df_gene_NG1 <- as.data.frame(gene_NG1)
getOption(timeout)
options(timeout=1000)

cols <- "SYMBOL"
ensids <- gene_list$gene_list

select(org.Hs.eg.db, keys='ABCA4', columns='ENSEMBL', keytype="SYMBOL")
gene_NG1 <- rownames(beta_NG1)

ensids <- gene_list$gene_list



write.table(df_gene_NG1, "gene_list.txt", sep = " ", row.names = FALSE)



M_NG1 <- beta2m(beta_NG1)
group <- factor(phNG1$label,levels=c(0,1))
id <- factor(phNG1$Sample_Name)
design <- model.matrix(~id + group)
design

fit.reduced <- lmFit(M_NG1,design)
fit.reduced <- eBayes(fit.reduced)

rownames(NG1$ph_filter) <- NG1$ph_filter$Sample_Name

beta_NG1 <-NG1$beta_filter 
ph_NG1 <- NG1$ph_filter

dds <- DESeqDataSetFromMatrix(counts = beta_NG1, 
                       colData = ph_NG1,
                       design = ~label )














#================================================================================================
DMP_NG1 <- champ.DMP(beta = NG1$beta_filter, pheno = NG1$ph_filter$label)
DMP_NG2 <- champ.DMP(beta = NG2$beta_filter, pheno = NG2$ph_filter$label)
DMP_NG3 <- champ.DMP(beta = NG3$beta_filter, pheno = NG3$ph_filter$label)
DMP_G2G3 <- champ.DMP(beta = G2G3$beta_filter, pheno = G2G3$ph_filter$label, adjPVal = 0.01)
#================================================================================================
DMR_NG1 <- champ.DMR(beta=NG1$beta_filter,pheno=NG1$ph_filter$label,method="DMRcate")
DMR_NG2 <- champ.DMR(beta=NG2$beta_filter,pheno=NG2$ph_filter$label,method="Bumphunter")
DMR_NG3 <- champ.DMR(beta=NG3$beta_filter %>% as.data.frame(),pheno=NG3$ph_filter$label,method="Bumphunter")
DMR_G2G3 <- champ.DMR(beta=G2G3$beta_filter %>% as.data.frame(),pheno=G2G3$ph_filter$label,method="Bumphunter")
#================================================================================================

GSEA_NG2 <- champ.ebGSEA(beta=beta_NG2, pheno=phNG2, minN=5, adjPval=0.05, arraytype="450K", cores=2)
DMRG2 <-  DMRcate_manual_run(beta_NG2, phNG2$label, "450K")

library(dplyr)
phNG2 <- NG2$ph_filter
phNG2$label <- case_when(
  phNG2$label == 0 ~ 'NT',
  phNG2$label == 2 ~ 'G2'
)


DMR_NG2 <- champ.DMR(beta=beta_NG2, minProbes=1, pheno=phNG2$label,method="Bumphunter")
beta_NG2 <- NG2$beta_filter


length(NG1$ph_filter$label)
length(NG1$beta_filter$cg00034556)


filtered_result <- filter(ph_label, c(0, 1), X_df)
filtered_ph <- filtered_result$ph_filter
filtered_beta <- filtered_result$beta_filter
identical_check <- filtered_result$identical_chec

champ.DMR(beta=myNorm,
          pheno=myLoad$pd$Sample_Group,
          compare.group=NULL,
          arraytype="450K",
          method = "Bumphunter",
          minProbes=1,
          adjPvalDmr=0.05,
          cores=3,


DMP <- champ.DMP(beta = myNorm,
          pheno = myLoad$pd$Sample_Group,
          compare.group = NULL,
          adjPVal = 0.05,
          adjust.method = "BH",
          arraytype = "450K")
champ.DMP(beta = myNorm,
          pheno = myLoad$pd$Sample_Group,
          compare.group = NULL,
          adjPVal = 0.05,
          adjust.method = "BH",
          arraytype = "450K")


myLoad <- champ.load(directory=system.file("extdata",package="ChAMPdata"))	
myNorm <- champ.norm()
myDMR <- champ.DMR(method = 'DMRcate')

champ.GSEA(beta=myNorm,
           DMP=myDMP[[1]],
           DMR=myDMR,
           CpGlist=NULL,
           Genelist=NULL,
           pheno=myLoad$pd$Sample_Group,
           method="fisher",
           arraytype="450K",
           Rplot=TRUE,
           adjPval=0.05,
           cores=1)
