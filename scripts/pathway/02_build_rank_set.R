library(readxl)
library(dplyr)
library(scales) 
library(org.Hs.eg.db)
library(clusterProfiler)

G1_hipo <- read.delim("D:/TA/KEGG_PATHWAY/rank_G1_hipo.rnk", header=FALSE)
G2_hipo <- read.delim("D:/TA/KEGG_PATHWAY/rank_G2_hipo.rnk", header=FALSE)
G3_hipo <- read.delim("D:/TA/KEGG_PATHWAY/rank_G3_hipo.rnk", header=FALSE)
G1_hiper <- read.delim("D:/TA/KEGG_PATHWAY/rank_G1_hiper.rnk", header=FALSE)
G2_hiper <- read.delim("D:/TA/KEGG_PATHWAY/rank_G2_hiper.rnk", header=FALSE)
G3_hiper <- read.delim("D:/TA/KEGG_PATHWAY/rank_G3_hiper.rnk", header=FALSE)

G1 <- rbind(G1_hiper, G1_hipo)
G2 <- rbind(G2_hiper, G2_hipo)
G3 <- rbind(G3_hiper, G3_hipo)

G1 <- G1[order(G1$V1), ]
rownames(G1) <- NULL
G2 <- G2[order(G2$V1), ]
rownames(G2) <- NULL
G3 <- G3[order(G3$V1), ]
rownames(G3) <- NULL

identical(G1$V1, G2$V1)
identical(G1$V1, G3$V1)
identical(G2$V1, G3$V1)

G1 <- rename(G1, SYMBOL = V1)
G1 <- rename(G1, G1 = V2)
G2 <- rename(G2, G2 = V2)
G3 <- rename(G3, G3 = V2)

G1G2 <- cbind(G1, G2$G2)
RankSet <- cbind(G1G2, G3$G3)
RankSet<- rename(RankSet, G2 = 'G2$G2')
RankSet<- rename(RankSet, G3 = 'G3$G3')

df <- RankSet
genes <- bitr(df$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
filtered_RankSet<- RankSet[RankSet$SYMBOL %in% genes$SYMBOL, ]
iden <- identical(filtered_RankSet$SYMBOL, genes$SYMBOL)
print(iden)
genes <- cbind(genes, filtered_RankSet[, 2:4])
genes <- rename(genes, GeneID = ENTREZID)
head(genes)

RankSet_New <- genes
head(RankSet_New)

normalize <- function(x) {
  if (all(x >= 0)) {
    return(rescale(x, to = c(0, 1)))
  } else if (all(x <= 0)) {
    return(rescale(x, to = c(-1, 0)))
  } else {
    return((x - min(x)) / (max(x) - min(x)) * 2 - 1)
  }
}

genes[, 3:5] <- t(apply(genes[, 3:5], 1, function(row) normalize(row)))


RankSet_Norm <- genes
head(RankSet_Norm)

RankSet_GeneID <- genes %>% select(-SYMBOL)
head(RankSet_GeneID)

write.table(RankSet_GeneID, file = "RankSet.txt", sep = "\t", row.names = FALSE, quote = FALSE)

total_enrich <- enrichKEGG(gene = RankSet_GeneID$GeneID, organism = 'hsa', pvalueCutoff = 0.01)
total_enrich_05 <- total_enrich %>% filter(pvalue < 0.05)
total_enrich_05 <- total_enrich_05@result
write.table(kegg_cm_df_05, file = "kegg_cm_05.txt", sep = "\t", row.names = TRUE, quote = FALSE)

RankSet_2<- RankSet
RankSet_2[, 2:4] <- t(apply(genes[, 2:4], 1, function(row) normalize(row)))
genes_2 <- genes
genes_2[, 3:5] <- t(apply(genes_2[, 3:5], 1, function(row) normalize(row)))


# Check the data types of columns G1, G2, and G3
str(RankSet_2[, 2:4])

# Ensure columns G1, G2, and G3 are numeric
RankSet_2[, 2:4] <- lapply(RankSet_2[, 2:4], as.numeric)

# Apply the normalization function
RankSet_2[, 2:4] <- t(apply(RankSet_2[, 2:4], 1, normalize))
write.table(RankSet_2, file = "RankSet_2.txt", sep = "\t", row.names = FALSE, quote = FALSE)


write.table(RankSet, file = "RankSet_ALL.txt", sep = "\t", row.names = FALSE, quote = FALSE)
