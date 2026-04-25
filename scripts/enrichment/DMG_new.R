#cek is all feat_lsvm 247 in feat psvm 370?

'output: 
feats_400
feats_350
feats_370'

X_df <- read.table("X_smt_df.txt", header = TRUE, check.names = FALSE)
X_df_t <- t(X_df)
colnames(X_df_t) <- y_df$Sample_Name
identical(y_df$Sample_Name, colnames(X_df_t))

beta_400 <- X_df_t
beta_400 <- beta_400[rownames(beta_400) %in% feats_400$gene, ]

all(names(bin_p_56_01) %in% rownames(beta_400))
identical(names(bin_p_56_01), rownames(beta_400))

gene2cpg <- function(beta, bin) {
  for (rowname in rownames(beta)) {
    rownames(beta)[rownames(beta) == rowname] <- bin[[rowname]][1]
  }
  return(beta)
}

gene2cpg <- function(bins, beta){
  first_cg_mapping <- sapply(bins, `[`, 1)
  # Replace row names in beta_NG1
  rownames(beta) <- first_cg_mapping[rownames(beta)]
  return(beta)
}

beta_400_cpg <- gene2cpg(beta_400, bin_p_56_01)


feats_lsvm247 <- head(feats_lsvm, 247)
feats_psvm370 <- head(feats_psvm, 370)
feats_psvmg2g3 <- head(feats_psvm, 350)
all(rownames(feats_lsvm247) %in% rownames(feats_psvm370))
setdiff(rownames(feats_lsvm247), rownames(feats_psvm370))
setdiff(rownames(feats_lsvm247), rownames(feats_psvmg2g3))

combine <- unique(c(rownames(feats_psvm370), rownames(feats_lsvm247)))
feats_400 <- data.frame(gene = combine)
feats_350 <- feats_psvmg2g3
all(rownames(feats_lsvm247) %in% feats_400$gene)
all(rownames(feats_psvm370) %in% feats_400$gene)

all(rownames(feats_psvm) %in% colnames(X_df))
all(rownames(feats_lsvm) %in% colnames(X_df))

setdiff(rownames(feats_psvm), colnames(X_df))
setdiff(rownames(feats_lsvm), colnames(X_df))


replace_label <- function(ph, label){
  ph_replace <- ph
  ph_replace$label <- ifelse(ph_replace$label == label, label, 'BG') 
  return(ph_replace)
}

ph_normal_else <- replace_label(y_df, 'NT')
ph_g1_else <- replace_label(y_df, 'G1')
ph_g2_else <- replace_label(y_df, 'G2')
ph_g3_else <- replace_label(y_df, 'G3')

#1. search for most hypermethylation or hypomethylation
'
input: beta_400_cpg, ph_g1_else
'
filter <- function(ph, listpair, beta) {
  # Step 1: Filter ph to create ph_filter
  ph_filter <- ph[ph$label %in% listpair, ]
  
  # Step 2: Map ph_filter$Sample_Name to the column names of beta and filter out unmapped columns
  sample_names <- ph_filter$Sample_Name
  beta_filter <- beta[, colnames(beta) %in% sample_names]
  
  # Ensure the column names of beta_filter are identical to the Sample_Name in ph_filter
  identical_check <- identical(colnames(beta_filter), sample_names)
  print(paste("Column names identical:", identical_check))
  
  return(list(ph_filter = ph_filter, beta_filter = beta_filter))
}

NT_G1_set <- filter(y_df, c('NT','G1'), beta_400_cpg)
NT_G2_set <- filter(y_df, c('NT','G2'), beta_400_cpg)
NT_G3_set <- filter(y_df, c('NT','G3'), beta_400_cpg)

DMP_NT_G1 <- champ.DMP(beta = NT_G1_set$beta_filter, pheno = NT_G1_set$ph_filter$label, adjPVal = 1)
DMP_NT_G2 <- champ.DMP(beta = NT_G2_set$beta_filter, pheno = NT_G2_set$ph_filter$label, adjPVal = 1)
DMP_NT_G3 <- champ.DMP(beta = NT_G3_set$beta_filter, pheno = NT_G3_set$ph_filter$label, adjPVal = 1)

sigs_NT_G1 <- DMP_NT_G1$NT_to_G1
sigs_NT_G2 <- DMP_NT_G2$NT_to_G2
sigs_NT_G3 <- DMP_NT_G3$NT_to_G3

sigs_NT_G1$log2FC <- ifelse(sigs_NT_G1$NT_AVG != 0, log2(sigs_NT_G1$G1_AVG / sigs_NT_G1$NT_AVG), NA)
sigs_NT_G2$log2FC <- ifelse(sigs_NT_G2$NT_AVG != 0, log2(sigs_NT_G2$G2_AVG / sigs_NT_G2$NT_AVG), NA)
sigs_NT_G3$log2FC <- ifelse(sigs_NT_G3$NT_AVG != 0, log2(sigs_NT_G3$G3_AVG / sigs_NT_G3$NT_AVG), NA)

ranking_df <- function(sigs_df){
  rank_df <- sigs_df
  rank_df$rank <- -log10(rank_df$P.Value) * sign(rank_df$log2FC)
  #rank_df <- subset(rank_df, select = c(gene, rank))
  rank_df <- rank_df[order(rank_df$rank, decreasing = TRUE), ]
  return(rank_df)
}

sigs_NT_G1 <-ranking_df(sigs_NT_G1)
sigs_NT_G2 <-ranking_df(sigs_NT_G2)
sigs_NT_G3 <-ranking_df(sigs_NT_G3)

any(is.na(sigs_NT_G1$log2FC))
any(is.na(sigs_NT_G2$log2FC))
any(is.na(sigs_NT_G3$log2FC))


any(sigs_NT_G1$log2FC == 0)
any(sigs_NT_G2$log2FC == 0)
any(sigs_NT_G3$log2FC == 0)



# Define the function to split the data frame
split_dataframe <- function(df) {
  # Filter for positive log2FC
  df_pos <- subset(df, rank > 0)
  # Filter for negative log2FC
  df_neg <- subset(df, rank < 0)
  # Return a list containing the two data frames
  return(list(df_pos = df_pos, df_neg = df_neg))
}


'
rnk_NT_G1 <- ranking(sigs_NT_G1)
rnk_NT_G2 <- ranking(sigs_NT_G2)
rnk_NT_G3 <- ranking(sigs_NT_G3)
'
sigs_NT_G1_split <- split_dataframe(sigs_NT_G1)
G1_pos <- sigs_NT_G1_split$df_pos
G1_neg <- sigs_NT_G1_split$df_neg

sigs_NT_G2_split <- split_dataframe(sigs_NT_G2)
G2_pos <- sigs_NT_G2_split$df_pos
G2_neg <- sigs_NT_G2_split$df_neg

sigs_NT_G3_split <- split_dataframe(sigs_NT_G3)
G3_pos <- sigs_NT_G3_split$df_pos
G3_neg <- sigs_NT_G3_split$df_neg

'
rnk_NT_G1 <- ranking(sigs_NT_G1)
rnk_NT_G2 <- ranking(sigs_NT_G2)
rnk_NT_G3 <- ranking(sigs_NT_G3)
'
setdiff(G3_pos$gene, G2_pos$gene)
setdiff(G3_pos$gene, G1_pos$gene)
setdiff(G3_neg$gene, G2_neg$gene)
setdiff(G3_neg$gene, G1_neg$gene)

rank_generator <- function(sigs, sigs_ref){
  new_sigs <- sigs[sigs$gene %in% sigs_ref$gene,]
  rank_df <- subset(new_sigs, select = c(gene, rank))
  return(rank_df)
}
rank_G1_hiper <- rank_generator(sigs_NT_G1, G3_pos)
rank_G1_hipo <- rank_generator(sigs_NT_G1, G3_neg)

rank_G2_hiper <- rank_generator(sigs_NT_G2, G3_pos)
rank_G2_hipo <- rank_generator(sigs_NT_G2, G3_neg)

rank_G3_hiper <- rank_generator(sigs_NT_G3, G3_pos)
rank_G3_hipo <- rank_generator(sigs_NT_G3, G3_neg)

rank_G1_hiper <- rnk_NT_G1[rnk_NT_G1$gene %in% G3_pos$gene, ]
rank_G2_hiper <- rnk_NT_G2[rnk_NT_G2$gene %in% G3_pos$gene, ]
rank_G3_hiper <- rnk_NT_G3[rnk_NT_G3$gene %in% G3_pos$gene, ]

rank_G1_hipo <- rnk_NT_g1[rnk_NT_g1$gene %in% G3_neg$gene, ]
rank_G2_hipo <- rnk_NT_g2[rnk_NT_g2$gene %in% G3_neg$gene, ]
rank_G3_hipo <- rnk_NT_g3[rnk_NT_g3$gene %in% G3_neg$gene, ]

write_rank <- function(rank_df, filename){
  write.table(rank_df, filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

write_rank(rank_G1_hiper, 'rank_G1_hiper.rnk')
write_rank(rank_G2_hiper, 'rank_G2_hiper.rnk')
write_rank(rank_G3_hiper, 'rank_G3_hiper.rnk')

write_rank(rank_G1_hipo, 'rank_G1_hipo.rnk')
write_rank(rank_G2_hipo, 'rank_G2_hipo.rnk')
write_rank(rank_G3_hipo, 'rank_G3_hipo.rnk')


all(rank_G3_hipo$gene %in% hipo_g3$Genes)

all(G3_neg$gene %in% rank_G3_hipo$gene)

all(G1_neg$gene %in% G2_neg$gene)

write.table(G3_neg$gene, 'hipo_g3.txt', quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(G3_pos$gene, 'hiper_g3.txt', quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

library(ChAMP)
DMP_G1_else <- champ.DMP(beta = beta_400_cpg, pheno = ph_g1_else$label, adjPVal = 1)
DMP_G1_else_01 <- champ.DMP(beta = beta_400_cpg, pheno = ph_g1_else$label, adjPVal = 0.01)
DMP_G2_else <- champ.DMP(beta = beta_400_cpg, pheno = ph_g2_else$label, adjPVal = 1)
DMP_G2_else_01 <- champ.DMP(beta = beta_400_cpg, pheno = ph_g2_else$label, adjPVal = 0.01)
DMP_G3_else <- champ.DMP(beta = beta_400_cpg, pheno = ph_g3_else$label, adjPVal = 1)
DMP_G3_else_01 <- champ.DMP(beta = beta_400_cpg, pheno = ph_g3_else$label, adjPVal = 0.01)
DMP_NT_else <- champ.DMP(beta = beta_400_cpg, pheno = ph_normal_else$label, adjPVal = 1)

sigs_g1_else <- DMP_G1_else$BG_to_G1
sigs_g1_else_01 <- DMP_G1_else_01$BG_to_G1
sigs_g2_else <- DMP_G2_else$BG_to_G2
sigs_g2_else_01 <- DMP_G2_else_01$BG_to_G2
sigs_g3_else <- DMP_G3_else$BG_to_G3
sigs_g3_else_01 <- DMP_G3_else_01$BG_to_G3
#2. search for DMP each class vs normal
#3. GSEA
