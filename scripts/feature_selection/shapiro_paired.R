library(dplyr)
"
input:
beta.global: beta with no NA gene
ph_label:
"


beta_global_433 <- beta_global_433[ph_label$Sample_Name, , drop = FALSE]
identical(rownames(beta_global_433), ph_label$Sample_Name)
beta_global_433 <- round(beta_global_433, 5)
rm(beta_global_t)
#==================================================================================================
#Shapiro_ function
norm_shap <- function(beta, label) {
  norm_feat <- c()
  unnorm_feat <- c()
  count_label <- length(unique(label))
  idx <- 0
  count_feat <- length(colnames(beta))
  for (feat in colnames(beta)) {
    idx <- idx + 1
    progress <- round(idx/count_feat, 2)
    normality <- by(beta[[feat]], label, shapiro.test)
    normal <- TRUE
    for (i in 1:count_label) {
      pval <- normality[[i]]$p.value
      if (pval > 0.05) {
        normal <- FALSE
        break
      }
    }
    if (normal) {
      norm_feat <- c(norm_feat, feat)
      print(paste(progress, '%', feat, 'normal'))
    } else {
      unnorm_feat <- c(unnorm_feat, feat)
      print(paste(progress, '%', feat, '~normal'))
    }
  }
  return(list(norm_feat = norm_feat, unnorm_feat = unnorm_feat))
}

normality_cpg <- norm_shap(beta_global_433, ph_label$label)

#==================================================================================================
#pairwise
pairwise_test <- function(beta, label, normal) {
  sigs <- data.frame(genes = character(),
                     `1-0` = numeric(),
                     `2-0` = numeric(),
                     `3-0` = numeric(),
                     `2-1` = numeric(),
                     `3-1` = numeric(),
                     `3-2` = numeric(),
                     stringsAsFactors = FALSE)
  colnames(sigs) <- c("genes", "1-0", "2-0", "3-0", "2-1", "3-1", "3-2")
  total_feats <- length(colnames(beta))
  idx <- 0
  for (feat in colnames(beta)) {
    idx <- idx + 1
    if (feat %in% normal) {
      sig <- aov(beta[[feat]] ~ label)
      sig.sum <- summary(sig)
      sig.pval <- sig.sum[[1]]$`Pr(>F)`[1]
      if (sig.pval < 0.05) {
        sig_pair <- TukeyHSD(sig)
        sig_pair.pval <- sig_pair$label[, 4]
        sig_pair_df <- data.frame(genes = feat,
                                  `1-0` = sig_pair.pval['1-0'],
                                  `2-0` = sig_pair.pval['2-0'],
                                  `3-0` = sig_pair.pval['3-0'],
                                  `2-1` = sig_pair.pval['2-1'],
                                  `3-1` = sig_pair.pval['3-1'],
                                  `3-2` = sig_pair.pval['3-2'],
                                  stringsAsFactors = FALSE)
        colnames(sig_pair_df) <- c("genes", "1-0", "2-0", "3-0", "2-1", "3-1", "3-2")
        sigs <- rbind(sigs, sig_pair_df)
      }
    } else {
      sig <- kruskal.test(beta[[feat]] ~ label)
      sig_pval <- sig[['p.value']]
      if (sig_pval < 0.05) {
        sig_pair <- pairwise.wilcox.test(beta[[feat]], label, p.adjust.method = "BH")
        sig_pair.pval <- sig_pair[["p.value"]]
        sig_pair_df <- data.frame(genes = feat,
                                  `1-0` = sig_pair.pval[1],
                                  `2-0` = sig_pair.pval[2],
                                  `3-0` = sig_pair.pval[3],
                                  `2-1` = sig_pair.pval[5],
                                  `3-1` = sig_pair.pval[6],
                                  `3-2` = sig_pair.pval[9],
                                  stringsAsFactors = FALSE)
        colnames(sig_pair_df) <- c("genes", "1-0", "2-0", "3-0", "2-1", "3-1", "3-2")
        sigs <- rbind(sigs, sig_pair_df)
      } 
    }
    percent_done <- round((idx / total_feats) * 100, 2)
    cat("Progress:", percent_done, "% ", feat, "\n")
  }
  return(sigs)
}

sigs_pairwise <- pairwise_test(beta_global_433, ph_label$label, normality_cpg$norm_feat)

#==================================================================================================
#filter promoter and non promoter
library(dplyr)
promoter.features <- c('TSS1500', 'TSS200', "5'UTR")
body.features <- 'Body'

cpg_05 <- cpg.global[rownames(cpg.global) %in% sigs_pairwise$genes,]
cpg.p <- cpg_05 %>% filter(feature %in% promoter.features)
cpg.b <- cpg_05 %>% filter(feature %in% body.features)
'
> length(cpg.p$gene)
[1] 103342
> length(cpg.b$gene)
[1] 103620
'
sigs_p <- sigs_pairwise[sigs_pairwise$genes %in% rownames(cpg.p),]
sigs_b <- sigs_pairwise[sigs_pairwise$genes %in% rownames(cpg.b),]
'
> length(sigs_p$genes)
[1] 103342
> length(sigs_b$genes)
[1] 103620
'

less_than_0.05_p <- sigs_p[, 2:7] < 0.05
less_than_0.01_p <- sigs_p[, 2:7] < 0.01


less_than_0.05_b <- sigs_b[, 2:7] < 0.05
less_than_0.01_b <- sigs_b[, 2:7] < 0.01

p05_p <- rowSums(less_than_0.05_p)
p01_p <- rowSums(less_than_0.01_p)
p05_b <- rowSums(less_than_0.05_b)
p01_b <- rowSums(less_than_0.01_b)

sigs_p <- cbind(sigs_p, p05_p, p01_p)
sigs_b <- cbind(sigs_b, p05_b, p01_b)

table(sigs_p$p05)
table(sigs_p$p01)
table(sigs_b$p05)
table(sigs_b$p01)
'
> table(sigs_p$p05)

    0     1     2     3     4     5     6 
 1443  8277 14816 65923  8192  4392   299 
> table(sigs_p$p01)

    0     1     2     3     4     5     6 
 9315 10367 17334 61934  3348  1025    19 
> table(sigs_b$p05)

    0     1     2     3     4     5     6 
 1359  8686 15187 64185  8720  5111   372 
> table(sigs_b$p01)

    0     1     2     3     4     5     6 
 9639 11387 17130 60291  3823  1320    30 
'



round(table(sigs_p$p05)*100/length(sigs_p$genes), 2)
round(table(sigs_p$p01)*100/length(sigs_p$genes), 2)
round(table(sigs_b$p05)*100/length(sigs_b$genes), 2)
round(table(sigs_b$p01)*100/length(sigs_b$genes), 2)

'
> round(table(sigs_p$p05)*100/length(sigs_p$genes), 2)

    0     1     2     3     4     5     6 
 1.40  8.01 14.34 63.79  7.93  4.25  0.29 
> round(table(sigs_p$p01)*100/length(sigs_p$genes), 2)

    0     1     2     3     4     5     6 
 9.01 10.03 16.77 59.93  3.24  0.99  0.02 
> round(table(sigs_b$p05)*100/length(sigs_b$genes), 2)

    0     1     2     3     4     5     6 
 1.31  8.38 14.66 61.94  8.42  4.93  0.36 
> round(table(sigs_b$p01)*100/length(sigs_b$genes), 2)

    0     1     2     3     4     5     6 
 9.30 10.99 16.53 58.18  3.69  1.27  0.03 
'




#============================================================================================
#filter 56
sigs_p_56 <- subset(sigs_p, sigs_p$p05_p %in% c(5, 6))
length(sigs_p_56$genes) #4691
cpg_p_56 <- cpg.p[rownames(cpg.p) %in% sigs_p_56$genes, ]
length(unique(cpg_p_56$gene)) #3556

sigs_p_56_01 <- subset(sigs_p, sigs_p$p01_p %in% c(5, 6))
length(sigs_p_56_01$genes) #1044
cpg_p_56_01 <- cpg.p[rownames(cpg.p) %in% sigs_p_56_01$genes, ]
length(unique(cpg_p_56_01$gene)) #928


sigs_b_56 <- subset(sigs_b, sigs_b$p05_b %in% c(5, 6))
length(sigs_b_56$genes) #[1] 5483
cpg_b_56 <- cpg.b[rownames(cpg.b) %in% sigs_b_56$genes, ]
length(unique(cpg_b_56$gene)) #[1] 3315

sigs_b_56_01 <- subset(sigs_b, sigs_b$p01_b %in% c(5, 6))
length(sigs_b_56_01$genes) #[1] 1350
cpg_b_56_01 <- cpg.b[rownames(cpg.b) %in% sigs_b_56_01$genes, ]
length(unique(cpg_b_56_01$gene)) #[1] 1042


count_pval <- function(sigs_df,pval) {
  # Create an empty dataframe to store the counts
  p_value_counts <- data.frame(column_name = colnames(sigs_df)[2:7], count = NA)
  # Loop through each column from column 2 to column 7 of sigs_p
  for (col in colnames(sigs_df)[2:7]) {
    # Count the number of rows where the p-value is less than 0.05 for the current column
    count <- sum(sigs_df[[col]] < pval, na.rm = TRUE)
    # Store the count in the dataframe
    p_value_counts[p_value_counts$column_name == col, "count"] <- count
  }
  # Return the dataframe with counts
  return(p_value_counts)
}
count_p <-  count_pval(sigs_p)
count_b <-  count_pval(sigs_b)
count_p56 <-  count_pval(sigs_p_56)
count_b56 <- count_pval(sigs_b_56)
count_p56_01 <- count_pval(sigs_p_56_01)
count_b56_01 <- count_pval(sigs_b_56_01)


#============================================================================================
#filter beta
filter_beta <- function(beta, cpg) {
  filtered_beta <- beta[, colnames(beta) %in% cpg]
}

beta_p_56 <- filter_beta(beta_global_433, rownames(cpg_p_56))
length(colnames(beta_p_56)) == length(cpg_p_56$gene) #TRUE
length(colnames(beta_p_56)) #4691

beta_b_56 <- filter_beta(beta_global_433, rownames(cpg_b_56))
length(colnames(beta_b_56)) == length(cpg_b_56$gene) #TRUE
length(colnames(beta_b_56)) #5483

beta_p_56_01 <- filter_beta(beta_global_433, rownames(cpg_p_56_01))
length(colnames(beta_p_56_01)) == length(cpg_p_56_01$gene) #TRUE
length(colnames(beta_p_56_01)) #1044

beta_b_56_01 <- filter_beta(beta_global_433, rownames(cpg_b_56_01))
length(colnames(beta_b_56_01)) == length(cpg_b_56_01$gene) #TRUE
length(colnames(beta_b_56_01)) #1350

#============================================================================================
#binning
gene_binning <- function(cpg_annotation) {
  # Convert row names into a column
  gene_list <- cpg_annotation
  gene_list$CpG <- rownames(cpg_annotation)
  # Select columns CpG and gene
  gene_list <- gene_list[, c("CpG", "gene")]
  # Arrange by gene
  gene_list <- gene_list[order(gene_list$gene), ]
  # Reset row names
  rownames(gene_list) <- NULL 
  # Split CpGs by gene
  genes_bins <- split(gene_list$CpG, gene_list$gene)
  # Remove empty elements
  genes_bins <- genes_bins[sapply(genes_bins, function(x) length(x) > 0)]
  return(genes_bins)
}

bin_p_56 <- gene_binning(cpg_p_56)
bin_b_56 <- gene_binning(cpg_b_56)
bin_p_56_01 <- gene_binning(cpg_p_56_01)
bin_b_56_01 <- gene_binning(cpg_b_56_01)

#===================================================================================================

averaging <- function(beta, bins) {
  df_beta <- data.frame(row.names = rownames(beta), stringsAsFactors = FALSE)
  # Iterate over the list
  for (i in seq_along(bins)) {
    cpgs <- bins[[i]]
    gene <- names(bins)[i]
    # Calculate row means for the selected CpGs
    means <- rowMeans(beta[cpgs])
    cat(i,'calc:', gene, '\n')
    # Add means as a new column with the gene name
    df_beta[, gene] <- means
  }
  colnames(df_beta) <- toupper(colnames(df_beta))
  df_beta <- round(df_beta, 5)
  return(df_beta)  # Moved outside the loop
}
beta_avg_p56 <- averaging(beta_p_56, bin_p_56)
beta_avg_b56 <- averaging(beta_b_56, bin_b_56)

averaging_noupper <- function(beta, bins) {
  df_beta <- data.frame(row.names = rownames(beta), stringsAsFactors = FALSE)
  # Iterate over the list
  for (i in seq_along(bins)) {
    cpgs <- bins[[i]]
    gene <- names(bins)[i]
    # Calculate row means for the selected CpGs
    means <- rowMeans(beta[cpgs])
    cat(i,'calc:', gene, '\n')
    # Add means as a new column with the gene name
    df_beta[, gene] <- means
  }
  df_beta <- round(df_beta, 5)
  return(df_beta)  # Moved outside the loop
}
beta_avg_noup_56_pr <- averaging_noupper(beta_p_56, bin_p_56)
beta_avg_noup_56_bd <- averaging_noupper(beta_b_56, bin_b_56)
beta_avg_p56_01 <- averaging_noupper(beta_p_56_01, bin_p_56_01)
beta_avg_b56_01 <- averaging_noupper(beta_b_56_01, bin_b_56_01)


length(colnames(beta_avg_p56)) == length(bin_p_56)
length(colnames(beta_avg_b56)) == length(bin_b_56)

#==================================================================================================
#check significance again genes
#Shapiro_ function
norm_shap <- function(beta, label) {
  norm_feat <- c()
  unnorm_feat <- c()
  count_label <- length(unique(label))
  idx <- 0
  count_feat <- length(colnames(beta))
  for (feat in colnames(beta)) {
    idx <- idx + 1
    progress <- round(idx/count_feat, 2)
    normality <- by(beta[[feat]], label, shapiro.test)
    normal <- TRUE
    for (i in 1:count_label) {
      pval <- normality[[i]]$p.value
      if (pval > 0.05) {
        normal <- FALSE
        break
      }
    }
    if (normal) {
      norm_feat <- c(norm_feat, feat)
      print(paste(progress, '%', feat, 'normal'))
    } else {
      unnorm_feat <- c(unnorm_feat, feat)
      print(paste(progress, '%', feat, '~normal'))
    }
  }
  return(list(norm_feat = norm_feat, unnorm_feat = unnorm_feat))
}
identical(rownames(beta_avg_noup_56_pr), ph_label$Sample_Name)
normality_56_pr <- norm_shap(beta_avg_noup_56_pr, ph_label$label)
normality_56_bd <- norm_shap(beta_avg_noup_56_bd, ph_label$label)
normality_p56_01 <- norm_shap(beta_avg_p56_01, ph_label$label)
normality_b56_01 <- norm_shap(beta_avg_b56_01, ph_label$label)
'
> length(normality_56_pr$norm_feat)
[1] 812
> length(normality_56_pr$unnorm_feat)
[1] 2744
'
normality_56_cg_pr <- norm_shap(beta_p_56,ph_label$label)
normality_56_cg_bd <- norm_shap(beta_b_56,ph_label$label)
normality_p56_01_cg <- norm_shap(beta_p_56_01, ph_label$label)
normality_b56_01_cg <- norm_shap(beta_b_56_01, ph_label$label)
'
> length(normality_56_cg_pr$norm_feat)
[1] 1108
> length(normality_56_cg_pr$unnorm_feat)
[1] 3583

'



sigs_beta_56_pr <- pairwise_test(beta_avg_noup_56_pr, ph_label$label, normality_56_pr$norm_feat)
length(sigs_beta_56_pr$genes)
count_sigs_56_pr <- count_pval(sigs_beta_56_pr)
length(colnames(beta_avg_noup_56_pr))

sigs_beta_56_bd <- pairwise_test(beta_avg_noup_56_bd, ph_label$label, normality_56_bd$norm_feat)
length(sigs_beta_56_bd$genes)
count_sigs_56_pr <- count_pval(sigs_beta_56_pr)
length(colnames(beta_avg_noup_56_bd))


sigs_beta_p56_01 <- pairwise_test(beta_avg_p56_01, ph_label$label, normality_p56_01$norm_feat)
sigs_beta_b56_01 <- pairwise_test(beta_avg_b56_01, ph_label$label, normality_b56_01$norm_feat)
length(sigs_beta_p56_01$genes)
count_sigs_p56_01 <- count_pval(sigs_beta_p56_01, 0.01)
count_sigs_b56_01 <- count_pval(sigs_beta_b56_01, 0.01)
length(colnames(beta_avg_noup_56_pr))







#=========================================================================================
'
output( promotoer non gnn):
1. beta avg noup 56 pr
2. ph_label
3. beta avg noup 56 bd
4. beta_avg_p56_01
5. beta_avg_b56_01
'


identical(rownames(beta_avg_noup_56_pr), ph_label$Sample_Name)
identical(rownames(beta_avg_p56_01), ph_label$Sample_Name)
identical(rownames(beta_avg_b56_01), ph_label$Sample_Name)
write.table(beta_avg_noup_56_pr, "beta_avg_56_pr.txt", sep = " ", row.names = FALSE)
write.table(ph_label, "ph_label_universal.txt", sep = " ", row.names = FALSE)
write.table(beta_avg_noup_56_bd, "beta_avg_56_bd.txt", sep = " ", row.names = FALSE)
write.table(beta_avg_p56_01, "beta_p56_01_noppi.txt", sep = " ", row.names = FALSE)
write.table(beta_avg_b56_01, "beta_b56_01_noppi.txt", sep = " ", row.names = FALSE)








#+========================================================================================
#+
library(igraph)
library(STRINGdb)


'
input:
beta_avg_p56_01
beta_avg_b56_01
'
beta_01_ppi_p <- beta_avg_p56_01 
colnames(beta_01_ppi_p ) <- toupper(colnames(beta_01_ppi_p))

beta_01_ppi_bd <- beta_avg_b56_01 
colnames(beta_01_ppi_bd ) <- toupper(colnames(beta_01_ppi_bd))

getOption('timeout')
options(timeout=3600)

string_db <- STRINGdb$new(species=9606, score_threshold=700)
human_graph <- string_db$get_graph() # 19488 elements
length(human_graph) #[1] 16201
na_edges <- E(human_graph)[is.na(E(human_graph)$combined_score)]
print(na_edges)
V_na <- any(is.na(V(human_graph)$name))
E_na <- any(is.na(E(human_graph)))
E_cs_na <- any(is.na(E(human_graph)$combined_score))
print(paste(V_na, E_na, E_cs_na))

generate_ppi <- function(beta) {
  string_db <- STRINGdb$new(species=9606, score_threshold=700)
  genes <- data.frame(gene=colnames(beta))
  mapped_genes <- string_db$map(genes, "gene", removeUnmappedRows = TRUE )
  mapped_genes <- mapped_genes[!duplicated(mapped_genes$STRING_id), ]
  mapped_genes <- mapped_genes[!duplicated(mapped_genes$gene), ]
  unmapped_genes <- setdiff(genes$gene, mapped_genes$gene)
  unmapped_genes <- data.frame(gene=unmapped_genes)
  print(paste('first mapping: ', length(genes$gene), '->', length(mapped_genes$gene)))
  print(paste('removed genes: ', length(unmapped_genes$gene)))
  
  print('string_db get graph...')
  #human_graph <- string_db$get_graph() # 19488 elements
  
  V_na <- any(is.na(V(human_graph)$name))
  E_na <- any(is.na(E(human_graph)))
  E_cs_na <- any(is.na(E(human_graph)$combined_score))
  if (V_na || E_na || E_cs_na) {
    print("There are missing values, breaking function")
    print(V_na)
    print(E_na)
    print(E_cs_na)
    return()
  }
  print("There are no missing values")
  
  print('mapping graph...')
  mapped_vertices <- V(human_graph)[name %in% mapped_genes$STRING_id]
  mapped_human_graph <- subgraph(human_graph, mapped_vertices)
  
  
  V(mapped_human_graph)$name <- mapped_genes$gene[match(V(mapped_human_graph)$name, mapped_genes$STRING_id)]
  
  # Find genes not in vertices
  mapped_notin_vertices <- setdiff(mapped_genes$gene, V(mapped_human_graph)$name)
  # Filter out rows from mapped_genes where gene is not in vertices
  mapped_genes <- mapped_genes[!(mapped_genes$gene %in% mapped_notin_vertices), ]
  
  V_len <- length(V(mapped_human_graph)$name)
  print(paste('length graph: ', V_len))
  print(paste('length mapped_genes ', length(mapped_genes$gene)))
  # Check if all elements in df$name are also listed in df2$name
  if (!all(V(mapped_human_graph)$name %in% mapped_genes$gene)) {
    print("Not all elements in graph are listed in mapped_gene")
    return()
  }
  # Check if all elements in df2$name are also listed in df$name
  if (!all(mapped_genes$gene %in% V(mapped_human_graph)$name)) {
    print("Not all elements in mapped gene are listed in graph")
    return()
  }
  print("graph genes same with mapped gene, vice versa")
  
  
  #filter beta
  genes_mapped <- data.frame(gene=mapped_genes$gene)
  beta_mapped <- beta[, colnames(beta) %in% genes_mapped$gene]
  print(paste('length beta: ',length(colnames(beta_mapped))))
  
  return(list(beta_mapped = beta_mapped, mapped_human_graph = mapped_human_graph, genes_mapped = genes_mapped, unmapped_genes = unmapped_genes))
}

ppiset_01_p <- generate_ppi(beta_01_ppi_p)
ppiset_01_bd <- generate_ppi(beta_01_ppi_bd)

beta_avg <- ppiset_01_p$beta_mapped
ppi_01 <- ppiset_01_p$mapped_human_graph
length(colnames(beta_01_ppi_merge))

beta_01_pp

#export beta, y, ppi
write.table(beta_avg_p56_01, "beta_p56_01_noppi.txt", sep = " ", row.names = FALSE)
write.table(beta_avg_b56_01, "beta_b56_01_noppi.txt", sep = " ", row.names = FALSE)



result_hcc <- result
result_g <- result













if (!is.null(result)) {
  beta_hcc_ppi <- result$beta_mapped
  ppi_hcc <- result$mapped_human_graph
  genes_hcc <- result$genes_mapped
  unmap_genes_hcc <- result$unmapped_genes
} else {
  # Handle the case where the function returns NULL
  print("Error in generating PPI. Check the function output.")
}
#==========================================================================================================
all(colnames(beta_g.p.avg.05) %in% colnames(beta_hcc.p.avg.05))


#==========================================================================================================
write.table(ph_g, "ph_g.txt", sep = " ", row.names = FALSE)
write.table(beta_g.p.avg.05, "beta_g.txt", sep = " ", row.names = FALSE)


sigs_bionfo <- sigs_g.np.05[sigs_g.np.05$kruskal_pval < 0.0002, ]
length(sigs_bionfo$genes)  

beta_bioinfo <- beta_filter(beta_g.np.avg, sigs_bionfo$genes) #8270
identical(rownames(beta_bioinfo), ph_g$Sample_Name)
rownames(beta_bioinfo) <- NULL
beta_bioinfo_df <- cbind(ph_g, beta_bioinfo)

library(dplyr)

# Assuming df is your dataframe

# Group by the label column and sample 40 rows from each group
sampled_df <- beta_bioinfo_df %>%
  group_by(label) %>%
  sample_n(size = 40, replace = FALSE) %>%
  ungroup()
sampled_df$label <- ifelse(sampled_df$label == 0, "G1",
                           ifelse(sampled_df$label == 1, "G2",
                                  ifelse(sampled_df$label == 2, "G3", sampled_df$label)))
write.csv(sampled_df, file = "beta_bionfo.csv", row.names = FALSE)









#==================================================================================================
#create ph for hcc and g
ph_hcc <- ph_label
ph_hcc$label <- ifelse(ph_hcc$label != 0, 1, ph_hcc$label)
table(ph_hcc$label)
#0   1 
#100 333 

ph_g <- ph_label
ph_g <- ph_g[ph_g$label != 0, ]
ph_g$label <- ifelse(ph_g$label == 1, 0,
                     ifelse(ph_g$label == 2, 1,
                            ifelse(ph_g$label == 3, 2, ph_g$label)))
table(ph_g$label)
#0   1   2 
#46 167 120

#==================================================================================================
#beta_hcc
beta_hcc <- beta_global_433
# Match the row names of beta_hcc to the Sample_Name column of ph_hcc and reorder beta_hcc
beta_hcc <- beta_hcc[match(ph_hcc$Sample_Name, rownames(beta_hcc)), , drop = FALSE]
all(ph_hcc$Sample_Name == rownames(beta_hcc))

#beta_g
beta_g <- beta_global_433
beta_g <- subset(beta_g, rownames(beta_g) %in% ph_g$Sample_Name)
beta_g <- beta_g[match(ph_g$Sample_Name, rownames(beta_g)), , drop = FALSE]
all(ph_g$Sample_Name == rownames(beta_g))



#==================================================================================================
#Kruskal function

kruskal <- function(beta, ph) {
  sigs <- data.frame(genes = character(), kruskal_pval = numeric(), stringsAsFactors = FALSE)
  i <- 0
  for (gene in colnames(beta)) {
    # Perform Kruskal-Wallis test for each gene
    sig <- kruskal.test(beta[[gene]] ~ ph$label)
    sigs <- rbind(sigs, data.frame(genes = gene, kruskal_pval = sig$p.value))
    i <- i + 1
    print(paste(i, "Gene:", gene))
  }
  return(sigs)
}

beta_hcc <- round(beta_hcc, 5)
beta_g <- round(beta_g, 5)
save.image('1_data_prep2.RData')
gc()
sigs_hcc <- kruskal(beta_hcc, ph_hcc)
sigs_g <- kruskal(beta_g, ph_g)


#==================================================================================================
#filter p<0.05 function
sigs_hcc_05 <- sigs_hcc[sigs_hcc$kruskal_pval < 0.05, ]
length(sigs_hcc_05$genes) #[1] 240801
hist(sigs_hcc_05$kruskal_pval, main = "sigs_hcc_05", xlab = "p-value", ylab = "Frequency")

sigs_g_05 <- sigs_g[sigs_g$kruskal_pval < 0.05, ]
length(sigs_g_05$genes) # 62845
hist(sigs_g_05$kruskal_pval, main = "sigs_g_05", xlab = "p-value", ylab = "Frequency")


#==================================================================================================
#inspect GG and CGI.features function

library(ChAMP)
CpG.GUI(CpG=sigs_g_05$genes, arraytype = '450K')
CpG.GUI(CpG=sigs_hcc_05$genes, arraytype = '450K')



#==================================================================================================
#filter cpg list promoter or non promoter (body)
promoter.features <- c('TSS1500', 'TSS200', "5'UTR")
body.features <- 'Body'
cpg.hcc <- cpg.global[rownames(cpg.global) %in% sigs_hcc_05$genes,]
cpg.hcc.p <- cpg.hcc %>% filter(feature %in% promoter.features)
cpg.hcc.np <- cpg.hcc %>% filter(feature %in% body.features)
length(unique(cpg.hcc.p$gene))
length(unique(cpg.hcc.np$gene))


cpg.g <- cpg.global[rownames(cpg.global) %in% sigs_g_05$genes,]
cpg.g.p <- cpg.g %>% filter(feature %in% promoter.features)
cpg.g.np <- cpg.g %>% filter(feature %in% body.features)
length(unique(cpg.g.p$gene))
length(unique(cpg.g.np$gene))
CpG.GUI(CpG=rownames(cpg.g.p), arraytype = '450K')
CpG.GUI(CpG=rownames(cpg.g.np), arraytype = '450K')

#==================================================================================================
#filter beta
beta_hcc.p <- beta_hcc[, colnames(beta_hcc) %in% rownames(cpg.hcc.p)]
beta_hcc.np <- beta_hcc[, colnames(beta_hcc) %in% rownames(cpg.hcc.np)]

beta_g.p <- beta_g[, colnames(beta_g) %in% rownames(cpg.g.p)]
beta_g.np <- beta_g[, colnames(beta_g) %in% rownames(cpg.g.np)]

#==================================================================================================
# Genes binning

library(dplyr)




bins_hcc_p <- gene_binning(cpg.hcc.p)
bins_hcc_np <- gene_binning(cpg.hcc.np)
bins_g_p <- gene_binning(cpg.g.p)
bins_g_np <- gene_binning(cpg.g.np)


#==================================================================================================
#Averaging


beta_g.p.avg <- averaging(beta_g.p, bins_g_p)
beta_g.np.avg <- averaging(beta_g.np, bins_g_np)



rownames(beta_hcc.p.avg) <- beta_hcc.p.avg$rownames
beta_hcc.p.avg <- beta_hcc.p.avg[, -which(names(beta_hcc.p.avg) == "rownames")]
beta_hcc.p.avg <- round(beta_hcc.p.avg, 5)
colnames(beta_hcc.p.avg ) <- toupper(colnames(beta_hcc.p.avg))


#==================================================================================================
#inspect again the significance
kruskal_filter <- function(beta, ph, pval) {
  sigs <- kruskal(beta,ph)
  sigs <- sigs[sigs$kruskal_pval < pval, ]
  length_diff <- length(rownames(sigs))-length(colnames(beta))
  print(paste('length: ', length(rownames(sigs))))
  print(paste("length difference:", length_diff))
  return(sigs)
}

sigs_g.p.05 <- kruskal_filter(beta_g.p.avg, ph_g, 0.05) #[1] 10097 [1] -1679
sigs_g.np.05 <- kruskal_filter(beta_g.np.avg, ph_g, 0.05) #"8270 length difference: -1161"
sigs_hcc.p.05 <- kruskal_filter(beta_hcc.p.avg, ph_hcc, 0.05) #length:  17185" "length difference: -1057"
sigs_hcc.np.05 <- kruskal_filter(beta_hcc.np.avg, ph_hcc, 0.05) #"length:  14864" "length difference: -961"

beta_filter <- function(beta, ref){
  beta <- beta[,colnames(beta) %in% ref]
  print(length(colnames(beta)))
  return(beta)
}

beta_g.p.avg.05 <- beta_filter(beta_g.p.avg, sigs_g.p.05$genes) #[1] 10097
beta_g.np.avg.05 <- beta_filter(beta_g.np.avg, sigs_g.np.05$genes) #8270
beta_hcc.p.avg.05 <- beta_filter(beta_hcc.p.avg, sigs_hcc.p.05$genes) #17185
beta_hcc.np.avg.05 <- beta_filter(beta_hcc.np.avg, sigs_hcc.np.05$genes) #14864


#==================================================================================================

