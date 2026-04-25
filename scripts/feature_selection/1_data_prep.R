library(dplyr)
"
input:
beta.global: beta with no NA gene
ph_label:

output:
beta_global_t
beta_global_433
ph_hcc
ph_g
beta_hcc
beta_g
sigs_hcc_05
sigs_g_05

ML non graph (promoter):
-beta_g.p.avg.05
-ph_g
-beta_hcc.p.avg.05
-ph_hcc

ML graph
-





beta_g.p.avg.o5| ph_hcc |
beta_g.np.avg.05
beta_hcc.p.avg.05| ph_g
beta_hcc.np.avg.05



"


#==================================================================================================
# Tranpose beta.global
beta_global_t <- beta.global %>% t() %>% as.data.frame() #[1] 288747
length(rownames(beta_global_t)) #[1] 446
length(colnames(beta_global_t)) #[1] 288747

#==================================================================================================
#filter based on ph_label Sample_Name
beta_global_433 <- subset(beta_global_t, rownames(beta_global_t) %in% ph_label$Sample_Name)
length(rownames(beta_global_433)) #[1] 433
length(colnames(beta_global_433)) #[1] 288747

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


bins_hcc_p <- gene_binning(cpg.hcc.p)
bins_hcc_np <- gene_binning(cpg.hcc.np)
bins_g_p <- gene_binning(cpg.g.p)
bins_g_np <- gene_binning(cpg.g.np)


#==================================================================================================
#Averaging

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



beta_hcc.p.avg <- averaging(beta_hcc.p, bins_hcc_p)
beta_hcc.np.avg <- averaging(beta_hcc.np, bins_hcc_np)
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

library(igraph)
library(STRINGdb)



generate_ppi <- function(beta) {

  
  genes <- data.frame(gene=colnames(beta))
  mapped_genes <- string_db$map(genes, "gene", removeUnmappedRows = TRUE )
  mapped_genes <- mapped_genes[!duplicated(mapped_genes$STRING_id), ]
  mapped_genes <- mapped_genes[!duplicated(mapped_genes$gene), ]
  unmapped_genes <- setdiff(genes$gene, mapped_genes$gene)
  print(paste('first mapping: ', length(genes$gene), '->', length(mapped_genes$gene)))
  print(paste('removed genes: ', length(unmapped_genes$gene)))
  
  print('string_db get graph...')
  human_graph <- string_db$get_graph() # 19488 elements
  
  V_na <- any(is.na(V(human_graph)$name))
  E_na <- any(is.na(E(human_graph)))
  E_cs_na <- any(is.na(E(human_graph)$combined_score))
  if (V_na || E_na || E_cs_na) {
    print("There are missing values, breaking function")
    return()
  }
  print("There are no missing values")
  
  print('mapping graph...')
  mapped_vertices <- V(human_graph)[name %in% mapped_genes$STRING_id]
  mapped_human_graph <- subgraph(human_graph, mapped_vertices)
  
  V(mapped_human_graph)$name <- mapped_genes$gene[match(V(mapped_human_graph)$name, mapped_genes$STRING_id)]
  
  V_len <- length(V(mapped_human_graph)$name)
  print(paste('length graph: ', V_len))
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
  
  return(beta_mapped, mapped_human_graph, genes_mapped)
  
}
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

result <- generate_ppi(beta_g.p.avg.05)
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
