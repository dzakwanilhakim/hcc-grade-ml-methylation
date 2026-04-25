library(readxl)
library(dplyr)
read_xlsx <- function(file, sheet){
  df <- read_excel(file, sheet = sheet)
  df[, 2] <- ifelse(df[, 2] < 0, -1, 1)
  df <- df[, 1:2]
  return(df)
}
combine_unique <- function(df1, df2) {
  combined_df <- rbind(df1, df2)
  unique_df <- unique(combined_df)
  return(unique_df)
}

pos_cm <- read_xlsx('CM.xlsx', 'POS_CM')
neg_cm <- read_xlsx('CM.xlsx', 'NEG_CM')

pos_diff <- read_xlsx('DF.xlsx', 'POS_DIFF')
neg_diff <- read_xlsx('DF.xlsx', 'NEG_DIFF')

pos_pro <- read_xlsx('PRO.xlsx', 'POS_PRO')
neg_pro <- read_xlsx('PRO.xlsx', 'NEG_PRO')

cm <- read_excel('CM.xlsx', 'CM')
diff <- read_excel('DF.xlsx', 'DIFF')
pro <- read_excel('PRO.xlsx', 'PRO')

normalize <- function(x) {
  if (all(x >= 0)) {
    return(rescale(x, to = c(0.25, 1)))
  } else if (all(x <= 0)) {
    return(rescale(x, to = c(-1, 0.25)))
  } else {
    stop("Mixed signs in the data")
    return(rescale(x, to = c(-1, 1)))
  }
}

# Apply the normalization to each row



entrez <- function(df){
  genes <- bitr(df$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  print('1')
  iden <- identical(df$SYMBOL, genes$SYMBOL)
  print(iden)
  genes <- cbind(genes, df[, 2:4])
  print('4')
  genes <- dplyr::rename(genes, GeneID = ENTREZID)
  print(genes)
  genes[, 3:5] <- t(apply(genes[, 3:5], 1, function(row) normalize(row)))
  print(genes)
  genes_entrez <- genes %>% select(-SYMBOL)
  return(list(genes, genes_entrez))
}

cm_set <- entrez(cm)
diff_set <- entrez(diff)
pro_set <- entrez(pro)

write.table(cm_set[[2]], file = "cm_all.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(diff_set[[2]], file = "diff_all.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(pro_set[[2]], file = "pro_all.txt", sep = "\t", row.names = FALSE, quote = FALSE)



kegg_cm <- enrichKEGG(gene = cm_set[[2]]$GeneID, organism = 'hsa', pvalueCutoff = 0.01)
kegg_cm_df <- kegg_cm@result
kegg_cm_df_05 <- kegg_cm_df %>% filter(pvalue < 0.05)
write.table(kegg_cm_df_05, file = "kegg_cm_05.txt", sep = "\t", row.names = TRUE, quote = FALSE)





genes <- bitr(df$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez_genes <- genes$ENTREZID
identical(genes$SYMBOL, df$SYMBOL)
genes <- genes[, -1]
names(genes)[names(c) == "SYMBOL"] <- "GeneID"

library(org.Hs.eg.db)
library(pathview)

# Function for KEGG pathway enrichment analysis
perform_kegg_enrichment <- function(df) {
  # Convert gene symbols to Entrez IDs
  genes <- bitr(df$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  entrez_genes <- genes$ENTREZID
  names(entrez_genes) <- genes$SYMBOL
  
  # Create a named vector of pval_signlogfc values
  pval_signlogfc <- df$G1
  names(pval_signlogfc) <- entrez_genes
  
  # Perform KEGG pathway enrichment analysis
  kegg_enrichment <- enrichKEGG(gene = entrez_genes, organism = 'hsa', pvalueCutoff = 0.05)
  
  return(kegg_enrichment)
}

# Function to visualize and export the pathway
visualize_and_export_pathway <- function(kegg_enrichment, gene_data, pathway_order, title) {
  # Prepare gene data for pathview
  names(gene_data) <- kegg_enrichment@result$GeneID
  
  # Extract the pathway ID based on the order
  pathway_id <- kegg_enrichment@result$ID[pathway_order]
  
  # Define a function to map values to colors
  value_to_color <- function(value) {
    if (value < 0) {
      return("blue")
    } else {
      return("red")
    }
  }
  
  # Apply the function to your pval_signlogfc values
  colors <- sapply(gene_data, value_to_color)
  
  # Visualize the pathway with custom colors
  pv <- pathview(gene.data = gene_data, pathway.id = pathway_id, species = "hsa", 
                 low = list(gene = "blue"), high = list(gene = "red"))
  
  # Export the pathway
  output_file <- paste0(title, ".png")
  png(output_file)
  plot(pv)
  dev.off()
  
  return(pv)
}

visualize_and_export_pathway <- function(kegg_enrichment, gene_data, pathway_order, title_prefix) {
  # Extract the pathway ID based on the order
  pathway_id <- kegg_enrichment@result$ID[pathway_order]
  
  # Define a function to map values to colors
  value_to_color <- function(value) {
    if (value < 0) {
      return("blue")
    } else {
      return("red")
    }
  }
  
  # Loop through each column (G1, G2, G3) to generate separate heatmaps
  for (gene_col in colnames(gene_data)[-1]) { # assuming the first column is GeneID
    # Prepare gene data for pathview
    gene_col_data <- gene_data[[gene_col]]
    names(gene_col_data) <- gene_data$GeneID
    
    # Apply the function to the gene column data
    colors <- sapply(gene_col_data, value_to_color)
    
    # Visualize the pathway with custom colors
    pv <- pathview(gene.data = gene_col_data, pathway.id = pathway_id, species = "hsa", 
                   low = list(gene = "blue"), high = list(gene = "red"))
    
    # Export the pathway
    output_file <- paste0(title_prefix, "_", gene_col, ".png")
    png(output_file)
    plot(pv)
    dev.off()
  }
  
  return(pv)
}

# Example data
df <- data.frame(
  GeneID = c(283, 652, 57053),
  G1 = c(33.989, 49.303, -13.29),
  G2 = c(40.798, 54.549, -4.515),
  G3 = c(83.766, 83.836, -4.69)
)

# Convert gene symbols to Entrez IDs if needed (assuming you have this conversion in your perform_kegg_enrichment function)
# For this example, we will assume the GeneID column is already Entrez IDs

# Perform KEGG enrichment
kegg_enrichment <- perform_kegg_enrichment(cm_entrez)

# Visualize and export pathways for each gene data column
visualize_and_export_pathway(kegg_enrichment, df, pathway_order = 1, title_prefix = "KEGG_Pathway")

enrich_cm <- perform_kegg_enrichment(cm)
head(enrich_cm)

enrich_cm <- perform_kegg_enrichment(cm)
head(enrich_cm)


[1]    "23603"  "11221"  "10457" "57531"  
[11]   "4223"  "57509"   "5270"  
