library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations
library(ggplot2)
library(dplyr)

vc_plot_no_label <- function(sigs_df, logfc, pval, title) {
  # Define diffexpressed based on logFC and adjusted p-value
  sigs_df$diffexpressed <- "NO"
  sigs_df$diffexpressed[sigs_df$logFC > logfc & sigs_df$adj.P.Val < pval] <- "UP"
  sigs_df$diffexpressed[sigs_df$logFC < -logfc & sigs_df$adj.P.Val < pval] <- "DOWN"
  
  # Initialize delabel as gene names
  sigs_df$delabel <- sigs_df$gene
  # Set delabel to NA where diffexpressed is 'NO'
  sigs_df$delabel[sigs_df$diffexpressed == 'NO'] <- NA
  
  sigs_df$area <- -log10(sigs_df$adj.P.Val) * sigs_df$logFC
  sigs_df$area[sigs_df$diffexpressed == 'NO'] <- NA
  
  sigs_df <- sigs_df %>% arrange(desc(area))
  
  theme_set(theme_classic(base_size = 20) +
              theme(
                axis.title.y = element_text(face = "bold", margin = margin(0,5,0,0), size = rel(1), color = 'black'),
                axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(5,0,0,0), size = rel(1), color = 'black'),
                plot.title = element_text(hjust = 0.5)
              ))
  
  # Plotting with corrected color mapping and labels
  plt <- ggplot(data = sigs_df, aes(x = logFC, y = -log10(adj.P.Val), col = diffexpressed)) +
    geom_vline(xintercept = c(-logfc, logfc), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(pval), col = "gray", linetype = 'dashed') +
    geom_point(size = 2) +
    scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red"),
                       labels = c("DOWN" = "Hypomethylated", "NO" = "Not significant", "UP" = "Hypermethylated")) +
    labs(color = "", x = expression("log"[2]*"FC"), y = expression("-log"[10]*"adjp-value")) +
    ggtitle(title)
  
  print(plt)  # Use print to display the plot
  
  return(sigs_df)
}

vc_plot <- function(sigs_df, logfc, pval, title ) {
  # Define diffexpressed based on logFC and adjusted p-value
  sigs_df$diffexpressed <- "NO"
  sigs_df$diffexpressed[sigs_df$logFC > logfc & sigs_df$adj.P.Val < pval] <- "UP"
  sigs_df$diffexpressed[sigs_df$logFC < -logfc & sigs_df$adj.P.Val < pval] <- "DOWN"
  
  # Initialize delabel as gene names
  sigs_df$delabel <- sigs_df$gene
  # Set delabel to NA where diffexpressed is 'NO'
  sigs_df$delabel[sigs_df$diffexpressed == 'NO'] <- NA
  
  sigs_df$area <- -log10(sigs_df$adj.P.Val) * sigs_df$logFC
  sigs_df$area[sigs_df$diffexpressed == 'NO'] <- NA
  
  sigs_df <- sigs_df %>% arrange(desc(area))
  
  theme_set(theme_classic(base_size = 20) +
              theme(
                axis.title.y = element_text(face = "bold", margin = margin(0,5,0,0), size = rel(1), color = 'black'),
                axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(5,0,0,0), size = rel(1), color = 'black'),
                plot.title = element_text(hjust = 0.5)
              ))
  
  # Plotting with corrected color mapping and labels
  plt <- ggplot(data = sigs_df, aes(x = logFC, y = -log10(adj.P.Val), col = diffexpressed, label = delabel)) +
    geom_vline(xintercept = c(-logfc, logfc), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(pval), col = "gray", linetype = 'dashed') +
    geom_point(size = 2) +
    scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red"),
                       labels = c("DOWN" = "Hypomethylated", "NO" = "Not significant", "UP" = "Hypermethylated")) +
    #coord_cartesian(ylim = c(0, 30), xlim = c(-0.5, 0.5)) +
    labs(color = "", #legend_title, 
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"adjp-value")) +
    ggtitle(title) +
    geom_text_repel(aes(label = delabel), size = 3, max.overlaps = Inf)
  
  print(plt)  # Use print to display the plot
  
  return(sigs_df)
}
vc_plot_pval <- function(sigs_df, logfc, pval, title ) {
  # Define diffexpressed based on logFC and adjusted p-value
  sigs_df$diffexpressed <- "NO"
  sigs_df$diffexpressed[sigs_df$logFC > logfc & sigs_df$P.Value < pval] <- "UP"
  sigs_df$diffexpressed[sigs_df$logFC < -logfc & sigs_df$P.Value < pval] <- "DOWN"
  
  # Initialize delabel as gene names
  sigs_df$delabel <- sigs_df$gene
  # Set delabel to NA where diffexpressed is 'NO'
  sigs_df$delabel[sigs_df$diffexpressed == 'NO'] <- NA
  
  sigs_df$area <- -log10(sigs_df$P.Value) * sigs_df$logFC
  sigs_df$area[sigs_df$diffexpressed == 'NO'] <- NA
  
  sigs_df <- sigs_df %>% arrange(desc(area))
  
  theme_set(theme_classic(base_size = 20) +
              theme(
                axis.title.y = element_text(face = "bold", margin = margin(0,5,0,0), size = rel(1), color = 'black'),
                axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(5,0,0,0), size = rel(1), color = 'black'),
                plot.title = element_text(hjust = 0.5)
              ))
  
  # Plotting with corrected color mapping and labels
  plt <- ggplot(data = sigs_df, aes(x = logFC, y = -log10(P.Value), col = diffexpressed, label = delabel)) +
    geom_vline(xintercept = c(-logfc, logfc), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(pval), col = "gray", linetype = 'dashed') +
    geom_point(size = 2) +
    scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red"),
                       labels = c("DOWN" = "Hypomethylated", "NO" = "Not significant", "UP" = "Hypermethylated")) +
    #coord_cartesian(ylim = c(0, 30), xlim = c(-0.5, 0.5)) +
    labs(color = "", #legend_title, 
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"adjp-value")) +
    ggtitle(title) +
    geom_text_repel(aes(label = delabel), size = 3, max.overlaps = Inf)
  
  print(plt)  # Use print to display the plot
  
  return(sigs_df)
}

sigs_g1_else_new <- vc_plot(sigs_g1_else, 0.1, 0.01, 'G1')
sigs_g2_else_new <- vc_plot(sigs_g2_else, 0.1, 0.01, 'G2')
sigs_g3_else_new <- vc_plot(sigs_g3_else, 0.1, 0.01, 'G3')

sigs_g1_else_new_pval <- vc_plot_pval(sigs_g1_else, 0.1, 0.01, 'G1')
sigs_g2_else_new_pval <- vc_plot_pval(sigs_g2_else, 0.1, 0.01, 'G2')
sigs_g3_else_new_pval <- vc_plot_pval(sigs_g3_else, 0.1, 0.01, 'G3')

sigs_NT_G1_NEW <- vc_plot(sigs_NT_G1, 0.25, 0.01, 'G1')
sigs_NT_G2_NEW <- vc_plot(sigs_NT_G2, 0.25, 0.01, 'G2')
sigs_NT_G3_NEW <- vc_plot(sigs_NT_G3, 0.25, 0.01, 'G3')

sigs_NT_G1_NEW <- vc_plot_no_label(sigs_NT_G1, 0.25, 0.01, 'G1')
sigs_NT_G2_NEW <- vc_plot_no_label(sigs_NT_G2, 0.25, 0.01, 'G2')
sigs_NT_G3_NEW <- vc_plot_no_label(sigs_NT_G3, 0.25, 0.01, 'G3')

dev.off()
ranking <- function(sigs_df){
  rank_df <- sigs_df
  rank_df$rank <- -log10(rank_df$P.Value) * sign(rank_df$logFC)
  rank_df <- subset(rank_df, select = c(gene, rank))
  rank_df <- rank_df[order(rank_df$rank, decreasing = TRUE), ]
  return(rank_df)
}
#sigs_NT_G1
rnk_NT_G1 <- ranking(sigs_NT_G1)
rnk_NT_G2 <- ranking(sigs_NT_G2)
rnk_NT_G3 <- ranking(sigs_NT_G3)

library(org.Hs.eg.db)
library(AnnotationDbi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)



# Get the annotation data
annotation_data <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# Extract CpG ID and gene symbol columns
cpg_gene_mapping <- annotation_data[, c("Name", "UCSC_RefGene_Name")]
# Display the first few rows
head(cpg_gene_mapping)
anno_df <- cpg_gene_mapping@listData

anno_table <- as.data.frame(anno_df)
filtered_anno_df <- anno_table[anno_table$Name %in% rownames(rnk_g1), ]
all(rownames(rnk_g1) %in% anno_df$Name)

# Define the columns and keys
columns <- c("SYMBOL")
keys <- keys(org.Hs.eg.db, keytype = "ENTREZID")
# Retrieve the gene symbols
gene_symbols <- select(org.Hs.eg.db, keys = keys, columns = columns, keytype = "ENTREZID")
# Display the first few rows
head(gene_symbols)

all(rnk_g1$gene %in% gene_symbols$SYMBOL)
setdif <- setdiff(rnk_g1$gene, gene_symbols$SYMBOL)

head(anno_df$Name)
rm(annotation_data)
rm(cpg_gene_mapping)



write_rank <- function(rank_df, filename){
  write.table(rank_df, filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

write_rank(rnk_NT_G1, 'rank_NT_G1.rnk')
write_rank(rnk_NT_G2, 'rank_NT_G2.rnk')
write_rank(rnk_NT_G3, 'rank_NT_G3.rnk')

