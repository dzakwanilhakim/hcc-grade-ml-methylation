library(dplyr)
# Transpose the data frame using dplyr and tibble
beta_t_promoter <- beta.promoter %>% t() %>% as.data.frame()

#ph_sort <- ph %>% arrange(desc(tumor_grade))
# Reset the row index
ph_sort <- ph
ph_sort <- ph_sort %>% mutate(label = tumor_grade)
ph_sort$label <- ifelse(ph_sort$label == 'Normal', 0, ph_sort$label)
ph_sort$label <- ifelse(ph_sort$label == 'G1', 1, ph_sort$label)
ph_sort$label <- ifelse(ph_sort$label == 'G2', 2, ph_sort$label)
ph_sort$label <- ifelse(ph_sort$label == 'G3', 3, ph_sort$label)
ph_sort$label <- ifelse(ph_sort$label == 'G4', 4, ph_sort$label)

ph_sort <- ph_sort %>% arrange(label)

rownames(ph_sort) <- NULL #reset index
ph_sort <- ph_sort %>% filter(label != 4)

beta_t_promoter <- beta_t_promoter[rownames(beta_t_promoter) %in% ph_sort$Sample_Name, , drop = FALSE]

# Identify common row names
#common_row_names <- intersect(rownames(beta_t), ph_sort$Sample_Name)
# Combine the data frames using cbind
#beta_t_label <- data.frame(label = ph_sort[ph_sort$Sample_Name %in% common_row_names, "label", drop = FALSE],
#                          beta_t[common_row_names, , drop = FALSE])

# Assuming ph_sort and beta_t are defined
#common_row_names <- intersect(rownames(ph_sort), rownames(beta_t))

ph_label <- ph_sort[, c("Sample_Name", "label")]
# Merge based on row names
beta_t_promoter <- merge(beta_t_promoter, ph_label, by.x = "row.names", by.y = "Sample_Name", all.x = TRUE)

# Rename the row names column to its original name
rownames(beta_t_promoter) <- beta_t_promoter$Row.names
beta_t_promoter$Row.names <- NULL

#df <- df[, c(column_to_move, setdiff(names(df), column_to_move))]
#beta_t_label$label
# Print the merged data frame
#print(beta_t_label)





beta_t_promoter <- beta_t_promoter %>% arrange(label)
beta_t_promoter$label

write.table(beta_t_label, 'beta_global.tsv', sep = "\t", row.names = TRUE)

#gene count
list_promoter <- cpg.promoter %>% rownames_to_column(var = 'CpG')
list_promoter <- list_promoter[, c("CpG", "gene")]
list_promoter <- list_promoter %>% arrange(gene)
rownames(list_promoter) <- NULL #reset index
genes_promoter <- split(list_promoter$CpG, list_promoter$gene)
genes_promoter <- genes_promoter[sapply(genes_promoter, function(x) length(x) > 0)]

#gene count global
list_global <- cpg.global %>% rownames_to_column(var = 'CpG')
list_global <- list_global[, c("CpG", "gene")]
list_global <- list_global %>% arrange(gene)
genes_global <- split(list_global$CpG, list_global$gene)
genes_global <- genes_global[sapply(genes_global, function(x) length(x) > 0)]

write.csv(df, file = "filename.csv", row.names = FALSE)






#gene_count <- table(cpg.global$gene)

