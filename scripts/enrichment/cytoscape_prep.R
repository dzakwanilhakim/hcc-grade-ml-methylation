#
beta_cytoscape <- data.frame(Name = rownames(beta_400))
identical(beta_cytoscape$Name, rownames(beta_400))

columns <- c("ENTREZID", "SYMBOL")
# Retrieve the mapping for the gene symbols in beta_cytoscape$Name
gene_mapping <- select(org.Hs.eg.db, keys = beta_cytoscape$Name, columns = columns, keytype = "SYMBOL")
names(gene_mapping)[names(gene_mapping) == "ENTREZID"] <- "geneid"
names(gene_mapping)[names(gene_mapping) == "SYMBOL"] <- "Name"

beta_df <- beta_400
rownames(beta_df) <- NULL
beta_cytoscape <- cbind(gene_mapping, beta_df)

write.table(beta_cytoscape, 'beta_cytoscape.txt', 
            quote = FALSE, sep = "\t", row.names = FALSE)


list_ph_g1 <- ph_g1_else$label
list_ph_g2 <- ph_g2_else$label
list_ph_g3 <- ph_g3_else$label
write.table(list_ph_g1, 'g1_class.cls', quote = FALSE, sep = "\t", row.names = FALSE)
write.table(list_ph_g2, 'g2_class.cls', quote = FALSE, sep = "\t", row.names = FALSE)
write.table(list_ph_g3, 'g3_class.cls', quote = FALSE, sep = "\t", row.names = FALSE)
