library(igraph)
#edge_list
edge_list <- as_edgelist(ppi_compmax)
combined_score <- edge_attr(ppi_compmax)
score <- data.frame(combined_score=combined_score[[1]])
length(E(ppi_compmax)) #367492
length(combined_score[[1]]) #367492
length(edge_list[,1]) #367492
ppi_edgelist <- cbind(edge_list, score)
colnames(ppi_edgelist) <- c("protein1", "protein2", "combined_score")
ppi_edgelist <- unique(ppi_edgelist)
ppi_edgelist_inv <- ppi_edgelist[, c("protein2", "protein1", "combined_score")]
colnames(ppi_edgelist_inv) <- c("protein1", "protein2", "combined_score")
combined_ppi_edgelist <- rbind(ppi_edgelist, ppi_edgelist_inv)
combined_ppi_edgelist <- combined_ppi_edgelist[order(combined_ppi_edgelist$protein1, 
                                                     combined_ppi_edgelist$protein2), ]
rownames(combined_ppi_edgelist)
combined_ppi_edgelist$combined_score <- combined_ppi_edgelist$combined_score/1000

all(ph_label$Sample_Name == ph_sort$Sample_Name) #TRUE
#So set order based on ph_label/ph_sort

# Use match function to find the corresponding indices in df
sorted_indices <- match(ph_label$Sample_Name, rownames(beta_700_compmax))
# Reorder df based on sorted_indices
beta_700_compmax_sort <- beta_700_compmax[sorted_indices, ]
# Sort column names alphabetically
sorted_colnames <- sort(colnames(beta_700_compmax_sort))
# Reorder the columns of the dataframe
beta_700_compmax_sort <- beta_700_compmax_sort[, sorted_colnames]
beta_round <- round(beta_700_compmax_sort, digits = 5)
all(rownames(beta_round) == beta_round$sample_id)
rownames(beta_round)<- NULL

dataframe <- cbind(ph_label, beta_round)
beta_round <- cbind(sample_id = rownames(beta_round),beta_round)

all(rownames(dataframe) == dataframe$Sample_Name)

rownames(dataframe) <- NULL

all(ph_sort$Sample_Name == rownames(beta_700_compmax_sort))
all(unique(combined_ppi_edgelist$protein1) == colnames(beta_700_compmax_sort))

#label
labels <- ph_label[, 2]
# Convert labels into a single-row dataframe
label_df <- as.data.frame(t(labels))
# Rename the column of the label dataframe if necessary
colnames(label_df) <- ph_label$Sample_Name
all(colnames(label_df) == rownames(beta_700_compmax_sort))


write.table(combined_ppi_edgelist, file = "ppi_edgelist.txt", sep = "\t", row.names = FALSE)
write.table(beta_round, file = "beta_433.txt", sep = "\t", row.names = FALSE)
write.table(label_df, file = "label.txt", sep = "\t", row.names = FALSE)
write.table(dataframe, file = "beta_label_433.txt", sep = "\t", row.names = FALSE)

#output
#1.beta_700_compmax_sort
#2.combined_ppi_edgelist
#3.label_df
