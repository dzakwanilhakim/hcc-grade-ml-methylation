empty_df1 <- data.frame(rownames = rownames(beta_t_promoter), stringsAsFactors = FALSE)
empty_df1 <- cbind(empty_df1, beta_t_promoter$cg12045430)

# Assuming beta_t_promoter and genes_global are defined
library(dplyr)
# Initialize an empty data frame for beta means
beta_t_promoter_num <- beta_t_promoter %>% select(-label)
betameans_promoter <- data.frame(rownames = rownames(beta_t_promoter_num), stringsAsFactors = FALSE)

# Iterate over the list
for (i in seq_along(genes_promoter)) {
  cpgs <- genes_promoter[[i]]
  gene <- names(genes_promoter)[i]
  
  # Calculate row means for the selected CpGs
  means <- rowMeans(beta_t_promoter_num[cpgs])
  cat(i,' Calculating gene:', gene, '\n')
  # Add means as a new column with the gene name
  betameans_promoter[, gene] <- means
}

rownames(betameans_promoter) <- betameans_promoter$rownames
betameans_promoter <- betameans_promoter %>% select(-rownames)

list_promoter_genes <- names(genes_promoter)

#selected_genes <- c('cg03817621','cg24411946')
#means <- rowMeans(beta_t_promoter_num[cpgs])

#means <- rowMeans(beta_t_promoter_num[c("cg03817621","cg24411946"), , drop = FALSE])

#beta_t_promoter_num[, "cg24411946"]





# Display the resulting data frame
print(betameans_promoter)


# Create a sample data frame
data <- data.frame(
  Gene1 = c(1, 2, 3, 4),
  Gene2 = c(5, 6, 7, 8),
  Gene3 = c(9, 10, 11, 12),
  Gene4 = c(1, 2, 3, 4),
  Gene5 = c(5, 6, 7, 8),
  Gene6 = c(9, 10, 11, 12),
  Gene7 = c(1, 2, 3, 4),
  Gene8 = c(5, 6, 7, 8),
  Gene9 = c(9, 10, 11, 12),
  Gene10 = c(1, 2, 3, 4),
  Gene11 = c(5, 6, 7, 8),
  Gene12  = c(9, 10, 11, 12)
)


# Assuming means is a vector
means_ex <- c(1, 2, 3, 4)

# Creating an empty data frame
empty_df1 <- data.frame(NULL)

# Adding means as a new column to empty_df1
empty_df1$new_column <- means_ex
