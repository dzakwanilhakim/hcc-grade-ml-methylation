library(igraph)
library(dplyr)
library(plotly)
#check graph_ppi connected
is_connected(graph_ppi)
# get the component
ppi_clusters <- components(graph_ppi)
count_components(graph_ppi)
components_gene <- groups(ppi_clusters)
comp_genes <- data.frame(gene=components_gene[[1]])
# Compmax have 13744 vertices

#mapped vertices final
compmax_vertices <- V(graph_ppi)[name %in% comp_genes$gene]
length(compmax_vertices) #14105
any(is.na(compmax_vertices))

# new PPI
ppi_compmax <- subgraph(graph_ppi, compmax_vertices)
length(ppi_compmax)
#[1] 13744
print(ppi_compmax)
length(E(ppi_compmax))
# 367492
any(is.na(V(ppi_compmax)))
any(is.na(E(ppi_compmax)))
any(is.na(E(ppi_compmax)$combined_score))
# no NA data

# ordered beta column
beta_700_compmax <- beta_ppi_compmax[, comp_genes$gene]
# Check if column names are in the desired order
all(colnames(beta_700_compmax) == comp_genes$gene) #TRUE
all(colnames(beta_700_compmax) == adj_ppi_compmax@Dimnames[[1]]) #TRUE

# adj_matrix
adj_matrix <- as.matrix(adj_ppi_compmax)

#create gene-edges count
genes_edges <- data.frame(genes = V(ppi_compmax)$name, edges = degree(ppi_compmax))
rownames(genes_edges) <- NULL

#plot
plot_ly(genes_edges, x = ~edges, type = "histogram") %>%
  layout(title = "Interactive Histogram of Edge Counts",
         xaxis = list(title = "Edges"),
         yaxis = list(title = "Frequency"),
         bargap = 0.05) # Adjust the gap between bars if needed

all(colnames(adj_matrix) == colnames(beta_700_compmax))
#[1] TRUE

# Assuming you have an adjacency matrix named 'adj_matrix'

# Create a heatmap from the adjacency matrix
heatmap(adj_matrix, 
        col = heat.colors(256),  # Choose the color palette for the heatmap
        scale = "column",        # Scale each column individually
        Rowv = NA,               # Disable row clustering
        Colv = NA,               # Disable column clustering
        xlab = "Nodes",          # X-axis label
        ylab = "Nodes",          # Y-axis label
        main = "Heatmap of Adjacency Matrix")  # Title of the heatmap



#export
write.csv(adj_matrix, file = "PPI_StringDB_700.csv", row.names = FALSE)
write.csv(beta_700_compmax, file = "beta_700.csv", row.names = TRUE)


#output for stringDB PPI
#1. beta_700_compmax : beta filtered 13744 gene -> csv
#2. comp_genes: gene 13744 list
#3. ppi_compmax: filtered graph 13744
#4. adj_beta_ppi_compmax: adjmatrix 13744
#5. adj_matrix: matrix 13744x13744 -> csv








