library(STRINGdb)
library(igraph)
library(biomaRt)
data(diff_exp_example1)

# 1. getSTRINGdb for human
string_db <- STRINGdb$new(species=9606)
getOption('timeout')
options(timeout=200)
getOption('timeout')
human_graph <- string_db$get_graph()

# 2. NA Handling
na_edges <- E(human_graph)[is.na(E(human_graph)$combined_score)]
print(na_edges)
# + 0/1858944 edges from 17baff0 (vertex names):
# Assuming E(human_graph)$combined_score is the edge attribute 'combined_score'
edge.scores <- E(human_graph)$combined_score
E(human_graph)[1619937]
# Find the indices of edges with non-NA combined scores
#valid_edges <- which(!is.na(edge.scores))
# Delete edges with NA combined scores
#human_graph <- delete_edges(human_graph, valid_edges)

human_graph_filtered <- delete_edges(human_graph, 1619937)
E(human_graph_filtered)[1619937]

#3. get edges with high confidence score
edge.scores <- E(human_graph_filtered)$combined_score
percentile <- quantile(edge.scores, 0.747)
thresh <- data.frame(name='74.7th_percentile', val=percentile)
#we use
human_graph_filtered <- subgraph.edges(human_graph_filtered,
                              E(human_graph_filtered)[combined_score >= percentile])

#4 Adj matrices
# 3. create adjacency matrix
adj_matrix <- as_adjacency_matrix(human_graph_filtered)



# 5. map gene ids to protein ids

### get gene/protein ids via Biomart
mart=useMart(host = 'https://grch37.ensembl.org/',
             biomart='ENSEMBL_MART_ENSEMBL',
             dataset='hsapiens_gene_ensembl')

### extract protein ids from the human network
protein_ids <- sapply(strsplit(rownames(adj_matrix), '\\.'),
                      function(x) x[2])

### get protein to gene id mappings
mart_results <- getBM(attributes = c("external_gene_name",
                                     "ensembl_peptide_id"),
                      filters = "ensembl_peptide_id", values = protein_ids,
                      mart = mart)

### replace protein ids with gene ids
ix <- match(protein_ids, mart_results$ensembl_peptide_id)
ix <- ix[!is.na(ix)]

newnames <- protein_ids
newnames[match(mart_results[ix,'ensembl_peptide_id'], newnames)] <-
  mart_results[ix, 'external_gene_name']
rownames(adj_matrix) <- newnames
colnames(adj_matrix) <- newnames

ppi <- adj_matrix[!duplicated(newnames), !duplicated(newnames)]
nullrows <- Matrix::rowSums(ppi)==0
ppi <- ppi[!nullrows,!nullrows] ## ppi is the network with gene ids

#2. Mapping
example1_mapped <- string_db$map( diff_exp_example1, "gene", removeUnmappedRows = TRUE )
hits <- example1_mapped$STRING_id[1:200]
string_db$plot_network( hits )
httr::set_config(httr::config(ssl_verifypeer=0L))
