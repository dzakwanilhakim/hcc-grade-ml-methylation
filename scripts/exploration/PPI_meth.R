library(STRINGdb)
library(igraph)
library(biomaRt)
library(dplyr)

#1. getSTRINGdb for human
string_db <- STRINGdb$new(species=9606)
getOption('timeout')
options(timeout=200)
getOption('timeout')
human_graph <- string_db$get_graph() # 19488 elements
length(human_graph)


#2. create mapping list
# gene list
df_genes <- data.frame(gene=list_promoter_genes)
#change all gene names to upper case in df_genes
df_genes$gene <- toupper(df_genes$gene)
colnames(betameans_promoter) <- toupper(colnames(betameans_promoter))

#check rows value between colnames and df gene
all(df_genes$gene == colnames(betameans_promoter))
# TRUE, colnames have same rows value

# mapping with gene list
mapped <- string_db$map(df_genes, "gene", removeUnmappedRows = TRUE )
dim(mapped)
#Warning:  we couldn't map to STRING 9% of your identifiers
# 18464 -mapped-> 16827 genes
length(V(human_graph)$name)

# is all mapped element also contained in df_genes?
mapped_vs_genes <- all(mapped$gene %in% df_genes$gene) #FALSE because case letter

# If not all values are included, identify the values in df1$gene that are not in df2$gene
if (!mapped_vs_genes) {
  values_not_included1 <- setdiff(mapped$gene, df_genes$gene)
  values_not_included2 <- setdiff(df_genes$gene, mapped$gene)
  print(values_not_included)
}

# mapped vertices in human_graph that intersect with mapped
any(is.na(mapped$STRING_id))
vertices <- data.frame(STRING_id = V(human_graph)$name)
any(is.na(vertices$STRING_id))
mapped_vs_vertices <- all(mapped$STRING_id %in% vertices$STRING_id) #FALSE
if (!mapped_vs_vertices) {
  mapped_notin_vertices <- setdiff(mapped$STRING_id, vertices$STRING_id )
  vertices_notin_mapped <- setdiff(vertices$STRING_id, mapped$STRING_id)
  print(mapped_notin_vertices)
}
#+ } there are 12 genes that under 700 threshold
#[1] "9606.ENSP00000433646" "9606.ENSP00000489684" "9606.ENSP00000500329" "9606.ENSP00000419502"
#[5] "9606.ENSP00000479346" "9606.ENSP00000383474" "9606.ENSP00000386193" "9606.ENSP00000367065"
#[9] "9606.ENSP00000499042" "9606.ENSP00000451945" "9606.ENSP00000246104" "9606.ENSP00000370192"

#delete the 12 genes that in mapped gene
mapped_graph <- mapped[!(mapped$STRING_id %in% mapped_notin_vertices), ] #there are 16815 graph
all(mapped_graph$STRING_id %in% vertices$STRING_id) #TRUE

mapped_vertices <- V(human_graph)[name %in% mapped_graph$STRING_id]
all(mapped_vertices$name %in% vertices$STRING_id)
all(mapped_graph$STRING_id %in% mapped_vertices$name)
#all(mapped_graph == colnames(betameans_promoter))

#check are there double value in mapepd graph$gene? 
duplicate_gene <- table(mapped_graph$gene[duplicated(mapped_graph$gene)]) + 1
duplicate_stringid <- table(mapped_graph$gene[duplicated(mapped_graph$STRING_id)])
duplicate_vertices <- table(mapped_vertices$name[duplicated(mapped_vertices$name)])
# there are duplicated value in gene name and in string ID

#how to deal with the duplicate value?
no_dupID <- mapped_graph[!duplicated(mapped_graph$STRING_id), ]
table(no_dupID$gene[duplicated(no_dupID$gene)]) + 1
table(no_dupID$STRING_id[duplicated(no_dupID$STRING_id)]) + 1

any(is.na(mapped_vertices))
length(mapped_vertices)

#3. NA Handling
na_edges <- E(human_graph)[is.na(E(human_graph)$combined_score)]
print(na_edges)
# + 0/1858944 edges from 17baff0 (vertex names):

#4. get edges with high confidence score
edge.scores <- E(human_graph)$combined_score
percentile <- quantile(edge.scores, 0.746)
thresh <- data.frame(name='74.6th_percentile', val=percentile)
#we use
human_graph_700 <- subgraph.edges(human_graph, E(human_graph)[combined_score >= percentile])
#19488 -700-> 16201
length(E(human_graph))
length(E(human_graph_700))
#edges from  1858944 to 473860


#5. filter filter subgraph based on mapped vertices
# Filter vertices based on the list of nodes to keep
mapped_vertices <- V(human_graph_700)[name %in% mapped$STRING_id]
any(is.na(mapped_vertices))
length(mapped_vertices)
# 14185/16201 vertices, named, from 9eb0897:

graph_mapped <- subgraph(human_graph_700, mapped_vertices)
length(graph_mapped)
#[1] 14185
print(graph_mapped)
length(E(graph_mapped))
# edges from 473860 -> 380318
any(is.na(V(graph_mapped)))
any(is.na(E(graph_mapped)))
any(is.na(E(graph_mapped)$combined_score))
# no NA data

# replace STRING_id into gene name
V(graph_mapped)$name <- mapped$gene[match(V(graph_mapped)$name, mapped$STRING_id)]

# intersect genes
intersect_genes <- data.frame(gene = V(graph_mapped)$name)
any(is.na(V(graph_mapped)$name))

# Display the filtered and renamed subgraph
print(graph_mapped)

##plot
edge.scores_mapped <- E(graph_mapped)$combined_score
hist(edge.scores_mapped)
plot(density(edge.scores_mapped))
plot(density(edge.scores))


#6. Adj matrices
adj_matrix <- as_adjacency_matrix(graph_mapped)

#7. filter the beta_means meth
library(dplyr)
beta_ppi <- betameans_promoter[, colnames(betameans_promoter) %in% intersect_genes$gene]
