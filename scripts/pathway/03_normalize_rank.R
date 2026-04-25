normalize_custom <- function(x) {
  if (all(x >= 0)) {
    return(rescale(x, to = c(0.25, 1)))
  } else if (all(x <= 0)) {
    return(rescale(x, to = c(-1, -0.25)))
  } else {
    return(rescale(x, to = c(-1, 1)))
  }
}
rankset_custom <- RankSet_New
rankset_custom[, 3:5] <- t(apply(rankset_custom[, 3:5], 1, function(row) normalize_custom(row)))

write.table(rankset_custom, file = "RankSet_with_SYMBOL.txt", sep = "\t", row.names = FALSE, quote = FALSE)

RankSet_GeneID <- rankset_custom %>% select(-SYMBOL)
write.table(RankSet_GeneID, file = "RankSet_GeneID.txt", sep = "\t", row.names = FALSE, quote = FALSE)

g1_geneid <- semi_join(rankset_custom, g1_importance, by = "SYMBOL")
g2_geneid <- semi_join(rankset_custom, g2_importance, by = "SYMBOL")
g3_geneid <- semi_join(rankset_custom, g3_importance, by = "SYMBOL")
g1_wo_sym <- g1_geneid %>% select(-SYMBOL)
g2_wo_sym <- g2_geneid %>% select(-SYMBOL)
g3_wo_sym <- g3_geneid %>% select(-SYMBOL)
write.table(g1_geneid, file = "g1_geneid.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(g2_geneid, file = "g2_geneid.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(g3_geneid, file = "g3_geneid.txt", sep = "\t", row.names = FALSE, quote = FALSE)

write.table(g1_wo_sym, file = "g1_wo_sym.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(g2_wo_sym, file = "g2_wo_sym.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(g3_wo_sym, file = "g3_wo_sym.txt", sep = "\t", row.names = FALSE, quote = FALSE)

enrich_kegg <- function(rankset_geneid, pval){
  enrich <- enrichKEGG(gene = rankset_geneid, organism = 'hsa', pvalueCutoff = pval)
  df_enrich <- enrich@result
  return(df_enrich)
} 

enrich_g1 <- enrich_kegg(g1_wo_sym$GeneID)
enrich_g2 <- enrich_kegg(g2_wo_sym$GeneID)
enrich_g3 <- enrich_kegg(g3_wo_sym$GeneID)
enrich_g1_no <- enrich_kegg(g1_wo_sym$GeneID, 0.001)
enrich_g2_no <- enrich_kegg(g2_wo_sym$GeneID, 1)
enrich_g3_no <- enrich_kegg(g3_wo_sym$GeneID, 1)
