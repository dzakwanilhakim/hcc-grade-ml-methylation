set.seed(10)
sim.cpd.data2 = matrix(sample(sim.cpd.data, 18000, replace = T), ncol = 6)
rownames(sim.cpd.data2) = names(sim.cpd.data)
colnames(sim.cpd.data2) = paste("exp", 1:6, sep = "")
head(sim.cpd.data2, 3)




# Perform KEGG pathway enrichment analysis
kegg_enrichment <- enrichKEGG(gene = cm_entrez$GeneID, organism = 'hsa', pvalueCutoff = 1)
head(kegg_enrichment)
#KEGG view
pv.out <- pathview(gene.data = cm_entrez, 
                   pathway.id = 'hsa04062',
                   species = "hsa", 
                   out.suffix = "hsa05206", 
                   keys.align = "y", 
                   kegg.native = T, 
                   match.data = F, 
                   multi.state = T, 
                   same.layer = T)

demo.paths$sel.paths

hsa05206




#KEGG view with data match
pv.out <- pathview(gene.data = gse16873.d[, 1:3],
                   cpd.data = sim.cpd.data2[, 1:2], 
                   pathway.id = demo.paths$sel.paths[i],
                   species = "hsa", 
                   out.suffix = "gse16873.cpd.3-2s.match",
                   keys.align = "y", 
                   kegg.native = T, 
                   match.data = T, 
                   multi.state = T,
                   same.layer = T)
#graphviz view
pv.out <- pathview(gene.data = gse16873.d[, 1:3],
                   cpd.data = sim.cpd.data2[, 1:2], 
                   pathway.id = demo.paths$sel.paths[i],
                   species = "hsa", 
                   out.suffix = "gse16873.cpd.3-2s", 
                   keys.align = "y",
                   kegg.native = F, 
                   match.data = F, 
                   multi.state = T, 
                   same.layer = T)



head(pv.out$plot.data.cpd)


head(gse16873.d[, 1:3])
head(sim.cpd.data2[, 1:2])

demo.paths