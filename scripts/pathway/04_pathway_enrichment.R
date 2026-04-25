require(org.Hs.eg.db)
data(gse16873.d)
head(gse16873.d)

gse16873.t <- apply(gse16873.d, 1, function(x) t.test(x,alternative = "two.sided")$p.value)
gse16873.t

sel.genes <- names(gse16873.t)[gse16873.t < 0.1]
head(sel.genes)

sim.cpd.data=sim.mol.data(mol.type="cpd", nmol=3000)
head(sim.cpd.data)

sel.cpds <- names(sim.cpd.data)[abs(sim.cpd.data) > 0.5]

i <- 1
pv.out <- pathview(gene.data = sel.genes, cpd.data = sel.cpds,
                   pathway.id = demo.paths$sel.paths[i], species = "hsa", 
                   out.suffix = "sel.genes.sel.cpd",
                   keys.align = "y", kegg.native = T, key.pos = demo.paths$kpos1[i],
                   limit = list(gene = 5, cpd = 2), bins = list(gene = 5, cpd = 2),
                   na.col = "gray", discrete = list(gene = T, cpd = T))

pv.out <- pathview(gene.data = sel.genes, cpd.data = sim.cpd.data,
                   pathway.id = demo.paths$sel.paths[i], species = "hsa", out.suffix = "sel.genes.cpd",
                   keys.align = "y", kegg.native = T, key.pos = demo.paths$kpos1[i],
                   limit = list(gene = 5, cpd = 1), bins = list(gene = 5, cpd = 10),
                   na.col = "gray", discrete = list(gene = T, cpd = F))
head(gse16873.d)