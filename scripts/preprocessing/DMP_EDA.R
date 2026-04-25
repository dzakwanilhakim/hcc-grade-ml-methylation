library(dplyr)
library(ggplot2)
library(ChAMP)
CpG.GUI(rownames(beta))
#extract rownames beta
cpg <- rownames(beta.global)
#filter in probe features
cpg.promoter <- probe.features
cpg.promoter <- cpg.promoter %>% filter(row.names(cpg.promoter) %in% as.character(cpg))
#filter out utr
features <- c('TSS1500', 'TSS200', "5'UTR")
cpg.promoter <- cpg.promoter %>% filter(feature %in% features)
CpG.GUI(rownames(cpg.features))
#intersect beta based on probe features, filter out
beta.promoter <- beta.global
beta.promoter <- beta.promoter[rownames(beta.promoter) %in% rownames(cpg.promoter), ]
dim(beta.promoter)
CpG.GUI(rownames(cpg.global))

#EDA features in beta
##. cgi-feature
# Create a table for heatmap
cgi_features <- table(cpg.features$cgi, cpg.features$feature)

#check NA in gene
any(is.na(cpg.features$gene))
#check NA in cpg
any(is.na(rownames(cpg.features)))
#check NA in cgi
any(is.na(cpg.features$cgi))
#check NA in feature
any(is.na(cpg.features$feature))
#check NA in chromosome
any(is.na(cpg.features$CHR))

#QC.plot
QC.GUI(beta = beta.promoter, pheno = ph$tumor_grade)

#plot gene's probe with count value with order
list_gene <- unique(cpg.features$gene)

#filter in cpg island
cpgisland_features <- cpg.features %>% filter(cgi %in% 'island')
listisland_gene <- unique(cpgisland_features$gene)
beta.island <- beta.promoter[rownames(beta.promoter) %in% rownames(cpgisland_features), ]
dim(beta.island)
QC.GUI(beta = beta.island, pheno = ph$tumor_grade)

#champ.DMP beta island
DMP.promoter <- champ.DMP(beta = beta.promoter,
                          pheno = ph$tumor_grade,
                          compare.group = NULL,
                          adjPVal = 0.05,
                          adjust.method = "BH",
                          arraytype = "450K")
DMP.island <- champ.DMP(beta = beta.island,
                          pheno = ph$tumor_grade,
                          compare.group = NULL,
                          adjPVal = 0.05,
                          adjust.method = "BH",
                          arraytype = "450K")
DMP.beta <- champ.DMP(beta = beta,
                      pheno = ph$tumor_grade,
                      compare.group = NULL,
                      adjPVal = 0.05,
                      adjust.method = "BH",
                      arraytype = "450K")

#2.Global
cpg.global_withna <- probe.features %>% filter(row.names(probe.features) %in% as.character(cpg))
cpg.global <- probe.features %>% filter(row.names(probe.features) %in% as.character(cpg))
#check NA in gene
any(cpg.global$gene == "")
#check NA in cpg
any(rownames(cpg.global) == "")
#check NA in cgi
any(cpg.global$cgi == "")
#check NA in feature
any(cpg.global$feature == "")
#check NA in chromosome
any(cpg.global$CHR == "")

#filter empty cells in gene
cpg.global <- cpg.global %>% filter(gene != "")
cpg.global.na <- cpg.global_withna %>% filter(gene == "")

# Create a table for heatmap
cgi.features.global <- table(cpg.global$cgi, cpg.global$feature)
norm.matrix <- scale(cgi.features.global)
# Plot heatmap using base R
heatmap(norm.matrix, 
        main = "cgi-features global")

listglobal_gene <- unique(cpg.global$gene)
listglobalwithna_gene <- unique(cpg.global_withna$gene)

#beta global
beta.global <- beta[rownames(beta) %in% rownames(cpg.global), ]
dim(beta.global)
beta.na <- beta[rownames(beta) %in% rownames(cpg.global.na), ]
dim(beta.na)

#DMP
DMP.beta.global <- champ.DMP(beta = beta.global,
                      pheno = ph$tumor_grade,
                      compare.group = NULL,
                      adjPVal = 0.05,
                      adjust.method = "BH",
                      arraytype = "450K")
DMP.beta.na <- champ.DMP(beta = beta.na,
                             pheno = ph$tumor_grade,
                             compare.group = NULL,
                             adjPVal = 0.05,
                             adjust.method = "BH",
                             arraytype = "450K")

#empty.gene <- cpg.global$gene == ""
#any(cpg.global$gene == "")
#print(empty.gene)
#plot cgi with count value with order
#plot features with count value with order
#plot chromosomes 
#plot cgi-feature with count value with order

#ggsave(file.path(resultDir, "cgi_features.png"), plot, width = 8, height = 6)
#clustering