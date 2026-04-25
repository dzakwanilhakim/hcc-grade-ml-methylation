# Use match function to find the corresponding indices in df
sorted_indices <- match(ph_label$Sample_Name, rownames(beta_ppi))
# Reorder df based on sorted_indices
beta_ppi_sort <- beta_ppi
beta_ppi_sort <- beta_ppi_sort[sorted_indices, ]
all(ph_label$Sample_Name == rownames(beta_ppi_sort))
beta_ppi_sort <- round(beta_ppi_sort, 5) 

#=========================================================================================
#Kruskal
sigs <- data.frame(genes = character(), kruskal_pval = numeric(), stringsAsFactors = FALSE)
i <- 0
for (gene in colnames(beta_ppi_sort)) {
  # Perform Kruskal-Wallis test for each gene
  sig <- kruskal.test(beta_ppi_sort[[gene]] ~ ph_label$label, data = beta_ppi_sort)
  sigs <- rbind(sigs, data.frame(genes = gene, kruskal_pval = sig$p.value))
  i <- i + 1
  print(paste(i, "Gene:", gene))
}

#==========================================================================================
ph_label_noct <- ph_label[ph_label$label != 0, ]
rownames(ph_label_noct) <- NULL
beta_ppi_sort_noct <- beta_ppi_sort
beta_ppi_sort_noct <- beta_ppi_sort_noct[rownames(beta_ppi_sort_noct) %in% ph_label_noct$Sample_Name, ]
all(ph_label_noct$Sample_Name == rownames(beta_ppi_sort_noct))

sigs_noct <- data.frame(genes = character(), kruskal_pval = numeric(), stringsAsFactors = FALSE)
i <- 0
for (gene in colnames(beta_ppi_sort_noct)) {
  # Perform Kruskal-Wallis test for each gene
  sig <- kruskal.test(beta_ppi_sort_noct[[gene]] ~ ph_label_noct$label, data = beta_ppi_sort_noct)
  sigs_noct <- rbind(sigs_noct, data.frame(genes = gene, kruskal_pval = sig$p.value))
  i <- i + 1
  print(paste(i, "Gene:", gene))
}

# Assuming you have already run the code to generate `sigs`
hist(sigs_noct$kruskal_pval, main = "Histogram of p-values", xlab = "p-value", ylab = "Frequency")

sigs_noct_1 <- sigs_noct[sigs_noct$kruskal_pval < 0.1, ]
length(sigs_noct_1$genes) #[1] 3771
hist(sigs_noct_1$kruskal_pval, main = "sign_noct_1", xlab = "p-value", ylab = "Frequency")

sigs_noct_05 <- sigs_noct[sigs_noct$kruskal_pval < 0.05, ]
length(sigs_noct_05$genes) #[1] 2609
hist(sigs_noct_05$kruskal_pval, main = "sign_noct_05", xlab = "p-value", ylab = "Frequency")


sigs_05 <- sigs[sigs$kruskal_pval < 0.05, ]
length(sigs_05$genes)
# Check the structure of the subset dataframe

sigs_01 <- sigs[sigs$kruskal_pval < 0.01, ]
length(sigs_01$genes)
# Check the structure of the subset dataframe

sigs_001 <- sigs[sigs$kruskal_pval < 0.001, ]
length(sigs_001$genes)
# Check the structure of the subset dataframe

sigs_0001 <- sigs[sigs$kruskal_pval < 0.0001, ]
length(sigs_0001$genes)


str(sigs_05)

#==========================================================================================

sigs_05 <- sigs[sigs$kruskal_pval < 0.05, ]
length(sigs_05$genes)
# Check the structure of the subset dataframe

sigs_01 <- sigs[sigs$kruskal_pval < 0.01, ]
length(sigs_01$genes)
# Check the structure of the subset dataframe

sigs_001 <- sigs[sigs$kruskal_pval < 0.001, ]
length(sigs_001$genes)
# Check the structure of the subset dataframe

sigs_0001 <- sigs[sigs$kruskal_pval < 0.0001, ]
length(sigs_0001$genes)


str(sigs_05)

#==========================================================================================










sigs <- data.frame(genes = character(), pval = numeric(), stringsAsFactors = FALSE)
i <- 0
for (gene in colnames(beta_700_compmax_sort)) {
  # Perform Kruskal-Wallis test for each gene
  i <- i + 1
  sig <- kruskal.test(gene ~ ph_label$label, data = beta_700_compmax_sort)
  sigs$genes.append(gene)
  sigs$pval.append(sig)
  print(i, gene)
}

sig <- kruskal.test(as.formula(paste(HIST1H1A, "~ ph_label$label")), data = beta_700_compmax_sort)
sig$p.value

sig <- kruskal.test(HIST1H1A ~ ph_label$label, data = beta_700_compmax_sort)
sig$data.name[1]
pairwise.wilcox.test(beta_ppi_sort$A2BP1, ph_label$label,
                     p.adjust.method = "BH")

