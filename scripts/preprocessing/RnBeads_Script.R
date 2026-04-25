library(RnBeads)
library(beepr)

#0. preparation
setwd("D:/TA/epigenomics2016_code/src")
dataDir <- file.path(getwd(), "data")
resultDir <- file.path(getwd(), "results")
#dataSource <- c(idatDir, sampleSheet)
datasetDir <- file.path(dataDir, "Ziller2011_PLoSGen_450K")
idatDir <- file.path(datasetDir, "dataset", "idat")
sampleSheet <- file.path(datasetDir, "dataset", "sample_annotation.csv")
reportDir <- file.path(resultDir, "report_Ziller2011_stepByStep")
samp <- file.path(idatDir, "sample_annotation.csv")
dataSource <- c(idatDir, sampleSheet)
gc()
#1. import
rnbs <- rnb.execute.import(dataSource, data.type="infinium.idat.dir")

#2. Filtering
rnb.set.unfiltered <- rnbs
nrow(meth(rnb.set.unfiltered)) # the number of sites in the unfiltered object

#.p-Value (optional)
any.bad.p.val <- apply(dpval(rnb.set.unfiltered)>0.01, 1, any)
rnb.set.unfiltered <- remove.sites(rnb.set.unfiltered, any.bad.p.val)
nsites(rnb.set.unfiltered)

# 2.1 Remove probes outside of CpG context
rnb.set.filtered <- rnb.execute.context.removal(rnb.set.unfiltered)$dataset
nrow(meth(rnb.set.filtered)) # the number of CpG sites in the unfiltered object

# 2.2 SNP filtering allowing no SNPs in the probe sequence
rnb.set.filtered <- rnb.execute.snp.removal(rnb.set.filtered, snp="any")$dataset
nrow(meth(rnb.set.filtered))

# 2.3 Remove CpGs on sex chromosomes
rnb.set.filtered <- rnb.execute.sex.removal(rnb.set.filtered)$dataset
nrow(meth(rnb.set.filtered))

# 2.4 Remove probes and samples based on a greedy approach
#greedy <- rnb.execute.greedycut(rnb.set.filtered)
#nrow(meth(rnb.set.filtered))

# 2.5 Remove probes containing NA for beta values
rnb.set.filtered <- rnb.execute.na.removal(rnb.set.filtered, 
                                           threshold = 0.05)$dataset
nrow(meth(rnb.set.filtered))

# 2.6 Remove probes for which the beta values have low standard deviation
rnb.set.filtered <- rnb.execute.variability.removal(rnb.set.filtered, 0.005)$dataset
nrow(meth(rnb.set.filtered))

# 2.7 Imputation
rnb.set.filtered.impute <- rnb.execute.imputation(rnb.set.filtered,
                                           method = "knn",
                                           update.ff = TRUE)
nrow(meth(rnb.set.filtered.impute))

rm(rnb.set.unfiltered)
beep(sound = 8, expr = NULL)

#3. Normalization
rnb.set.unnorm <- rnb.set.filtered.impute
rnb_norm <- rnb.execute.normalization(rnb.set.unnorm, method="bmiq",
                                      bgcorr.method="methylumi.noob")
beep(sound = 8, expr = NULL)


rm(rnb.set.filtered)
rm(rnbs)
rm(rnb.set.filtered.impute)
rm(rnb.set.unnorm)

