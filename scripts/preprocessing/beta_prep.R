library(RnBeads.hg19)
library(RnBeads)
library(ChAMP)
idat <- "C:/Users/lenovo/AppData/Local/R/win-library/4.3/ChAMPdata/extdata"
samp <- "C:/Users/lenovo/AppData/Local/R/win-library/4.3/ChAMPdata/lung_test_set.csv"
datasaus <- c(idat, samp)
rnbs <- rnb.execute.import(datasaus, data.type="infinium.idat.dir")
rnbs <- rnb.execute.sex.removal(rnbs)$dataset
rnbs <- rnb.execute.na.removal(rnbs,threshold = 0)$dataset
#create cpg annotation


ann.prom<-annotation(rnbs, add.names=TRUE)
#data.champ <- champ.load(idat)
#beta.champ <- data.champ$beta
#pd.champ <- data.champ$pd

#retrieve sample_ID
ph <- pheno(rnbs)
beta <- meth(rnbs)
Mval <- mval(rnbs)
col <- ph$Sample_Name
row <- ann.prom$ID
CpG.GUI(row)
#beta_promoter <- meth(rnbs, type = "promoters")
#row_promoter <- ann.prom$ID

#add column
colnames(beta) <- col
rownames(beta) <- row
colnames(Mval) <- col
rownames(Mval) <- row
#colnames(beta_promoter) <- col
#rownames(beta_promoter) <- row
#SVD
#?champ.SVD
svd <-  champ.SVD(beta = Mval %>% as.data.frame(),
                  rgSet=NULL,
                  pd=ph,
                  RGEffect=FALSE,
                  PDFplot=TRUE,
                  Rplot=TRUE,
                  resultsDir=resultDir)

#Batch effect Correction
beta.correction <- champ.runCombat(beta=Mval,pd=ph,
                                         variablename = "Sample_Group",
                                         batchname=c("Cell_Line"))

#DMP
DMP <- champ.DMP(beta = beta,
                 pheno = ph$Sample_Group,
                 compare.group = NULL,
                 adjPVal = 0.1,
                 adjust.method = "BH",
                 arraytype = "450K")

DMP.mval <- champ.DMP(beta = Mval,
                 pheno = ph$Sample_Group,
                 compare.group = NULL,
                 adjPVal = 0.1,
                 adjust.method = "BH",
                 arraytype = "450K")

DMP_promoter <- champ.DMP(beta = beta_promoter,
                          pheno = ph$Sample_Group,
                          compare.group = NULL,
                          adjPVal = 0.1,
                          adjust.method = "BH",
                          arraytype = "450K")

DMP.champ <- champ.DMP(beta = beta.champ,
                 pheno = data.champ$pd$Sample_Group,
                 compare.group = NULL,
                 adjPVal = 0.1,
                 adjust.method = "BH",
                 arraytype = "450K")

?QC.GUI
QC.GUI(beta=Mval,
       pheno=ph$Sample_Group,
       arraytype="450K")

library(RnBeads.hg19)
data(small.example.object)
logger.start(fname=NA)
dm <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"))
get.comparisons(dm)