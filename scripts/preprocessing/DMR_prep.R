library(dplyr)
ph_normal <- ph %>% filter(tumor_grade_merge1 == "Normal")
ph_G1 <- ph %>% filter(tumor_grade_merge1 == "G1")
ph_G2 <- ph %>% filter(tumor_grade_merge1 == "G2")
ph_G3G4 <- ph %>% filter(tumor_grade_merge1 == "G3G4")

ph_normal <- rbind(df1, df2)

beta_normal <- beta[, colnames(beta) %in% ph_normal$Sample_Name]
beta_G1 <- beta[, colnames(beta) %in% ph_G1$Sample_Name]
beta_G2 <- beta[, colnames(beta) %in% ph_G2$Sample_Name]
beta_G3G4 <- beta[, colnames(beta) %in% ph_G3G4$Sample_Name]

DMR_myLoad <- champ.DMR(beta=myLoad$beta %>% as.data.frame(),
                       pheno=myLoad$pd$Sample_Group,
                       compare.group=NULL,
                       arraytype="450K",
                       method = "DMRcate",
                       minProbes=7,
                       adjPvalDmr=0.05,
                       cores=2,
                       rmSNPCH=T,
                       fdr=0.05,
                       dist=2,
                       mafcut=0.05,
                       lambda=1000,
                       C=2)

