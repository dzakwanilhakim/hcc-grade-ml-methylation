na_count <- sum(is.na(ph_selected))
na_count
library(dplyr)
columns_with_na <- colSums(is.na(ph_na)) > 0
print(columns_with_na)

ph_na <- ph
ph_na$race <- ifelse(is.na(ph_na$race), "[Not Available]", ph_na$race)
ph_na$ajcc_tumor_pathologic_pt <- ifelse(is.na(ph_na$ajcc_tumor_pathologic_pt), "[Not Available]", ph_na$ajcc_tumor_pathologic_pt)
ph_na$vascular_invasion <- ifelse(is.na(ph_na$vascular_invasion), "[Not Available]", ph_na$vascular_invasion)
ph_na$pt_merge1 <- ifelse(is.na(ph_na$pt_merge1), "[Not Available]", ph_na$pt_merge1)
ph_na$race <- ifelse(is.na(ph_na$race), "[Not Available]", ph_na$race)

ph_na <- select(ph_na, -age_at_diagnosis)
ph_na <- select(ph_na, -barcode)
ph_na <- select(ph_na, -Sample_Well)
ph_na <- select(ph_na, -Pool_ID)
ph_na <- select(ph_na, -Sample_Plate)
ph_selected <- ph_na %>% select(Sample_Name, gender)
ph3 <- unique(ph_na)

[1] "Sample_Name"              "Case.ID"                 
[3] "Sentrix_ID"               "Sentrix_Position"        
[5] "Sample.Type"              "Project.ID"              
[7] "tissue_source_site"       "gender"                  
[9] "race"                     "history_other_malignancy"
[11] "tumor_grade"              "ajcc_tumor_pathologic_pt"
[13] "vascular_invasion"        "tumor_grade_merge1"      
[15] "tumor_grade2"             "pt_merge1"               
[17] "icd_o_3_histology"   