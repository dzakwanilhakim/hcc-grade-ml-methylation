gdc_sample_HCC750$File.Name <- gsub("_.+", "", gdc_sample_HCC750$File.Name)
gdc_sample_HCC750 <- gdc_sample_HCC750 %>% rename(Sentrix_ID = File.Name)
gdc_sample_HCC750 <- gdc_sample_HCC750%>% select(-File.ID)
gdc_sample_HCC750 <- unique(gdc_sample_HCC750)

#clinical
library(dplyr)
clinical_sheet_HCC750<- clinical_sheet_HCC750 %>% select(-treatment_or_therapy, -treatment_type)
clinical_sheet_HCC750 <- unique(clinical_sheet_HCC750)
clinical_sheet_HCC750 <- clinical_sheet_HCC750 %>% rename(Case.ID = case_submitter_id)


#combined
clinical_hcc_merge <- merge(gdc_sample_HCC750, clinical_sheet_HCC750, by = "Case.ID")
#clinical_hcc_merge <- merge(id_hcc, clinical_hcc, by = "Case.ID")

#NA handling
NA_candidate <- c("'--", "not reported", "TX","NX","MX","Not Reported")
clinical_hcc_merge <- clinical_hcc_merge %>% mutate_all(~replace(., . %in%  NA_candidate, NA))
na_percent <- colMeans(is.na(clinical_hcc_merge))
na_percent
# Drop columns with more than 50% NA values
threshold <- 0.5 # 50%
clinical_hcc_merge <- clinical_hcc_merge[, na_percent <= threshold]
colMeans(is.na(clinical_hcc_merge))

#unique value
unique_counts_list <- list()

for (col in names(clinical_hcc_merge)) {
  unique_counts <- table(clinical_hcc_merge[[col]])
  unique_counts_list[[col]] <- as.data.frame(unique_counts)
}

#final
clinical_hcc_merge <- clinical_hcc_merge %>% select(-Data.Category, 
                                               -Data.Type, 
                                               -project_id,
                                               -days_to_diagnosis,
                                               -icd_10_code,
                                               -prior_treatment,
                                               -site_of_resection_or_biopsy,
                                               -synchronous_malignancy,
                                               -tissue_or_organ_of_origin)


#samplesheet preparation
clinical_hcc_merge$Sentrix_Position <- "noid"
clinical_hcc_merge$Sample_Plate <- NA
clinical_hcc_merge$Pool_ID <- NA
clinical_hcc_merge$Sample_Well <- NA

write.csv(clinical_hcc_merge, file = "samplesheetHCC350_RAW.csv", row.names = FALSE)

#final filter
clinical_hcc_merge <- clinical_hcc_merge %>% filter(ajcc_pathologic_t != "NA")


count_stage <- data_filtered %>%count(Sample_Name)
print(count_stage)

data_filtered$Sample_Name <- data_filtered$cancer_stage

#GEO Sample preparation
geoDir <- "D:/TA/GSE61278_RAW_normal_methylation"
geofilenames <- list.files(geoDir, pattern = "\\Grn.idat$", full.names = TRUE)
idat_geofilenames <- basename(geofilenames)
print(idat_geofilenames)

df_geo <- data.frame(Sentrix_ID = idat_geofilenames)
df_geo$Sentrix_Position <- df_geo$Sentrix_ID
df_geo$Sentrix_ID <- sub("^(.*?)_(.*?)_(.*)$", "\\1_\\2", df_geo$Sentrix_ID)

df_geo$Sentrix_Position <- idat_geofilenames
# Split the strings in col1 based on underscore "_"
split_col <- strsplit(as.character(df_geo$Sentrix_Position), "_")

# Extract the third element from each split segment
df_geo$Sentrix_Position <- sapply(split_col, function(x) ifelse(length(x) >= 3, x[[3]], ""))


# Print the updated data frame
print(df_geo$Sentrix_Position)
write.csv(df_geo, file = "presampsheet_geo.csv")


#combine
TCGA_samplesheet <- data_filtered
samplesheet_merge<- rbind(TCGA_samplesheet, GSE61278_samplesheet)
samplesheet_merge$Sample_Name <- samplesheet_merge$cancer_stage
samplesheet_merge <- samplesheet_merge %>% rename(Sample_Group = cancer_stage)
count_group <- samplesheet_merge %>%count(Sample_Group)
print(count_group)


#set sample_name
write.csv(samplesheet_merge, file = "samplesheet_final.csv", row.names = FALSE)

