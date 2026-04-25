#id
id_normal$File.Name <- gsub("_.+", "", id_normal$File.Name)
id_hcc$File.Name <- gsub("_.+", "", id_hcc$File.Name)
id_normal <- id_normal %>% rename(Sentrix_ID = File.Name)
id_hcc <- id_hcc %>% rename(Sentrix_ID = File.Name)
id_normal <- id_normal%>% select(-File.ID)
id_hcc <- id_hcc%>% select(-File.ID)
id_normal <- unique(id_normal)
id_hcc <- unique(id_hcc)

#clinical
library(dplyr)
clinical_hcc <- clinical_hcc %>% select(-treatment_or_therapy, -treatment_type)
clinical_hcc <- unique(clinical_hcc)
clinical_normal <- clinical_normal %>% select(-treatment_or_therapy, -treatment_type)
clinical_normal <- unique(clinical_normal)
clinical_hcc <- clinical_hcc %>% rename(Case.ID = case_submitter_id)
clinical_normal <- clinical_normal %>% rename(Case.ID = case_submitter_id)

#combined
clinical_normal_merge <- merge(id_normal, clinical_normal, by = "Case.ID")
clinical_hcc_merge <- merge(id_hcc, clinical_hcc, by = "Case.ID")
clinical_hcc_merge <- merge(id_hcc, clinical_hcc, by = "Case.ID")

#cancer_stage
clinical_normal_merge$cancer_stage <- "CT"
clinical_hcc_merge$cancer_stage <- clinical_hcc_merge$ajcc_pathologic_stage#
clinical_hcc_merge$cancer_stage <- sub("Stage ", "", clinical_hcc_merge$cancer_stage)

convert_roman <- function(x) {
  # Define a mapping of Roman numerals to their standardized forms
  roman_map <- c('I'='I', 'II'='II', 'III'='III', 'IV'='IV','IIIA' = 'III', 'IIIB' = 'III', 'IIIC'='III', 'IVA'= 'IV', 'IVB'='IV')
  # Apply the mapping to the input x
  sapply(x, function(y) roman_map[y])
}
clinical_hcc_merge$cancer_stage <- convert_roman(clinical_hcc_merge$cancer_stage)
unique(clinical_hcc_merge$cancer_stage)

#rbind
clinical_combined <- rbind(clinical_hcc_merge, clinical_normal_merge)

#NA handling
na_percent <- colMeans(is.na(clinical_combined))
# Drop columns with more than 50% NA values
threshold <- 0.5 # 50%
data_filtered <- clinical_combined[, na_percent <= threshold]
colMeans(is.na(data_filtered))

#final
data_filtered <- data_filtered %>% select(-Data.Category, -Data.Type, -project_id)

#samplesheet preparation
data_filtered$Sentrix_Position <- "noid"
data_filtered$Sample_Plate <- NA
data_filtered$Pool_ID <- NA
data_filtered$Sample_Well <- NA

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

