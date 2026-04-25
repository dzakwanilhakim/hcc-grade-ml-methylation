library(TCGAbiolinks)
library(dplyr)

query <- GDCquery(
  project = "TCGA-LIHC", 
  data.category = "Clinical",
  data.type = "Clinical Supplement", 
  data.format = "BCR Biotab"
)
res <- getResults(query)

GDCdownload(query)
clinical_tab_all <- GDCprepare(query)

# All available tables
names(clinical_tab_all)

# columns from clinical_patient
glimpse_patient <- dplyr::glimpse(clinical_tab_all$clinical_patient_lihc)

#rename Patient id columns
df_patient <- df_patient %>% rename(Case.ID = Sample.ID)


# delete rows
cols_df_patient <- c("Case.ID",                
          "gender",                  
          "height_cm_at_diagnosis",                
          "weight_kg_at_diagnosis",                
          "race",                                  
          "ethnicity",                             
          "history_other_malignancy",              
          "history_neoadjuvant_treatment",         
          "tumor_status",                          
          "vital_status",                          
          "family_history_cancer_indicator",       
          "history_hepato_carcinoma_risk_factors", 
          "histologic_diagnosis",                  
          "tumor_grade", 
          "ajcc_staging_edition",
          "ajcc_tumor_pathologic_pt",              
          "vascular_invasion",                     
          "days_to_initial_pathologic_diagnosis",     
          "alpha_fetoprotien_at_procurement",      
          "alpha_fetoprotien_norm_range_lower",   
          "alpha_fetoprotien_norm_range_upper",    
          "platelet_count_preresection",           
          "platelet_norm_range_lower",             
          "platelet_norm_range_upper",             
          "prothrombin_time_INR_at_procurement",   
          "prothrom_time_INR_norm_range_lower",    
          "prothrombin_time_norm_range_upper",     
          "serum_albumin_preresection",            
          "serum_albumin_norm_range_lower",        
          "serum_albumin_norm_range_upper",        
          "bilirubin_total_norm_range_upper",      
          "bilirubin_total_norm_range_lower",      
          "bilirubin_total",                       
          "creatinine_level_preresection",         
          "creatinine_norm_range_lower",           
          "creatinine_norm_range_upper",           
          "ishak_fibrosis_score", 
          "hepatic_inflammation_adj_tissue",       
          "new_tumor_event_dx_indicator",          
          "age_at_diagnosis",                      
          "birth_days_to",
          "year_of_initial_pathologic_diagnosis", 
          "death_days_to",                         
          "icd_o_3_histology",                     
          "informed_consent_verified",             
          "tissue_source_site")
filtered_df_patients <- df_patient %>%
  select(all_of(cols_df_patient))

cols_sheet <- c("Case.ID",
                "Sentrix_ID",
                "Sample.ID",
                "Sentrix_Position",
                "prior_treatment",
                "Sample.Type")
filtered_sheet <- sheet  %>%
  select(all_of(cols_sheet))


#merge
df_clinical_hcc <- merge(filtered_df_patients, filtered_sheet, by = "Case.ID")
df_clinical_hcc<- df_clinical_hcc %>% 
  filter(grepl("Primary Tumor", Sample.Type, fixed = TRUE))


#impute
df_clinical_hcc <- df_clinical_hcc %>% mutate(tumor_grade = ifelse(Case.ID == "TCGA-FV-A23B", "G1", tumor_grade))

df_clinical_hcc %>% filter(tumor_grade == "[Not Available]")

#NA convert
# Define the vector of values to be converted to NA
na <- c("[Not Available]", "[Not Evaluated]", "[Unknown]")
# Convert matching values to NA using mutate(across())
df_clinical_hcc_na <- df_clinical_hcc %>%
  mutate_all(~ replace(., . %in% na, NA))

df_tumor_grade_hcc <- df_clinical_hcc_na %>%
  filter(!is.na(!!sym("tumor_grade")))

#normal sheet
normal_sheet <- filtered_sheet %>% filter(Sample.Type != "Primary Tumor")
normal_sheet <- merge(normal_sheet, filtered_df_patients, by = "Case.ID")
normal_sheet <- normal_sheet %>%
  mutate_all(~ replace(., . %in% na, NA))

# Create a histogram
df_tumor_grade_hcc$age_at_diagnosis <- as.numeric(df_tumor_grade_hcc$age_at_diagnosis)

# Remove rows with NA values in the "age_at_diagnosis" column
df_tumor_grade_hcc <- na.omit(df_tumor_grade_hcc)

# Create a histogram without ggplot2
hist(df_tumor_grade_hcc$age_at_diagnosis, 
     main = "Age Distribution (excluding NA)",
     xlab = "Age",
     ylab = "Frequency",
     col = "blue",
     border = "black",
     breaks = 5)
write.csv(df_tumor_grade_hcc, file = "Tumor_grade_HCC349.csv", row.names = FALSE)
write.csv(normal_sheet, file = "normal_sheet_tcga.csv", row.names = FALSE)
