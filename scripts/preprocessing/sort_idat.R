dir <- getwd()

destDir <- file.path(dir, "data", "new_data_sorted")
sourceDir <- file.path(dir, "data", "new_data")

# List all ".idat" files recursively in the source directory
idat_files <-list.files(path = sourceDir, pattern = "\\.idat$", recursive = TRUE, full.names = TRUE)
idatidat_files <-sub("_.*", "", idat_files)
idat_files <- unique(idat_files)
# Move all ".idat" files to the destination directory
for (file_path in idat_files) {
  file.copy(file_path, file.path(destDir, basename(file_path)), overwrite = TRUE)
  cat(file_path, "has been copied to", destDir, "\n")
}

Sentrix_ID <- samplesheet_HCCNormal451$Sentrix_ID

# Find files in vector2 that are not in vector1
files_not_in_directory <- setdiff(Sentrix_ID, idat_files)

selected_elements <- files_not_in_directory[1:26]

red <- paste(selected_elements, "_noid_Red.idat", sep = "")
green <- paste(selected_elements, "_noid_Grn.idat", sep = "")
red_green <- c(red, green)
red_green

filtered_manifest <- subset(gdc_manifest_20240121_060937, filename %in% red_green)

write.table(filtered_manifest, file = "filtered_manifest.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(df, file = file_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# Print the result# Print the resultidat_files
print(files_not_samplesheet)
