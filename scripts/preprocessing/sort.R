library(fs)

# Specify the directory path
rawpdfDir <- "D:/TA/Clinical_Data_HCC_TCGA"

# List all PDF files recursively in the specified directory
pdf_files <- fs::dir_ls(rawpdfDir, regexp = "\\.pdf$", recurse = TRUE)


# Print the list of PDF files
print(pdf_files)
