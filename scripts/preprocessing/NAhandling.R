svd.Mval <-  champ.SVD(beta = Mval %>% as.data.frame(),
                       rgSet=NULL,
                       pd=ph,
                       RGEffect=FALSE,
                       PDFplot=TRUE,
                       Rplot=TRUE,
                       resultsDir=resultDir)
#unique value
unique_values <- unique(ph$)

#check NA
Mval <- mval(rnb.filtering)
if (any(is.na(Mval))) {
  print("There are NA values in the dataframe.")
} else {
  print("There are no NA values in the dataframe.")
}

# Identify rows with at least one NA value
rows_with_na <- apply(Mval, 1, function(row) any(is.na(row)))

#not reported replacing
ph <- replace(ph, is.na(ph), 'Not Reported')


# Calculate the percentage of NA values in each row with at least one NA
percent_na <- apply(Mval[rows_with_na, ], 1, function(row) {
  na_count <- sum(is.na(row))
  total_columns <- length(row)
  percent_na <- (na_count / total_columns) * 100
  return(percent_na)
})

# Print the result
print(percent_na)

na_count <- sum(is.na(Mval))
na_count
