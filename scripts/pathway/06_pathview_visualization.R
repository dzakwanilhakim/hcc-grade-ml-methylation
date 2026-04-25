library(dplyr)
library(scales)

# Sample data
df <- data.frame(
  SYMBOL = c("BMP4", "CCL28", "CORO1C", "CXCL16"),
  ENTREZID = c(652, 56477, 23603, 58191),
  G1 = c(49.303, -49.353, -23.015, -15.061),
  G2 = c(54.549, -39.251, -34.772, -18.548),
  G3 = c(83.836, -64.201, -43.159, -24.108)
)
head(df)
# Normalization function
normalize <- function(x) {
  if (all(x >= 0)) {
    return(rescale(x, to = c(0, 1)))
  } else if (all(x <= 0)) {
    return(rescale(x, to = c(-1, 0)))
  } else {
    stop("Mixed signs in the data")
  }
}
head(df)
# Apply the normalization to each row
df[, 3:5] <- t(apply(df[, 3:5], 1, function(row) normalize(row)))

print(df)
