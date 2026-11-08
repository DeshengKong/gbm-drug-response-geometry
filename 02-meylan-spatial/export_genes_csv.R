#!/usr/bin/env Rscript
# Export the 22-gene counts as a gzipped CSV (one row per cell, same order
# as meta.data) so the Python analysis can read it without R.
suppressMessages(library(Matrix))
m <- readRDS("/content/meylan_genes.rds")
cat("dim:", paste(dim(m), collapse=" x "), "\n")
# transpose to cells x genes, convert to dense data.frame (22 cols is small)
df <- as.data.frame(as.matrix(t(m)))
df$cell <- rownames(df)
# reorder columns: cell first
df <- df[, c("cell", setdiff(colnames(df), "cell"))]
# align row order to meta.data (colnames of m are already meta order, but verify)
meta <- read.csv(gzfile("/content/meylan_meta.csv.gz"), row.names = 1, nrows = 5)
cat("first meta rownames:", paste(rownames(meta), collapse=","), "\n")
cat("first gene-df cells:", paste(head(df$cell, 5), collapse=","), "\n")
write.csv(df, gzfile("/content/meylan_genes.csv.gz"), row.names = FALSE)
cat("saved meylan_genes.csv.gz\n")
