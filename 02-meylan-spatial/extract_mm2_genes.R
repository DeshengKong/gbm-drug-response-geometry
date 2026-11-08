#!/usr/bin/env Rscript
# Extract the 6 non-VEGFA MM2 (hypoxia) genes per-cell counts from the Meylan
# Seurat object, export as a gzipped long-format CSV (cell, gene, count)
# for nonzero entries only. Runs on Colab (gbm-meylan session).
suppressMessages(library(SeuratObject))
suppressMessages(library(Matrix))

obj <- readRDS("/content/GSE296577_GBM_oHSV_seurat_Full_20260107.rds")
cat("cells:", ncol(obj), " features:", nrow(obj[["RNA"]]), "\n")

genes <- c("MET", "CAV1", "LOX", "ICAM1", "CXCL10", "CD44")
found <- intersect(genes, rownames(obj[["RNA"]]))
cat("genes found:", paste(found, collapse=","), "\n")
cat("genes missing:", paste(setdiff(genes, found), collapse=","), "\n")

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")[found, , drop = FALSE]
cat("counts dim:", paste(dim(counts), collapse=" x "), "\n")
cat("nnz:", length(counts@x), "\n")

ct <- as(counts, "TsparseMatrix")
df <- data.frame(
  cell  = colnames(counts)[ct@j + 1],
  gene  = rownames(counts)[ct@i + 1],
  count = ct@x,
  stringsAsFactors = FALSE
)
cat("long rows:", nrow(df), "\n")
write.csv(df, gzfile("/content/meylan_analysis/mm2_genes_long.csv.gz"),
          row.names = FALSE)
cat("saved /content/meylan_analysis/mm2_genes_long.csv.gz\n")
cat("DONE\n")
