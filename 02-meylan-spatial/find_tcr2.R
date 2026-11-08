#!/usr/bin/env Rscript
# Where are the 43 TCR-P28/P34 clonotype IDs? Check commands, assay
# meta.features, FOV misc, neighbors, and any data.frame nested in the
# object that the earlier walk may have skipped (Assay / DimReduc).
suppressMessages(library(SeuratObject))
obj <- readRDS("/content/GSE296577_GBM_oHSV_seurat_Full_20260107.rds")

cat("=== commands (names + call strings) ===\n")
for (nm in names(obj@commands)) {
  cmd <- obj@commands[[nm]]
  cs <- tryCatch(slot(cmd, "call.string"), error = function(e) "")
  cat("[", nm, "] ", substr(cs, 1, 200), "\n", sep = "")
}
cat("\n=== assay meta.features columns (RNA) ===\n")
mf <- obj[["RNA"]]@meta.features
cat("cols:", paste(colnames(mf), collapse=", "), " rows:", nrow(mf), "\n")
for (cc in colnames(mf)) {
  v <- mf[[cc]]
  if (is.character(v) || is.factor(v)) {
    hit <- grep("TCR", as.character(v), value = TRUE)
    if (length(hit) > 0) cat("col", cc, "TCR hits:", length(hit),
                             "sample:", paste(head(unique(hit),3), collapse="; "), "\n")
  }
}
cat("\n=== TCR / TR genes in RNA features ===\n")
print(grep("TR[AB][VDJ]|TCR", rownames(obj[["RNA"]]), value = TRUE))
cat("\n=== FOV misc ===\n")
for (fv in Images(obj)) {
  m <- obj[[fv]]@misc
  if (length(m) > 0) cat(fv, "misc names:", paste(names(m), collapse=","), "\n")
}
cat("\n=== neighbors slot ===\n")
print(names(obj@neighbors))
cat("\n=== active.ident levels ===\n")
print(head(levels(Idents(obj)), 20))
cat("\n=== seurat_clusters values ===\n")
print(table(obj@meta.data$seurat_clusters))
