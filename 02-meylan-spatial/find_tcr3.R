#!/usr/bin/env Rscript
# Final targeted TCR hunt: descend into Assay (meta.data/misc) and
# DimReduc (misc) which the earlier walk skipped, plus neighbors.
suppressMessages(library(SeuratObject))
obj <- readRDS("/content/GSE296577_GBM_oHSV_seurat_Full_20260107.rds")

check_char <- function(x, path) {
  if (is.character(x)) {
    hit <- grep("TCR-P", x, value = TRUE)
    if (length(hit) > 0) cat(path, "| character n=", length(x),
        "| sample:", paste(head(unique(hit), 3), collapse="; "), "\n")
  } else if (is.factor(x)) {
    hit <- grep("TCR-P", levels(x), value = TRUE)
    if (length(hit) > 0) cat(path, "| factor levels hit:",
        paste(head(unique(hit), 3), collapse="; "), "\n")
  } else if (is.data.frame(x)) {
    for (cc in colnames(x)) {
      col <- x[[cc]]
      if (is.character(col) || is.factor(col)) {
        vals <- if (is.factor(col)) levels(col) else col
        hit <- grep("TCR-P", vals, value = TRUE)
        if (length(hit) > 0) cat(path, "$", cc, "| df nrow=", nrow(x),
            "| sample:", paste(head(unique(hit), 3), collapse="; "), "\n")
      }
    }
  }
}

cat("=== Assay5 slots ===\n")
for (an in Assays(obj)) {
  a <- obj[[an]]
  cat(an, "class:", class(a), "| slots:", paste(slotNames(a), collapse=","), "\n")
  for (sl in slotNames(a)) {
    v <- tryCatch(slot(a, sl), error = function(e) NULL)
    if (is.null(v)) next
    if (sl %in% c("counts", "data", "scale.data")) next  # big matrices
    check_char(v, paste0(an, "@", sl))
    if (is.list(v)) for (nm in names(v)) check_char(v[[nm]], paste0(an, "@", sl, "$", nm))
  }
}
cat("\n=== DimReduc misc ===\n")
for (rn in Reductions(obj)) {
  r <- obj[[rn]]
  m <- tryCatch(r@misc, error = function(e) NULL)
  cat(rn, "misc length:", length(m), "names:", paste(names(m), collapse=","), "\n")
  if (length(m) > 0) for (nm in names(m)) check_char(m[[nm]], paste0(rn, "@misc$", nm))
}
cat("\n=== neighbors ===\n")
for (nn in names(obj@neighbors)) {
  nb <- obj@neighbors[[nn]]
  cat(nn, "class:", class(nb), "| slots:", paste(slotNames(nb), collapse=","), "\n")
}
cat("\n=== RNA feature names containing TR / TCR ===\n")
print(grep("TR[AB]", rownames(obj[["RNA"]]), value = TRUE))
cat("\nDONE\n")
