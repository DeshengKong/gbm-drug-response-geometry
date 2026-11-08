#!/usr/bin/env Rscript
# Find where the TCR clonotype mapping lives in the Meylan Seurat object.
# Defensive recursive walk: only character / data.frame / list / S4 nodes
# that could contain TCR strings; skip numeric matrices and DimReduc.
suppressMessages(library(SeuratObject))
obj <- readRDS("/content/GSE296577_GBM_oHSV_seurat_Full_20260107.rds")

found <- list()
walk <- function(x, path, depth = 0) {
  if (is.null(x) || depth > 6) return(invisible(NULL))
  # skip big numeric-only objects
  if (inherits(x, c("DimReduc", "Assay", "dgCMatrix", "dgTMatrix",
                    "matrix", "array"))) return(invisible(NULL))
  if (is.character(x)) {
    hit <- grep("TCR", x, value = TRUE)
    if (length(hit) > 0)
      found[[path]] <<- list(kind = "character", n = length(x),
                             sample = head(unique(hit), 5))
    return(invisible(NULL))
  }
  if (is.factor(x)) {
    hit <- grep("TCR", levels(x), value = TRUE)
    if (length(hit) > 0)
      found[[path]] <<- list(kind = "factor", n = length(x),
                             sample = head(unique(hit), 5))
    return(invisible(NULL))
  }
  if (is.data.frame(x)) {
    for (cc in colnames(x)) {
      col <- x[[cc]]
      if (is.character(col) || is.factor(col)) {
        vals <- if (is.factor(col)) levels(col) else col
        hit <- grep("TCR", vals, value = TRUE)
        if (length(hit) > 0)
          found[[paste0(path, "$", cc)]] <<- list(kind = "df-col",
              n = nrow(x), sample = head(unique(hit), 5))
      }
    }
    return(invisible(NULL))
  }
  if (isS4(x)) {
    nms <- tryCatch(slotNames(x), error = function(e) NULL)
    if (!is.null(nms)) for (nm in nms) {
      v <- tryCatch(slot(x, nm), error = function(e) NULL)
      walk(v, paste0(path, "@", nm), depth + 1)
    }
    return(invisible(NULL))
  }
  if (is.list(x) || is.environment(x)) {
    nms <- tryCatch(names(x), error = function(e) NULL)
    if (!is.null(nms)) {
      for (nm in nms) {
        v <- tryCatch(x[[nm]], error = function(e) NULL)
        walk(v, paste0(path, "$", nm), depth + 1)
      }
    } else {
      for (i in seq_along(x)) {
        v <- tryCatch(x[[i]], error = function(e) NULL)
        walk(v, paste0(path, "[[", i, "]]"), depth + 1)
      }
    }
    return(invisible(NULL))
  }
  invisible(NULL)
}
walk(obj, "obj")
cat("=== TCR-containing nodes ===\n")
if (length(found) == 0) cat("NONE FOUND\n") else {
  for (p in names(found)) {
    cat(p, " | kind:", found[[p]]$kind, "| n:", found[[p]]$n,
        "| sample:", paste(found[[p]]$sample, collapse="; "), "\n")
  }
}
cat("\n=== slot names ===\n"); print(slotNames(obj))
cat("\n=== tools ===\n"); print(names(obj@tools))
cat("\n=== graphs ===\n"); print(names(obj@graphs))
cat("\n=== images slot names ===\n")
for (fv in Images(obj)) {
  f <- obj[[fv]]
  cat(fv, "class:", class(f), "| slots:", paste(slotNames(f), collapse=","), "\n")
}
