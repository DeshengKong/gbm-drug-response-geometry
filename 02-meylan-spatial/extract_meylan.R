#!/usr/bin/env Rscript
# Extract lightweight analysis tables from the Meylan Seurat object.
# Runs on Colab (gbm-meylan session) after SeuratObject is installed.
# Outputs (all in /content):
#   meylan_meta.csv.gz      - meta.data (all 16 cols, ~2.4M rows)
#   meylan_centroids.csv.gz - per-cell FOV + x/y from the images slot
#   meylan_genes.rds        - sparse counts subset for genes of interest
#   meylan_structure.txt    - object structure + TCR search
RDS <- "/content/GSE296577_GBM_oHSV_seurat_Full_20260107.rds"
suppressMessages(library(SeuratObject))

cat("reading RDS ...\n")
obj <- readRDS(RDS)
cat("class:", class(obj), "\n")
sink("meylan_structure.txt")
cat("class:", class(obj), "\n")
cat("assays:", paste(Assays(obj), collapse=", "), "\n")
cat("reductions:", paste(Reductions(obj), collapse=", "), "\n")
cat("images/FOVs:", paste(Images(obj), collapse=", "), "\n")
cat("meta.data columns:", paste(colnames(obj@meta.data), collapse=", "), "\n")
cat("meta.data rows:", nrow(obj@meta.data), "\n")
cat("features (RNA):", nrow(obj[["RNA"]]), "\n")
cat("layers (RNA):", paste(Layers(obj[["RNA"]]), collapse=", "), "\n")
cat("\n--- misc slot names ---\n")
print(names(obj@misc))
cat("\n--- searching for TCR/clonotype/CDR3 ---\n")
for (nm in c("misc", "tools")) {
  s <- capture.output(str(slot(obj, nm), max.level = 1))
  hit <- grep("TCR|clonotype|CDR3|tcr", s, ignore.case = TRUE, value = TRUE)
  cat(nm, "hits:", length(hit), "\n"); print(head(hit, 10))
}
cat("\n--- meta.data head ---\n")
print(head(obj@meta.data))
cat("\n--- tier1 counts ---\n"); print(table(obj@meta.data$tier1, useNA = "ifany"))
cat("\n--- tier2 counts ---\n"); print(table(obj@meta.data$tier2, useNA = "ifany"))
cat("\n--- tier1_T counts ---\n"); print(table(obj@meta.data$tier1_T, useNA = "ifany"))
cat("\n--- timepoint x orig.ident ---\n"); print(table(obj@meta.data$timepoint, obj@meta.data$orig.ident))
sink()
cat("structure written\n")

# meta.data
cat("saving meta.data ...\n")
write.csv(obj@meta.data, gzfile("meylan_meta.csv.gz"), row.names = TRUE)

# centroids from images slot (Seurat v5 FOV objects)
cat("extracting centroids ...\n")
fovs <- Images(obj)
cent_list <- list()
for (fv in fovs) {
  fov_obj <- obj[[fv]]
  # FOV class: centroids stored in the `cells` slot or via GetTissueCoordinates
  coords <- tryCatch({
    tc <- GetTissueCoordinates(fov_obj)
    as.data.frame(tc)
  }, error = function(e) NULL)
  if (is.null(coords)) {
    # fallback: try the cells slot directly
    coords <- tryCatch(as.data.frame(fov_obj@cells), error = function(e) NULL)
  }
  if (!is.null(coords) && nrow(coords) > 0) {
    coords$cell <- rownames(coords)
    coords$FOV <- fv
    cent_list[[fv]] <- coords
    cat(fv, ":", nrow(coords), "cells; cols:", paste(colnames(coords), collapse=","), "\n")
  } else {
    cat(fv, ": NO coordinates extracted\n")
  }
}
if (length(cent_list) > 0) {
  cent <- do.call(rbind, cent_list)
  write.csv(cent, gzfile("meylan_centroids.csv.gz"), row.names = FALSE)
  cat("saved centroids:", nrow(cent), "rows\n")
}

# counts subset for genes of interest
genes <- c("VEGFA","VCAN","ANXA1","CD3D","CD3E","CD8A","CD4","GZMB","PDCD1",
           "CD68","CD163","MS4A1","NCAM1","KLRD1","FOXP3","IFNG","CXCL10",
           "CXCL9","CCL5","HLA-DRA","CD69","NR4A1","MKI67","TOP2A",
           "HSV1-gB","HSV1-gD","UL27","UL44")
genes <- intersect(genes, rownames(obj[["RNA"]]))
cat("genes found:", length(genes), "of requested\n")
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")[genes, , drop = FALSE]
saveRDS(counts, "meylan_genes.rds")
cat("saved meylan_genes.rds:", nrow(counts), "x", ncol(counts), "\n")
cat("DONE\n")
