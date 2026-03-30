.libPaths(c("C:/Users/win10/R/win-library/4.4", .libPaths()))
devtools::load_all("C:/Users/win10/SpatialDataR", quiet = TRUE)

cat("=== Bioconductor Bridge ===\n")
sd <- readSpatialData("C:/Users/win10/xenium_breast/data.zarr")

t0 <- proc.time()
sce <- toSingleCellExperiment(sd)
cat("toSingleCellExperiment():", round((proc.time()-t0)[3],1), "s\n")
cat("Class:", class(sce)[1], "\n")
cat("Dims:", nrow(sce), "genes x", ncol(sce), "cells\n")
cat("First 5 genes:", paste(head(rownames(sce), 5), collapse=", "), "\n")
cat("colData cols:", paste(names(SummarizedExperiment::colData(sce)), collapse=", "), "\n")

cat("\n=== Moran's I ===\n")
coords <- geometryCentroids(shapes(sd)[["cell_circles"]][["geometry"]])
set.seed(42); idx <- sample(ncol(sce), 2000)
sub_expr <- t(as.matrix(SummarizedExperiment::assay(sce[, idx], "counts")))
sub_coords <- coords[idx, ]

t0 <- proc.time()
morans <- spatialAutocorrelation(sub_expr, sub_coords)
cat("spatialAutocorrelation():", round((proc.time()-t0)[3],1), "s\n")
cat("Genes tested:", nrow(morans), "\n")

## Top spatially variable genes
morans_sig <- morans[order(morans$p.value), ]
cat("\nTop 10 spatially variable genes:\n")
print(head(morans_sig[, c("gene","observed","p.value","adjusted.p")], 10))

## Check known markers
for (g in c("EPCAM","LUM","PTPRC","PECAM1")) {
    row <- morans[morans$gene == g, ]
    if (nrow(row) > 0)
        cat(g, ": I=", round(row$observed, 3), " p=", format(row$p.value, digits=3), "\n")
}

cat("\n=== spatstat ===\n")
pts <- spatialPoints(sd)[["transcripts"]]
set.seed(42); sub_pts <- pts[sample(nrow(pts), 50000), ]
pp <- toPointPattern(sub_pts, marks = "feature_name")
cat("ppp class:", class(pp)[1], "\n")
cat("ppp npoints:", pp$n, "\n")
cat("ppp marks levels:", length(levels(pp$marks)), "\n")

cat("\n=== spatialDiffExpression ===\n")
## Epithelial vs Stromal
epcam <- sub_expr[, "EPCAM"]
epithelial <- epcam > median(epcam)
de <- spatialDiffExpression(sub_expr, epithelial, !epithelial)
de_sig <- de[order(de$p.value), ]
cat("Top 5 DE genes (epithelial vs rest):\n")
print(head(de_sig[, c("gene","log2FC","p.value","adjusted.p")], 5))
