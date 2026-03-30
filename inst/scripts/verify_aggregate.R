.libPaths(c("C:/Users/win10/R/win-library/4.4", .libPaths()))
devtools::load_all("C:/Users/win10/SpatialDataR", quiet = TRUE)

cat("=== Loading Xenium ===\n")
sd <- readSpatialData("C:/Users/win10/xenium_breast/data.zarr")

pts <- spatialPoints(sd)[["transcripts"]]
cat("Points class:", class(pts)[1], "\n")
cat("Points cols:", paste(names(pts), collapse = ", "), "\n")
cat("Points rows:", nrow(pts), "\n")
cat("cell_id range:", range(pts$cell_id), "\n")
cat("cell_id==0 (unassigned):", sum(pts$cell_id == 0), "\n\n")

cells <- shapes(sd)[["cell_circles"]]
cat("Cells class:", class(cells)[1], "\n")
cat("Cells cols:", paste(names(cells), collapse = ", "), "\n")
cat("Cells rows:", nrow(cells), "\n")
cat("cell_id range:", range(cells$cell_id), "\n\n")

## The feature column is "feature_name" not "gene"
cat("Feature column check: 'gene' in pts?", "gene" %in% names(pts), "\n")
cat("Feature column check: 'feature_name' in pts?", "feature_name" %in% names(pts), "\n\n")

## Try aggregatePoints with correct column names
cat("=== Running aggregatePoints() ===\n")
## Filter to QV >= 20 and no controls
pts_filt <- pts[pts$qv >= 20 & !grepl("^NegControl|^BLANK", as.character(pts$feature_name)), ]
## Remove unassigned (cell_id == 0)
pts_filt <- pts_filt[pts_filt$cell_id > 0, ]
cat("Filtered points:", format(nrow(pts_filt), big.mark = ","), "\n")

t0 <- proc.time()
counts <- aggregatePoints(
    pts_filt,
    cells,
    feature_col = "feature_name",
    region_col = "cell_id")
elapsed <- (proc.time() - t0)[3]

cat("\n=== Result ===\n")
cat("Class:", class(counts)[1], "\n")
cat("Dimensions:", nrow(counts), "x", ncol(counts), "\n")
cat("Time:", round(elapsed, 1), "seconds\n")
cat("First 3 rows, first 5 cols:\n")
print(head(as.data.frame(counts[, 1:min(5, ncol(counts))]), 3))
cat("\nColumn names (first 10):", paste(head(names(counts), 10), collapse = ", "), "\n")
cat("Total counts per cell (first 5):", head(rowSums(as.matrix(counts[, -which(names(counts) == "cell_id")])), 5), "\n")
