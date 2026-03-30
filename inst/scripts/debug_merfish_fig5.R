.libPaths(c("C:/Users/win10/R/win-library/4.4", .libPaths()))
devtools::load_all("C:/Users/win10/SpatialDataR", quiet = TRUE)

sd <- readSpatialData("C:/Users/win10/merfish_scverse/data.zarr")
pts <- as.data.frame(spatialPoints(sd)[["single_molecule"]])

cat("Columns:", paste(names(pts), collapse = ", "), "\n")
cat("Rows:", nrow(pts), "\n")
cat("\ncell_type distribution:\n")
print(table(pts$cell_type))
cat("\nx range:", range(pts$x), "\n")
cat("y range:", range(pts$y), "\n")

## The key question: does this dataset have gene names?
## cell_type column contains REGION labels, not gene names
## So this is a region-annotated molecule dataset, not a gene-expression one
cat("\nFirst 10 rows:\n")
print(head(pts, 10))
