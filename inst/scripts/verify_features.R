.libPaths(c("C:/Users/win10/R/win-library/4.4", .libPaths()))
devtools::load_all("C:/Users/win10/SpatialDataR", quiet = TRUE)

cat("=== parseGeometry() ===\n")
sd <- readSpatialData("C:/Users/win10/xenium_breast/data.zarr")
cells <- shapes(sd)[["cell_circles"]]
coords <- geometryCentroids(cells[["geometry"]])
cat("Xenium centroids:", nrow(coords), "x", ncol(coords), "\n")
cat("First 3:\n"); print(head(coords, 3))

cat("\n=== spatialJoin() ===\n")
sd2 <- readSpatialData("C:/Users/win10/merfish_scverse/data.zarr")
pts <- spatialPoints(sd2)[["single_molecule"]]
anat <- shapes(sd2)[["anatomical"]]
t0 <- proc.time()
assignments <- spatialJoin(pts, anat,
    region_names = c("VISp_I","VISp_II/III","VISp_IV","VISp_V","VISp_VI","VISp_wm"))
elapsed <- (proc.time() - t0)[3]
cat("Time:", round(elapsed, 1), "seconds\n")
cat("Assignments:\n"); print(table(assignments, useNA = "always"))

cat("\n=== plotSpatialData() ===\n")
cat("Function exists:", exists("plotSpatialData"), "\n")
