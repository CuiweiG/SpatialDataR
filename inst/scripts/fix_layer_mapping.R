.libPaths(c("C:/Users/win10/R/win-library/4.4", .libPaths()))
devtools::load_all("C:/Users/win10/SpatialDataR", quiet = TRUE)

sd2 <- readSpatialData("C:/Users/win10/merfish_scverse/data.zarr")
pts <- spatialPoints(sd2)[["single_molecule"]]
anat <- shapes(sd2)[["anatomical"]]

cat("=== Strategy: use transcript cell_type annotations ===\n")
cat("Each transcript has a cell_type label (VISp_I, VISp_II/III, etc.)\n")
cat("For each polygon, count which cell_type is most common inside it.\n\n")

## Step 1: spatialJoin to assign each transcript to a polygon
## (not to a layer name — just to polygon index 1-6)
t0 <- proc.time()
poly_assignment <- spatialJoin(pts, anat,
    region_names = paste0("poly", 1:6))
elapsed <- (proc.time() - t0)[3]
cat("spatialJoin:", round(elapsed, 1), "s\n")

## Step 2: for each polygon, find dominant cell_type from transcripts
pts_df <- as.data.frame(pts)
pts_df$polygon <- poly_assignment

for (pid in 1:6) {
    in_poly <- pts_df[pts_df$polygon == paste0("poly", pid) & !is.na(pts_df$polygon), ]
    cat("\nPolygon", pid, ":", nrow(in_poly), "transcripts\n")
    if (nrow(in_poly) > 0) {
        ct_table <- sort(table(in_poly$cell_type), decreasing = TRUE)
        cat("  Top cell_types:\n")
        print(head(ct_table, 3))
        cat("  Dominant:", names(ct_table)[1], "\n")
    }
}
