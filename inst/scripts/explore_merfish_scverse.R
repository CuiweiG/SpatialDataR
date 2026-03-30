.libPaths(c("C:/Users/win10/R/win-library/4.4", .libPaths()))
devtools::load_all("C:/Users/win10/SpatialDataR", quiet = TRUE)

cat("=== Reading MERFISH scverse ===\n")
store <- "C:/Users/win10/merfish_scverse/data.zarr"

tryCatch({
    sd <- readSpatialData(store)
    print(sd)
}, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
})

cat("\n=== Manual check ===\n")
cat("Points dir:", list.dirs(file.path(store, "points"), recursive=FALSE, full.names=FALSE), "\n")
cat("Shapes dir:", list.dirs(file.path(store, "shapes"), recursive=FALSE, full.names=FALSE), "\n")

# Check points structure
pts_dir <- file.path(store, "points/single_molecule")
cat("Points contents:", list.files(pts_dir, recursive=FALSE), "\n")

# Try reading points
tryCatch({
    pts <- readParquetPoints(file.path(pts_dir, "points.parquet"))
    cat("Points rows:", nrow(pts), "\n")
    cat("Points cols:", paste(names(pts), collapse=", "), "\n")
    cat("First 3:\n")
    print(head(as.data.frame(pts), 3))
}, error = function(e) {
    cat("Points error:", conditionMessage(e), "\n")
})

# Check shapes
for (nm in c("cells", "anatomical")) {
    shp_dir <- file.path(store, "shapes", nm)
    if (dir.exists(shp_dir)) {
        cat("\nShape", nm, "contents:", list.files(shp_dir, recursive=TRUE), "\n")
        tryCatch({
            df <- as.data.frame(arrow::read_parquet(
                file.path(shp_dir, "shapes.parquet")))
            cat("  Rows:", nrow(df), " Cols:", paste(names(df), collapse=", "), "\n")
        }, error = function(e) cat("  Error:", e$message, "\n"))
    }
}

# Check tables
tbl_dir <- file.path(store, "tables/table")
cat("\nTable contents:", list.files(tbl_dir, recursive=FALSE), "\n")
# obs
obs_dir <- file.path(tbl_dir, "obs")
cat("obs zarr.json exists:", file.exists(file.path(obs_dir, "zarr.json")), "\n")
cat("obs .zattrs exists:", file.exists(file.path(obs_dir, ".zattrs")), "\n")
