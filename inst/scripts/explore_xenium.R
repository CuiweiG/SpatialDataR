## Explore Xenium breast cancer dataset structure
.libPaths(c("C:/Users/win10/R/win-library/4.4",
            "C:/Users/win10/AppData/Local/R/win-library/4.4",
            .libPaths()))

library(SpatialDataR)
library(arrow)

store <- "C:/Users/win10/xenium_breast/data.zarr"

cat("=== Attempting readSpatialData() ===\n")
tryCatch({
    sd <- readSpatialData(store)
    cat("Success!\n")
    print(sd)
}, error = function(e) {
    cat("Error:", conditionMessage(e), "\n\n")
    cat("=== Manual exploration ===\n")
})

## Manual exploration of zarr structure
cat("\n=== zarr.json ===\n")
zj <- file.path(store, "zarr.json")
if (file.exists(zj)) cat(readLines(zj, warn=FALSE), sep="\n")

## Check points structure
cat("\n=== Points structure ===\n")
pts_dir <- file.path(store, "points/transcripts")
cat("Contents:", paste(list.files(pts_dir, recursive=FALSE), collapse=", "), "\n")

## Check if it's parquet
pq_files <- list.files(pts_dir, pattern="\\.parquet$", recursive=TRUE)
cat("Parquet files:", paste(pq_files, collapse=", "), "\n")

## Try reading parquet directly
if (length(pq_files) > 0) {
    pq_path <- file.path(pts_dir, pq_files[1])
    cat("\nReading first parquet file:", pq_path, "\n")
    tryCatch({
        df <- arrow::read_parquet(pq_path)
        cat("Rows:", nrow(df), "\n")
        cat("Cols:", paste(names(df), collapse=", "), "\n")
        cat("First few rows:\n")
        print(head(df, 3))
    }, error = function(e) cat("Error reading parquet:", e$message, "\n"))
}

## Check shapes
cat("\n=== Shapes structure ===\n")
shp_dir <- file.path(store, "shapes/cell_circles")
if (dir.exists(shp_dir)) {
    cat("cell_circles contents:", paste(list.files(shp_dir, recursive=TRUE), collapse=", "), "\n")
    pq_s <- list.files(shp_dir, pattern="\\.parquet$", recursive=TRUE)
    if (length(pq_s) > 0) {
        tryCatch({
            df_s <- arrow::read_parquet(file.path(shp_dir, pq_s[1]))
            cat("Shape rows:", nrow(df_s), "\n")
            cat("Shape cols:", paste(names(df_s), collapse=", "), "\n")
            print(head(df_s, 3))
        }, error = function(e) cat("Error:", e$message, "\n"))
    }
}

## Check tables
cat("\n=== Tables structure ===\n")
tbl_dir <- file.path(store, "tables/table")
if (dir.exists(tbl_dir)) {
    cat("Table contents:", paste(list.files(tbl_dir, recursive=FALSE), collapse=", "), "\n")
    obs_dir <- file.path(tbl_dir, "obs")
    if (dir.exists(obs_dir)) {
        cat("obs contents:", paste(list.files(obs_dir), collapse=", "), "\n")
    }
    var_dir <- file.path(tbl_dir, "var")
    if (dir.exists(var_dir)) {
        cat("var contents:", paste(list.files(var_dir), collapse=", "), "\n")
    }
    # Check X matrix
    x_dir <- file.path(tbl_dir, "X")
    if (dir.exists(x_dir)) {
        cat("X contents:", paste(list.files(x_dir), collapse=", "), "\n")
    }
}

## Check images
cat("\n=== Images structure ===\n")
img_dir <- file.path(store, "images/morphology_focus")
if (dir.exists(img_dir)) {
    cat("morphology_focus contents:", paste(list.files(img_dir, recursive=FALSE), collapse=", "), "\n")
    ## Check for multiscale
    for (sc in list.files(img_dir)) {
        sc_path <- file.path(img_dir, sc)
        if (dir.exists(sc_path)) {
            zattrs <- file.path(sc_path, ".zarray")
            if (file.exists(zattrs)) {
                cat("  Scale", sc, ":", readLines(zattrs, warn=FALSE)[1:3], "\n")
            }
        }
    }
}
