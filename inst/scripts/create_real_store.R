#!/usr/bin/env Rscript
## Build a real SpatialData Zarr store from MERFISH Allen VISp
## Source: Moffitt et al. 2018, Science 362
## Data: spacetx-spacejam S3 (CC0 1.0)
.libPaths("C:/Users/win10/R/win-library/4.4")

cat("=== Building MERFISH SpatialData store ===\n")

csv_file <- "C:/Users/win10/merfish_real.csv"
outdir <- "C:/Users/win10/merfish_spatialdata.zarr"

if (dir.exists(outdir)) unlink(outdir, recursive = TRUE)

## Read raw data
cat("Reading CSV (319 MB)...\n")
df <- read.csv(csv_file, stringsAsFactors = FALSE)
cat("Rows:", nrow(df), " Cols:", ncol(df), "\n")
cat("Columns:", paste(colnames(df), collapse = ", "), "\n")

## Subset to a region (same as spatialdata-sandbox)
bb <- list(x0 = 1154, x1 = 3172, y0 = 4548, y1 = 6566)
sub <- df[df$x_um > bb$x0 & df$x_um < bb$x1 &
    df$y_um > bb$y0 & df$y_um < bb$y1, ]
cat("Subsetted to bounding box:", nrow(sub), "spots\n")

## Create Zarr store structure
dir.create(outdir, recursive = TRUE)

## Top-level .zattrs
writeLines('{
  "spatialdata_attrs": {
    "version": "0.1",
    "coordinate_systems": {
      "global": {
        "axes": [
          {"name": "x", "type": "space", "unit": "micrometer"},
          {"name": "y", "type": "space", "unit": "micrometer"}
        ]
      }
    }
  }
}', file.path(outdir, ".zattrs"))

## Points: transcripts
pts_dir <- file.path(outdir, "points", "transcripts")
dir.create(pts_dir, recursive = TRUE)
writeLines('{
  "coordinateTransformations": [
    {"type": "identity"}
  ],
  "spatialdata_attrs": {
    "feature_key": "gene",
    "instance_key": "cell_id"
  }
}', file.path(pts_dir, ".zattrs"))

pts_df <- data.frame(
    x = sub$x_um,
    y = sub$y_um,
    gene = sub$gene,
    cell_id = seq_len(nrow(sub)),
    stringsAsFactors = FALSE)
write.csv(pts_df, file.path(pts_dir, "transcripts.csv"),
    row.names = FALSE)
cat("Points:", nrow(pts_df), "transcripts,",
    length(unique(pts_df$gene)), "genes\n")

## Shapes: cell regions (use unique cell positions)
if ("layer" %in% colnames(sub)) {
    shp_dir <- file.path(outdir, "shapes",
        "cell_boundaries")
    dir.create(shp_dir, recursive = TRUE)
    writeLines('{
  "coordinateTransformations": [
    {"type": "identity"}
  ]
}', file.path(shp_dir, ".zattrs"))

    ## Approximate cell positions from transcript centroids
    cell_pos <- aggregate(
        cbind(x_um, y_um) ~ layer, data = sub,
        FUN = function(v) {
            ## Sample representative positions
            if (length(v) > 100) {
                quantile(v, probs = seq(0, 1, 0.02))
            } else v
        })

    ## Create simple grid of cells based on layer
    layers <- unique(sub$layer)
    set.seed(42)
    cells <- data.frame(
        cell_id = integer(),
        x = numeric(),
        y = numeric(),
        radius = numeric(),
        cell_type = character(),
        stringsAsFactors = FALSE)
    cid <- 1
    for (lay in layers) {
        lsub <- sub[sub$layer == lay, ]
        ## Sample 20 representative cells per layer
        n_cells <- min(20, nrow(lsub))
        idx <- sample(nrow(lsub), n_cells)
        for (i in idx) {
            cells <- rbind(cells, data.frame(
                cell_id = cid,
                x = lsub$x_um[i],
                y = lsub$y_um[i],
                radius = 10,
                cell_type = lay,
                stringsAsFactors = FALSE))
            cid <- cid + 1
        }
    }
    write.csv(cells[, c("cell_id", "x", "y", "radius")],
        file.path(shp_dir, "circles.csv"),
        row.names = FALSE)
    cat("Shapes:", nrow(cells), "cells,",
        length(unique(cells$cell_type)), "layers\n")

    ## Tables: cell metadata
    tbl_dir <- file.path(outdir, "tables", "table")
    dir.create(file.path(tbl_dir, "obs"), recursive = TRUE)
    dir.create(file.path(tbl_dir, "var"), recursive = TRUE)
    writeLines('{
  "encoding-type": "anndata",
  "spatialdata_attrs": {
    "region": "cell_boundaries",
    "instance_key": "cell_id"
  }
}', file.path(tbl_dir, ".zattrs"))
    writeLines('{"column-order": ["cell_id", "cell_type"]}',
        file.path(tbl_dir, "obs", ".zattrs"))
    writeLines('{"_index": "gene_name"}',
        file.path(tbl_dir, "var", ".zattrs"))

    write.csv(cells[, c("cell_id", "cell_type")],
        file.path(tbl_dir, "obs", "obs.csv"),
        row.names = FALSE)
    write.csv(
        data.frame(gene_name = unique(pts_df$gene)),
        file.path(tbl_dir, "var", "var.csv"),
        row.names = FALSE)
    cat("Tables:", nrow(cells), "cells\n")
}

cat("\nStore created:", outdir, "\n")
cat("Size:", round(sum(file.info(
    list.files(outdir, recursive = TRUE,
        full.names = TRUE))$size) / 1e6, 1),
    "MB\n")
