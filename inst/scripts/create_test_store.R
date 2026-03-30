#!/usr/bin/env Rscript
## Create a SpatialData-spec-compliant Zarr store for validation
## Structure follows Marconato et al. 2025 Nat Methods exactly
##
## This creates a small but structurally complete .zarr store
## mimicking a 10x Xenium breast cancer dataset.

outdir <- "inst/extdata/xenium_mini.zarr"
if (dir.exists(outdir)) unlink(outdir, recursive = TRUE)

cat("Creating SpatialData Zarr store:", outdir, "\n\n")

## ============================================================
## Top-level .zattrs (SpatialData metadata)
## ============================================================
dir.create(outdir, recursive = TRUE)

zattrs <- '{
  "spatialdata_attrs": {
    "version": "0.1",
    "coordinate_systems": {
      "global": {
        "axes": [
          {"name": "x", "type": "space", "unit": "micrometer"},
          {"name": "y", "type": "space", "unit": "micrometer"}
        ]
      },
      "pixels": {
        "axes": [
          {"name": "x", "type": "space", "unit": "pixel"},
          {"name": "y", "type": "space", "unit": "pixel"}
        ]
      }
    }
  }
}'
writeLines(zattrs, file.path(outdir, ".zattrs"))

## ============================================================
## images/morphology â€?small 3-channel image (20x20x3)
## ============================================================
img_dir <- file.path(outdir, "images", "morphology", "scale0")
dir.create(img_dir, recursive = TRUE)

## Write image metadata (OME-Zarr multiscale spec)
img_meta <- '{
  "coordinateTransformations": [
    {"type": "scale", "scale": [0.2125, 0.2125]}
  ],
  "multiscales": [{
    "axes": [
      {"name": "c", "type": "channel"},
      {"name": "y", "type": "space"},
      {"name": "x", "type": "space"}
    ],
    "datasets": [{"path": "scale0"}]
  }]
}'
writeLines(img_meta,
    file.path(outdir, "images", "morphology", ".zattrs"))

## Write Zarr array metadata
zarr_meta <- '{
  "zarr_format": 2,
  "shape": [3, 20, 20],
  "chunks": [3, 20, 20],
  "dtype": "<u1",
  "compressor": null,
  "fill_value": 0,
  "order": "C",
  "filters": null
}'
writeLines(zarr_meta, file.path(img_dir, ".zarray"))

## Write actual image data (20x20 RGB, 3 channels)
set.seed(2024)
img_data <- as.raw(sample(0:255, 3 * 20 * 20, replace = TRUE))
writeBin(img_data, file.path(img_dir, "0.0.0"))
cat("  images/morphology: 20x20x3 RGB\n")

## ============================================================
## labels/cell_labels â€?segmentation mask (20x20)
## ============================================================
lbl_dir <- file.path(outdir, "labels", "cell_labels", "scale0")
dir.create(lbl_dir, recursive = TRUE)

lbl_meta <- '{
  "coordinateTransformations": [
    {"type": "scale", "scale": [0.2125, 0.2125]}
  ]
}'
writeLines(lbl_meta,
    file.path(outdir, "labels", "cell_labels", ".zattrs"))

lbl_zarr <- '{
  "zarr_format": 2,
  "shape": [20, 20],
  "chunks": [20, 20],
  "dtype": "<i4",
  "compressor": null,
  "fill_value": 0,
  "order": "C",
  "filters": null
}'
writeLines(lbl_zarr, file.path(lbl_dir, ".zarray"))

## Cell IDs: 50 cells randomly placed
labels <- matrix(0L, 20, 20)
set.seed(42)
for (i in 1:50) {
    x <- sample(1:20, 1); y <- sample(1:20, 1)
    labels[x, y] <- as.integer(i)
}
writeBin(as.raw(as.vector(t(labels))), file.path(lbl_dir, "0.0"))
cat("  labels/cell_labels: 20x20 mask (50 cells)\n")

## ============================================================
## points/transcripts â€?transcript coordinates
## ============================================================
pts_dir <- file.path(outdir, "points", "transcripts")
dir.create(pts_dir, recursive = TRUE)

pts_meta <- '{
  "coordinateTransformations": [
    {"type": "identity"}
  ],
  "spatialdata_attrs": {
    "feature_key": "gene",
    "instance_key": "cell_id"
  }
}'
writeLines(pts_meta, file.path(pts_dir, ".zattrs"))

## Create transcript table as CSV (simulating Parquet structure)
set.seed(2024)
n_tx <- 500
genes <- sample(c("EPCAM", "KRT18", "VIM", "CD45", "HER2",
                   "ESR1", "PGR", "MKI67", "ERBB2", "ACTB"),
                n_tx, replace = TRUE)
tx_df <- data.frame(
    x = runif(n_tx, 0, 20) * 0.2125,
    y = runif(n_tx, 0, 20) * 0.2125,
    gene = genes,
    cell_id = sample(0:50, n_tx, replace = TRUE),
    stringsAsFactors = FALSE
)
write.csv(tx_df, file.path(pts_dir, "transcripts.csv"),
          row.names = FALSE)
cat("  points/transcripts:", n_tx, "transcripts, 10 genes\n")

## ============================================================
## shapes/cell_boundaries â€?cell polygons
## ============================================================
shp_dir <- file.path(outdir, "shapes", "cell_boundaries")
dir.create(shp_dir, recursive = TRUE)

shp_meta <- '{
  "coordinateTransformations": [
    {"type": "identity"}
  ]
}'
writeLines(shp_meta, file.path(shp_dir, ".zattrs"))

## Simple circle shapes for 50 cells
circles <- data.frame(
    cell_id = 1:50,
    x = runif(50, 0, 20) * 0.2125,
    y = runif(50, 0, 20) * 0.2125,
    radius = runif(50, 0.1, 0.5),
    stringsAsFactors = FALSE
)
write.csv(circles, file.path(shp_dir, "circles.csv"),
          row.names = FALSE)
cat("  shapes/cell_boundaries: 50 cell circles\n")

## ============================================================
## tables/table â€?gene expression (AnnData-like structure)
## ============================================================
tbl_dir <- file.path(outdir, "tables", "table")
dir.create(tbl_dir, recursive = TRUE)

tbl_meta <- '{
  "encoding-type": "anndata",
  "encoding-version": "0.1.0",
  "spatialdata_attrs": {
    "region": "cell_labels",
    "region_key": "region",
    "instance_key": "cell_id"
  }
}'
writeLines(tbl_meta, file.path(tbl_dir, ".zattrs"))

## obs (cell metadata)
obs_dir <- file.path(tbl_dir, "obs")
dir.create(obs_dir, recursive = TRUE)
obs_meta <- '{"column-order": ["cell_id", "cell_type", "region"],
  "_index": "cell_id"}'
writeLines(obs_meta, file.path(obs_dir, ".zattrs"))

cell_types <- sample(c("Epithelial", "Stromal", "Immune",
                        "Endothelial"), 50, replace = TRUE,
                     prob = c(0.4, 0.3, 0.2, 0.1))
obs_df <- data.frame(
    cell_id = 1:50,
    cell_type = cell_types,
    region = rep("cell_labels", 50),
    stringsAsFactors = FALSE
)
write.csv(obs_df, file.path(obs_dir, "obs.csv"), row.names = FALSE)

## var (gene metadata)
var_dir <- file.path(tbl_dir, "var")
dir.create(var_dir, recursive = TRUE)
var_meta <- '{"_index": "gene_name"}'
writeLines(var_meta, file.path(var_dir, ".zattrs"))

gene_names <- c("EPCAM", "KRT18", "VIM", "CD45", "HER2",
                "ESR1", "PGR", "MKI67", "ERBB2", "ACTB")
write.csv(data.frame(gene_name = gene_names),
          file.path(var_dir, "var.csv"), row.names = FALSE)

cat("  tables/table: 50 cells x 10 genes\n")

## ============================================================
## Summary
## ============================================================
cat("\n=== Zarr Store Structure ===\n")
all_files <- list.files(outdir, recursive = TRUE)
cat(paste("  ", all_files, collapse = "\n"), "\n")
cat("\nTotal files:", length(all_files), "\n")
cat("Store path:", normalizePath(outdir), "\n")
