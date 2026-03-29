# SpatialDataR 0.99.0

## New features

* `SpatialData` S4 class with validity checking for multi-modal
  spatial omics data.
* `readSpatialData()` reads SpatialData Zarr stores natively
  without Python dependencies. Images and labels stored as path
  references; points and shapes eagerly loaded from CSV/Parquet.
* `readZarrArray()` for array access via Rarr/pizzarr backends.
* `readParquetPoints()` for Parquet-backed transcript coordinates.
* `readCSVElement()` for CSV-backed points and shapes.
* `readSpatialTable()` converts AnnData-format tables to
  SpatialExperiment, with CSV fallback for obs/var.
* `CoordinateTransform` class with affine transformation support
  and validity checking.
* `transformCoords()` generic with methods for `DataFrame` and
  `matrix` inputs.
* Accessors: `images()`, `spatialLabels()`, `spatialPoints()`,
  `shapes()`, `tables()`, `coordinateSystems()`.
* Mini Zarr store (`xenium_mini.zarr`) shipped in `inst/extdata`
  for testing and examples.
