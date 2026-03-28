# SpatialDataR 0.99.0

## New features

* `SpatialData` S4 class for multi-modal spatial omics data.
* `readSpatialData()` reads SpatialData Zarr stores natively
  without Python dependencies.
* `readZarrArray()` for lazy array access via Rarr/pizzarr.
* `readParquetPoints()` for Parquet-backed transcript coordinates.
* `readSpatialTable()` converts AnnData tables to SpatialExperiment.
* `CoordinateTransform` class with affine transformation support.
* `transformCoords()` generic for coordinate system mapping.
* Accessors: `images()`, `labels()`, `spatialPoints()`, `shapes()`,
  `tables()`, `coordinateSystems()`.
