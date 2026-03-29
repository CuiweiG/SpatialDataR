# SpatialDataR 0.99.0

## New features

* `SpatialData` S4 class with validity checking for multi-modal
  spatial omics data, following the SpatialData specification
  (Marconato et al. 2024 *Nat Methods*).
* `readSpatialData()` reads SpatialData Zarr stores natively
  without Python dependencies. Points and shapes eagerly loaded
  from CSV/Parquet; images and labels stored as path references.
* `readZarrArray()` for array access via Rarr/pizzarr backends.
* `readParquetPoints()` for Parquet-backed transcript coordinates.
* `readCSVElement()` for CSV-backed points and shapes.
* `readSpatialTable()` converts AnnData-format tables to
  SpatialExperiment, with CSV fallback for obs/var.
* `CoordinateTransform` class with full affine transformation
  support:
    - 2D (3x3) and 3D (4x4) affine matrices
    - `composeTransforms()` for chaining transform sequences
    - `invertTransform()` for computing inverse transforms
    - `transformCoords()` methods for `DataFrame` and `matrix`
    - `.parseTransform()` supports identity, affine, scale,
      translation, and sequence types from OME-NGFF spec
* `bboxQuery()` for bounding box spatial queries on `DataFrame`
  and `SpatialData` objects (mirrors Python
  `spatialdata.bounding_box_query()`).
* `elementSummary()` for listing all elements with type, name,
  and row counts.
* `elementTransform()` for extracting coordinate transforms from
  element metadata.
* `coordinateSystemElements()` for mapping coordinate systems to
  their associated elements.
* `[` operator for subsetting `SpatialData` by element name.
* Accessors: `images()`, `spatialLabels()`, `spatialPoints()`,
  `shapes()`, `tables()`, `coordinateSystems()`.
* Mini Zarr store (`xenium_mini.zarr`) shipped in `inst/extdata`
  for testing and examples.
