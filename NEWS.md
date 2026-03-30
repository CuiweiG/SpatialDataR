# SpatialDataR 0.99.0

## New features

* `SpatialData` S4 class with validity checking for
  multi-modal spatial omics data, following the SpatialData
  specification (Marconato et al. 2025 *Nat Methods*).

### Store I/O
* `readSpatialData()` reads SpatialData Zarr stores
  natively without Python dependencies.
* `readZarrArray()` for Zarr arrays via Rarr/pizzarr.
* `readParquetPoints()` for Parquet transcript tables.
* `readCSVElement()` for CSV-backed points and shapes.
* `readSpatialTable()` converts AnnData tables to
  SpatialExperiment with CSV fallback.

### Coordinate transforms (OME-NGFF compliant)
* `CoordinateTransform` S4 class: 2D (3x3) and 3D (4x4)
  affine matrices with validity checking.
* `transformCoords()` methods for `DataFrame` and `matrix`.
* `composeTransforms()` chains transform sequences.
* `invertTransform()` computes inverse transforms.
* `.parseTransform()` supports identity, affine, scale,
  translation, and sequence types from OME-NGFF spec.
* `elementTransform()` extracts transforms from metadata.

### Spatial operations
* `bboxQuery()` for bounding box spatial queries on
  `DataFrame` and `SpatialData` objects (mirrors Python
  `spatialdata.bounding_box_query()`).
* `aggregatePoints()` creates cell-by-gene count matrices
  from molecule coordinates + region assignments (mirrors
  Python `spatialdata.aggregate()`).

### Multi-sample
* `combineSpatialData()` merges multiple stores with
  automatic name disambiguation and CS prefixing.
* `filterSample()` extracts elements for a single sample.

### Validation
* `validateSpatialData()` checks Zarr stores against the
  SpatialData on-disk specification. Reports structural
  issues, missing metadata, and transform parseability.

### Interoperability
* `elementSummary()` tabulates all elements.
* `coordinateSystemElements()` maps CS to elements.
* `length()`, `names()`, `[` methods for `SpatialData`.
* Accessors: `images()`, `spatialLabels()`,
  `spatialPoints()`, `shapes()`, `tables()`,
  `coordinateSystems()`.

### Data
* Mini Zarr store (`xenium_mini.zarr`) in `inst/extdata`.
* Validation scripts for bundled and public datasets in
  `inst/scripts/`.
