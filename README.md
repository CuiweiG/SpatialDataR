<div align="center">

# SpatialDataR

*Native R Infrastructure for Reading SpatialData Zarr Stores*

[![R-CMD-check](https://github.com/CuiweiG/SpatialDataR/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/CuiweiG/SpatialDataR/actions/workflows/R-CMD-check.yml)
[![License: Artistic-2.0](https://img.shields.io/badge/License-Artistic--2.0-blue.svg)](https://opensource.org/licenses/Artistic-2.0)
[![BioC status](https://img.shields.io/badge/BioC-0.99.0-orange.svg)](https://bioconductor.org/)

</div>

---

## The problem

[SpatialData](https://spatialdata.scverse.org/) (Marconato et al.
2024 *Nat Methods*) defines a universal Zarr-based format for
spatial omics. R-side support has remained limited, with no mature
native Bioconductor-oriented reader providing direct, zero-Python
access to SpatialData stores.

## What SpatialDataR provides

`SpatialDataR` reads SpatialData `.zarr` stores directly in R:

- **Store discovery**: `readSpatialData()` parses `.zattrs`
  metadata and discovers images, labels, points, shapes, tables
- **Element readers**: `readZarrArray()` (via Rarr/pizzarr),
  `readParquetPoints()` (via arrow), `readCSVElement()`,
  `readSpatialTable()`
- **Coordinate transforms**: `CoordinateTransform` class with
  affine transformation support for `DataFrame` and `matrix`
- **Eager loading**: points and shapes stored as CSV/Parquet are
  loaded as `DataFrame` objects on store read; images and labels
  are stored as path references for on-demand loading via
  `readZarrArray()`

```r
library(SpatialDataR)
sd <- readSpatialData("xenium_breast.zarr")
sd
#> SpatialData object
#>   path: /data/xenium_breast.zarr
#>   images(1): morphology
#>   spatialLabels(1): cell_labels
#>   spatialPoints(1): transcripts
#>   shapes(1): cell_boundaries
#>   tables(1): table
#>   coordinate_systems: global, pixels
```

---

<div align="center">
<img src="man/figures/fig1_spatial_overview.png" width="700"
  alt="Spatial overview"/>
</div>

> **Figure 1 | Multi-modal spatial data read from a
> SpatialData-spec Zarr store.** Demo data distributed with
> the package (50 cells, 10 genes, 500 transcripts).
> (**a**) Transcript spatial map. (**b**) Cell type
> composition. (**c**) Cell boundaries colored by type.

---

## Functions

### Store-level

| Function | Description |
|----------|-------------|
| `readSpatialData()` | Discover and read all elements from a .zarr store |
| `images()` | Access image element references |
| `spatialLabels()` | Access segmentation label references |
| `spatialPoints()` | Access point DataFrames |
| `shapes()` | Access shape DataFrames |
| `tables()` | Access annotation tables |
| `coordinateSystems()` | Access coordinate system metadata |

### Element-level readers

| Function | Description |
|----------|-------------|
| `readZarrArray()` | Read a Zarr array via Rarr/pizzarr backend |
| `readParquetPoints()` | Read Parquet-backed transcript coordinates |
| `readCSVElement()` | Read CSV-backed points or shapes |
| `readSpatialTable()` | Read AnnData-format table (obs/var/X) |

### Coordinate transforms

| Function | Description |
|----------|-------------|
| `CoordinateTransform()` | Construct affine/identity transformation |
| `transformCoords()` | Apply transformation (DataFrame or matrix) |

---

## Installation

```r
# Development version:
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("CuiweiG/SpatialDataR")
```

## References

- Marconato L et al. (2024). SpatialData: an open and universal
  data framework for spatial omics. *Nat Methods* 21:2196.
  doi:10.1038/s41592-024-02212-x
