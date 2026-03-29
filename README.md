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
spatial omics, building on the OME-NGFF specification (Moore et al.
2023). R-side support has remained limited, with no mature native
Bioconductor-oriented reader providing direct, zero-Python access
to SpatialData stores.

## What SpatialDataR provides

`SpatialDataR` reads SpatialData `.zarr` stores directly in R:

- **Store discovery**: `readSpatialData()` parses `.zattrs`
  metadata and discovers all element types
- **Element readers**: `readZarrArray()` (Rarr/pizzarr),
  `readParquetPoints()` (arrow), `readCSVElement()`,
  `readSpatialTable()` (AnnData → SpatialExperiment)
- **Coordinate transforms**: full OME-NGFF support — identity,
  scale, translation, affine, sequence, with `composeTransforms()`
  and `invertTransform()`; 2D and 3D
- **Spatial queries**: `bboxQuery()` for bounding box subsetting
  (mirrors Python `spatialdata.bounding_box_query()`)
- **Interoperability**: `elementSummary()`,
  `coordinateSystemElements()`, `[` subsetting, SpatialExperiment
  coercion

```r
library(SpatialDataR)
sd <- readSpatialData("xenium_breast.zarr")
sd
#> SpatialData object
#>   path: /data/xenium_breast.zarr
#>   images(1): morphology
#>   spatialLabels(1): cell_labels
#>   spatialPoints(1): transcripts [500 rows]
#>   shapes(1): cell_boundaries [50 rows]
#>   tables(1): table
#>   coordinate_systems: global, pixels

# Spatial query
sub <- bboxQuery(sd, xmin = 0, xmax = 100, ymin = 0, ymax = 100)

# Transform coordinates
ct <- elementTransform(images(sd)[["morphology"]])
inv <- invertTransform(ct)
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
| `elementSummary()` | Tabulate all elements with metadata |
| `coordinateSystemElements()` | Map coordinate systems → elements |
| `[` | Subset by element name |

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
| `CoordinateTransform()` | Construct affine/identity transform (2D/3D) |
| `transformCoords()` | Apply transform (DataFrame or matrix) |
| `composeTransforms()` | Chain two transforms |
| `invertTransform()` | Compute inverse transform |
| `elementTransform()` | Extract transform from element metadata |

### Spatial queries

| Function | Description |
|----------|-------------|
| `bboxQuery()` | Bounding box spatial query |

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
- Moore J et al. (2023). OME-Zarr: a cloud-optimized bioimaging
  file format. *Histochem Cell Biol* 160:223.
  doi:10.1007/s00418-023-02209-1
- Righelli D et al. (2022). SpatialExperiment: infrastructure for
  spatially-resolved transcriptomics data in R. *Bioinformatics*
  38:3128. doi:10.1093/bioinformatics/btac299
