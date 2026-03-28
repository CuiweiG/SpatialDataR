<div align="center">

# SpatialDataR

*Native R Interface to the SpatialData Zarr Format for Spatial Omics*

[![R-CMD-check](https://github.com/CuiweiG/SpatialDataR/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/CuiweiG/SpatialDataR/actions/workflows/R-CMD-check.yml)
[![License: Artistic-2.0](https://img.shields.io/badge/License-Artistic--2.0-blue.svg)](https://opensource.org/licenses/Artistic-2.0)
[![BioC status](https://img.shields.io/badge/BioC-0.99.0-orange.svg)](https://bioconductor.org/)

</div>

---

## The problem

[SpatialData](https://spatialdata.scverse.org/) (Marconato et al.
2024 *Nat Methods*) is the emerging universal format for spatial
omics. But it requires Python — R/Bioconductor users must use
`reticulate` bridges, losing native lazy loading, type safety,
and integration with `SpatialExperiment`.

## The solution

`SpatialDataR` reads SpatialData `.zarr` stores **directly in R**
with zero Python dependencies:

```r
library(SpatialDataR)

# Read a Xenium dataset
sd <- readSpatialData("xenium_breast.zarr")
sd
#> SpatialData object
#>   images(2): morphology, he_image
#>   labels(1): cell_labels
#>   points(1): transcripts
#>   tables(1): table

# Access elements as Bioconductor objects
spe <- tables(sd)[["table"]]       # SpatialExperiment
img <- readZarrArray("xenium_breast.zarr/images/morphology/scale0")
pts <- readParquetPoints("xenium_breast.zarr/points/transcripts")

# Coordinate transformation
ct <- CoordinateTransform("affine",
    affine = matrix(c(0.5, 0, 10, 0, 0.5, 20, 0, 0, 1),
                    nrow = 3, byrow = TRUE),
    input_cs = "pixels", output_cs = "microns")
pts_um <- transformCoords(pts, ct)
```

## Architecture

```
SpatialData .zarr store          SpatialDataR objects
========================         ========================
images/   (Zarr arrays)   --->   DelayedArray (lazy)
labels/   (Zarr arrays)   --->   DelayedArray (lazy)
points/   (Parquet)        --->   DataFrame
shapes/   (GeoParquet)     --->   DataFrame
tables/   (AnnData/Zarr)   --->   SpatialExperiment
.zattrs   (coord systems)  --->   CoordinateTransform
```

## Installation

```r
# Development version:
BiocManager::install("CuiweiG/SpatialDataR")
```

## Key design principles

1. **Zero Python.** Direct Zarr/Parquet reading via `Rarr` and
   `arrow` — no `reticulate`, no conda environments.
2. **Lazy by default.** Images stay on disk as `DelayedArray`
   until computation requires them.
3. **Bioconductor-native.** Tables become `SpatialExperiment`,
   not ad hoc data frames.
4. **Spec-compliant.** Follows the SpatialData Zarr spec
   including coordinate transformations.

## References

- Marconato L et al. (2024). SpatialData: an open and universal
  data framework for spatial omics. *Nat Methods* 21:2196.
  doi:10.1038/s41592-024-02212-x
