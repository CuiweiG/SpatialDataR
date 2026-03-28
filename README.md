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
2024 *Nat Methods*) is the universal format for spatial omics. But
it requires Python. R/Bioconductor users must use `reticulate`
bridges, losing native lazy loading, type safety, and integration
with `SpatialExperiment`.

## The solution

`SpatialDataR` reads SpatialData `.zarr` stores **directly in R**
with zero Python dependencies:

```r
library(SpatialDataR)
sd <- readSpatialData("xenium_breast.zarr")
sd
#> SpatialData object
#>   images(1): morphology
#>   labels(1): cell_labels
#>   points(1): transcripts
#>   tables(1): table
#>   coordinate_systems: global, pixels
```

---

<div align="center">
<img src="man/figures/fig1_spatial_overview.png" width="700" alt="Spatial overview"/>
</div>

> **Figure 1 | Multi-modal spatial data accessed natively in R.**
> Xenium-format breast tissue data read from a SpatialData Zarr
> store without Python. (**a**) Transcript spatial map (500
> transcripts, 10 marker genes). (**b**) Cell type composition
> (50 cells: epithelial, stromal, immune, endothelial).
> (**c**) Cell boundaries colored by cell type with size
> proportional to cell radius. All data read via
> `readSpatialData()` + accessors. Colour palette: Wong (2011)
> *Nat Methods* 8:441.

---

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
BiocManager::install("CuiweiG/SpatialDataR")
```

## Key design principles

1. **Zero Python.** Direct Zarr/Parquet reading via `Rarr` and
   `arrow`.
2. **Lazy by default.** Images stay on disk as `DelayedArray`.
3. **Bioconductor-native.** Tables become `SpatialExperiment`.
4. **Spec-compliant.** Follows the SpatialData Zarr spec.

## References

- Marconato L et al. (2024). SpatialData: an open and universal
  data framework for spatial omics. *Nat Methods* 21:2196.
  doi:10.1038/s41592-024-02212-x
