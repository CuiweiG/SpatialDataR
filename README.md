<div align="center">

# SpatialDataR

*Native R Interface to the SpatialData Zarr Format
for Spatial Omics*

[![R-CMD-check](https://github.com/CuiweiG/SpatialDataR/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/CuiweiG/SpatialDataR/actions/workflows/R-CMD-check.yml)
[![License: Artistic-2.0](https://img.shields.io/badge/License-Artistic--2.0-blue.svg)](https://opensource.org/licenses/Artistic-2.0)
[![BioC status](https://img.shields.io/badge/BioC-0.99.0-orange.svg)](https://bioconductor.org/)

</div>

---

## Motivation

SpatialData (Marconato et al. 2024) established a
universal Zarr-based format for spatial omics, adopted
by the scverse ecosystem and supported by 10x Genomics,
Vizgen, and NanoString platforms. However, R/Bioconductor
users currently require Python (via reticulate) to access
these stores — creating friction in analysis workflows
that otherwise run entirely in R.

**SpatialDataR** provides the first complete, zero-Python
R implementation of the SpatialData data model: store
reading, coordinate transforms, spatial queries, region
aggregation, and multi-sample operations.

---

## 1. Native Zarr Store Reading

> **No other R package can read SpatialData `.zarr`
> stores.** `readSpatialData()` discovers all five
> element types (images, labels, points, shapes, tables)
> and two coordinate systems from a single function call,
> returning Bioconductor-native objects.

<div align="center">
<img src="man/figures/fig1_store_reading.png"
    width="480" alt="Store reading"/>
</div>

> Multi-modal spatial map read from a SpatialData Zarr
> store. Transcripts (dots) overlaid on cell boundaries
> (colored outlines) with cell type annotations, loaded
> via `readSpatialData()`. Data: simulated Xenium breast
> cancer (50 cells, 10 genes, 500 transcripts).

```r
library(SpatialDataR)
sd <- readSpatialData("experiment.zarr")
sd
#> SpatialData object
#>   images(1): morphology
#>   spatialPoints(1): transcripts [500 rows]
#>   shapes(1): cell_boundaries [50 rows]
#>   tables(1): table
#>   coordinate_systems: global, pixels
```

**Comparison:**
[anndataR](https://bioconductor.org/packages/anndataR)
reads AnnData h5ad/zarr but has no spatial element
discovery or coordinate system support.
[SpatialExperiment](https://bioconductor.org/packages/SpatialExperiment)
(Righelli et al. 2022) stores spatial data but cannot
read SpatialData-format Zarr stores.

---

## 2. Bounding Box Spatial Query

> **R equivalent of Python
> `spatialdata.bounding_box_query()`.**
> `bboxQuery()` subsets all spatial elements to a
> rectangular region of interest — essential for
> analyzing tissue microenvironments.

<div align="center">
<img src="man/figures/fig2_spatial_query.png"
    width="480" alt="Spatial query"/>
</div>

> Bounding box query selecting 114/500 transcripts
> within the ROI [1, 3] × [1, 3] µm (pink box). Gray:
> excluded points. Colored: top 4 genes within ROI.

```r
sub <- bboxQuery(sd,
    xmin = 1, xmax = 3, ymin = 1, ymax = 3)
# Also works directly on DataFrames:
pts_roi <- bboxQuery(pts, 1, 3, 1, 3)
```

**Comparison:**
[Voyager](https://bioconductor.org/packages/Voyager)
(Moses & Pachter 2023) provides spatial autocorrelation
statistics but no bounding box query on SpatialData
stores.

---

## 3. Region Aggregation

> **Molecules → cell × gene matrix in one call.**
> `aggregatePoints()` mirrors Python
> `spatialdata.aggregate()`, converting transcript
> coordinates into quantification matrices grouped by
> spatial regions.

<div align="center">
<img src="man/figures/fig3_aggregation.png"
    width="480" alt="Aggregation"/>
</div>

> Cell × gene expression heatmap (row-normalized
> fractions) produced by `aggregatePoints()`, grouped
> by cell type. 50 cells × 10 genes. This is the
> critical bridge from molecule-level to cell-level
> analysis.

```r
counts <- aggregatePoints(
    spatialPoints(sd)[["transcripts"]],
    shapes(sd)[["cell_boundaries"]])
# Returns 50 x 10 DataFrame (cell × gene)
```

**Comparison:**
[MoleculeExperiment](https://bioconductor.org/packages/MoleculeExperiment)
(Parker et al. 2023) stores molecule-level data but
does not aggregate by arbitrary region DataFrames
from SpatialData stores.

---

## 4. Coordinate Transform Composition

> **Full OME-NGFF transform support: identity, scale,
> translation, affine, sequence.**
> `composeTransforms()` chains transforms;
> `invertTransform()` computes the inverse. Essential
> for aligning multi-modal spatial data (e.g., imaging
> + transcriptomics at different resolutions).

<div align="center">
<img src="man/figures/fig4_transforms.png"
    width="480" alt="Transforms"/>
</div>

> Five landmark points transformed from pixel
> coordinates (gray ×) to global coordinates
> (orange ●) via a composed scale + translation.
> The composed transform is a single affine matrix
> computed by `composeTransforms()`.

```r
ct <- composeTransforms(
    CoordinateTransform("affine",
        affine = diag(c(0.2125, 0.2125, 1))),
    CoordinateTransform("affine",
        affine = matrix(c(1,0,1, 0,1,0.5, 0,0,1),
            3, byrow = TRUE)))
inv <- invertTransform(ct)  # microns → pixels
```

**Comparison:**
No existing R/Bioconductor package implements OME-NGFF
coordinate transforms. Users currently construct ad hoc
affine matrices manually.

---

## Additional Features

| Function | Description |
|----------|-------------|
| `validateSpatialData()` | Spec compliance checker |
| `combineSpatialData()` | Multi-sample merge |
| `filterSample()` | Extract one sample |
| `elementSummary()` | Element overview table |
| `elementTransform()` | Extract transforms |
| `names()` / `length()` / `[` | R idioms |

---

## Installation

```r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("CuiweiG/SpatialDataR")
```

## References

1. Marconato L et al. (2024). SpatialData: an open
   and universal data framework for spatial omics.
   *Nat Methods* 21:2196--2209.
   doi:[10.1038/s41592-024-02212-x](https://doi.org/10.1038/s41592-024-02212-x)

2. Moore J et al. (2023). OME-Zarr: a cloud-optimized
   bioimaging file format. *Histochem Cell Biol*
   160:223--251.
   doi:[10.1007/s00418-023-02209-1](https://doi.org/10.1007/s00418-023-02209-1)

3. Righelli D et al. (2022). SpatialExperiment:
   infrastructure for spatially-resolved
   transcriptomics data in R. *Bioinformatics*
   38:3128--3131.
   doi:[10.1093/bioinformatics/btac299](https://doi.org/10.1093/bioinformatics/btac299)

4. Moses L & Pachter L (2023). Voyager: exploratory
   single-cell genomics data analysis with geospatial
   statistics. *Nat Methods* 20:1431--1441.
   doi:[10.1038/s41592-023-01920-2](https://doi.org/10.1038/s41592-023-01920-2)

5. Parker et al. (2023). MoleculeExperiment enables
   consistent infrastructure for molecule-resolved
   spatial omics. *Bioinformatics* 39:btad550.
   doi:[10.1093/bioinformatics/btad550](https://doi.org/10.1093/bioinformatics/btad550)
