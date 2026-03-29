<div align="center">

# SpatialDataR

*Native R/Bioconductor Interface to the SpatialData Zarr Format for Spatial Omics*

[![R-CMD-check](https://github.com/CuiweiG/SpatialDataR/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/CuiweiG/SpatialDataR/actions/workflows/R-CMD-check.yml)
[![License: Artistic-2.0](https://img.shields.io/badge/License-Artistic--2.0-blue.svg)](https://opensource.org/licenses/Artistic-2.0)
[![BioC status](https://img.shields.io/badge/BioC-0.99.0-orange.svg)](https://bioconductor.org/)

</div>

---

## Motivation

SpatialData (Marconato et al. 2024, *Nat Methods*) established a
universal Zarr-based on-disk format for spatial omics, adopted by the
scverse ecosystem and supported by 10x Genomics Xenium, Vizgen MERFISH,
and NanoString CosMx platforms. However, R/Bioconductor users currently
require Python (via `reticulate`) to access these stores, creating
friction in analysis workflows that otherwise run entirely in R.

**SpatialDataR** provides a native R interface for reading, querying,
aggregating, transforming, and writing SpatialData-formatted Zarr
stores, exposing elements through Bioconductor-standard S4 classes:

- **Points and shapes** are loaded as `DataFrame` objects
  (CSV, Parquet, or GeoParquet backends)
- **Images and labels** are represented as path references, loadable
  as in-memory arrays via `readZarrArray()` or as out-of-memory
  `DelayedArray` objects via `readZarrDelayed()`
- **Tables** are parsed from AnnData-style obs/var Zarr groups, with
  optional coercion to `SpatialExperiment` when both the package and
  an expression matrix are available

All coordinate transforms follow the OME-NGFF specification (Moore et
al. 2023), supporting identity, scale, translation, affine, and
sequence types with 2D/3D composition and inversion.

## Validation dataset

All figures below use the **MERFISH mouse primary visual cortex
(VISp)** dataset (Moffitt et al. 2018, *Science*): **3,714,642
transcripts** across **268 genes** in **8 cortical layers** (VISp I--VI,
white matter, and surrounding tissue). This dataset was chosen because:

1. It is a **canonical spatial transcriptomics benchmark** used by the
   scverse spatialdata-sandbox for format validation
2. It contains **multiple spatial element types** (transcript points,
   cell boundary shapes, annotation tables) that exercise the full
   SpatialDataR API
3. The cortical layer architecture provides **ground-truth spatial
   organization** for validating bounding box queries and region-based
   aggregation
4. It is **publicly available** under CC0 1.0 from the SpaceTx
   consortium (reproducible via `inst/scripts/create_real_store.R`)

---

## 1. Native Zarr Store Reading

> `readSpatialData()` discovers all five element types (images, labels,
> points, shapes, tables) and coordinate systems from a single function
> call. Points and shapes are eagerly loaded as `DataFrame` when CSV
> or Parquet data files are present; images and labels remain as
> lightweight path references until explicitly loaded.

<div align="center">
<img src="man/figures/fig1_store_reading.png"
    width="700" alt="Spatial transcript map of mouse VISp"/>
</div>

> **Fig. 1.** Spatial transcript map of mouse primary visual cortex read
> from a SpatialData Zarr store via `readSpatialData()`. 3,714,642
> transcripts, 268 genes, 8 cortical layers. Top 6 genes by frequency
> are coloured; remaining 262 genes in grey. Scale bar: 500 µm. Data:
> MERFISH (Moffitt et al. 2018; CC0 1.0).

```r
library(SpatialDataR)
sd <- readSpatialData("merfish_visp.zarr")
sd
#> SpatialData object
#>   path: merfish_visp.zarr
#>   images(0):
#>   spatialLabels(0):
#>   spatialPoints(1): transcripts [3714642 rows]
#>   shapes(1): cell_boundaries [160 rows]
#>   tables(1): table
#>   coordinate_systems: global
```

**Comparison.**
[anndataR](https://bioconductor.org/packages/anndataR) reads AnnData
h5ad/zarr but has no spatial element discovery or coordinate system
support.
[SpatialExperiment](https://bioconductor.org/packages/SpatialExperiment)
(Righelli et al. 2022) stores spatial data but cannot read
SpatialData-format Zarr stores.

---

## 2. Bounding Box Spatial Query

> `bboxQuery()` subsets a `SpatialData` object or individual `DataFrame`
> elements to a rectangular region of interest, analogous to Python
> `spatialdata.bounding_box_query()`. For `SpatialData` objects, all
> element types are queried simultaneously; image and label elements
> receive bounding box metadata for downstream cropping via
> `cropImage()`.

<div align="center">
<img src="man/figures/fig2_spatial_query.png"
    width="700" alt="Bounding box spatial query"/>
</div>

> **Fig. 2.** Bounding box query on real MERFISH data. (**a**) Full
> dataset overview with 400×400 µm ROI (orange dashed box).
> (**b**) Zoomed view of the queried region with transcript gene
> identity revealed. Scale bar: 100 µm.

```r
# Query a 400 x 400 um ROI in mouse cortex
sub <- bboxQuery(sd,
    xmin = 2000, xmax = 2400,
    ymin = 5200, ymax = 5600)
spatialPoints(sub)[["transcripts"]]  # DataFrame subset
```

**Comparison.**
[Voyager](https://bioconductor.org/packages/Voyager) (Moses & Pachter
2023) provides spatial autocorrelation statistics but no bounding box
query on multi-element SpatialData stores.

---

## 3. Region-Based Aggregation

> `aggregatePoints()` converts molecule-level transcript coordinates
> into cell-by-gene count matrices grouped by spatial regions, analogous
> to Python `spatialdata.aggregate()`. Supports count, sum, and mean
> aggregation functions.

<div align="center">
<img src="man/figures/fig3_aggregation.png"
    width="700" alt="Gene enrichment dot plot"/>
</div>

> **Fig. 3.** Gene enrichment across 6 cortical layers of mouse VISp,
> produced by `aggregatePoints()` on 3.7M real MERFISH transcripts.
> Dot size encodes the percentage of layer transcripts; colour encodes
> the fraction. Top 12 genes by total frequency are shown. Layer-specific
> expression patterns are consistent with known cortical marker genes.

```r
counts <- aggregatePoints(
    spatialPoints(sd)[["transcripts"]],
    shapes(sd)[["cell_boundaries"]],
    feature_col = "gene",
    region_col = "cell_id")
# Returns 160-cell x 268-gene DataFrame
dim(counts)
#> [1] 160 269
```

**Comparison.**
[MoleculeExperiment](https://bioconductor.org/packages/MoleculeExperiment)
(Parker et al. 2023) stores molecule-level data but does not aggregate
by arbitrary region `DataFrame` objects from SpatialData stores.

---

## 4. Coordinate Transform Composition

> `composeTransforms()` chains two affine transforms (second %*% first);
> `invertTransform()` computes the matrix inverse. The parser
> (`readSpatialData()`) automatically extracts and composes
> OME-NGFF-compliant transforms from element `.zattrs` metadata,
> supporting identity, scale, translation, affine, and recursive
> sequence types in both 2D (3×3) and 3D (4×4).

<div align="center">
<img src="man/figures/fig4_transforms.png"
    width="700" alt="Coordinate transform composition"/>
</div>

> **Fig. 4.** Five cell landmark coordinates in real MERFISH tissue
> space transformed from pixel (grey ×) to global (orange ●) via a
> composed scale(0.5) + translate(500, 2000) affine. Roundtrip error
> (forward then inverse) is at machine precision (~10⁻¹³).

```r
# Compose pixel -> micron -> global
ct <- composeTransforms(
    CoordinateTransform("affine",
        affine = diag(c(0.5, 0.5, 1)),
        input_cs = "pixels", output_cs = "scaled"),
    CoordinateTransform("affine",
        affine = matrix(c(1,0,500, 0,1,2000, 0,0,1),
            3, byrow = TRUE),
        input_cs = "scaled", output_cs = "global"))
inv <- invertTransform(ct)  # global -> pixels

# Apply to coordinates
pts_global <- transformCoords(pts_df, ct)
pts_back   <- transformCoords(pts_global, inv)
max(abs(pts_df$x - pts_back$x))  # ~1e-13
```

**Comparison.** No existing R/Bioconductor package provides OME-NGFF
coordinate transform parsing, composition, or inversion. Users
currently construct ad hoc affine matrices manually.

---

## 5. Roundtrip: Read → Query → Write → Verify

> `writeSpatialData()` produces SpatialData-formatted Zarr stores
> readable by Python `spatialdata`, enabling R analysis branches
> within mixed Python/R pipelines without lossy format conversion.

<div align="center">
<img src="man/figures/fig5_roundtrip.png"
    width="700" alt="Read-write roundtrip verification"/>
</div>

> **Fig. 5.** Full roundtrip on real MERFISH data. (**a**) Read 3.7M
> transcripts from Zarr store. (**b**) Spatial query selects 648,954
> transcripts in a 600×600 µm ROI; write to new `.zarr` via
> `writeSpatialData()`. (**c**) Read back and verify identical
> transcript count.

```r
sub <- bboxQuery(sd,
    xmin = 2039, xmax = 2639,
    ymin = 5091, ymax = 5691)
writeSpatialData(sub, "subset.zarr")
sd2 <- readSpatialData("subset.zarr")
# 648,954 transcripts preserved
```

**Comparison.** No existing R/Bioconductor package can write
SpatialData-formatted Zarr stores. Users must export to CSV and
convert in Python. `writeSpatialData()` eliminates this step.

---

## Additional Features

| Function | Description |
|---|---|
| `validateSpatialData()` | Checks Zarr store against the SpatialData on-disk specification |
| `combineSpatialData()` | Merges multiple stores with automatic element-name prefixing |
| `filterSample()` | Extracts a single sample from a combined multi-sample object |
| `cropImage()` | Crops a Zarr image array by pixel bounding box |
| `readZarrDelayed()` | Loads Zarr arrays as `DelayedArray` for out-of-memory access |
| `elementSummary()` | Tabulates all elements with type, name, and row counts |
| `elementTransform()` | Extracts `CoordinateTransform` from element metadata |
| `coordinateSystemElements()` | Maps coordinate system names to their member elements |
| `assignToRegions()` | Assigns points to nearest spatial regions |
| `loadElement()` | Eagerly loads a lazy element reference into memory |
| `names()` / `length()` / `[` | Standard R accessors for element introspection |

---

## Installation

```r
# From GitHub (development version)
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("CuiweiG/SpatialDataR")

# Optional backends
BiocManager::install("Rarr")        # Zarr array reading
install.packages("arrow")           # Parquet support
BiocManager::install("SpatialExperiment")  # Table coercion
```

## References

1. Marconato L et al. (2024). SpatialData: an open and universal data
   framework for spatial omics. *Nat Methods* 21:2196--2209.
   doi:[10.1038/s41592-024-02212-x](https://doi.org/10.1038/s41592-024-02212-x)

2. Moore J et al. (2023). OME-Zarr: a cloud-optimized bioimaging file
   format with a draft specification. *Histochem Cell Biol*
   160:223--251.
   doi:[10.1007/s00418-023-02209-1](https://doi.org/10.1007/s00418-023-02209-1)

3. Righelli D et al. (2022). SpatialExperiment: infrastructure for
   spatially-resolved transcriptomics data in R using Bioconductor.
   *Bioinformatics* 38:3128--3131.
   doi:[10.1093/bioinformatics/btac299](https://doi.org/10.1093/bioinformatics/btac299)

4. Moses L & Pachter L (2023). Voyager: exploratory single-cell
   genomics data analysis with geospatial statistics. *Nat Methods*
   20:1431--1441.
   doi:[10.1038/s41592-023-01920-2](https://doi.org/10.1038/s41592-023-01920-2)

5. Parker TJ et al. (2023). MoleculeExperiment enables consistent
   infrastructure for molecule-resolved spatial omics data in
   Bioconductor. *Bioinformatics* 39:btad550.
   doi:[10.1093/bioinformatics/btad550](https://doi.org/10.1093/bioinformatics/btad550)

6. Moffitt JR et al. (2018). Molecular, spatial, and functional
   single-cell profiling of the hypothalamic preoptic region. *Science*
   362:eaau5324.
   doi:[10.1126/science.aau5324](https://doi.org/10.1126/science.aau5324)
