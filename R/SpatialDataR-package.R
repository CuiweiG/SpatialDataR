#' SpatialDataR: Native R Interface to the SpatialData Zarr Format
#'
#' Provides native R access to SpatialData Zarr stores without
#' Python dependencies. Reads images, labels, points, shapes,
#' and annotation tables into Bioconductor-native objects.
#' Supports affine coordinate transforms with composition and
#' inversion, bounding box spatial queries, and interoperability
#' with \code{SpatialExperiment}.
#'
#' @section Why SpatialDataR:
#' The Python \code{spatialdata} library (Marconato et al. 2025)
#' established a universal data framework for spatial omics based
#' on OME-NGFF (Moore et al. 2023). R/Bioconductor users currently
#' need Python bridges (reticulate) to access SpatialData stores.
#' This package provides direct Zarr reading, eliminating the
#' Python dependency and enabling native integration with
#' \code{SpatialExperiment} and the Bioconductor ecosystem.
#'
#' @section Key functions:
#' \describe{
#'   \item{Store-level}{\code{\link{readSpatialData}},
#'     \code{\link{elementSummary}}}
#'   \item{Element readers}{\code{\link{readZarrArray}},
#'     \code{\link{readParquetPoints}},
#'     \code{\link{readCSVElement}},
#'     \code{\link{readSpatialTable}}}
#'   \item{Accessors}{\code{\link{images}},
#'     \code{\link{spatialLabels}},
#'     \code{\link{spatialPoints}},
#'     \code{\link{shapes}},
#'     \code{\link{tables}},
#'     \code{\link{coordinateSystems}}}
#'   \item{Transforms}{\code{\link{CoordinateTransform}},
#'     \code{\link{transformCoords}},
#'     \code{\link{composeTransforms}},
#'     \code{\link{invertTransform}},
#'     \code{\link{elementTransform}}}
#'   \item{Spatial queries}{\code{\link{bboxQuery}}}
#'   \item{Aggregation}{\code{\link{aggregatePoints}}}
#'   \item{Multi-sample}{\code{\link{combineSpatialData}},
#'     \code{\link{filterSample}}}
#'   \item{Writing}{\code{\link{writeSpatialData}}}
#'   \item{Lazy access}{\code{\link{readZarrDelayed}},
#'     \code{\link{loadElement}}}
#'   \item{Validation}{\code{\link{validateSpatialData}}}
#'   \item{Interoperability}{\code{\link{elementSummary}},
#'     \code{\link{coordinateSystemElements}}}
#' }
#'
#' @references
#' Marconato L et al. (2025). SpatialData: an open and universal
#' data framework for spatial omics. \emph{Nat Methods} 22:58-62.
#' \doi{10.1038/s41592-024-02212-x}
#'
#' Moore J et al. (2023). OME-Zarr: a cloud-optimized bioimaging
#' file format with international community support.
#' \emph{Histochem Cell Biol} 160:223-251.
#' \doi{10.1007/s00418-023-02209-1}
#'
#' Righelli D et al. (2022). SpatialExperiment: infrastructure for
#' spatially-resolved transcriptomics data in R using Bioconductor.
#' \emph{Bioinformatics} 38:3128-3131.
#' \doi{10.1093/bioinformatics/btac299}
#'
#' @examples
#' store <- system.file("extdata", "xenium_mini.zarr",
#'     package = "SpatialDataR")
#' sd <- readSpatialData(store)
#' sd
#'
#' # Element summary
#' elementSummary(sd)
#'
#' # Bounding box query
#' library(S4Vectors)
#' pts <- spatialPoints(sd)[["transcripts"]]
#' sub <- bboxQuery(pts, xmin = 0, xmax = 2, ymin = 0, ymax = 2)
#' nrow(sub)
#'
#' @docType package
#' @name SpatialDataR-package
#' @aliases SpatialDataR-package
#' @keywords package
"_PACKAGE"
