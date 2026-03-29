#' SpatialDataR: Native R Interface to the SpatialData Zarr Format
#'
#' Provides native R access to SpatialData Zarr stores without
#' Python dependencies. Reads images, labels, points, shapes,
#' and annotation tables into Bioconductor-native objects.
#'
#' @section Why SpatialDataR:
#' The Python \code{spatialdata} library (Marconato et al. 2024)
#' established a universal data framework for spatial omics. However,
#' R/Bioconductor users currently need Python bridges (reticulate)
#' to access SpatialData stores. This package provides direct Zarr
#' reading, eliminating the Python dependency and enabling native
#' integration with \code{SpatialExperiment} and the Bioconductor
#' ecosystem.
#'
#' @section Key functions:
#' \describe{
#'   \item{Store-level}{\code{\link{readSpatialData}}}
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
#'     \code{\link{transformCoords}}}
#' }
#'
#' @references
#' Marconato L et al. (2024). SpatialData: an open and universal
#' data framework for spatial omics. \emph{Nat Methods} 21:2196-2209.
#' \doi{10.1038/s41592-024-02212-x}
#'
#' @examples
#' store <- system.file("extdata", "xenium_mini.zarr",
#'     package = "SpatialDataR")
#' sd <- readSpatialData(store)
#' sd
#'
#' # Access elements
#' images(sd)
#' spatialPoints(sd)
#'
#' @docType package
#' @name SpatialDataR-package
#' @aliases SpatialDataR-package
#' @keywords package
"_PACKAGE"
