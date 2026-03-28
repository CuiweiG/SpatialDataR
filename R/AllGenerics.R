# R/AllGenerics.R
# Generic function definitions �?extension points for ecosystem

#' @include AllClasses.R
#' @importFrom methods setGeneric
NULL

# readSpatialData is a plain function (not generic) — see read-zarr.R

#' Access images from a SpatialData object
#'
#' @param x A \code{\linkS4class{SpatialData}} object.
#' @return A \code{SimpleList} of image references.
#' @export
#' @rdname SpatialData-accessors
#' @examples
#' sd <- new("SpatialData")
#' images(sd)
#' spatialLabels(sd)
#' spatialPoints(sd)
#' shapes(sd)
#' tables(sd)
setGeneric("images", function(x) standardGeneric("images"))

#' Access labels/segmentation masks
#'
#' @param x A \code{\linkS4class{SpatialData}} object.
#' @return A \code{SimpleList} of label references.
#' @export
#' @rdname SpatialData-accessors
setGeneric("spatialLabels", function(x) standardGeneric("spatialLabels"))

#' Access point coordinates
#'
#' @param x A \code{\linkS4class{SpatialData}} object.
#' @return A \code{SimpleList} of point DataFrames.
#' @export
#' @rdname SpatialData-accessors
setGeneric("spatialPoints", function(x) standardGeneric("spatialPoints"))

#' Access shape geometries
#'
#' @param x A \code{\linkS4class{SpatialData}} object.
#' @return A \code{SimpleList} of shape DataFrames.
#' @export
#' @rdname SpatialData-accessors
setGeneric("shapes", function(x) standardGeneric("shapes"))

#' Access annotation tables
#'
#' @param x A \code{\linkS4class{SpatialData}} object.
#' @return A \code{SimpleList} of \code{SpatialExperiment} objects.
#' @export
#' @rdname SpatialData-accessors
setGeneric("tables", function(x) standardGeneric("tables"))

#' Access coordinate systems
#'
#' @param x A \code{\linkS4class{SpatialData}} object.
#' @return A named list of coordinate system definitions.
#' @export
#' @rdname SpatialData-accessors
setGeneric("coordinateSystems", function(x)
    standardGeneric("coordinateSystems"))

#' Apply coordinate transformation
#'
#' @param x Data to transform (points DataFrame, shapes, etc.).
#' @param transform A \code{\link{CoordinateTransform}} object.
#' @param ... Additional arguments.
#' @return Transformed data (same class as input).
#' @export
#' @rdname transformCoords
setGeneric("transformCoords", function(x, transform, ...)
    standardGeneric("transformCoords"))
