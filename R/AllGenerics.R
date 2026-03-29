# R/AllGenerics.R
# Generic function definitions — extension points for ecosystem

#' @include AllClasses.R
#' @importFrom methods setGeneric
NULL

# readSpatialData is a plain function (not generic) — see read-zarr.R

#' Access images from a SpatialData object
#'
#' @param x A \code{\linkS4class{SpatialData}} object.
#' @return A \code{SimpleList} of image element descriptors.
#' @export
#' @rdname SpatialData-accessors
#' @examples
#' sd <- new("SpatialData")
#' images(sd)
#' spatialLabels(sd)
#' spatialPoints(sd)
#' shapes(sd)
#' tables(sd)
#' coordinateSystems(sd)
setGeneric("images", function(x) standardGeneric("images"))

#' Access labels/segmentation masks
#'
#' @param x A \code{\linkS4class{SpatialData}} object.
#' @return A \code{SimpleList} of label element descriptors.
#' @export
#' @rdname SpatialData-accessors
setGeneric("spatialLabels",
    function(x) standardGeneric("spatialLabels"))

#' Access point coordinates
#'
#' @param x A \code{\linkS4class{SpatialData}} object.
#' @return A \code{SimpleList} of point element descriptors or
#'   \code{DataFrame} objects.
#' @export
#' @rdname SpatialData-accessors
setGeneric("spatialPoints",
    function(x) standardGeneric("spatialPoints"))

#' Access shape geometries
#'
#' @param x A \code{\linkS4class{SpatialData}} object.
#' @return A \code{SimpleList} of shape element descriptors or
#'   \code{DataFrame} objects.
#' @export
#' @rdname SpatialData-accessors
setGeneric("shapes", function(x) standardGeneric("shapes"))

#' Access annotation tables
#'
#' @param x A \code{\linkS4class{SpatialData}} object.
#' @return A \code{SimpleList} of table element descriptors or
#'   \code{SpatialExperiment} objects.
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
#' Transforms coordinates using an affine matrix. Supported inputs
#' include \code{DataFrame} (with \code{x} and \code{y} columns)
#' and plain numeric matrices (Nx2 for 2D, Nx3 for 3D).
#'
#' @param x Data to transform (points DataFrame, matrix, etc.).
#' @param transform A \code{\link{CoordinateTransform}} object.
#' @param ... Additional arguments.
#' @return Transformed data (same class as input).
#' @export
#' @rdname transformCoords
setGeneric("transformCoords", function(x, transform, ...)
    standardGeneric("transformCoords"))
