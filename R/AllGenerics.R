# R/AllGenerics.R
# Generic function definitions -- extension points for ecosystem

#' @include AllClasses.R
#' @importFrom methods setGeneric
NULL

# readSpatialData is a plain function (not generic) -- see read-zarr.R

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
#' @examples
#' library(S4Vectors)
#' pts <- DataFrame(x = c(100, 200), y = c(50, 150))
#' ct <- CoordinateTransform("affine",
#'     affine = diag(c(0.5, 0.5, 1)))
#' transformCoords(pts, ct)
setGeneric("transformCoords", function(x, transform, ...)
    standardGeneric("transformCoords"))

#' Bounding box spatial query
#'
#' Subsets a \code{SpatialData} object or individual elements
#' to a rectangular region of interest. For point/shape
#' \code{DataFrame} elements, rows outside the bounding box are
#' removed. For image/label path references, the bounding box is
#' stored as metadata (load and crop manually with
#' \code{readZarrArray}).
#'
#' This mirrors the \code{bounding_box_query()} function from
#' the Python \code{spatialdata} library.
#'
#' @param x A \code{\linkS4class{SpatialData}} object, a
#'   \code{DataFrame} with \code{x} and \code{y} columns, or a
#'   \code{SimpleList} of such elements.
#' @param xmin Numeric. Minimum x coordinate.
#' @param xmax Numeric. Maximum x coordinate.
#' @param ymin Numeric. Minimum y coordinate.
#' @param ymax Numeric. Maximum y coordinate.
#'
#' @return Same class as input, subsetted to the bounding box.
#'
#' @references
#' Marconato L et al. (2025). SpatialData: an open and universal
#' data framework for spatial omics. \emph{Nat Methods} 22:58-62.
#' \doi{10.1038/s41592-024-02212-x}
#'
#' @export
#' @rdname bboxQuery
#' @examples
#' library(S4Vectors)
#' pts <- DataFrame(
#'     x = c(1, 2, 3, 4, 5),
#'     y = c(5, 4, 3, 2, 1),
#'     gene = c("A", "B", "C", "D", "E")
#' )
#' # Query points in [2,4] x [2,4]
#' sub <- bboxQuery(pts, xmin = 2, xmax = 4, ymin = 2, ymax = 4)
#' sub
setGeneric("bboxQuery",
    function(x, xmin, xmax, ymin, ymax)
    standardGeneric("bboxQuery"))
