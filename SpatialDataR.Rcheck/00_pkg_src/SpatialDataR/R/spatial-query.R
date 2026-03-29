# R/spatial-query.R
# Bounding box spatial queries on SpatialData elements
# Mirrors Python spatialdata.bounding_box_query()

#' @include AllClasses.R
#' @include AllGenerics.R
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom methods slot
NULL

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
#' Marconato L et al. (2024). SpatialData: an open and universal
#' data framework for spatial omics. \emph{Nat Methods} 21:2196-2209.
#' \doi{10.1038/s41592-024-02212-x}
#'
#' @export
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

#' @rdname bboxQuery
#' @export
setMethod("bboxQuery", "DataFrame",
    function(x, xmin, xmax, ymin, ymax) {
    if (!all(c("x", "y") %in% colnames(x))) {
        stop("DataFrame must have 'x' and 'y' columns",
            call. = FALSE)
    }
    keep <- x$x >= xmin & x$x <= xmax &
            x$y >= ymin & x$y <= ymax
    x[keep, , drop = FALSE]
})

#' @rdname bboxQuery
#' @export
setMethod("bboxQuery", "SpatialData",
    function(x, xmin, xmax, ymin, ymax) {
    ## Query points
    pts <- spatialPoints(x)
    pts_q <- .querySimpleList(pts, xmin, xmax, ymin, ymax)

    ## Query shapes
    shps <- shapes(x)
    shps_q <- .querySimpleList(shps, xmin, xmax, ymin, ymax)

    ## Images/labels: attach bbox metadata
    imgs <- .annotateBbox(images(x), xmin, xmax, ymin, ymax)
    lbls <- .annotateBbox(spatialLabels(x), xmin, xmax, ymin, ymax)

    new("SpatialData",
        images = imgs,
        labels = lbls,
        points = pts_q,
        shapes = shps_q,
        tables = tables(x),
        coordinate_systems = coordinateSystems(x),
        metadata = slot(x, "metadata"),
        path = slot(x, "path"))
})

#' Query a SimpleList of mixed elements
#' @param sl SimpleList of DataFrames or path refs.
#' @param xmin,xmax,ymin,ymax Bounding box.
#' @return Queried SimpleList.
#' @keywords internal
.querySimpleList <- function(sl, xmin, xmax, ymin, ymax) {
    if (length(sl) == 0L) return(sl)
    SimpleList(lapply(sl, function(elem) {
        if (is(elem, "DataFrame") &&
            all(c("x", "y") %in% colnames(elem))) {
            bboxQuery(elem, xmin, xmax, ymin, ymax)
        } else {
            elem
        }
    }))
}

#' Annotate path refs with bbox metadata
#' @param sl SimpleList of element refs.
#' @param xmin,xmax,ymin,ymax Bounding box.
#' @return Annotated SimpleList.
#' @keywords internal
.annotateBbox <- function(sl, xmin, xmax, ymin, ymax) {
    if (length(sl) == 0L) return(sl)
    SimpleList(lapply(sl, function(elem) {
        if (is.list(elem) && "path" %in% names(elem)) {
            elem$bbox <- list(
                xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax)
        }
        elem
    }))
}
