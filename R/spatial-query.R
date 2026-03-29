# R/spatial-query.R
# Bounding box spatial queries on SpatialData elements
# Mirrors Python spatialdata.bounding_box_query()

#' @include AllClasses.R
#' @include AllGenerics.R
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom methods slot
NULL

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
