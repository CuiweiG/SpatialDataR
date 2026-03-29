# R/coerce.R
# Coercion methods for interoperability with Bioconductor classes
# Bridges SpatialDataR ↔ SpatialExperiment / data.frame ecosystem

#' @include AllClasses.R
#' @include AllGenerics.R
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom methods is slot
NULL

#' Coerce SpatialData to a summary data.frame
#'
#' Creates a summary \code{data.frame} listing all elements in
#' the \code{SpatialData} object with their type, name, and
#' whether data has been loaded.
#'
#' @param x A \code{\linkS4class{SpatialData}} object.
#' @return A \code{data.frame} with columns \code{type},
#'   \code{name}, \code{loaded}, and \code{nrow} (for loaded
#'   tabular elements).
#'
#' @export
#' @examples
#' store <- system.file("extdata", "xenium_mini.zarr",
#'     package = "SpatialDataR")
#' sd <- readSpatialData(store)
#' elementSummary(sd)
elementSummary <- function(x) {
    .summarizeSlot <- function(sl, type_label) {
        if (length(sl) == 0L) return(list())
        lapply(names(sl), function(nm) {
            elem <- sl[[nm]]
            loaded <- is(elem, "DataFrame") || is(elem, "list")
            nr <- if (is(elem, "DataFrame")) nrow(elem) else NA
            data.frame(type = type_label, name = nm,
                loaded = loaded, nrow = nr,
                stringsAsFactors = FALSE)
        })
    }
    rows <- c(
        .summarizeSlot(images(x), "image"),
        .summarizeSlot(spatialLabels(x), "label"),
        .summarizeSlot(spatialPoints(x), "points"),
        .summarizeSlot(shapes(x), "shapes"),
        .summarizeSlot(tables(x), "table")
    )
    if (length(rows) == 0L) {
        return(data.frame(type = character(), name = character(),
            loaded = logical(), nrow = integer(),
            stringsAsFactors = FALSE))
    }
    do.call(rbind, rows)
}

#' Extract element transform from metadata
#'
#' Parses the coordinate transformation stored in an element's
#' \code{.zattrs} metadata. Works on element descriptors
#' (list with \code{metadata}) or \code{DataFrame} objects with
#' a \code{"spatialdata_metadata"} attribute.
#'
#' @param element An element from a \code{SpatialData} accessor
#'   (e.g., \code{images(sd)[["morphology"]]}).
#' @return A \code{\link{CoordinateTransform}} or \code{NULL} if
#'   no transform is defined.
#'
#' @export
#' @examples
#' store <- system.file("extdata", "xenium_mini.zarr",
#'     package = "SpatialDataR")
#' sd <- readSpatialData(store)
#' ct <- elementTransform(images(sd)[["morphology"]])
#' ct
elementTransform <- function(element) {
    meta <- NULL
    if (is.list(element) && "metadata" %in% names(element)) {
        meta <- element$metadata
    } else if (is(element, "DataFrame")) {
        meta <- attr(element, "spatialdata_metadata")
    }
    if (is.null(meta)) return(NULL)
    .parseTransform(meta)
}

#' List coordinate systems with their elements
#'
#' Returns a summary of which elements belong to each coordinate
#' system, by examining the
#' \code{coordinateTransformations} metadata of each element.
#'
#' @param x A \code{\linkS4class{SpatialData}} object.
#' @return A named list where each entry is a character vector
#'   of element names associated with that coordinate system.
#'
#' @export
#' @examples
#' store <- system.file("extdata", "xenium_mini.zarr",
#'     package = "SpatialDataR")
#' sd <- readSpatialData(store)
#' coordinateSystemElements(sd)
coordinateSystemElements <- function(x) {
    .scanSlot <- function(sl, type_label) {
        result <- list()
        for (nm in names(sl)) {
            elem <- sl[[nm]]
            meta <- if (is.list(elem) &&
                "metadata" %in% names(elem)) {
                elem$metadata
            } else if (is(elem, "DataFrame")) {
                attr(elem, "spatialdata_metadata")
            } else {
                NULL
            }
            if (is.null(meta)) next
            trs <- meta[["coordinateTransformations"]]
            if (is.null(trs)) next
            for (tr in trs) {
                ocs <- tr[["output"]]
                if (is.null(ocs)) ocs <- "global"
                label <- paste0(type_label, ":", nm)
                result[[ocs]] <- c(result[[ocs]], label)
            }
        }
        result
    }
    parts <- list(
        .scanSlot(images(x), "image"),
        .scanSlot(spatialLabels(x), "label"),
        .scanSlot(spatialPoints(x), "points"),
        .scanSlot(shapes(x), "shapes")
    )
    ## Merge all parts
    cs_map <- list()
    for (part in parts) {
        for (nm in names(part)) {
            cs_map[[nm]] <- c(cs_map[[nm]], part[[nm]])
        }
    }
    cs_map
}

#' Subset operator for SpatialData
#'
#' Select specific elements by name using bracket notation.
#' Returns a new \code{SpatialData} object containing only the
#' specified elements.
#'
#' @param x A \code{\linkS4class{SpatialData}} object.
#' @param i Character vector of element names to keep.
#' @param j Not used.
#' @param ... Not used.
#' @param drop Not used.
#' @return A \code{\linkS4class{SpatialData}} object.
#'
#' @export
#' @rdname SpatialData-class
#' @examples
#' store <- system.file("extdata", "xenium_mini.zarr",
#'     package = "SpatialDataR")
#' sd <- readSpatialData(store)
#' sd_sub <- sd["transcripts"]
#' sd_sub
setMethod("[", c("SpatialData", "character"),
    function(x, i, j, ..., drop = FALSE) {
    .filterSimpleList <- function(sl, names) {
        keep <- intersect(names(sl), names)
        if (length(keep) == 0L) return(SimpleList())
        sl[keep]
    }
    new("SpatialData",
        images = .filterSimpleList(images(x), i),
        labels = .filterSimpleList(spatialLabels(x), i),
        points = .filterSimpleList(spatialPoints(x), i),
        shapes = .filterSimpleList(shapes(x), i),
        tables = .filterSimpleList(tables(x), i),
        coordinate_systems = coordinateSystems(x),
        metadata = slot(x, "metadata"),
        path = slot(x, "path"))
})
