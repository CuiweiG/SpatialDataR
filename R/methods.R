# R/methods.R
# Accessor methods and show

#' @include AllClasses.R
#' @include AllGenerics.R
#' @importFrom methods setMethod show slot is
#' @importFrom S4Vectors SimpleList
NULL

#' @rdname SpatialData-accessors
#' @export
setMethod("images", "SpatialData",
    function(x) slot(x, "images"))

#' @rdname SpatialData-accessors
#' @export
setMethod("spatialLabels", "SpatialData",
    function(x) slot(x, "labels"))

#' @rdname SpatialData-accessors
#' @export
setMethod("spatialPoints", "SpatialData",
    function(x) slot(x, "points"))

#' @rdname SpatialData-accessors
#' @export
setMethod("shapes", "SpatialData",
    function(x) slot(x, "shapes"))

#' @rdname SpatialData-accessors
#' @export
setMethod("tables", "SpatialData",
    function(x) slot(x, "tables"))

#' @rdname SpatialData-accessors
#' @export
setMethod("coordinateSystems", "SpatialData",
    function(x) slot(x, "coordinate_systems"))

#' Count total elements in a SpatialData object
#' @param x A \code{SpatialData} object.
#' @return Integer. Total number of elements across all types.
#' @rdname SpatialData-class
#' @export
#' @examples
#' sd <- new("SpatialData")
#' length(sd)
setMethod("length", "SpatialData", function(x) {
    length(images(x)) + length(spatialLabels(x)) +
        length(spatialPoints(x)) + length(shapes(x)) +
        length(tables(x))
})

#' List all element names in a SpatialData object
#' @param x A \code{SpatialData} object.
#' @return Character vector of element names.
#' @rdname SpatialData-class
#' @export
#' @examples
#' sd <- new("SpatialData")
#' names(sd)
setMethod("names", "SpatialData", function(x) {
    c(names(images(x)), names(spatialLabels(x)),
        names(spatialPoints(x)), names(shapes(x)),
        names(tables(x)))
})

#' @rdname SpatialData-class
#' @param object A \code{SpatialData} object.
#' @export
setMethod("show", "SpatialData", function(object) {
    cat("SpatialData object\n")
    cat("  path:", slot(object, "path"), "\n")
    .showSlot <- function(label, sl) {
        n <- length(sl)
        nms <- if (n > 0L) {
            paste(names(sl), collapse = ", ")
        } else {
            ""
        }
        ## Show row counts for loaded DataFrames
        extra <- ""
        if (n > 0L) {
            counts <- vapply(sl, function(elem) {
                if (is(elem, "DataFrame")) nrow(elem)
                else NA_integer_
            }, integer(1))
            loaded <- sum(!is.na(counts))
            if (loaded > 0L) {
                extra <- paste0(
                    " [",
                    paste(counts[!is.na(counts)],
                        "rows", collapse = ", "),
                    "]")
            }
        }
        cat("  ", label, "(", n, "): ", nms, extra,
            "\n", sep = "")
    }
    .showSlot("images", images(object))
    .showSlot("spatialLabels", spatialLabels(object))
    .showSlot("spatialPoints", spatialPoints(object))
    .showSlot("shapes", shapes(object))
    .showSlot("tables", tables(object))
    cs <- coordinateSystems(object)
    if (length(cs) > 0L) {
        cat("  coordinate_systems: ",
            paste(names(cs), collapse = ", "), "\n")
    }
})

#' @rdname CoordinateTransform-class
#' @param object A \code{CoordinateTransform} object.
#' @export
setMethod("show", "CoordinateTransform", function(object) {
    aff <- slot(object, "affine")
    dims <- if (nrow(aff) == 4L) "3D" else "2D"
    cat("CoordinateTransform (", dims, ")\n", sep = "")
    cat("  type:", slot(object, "type"), "\n")
    cat("  ", slot(object, "input_cs"), " -> ",
        slot(object, "output_cs"), "\n", sep = "")
    if (slot(object, "type") == "affine") {
        cat("  affine:\n")
        print(aff)
    }
})
