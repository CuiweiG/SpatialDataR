# R/methods.R
# Accessor methods and show

#' @include AllClasses.R
#' @include AllGenerics.R
#' @importFrom methods setMethod show slot
#' @importFrom S4Vectors SimpleList
NULL

#' @rdname SpatialData-accessors
#' @export
setMethod("images", "SpatialData", function(x) slot(x, "images"))

#' @rdname SpatialData-accessors
#' @export
setMethod("spatialLabels", "SpatialData", function(x) slot(x, "labels"))

#' @rdname SpatialData-accessors
#' @export
setMethod("spatialPoints", "SpatialData",
    function(x) slot(x, "points"))

#' @rdname SpatialData-accessors
#' @export
setMethod("shapes", "SpatialData", function(x) slot(x, "shapes"))

#' @rdname SpatialData-accessors
#' @export
setMethod("tables", "SpatialData", function(x) slot(x, "tables"))

#' @rdname SpatialData-accessors
#' @export
setMethod("coordinateSystems", "SpatialData",
    function(x) slot(x, "coordinate_systems"))

#' @rdname SpatialData-class
#' @param object A \code{SpatialData} object.
#' @export
setMethod("show", "SpatialData", function(object) {
    cat("SpatialData object\n")
    cat("  path:", slot(object, "path"), "\n")
    cat("  images(", length(images(object)), "): ",
        paste(names(images(object)), collapse = ", "), "\n", sep = "")
    cat("  labels(", length(spatialLabels(object)), "): ",
        paste(names(labels(object)), collapse = ", "), "\n", sep = "")
    cat("  points(", length(spatialPoints(object)), "): ",
        paste(names(spatialPoints(object)), collapse = ", "),
        "\n", sep = "")
    cat("  shapes(", length(shapes(object)), "): ",
        paste(names(shapes(object)), collapse = ", "), "\n", sep = "")
    cat("  tables(", length(tables(object)), "): ",
        paste(names(tables(object)), collapse = ", "), "\n", sep = "")
    cs <- coordinateSystems(object)
    if (length(cs) > 0) {
        cat("  coordinate_systems: ",
            paste(names(cs), collapse = ", "), "\n")
    }
})

#' @rdname CoordinateTransform-class
#' @param object A \code{CoordinateTransform} object.
#' @export
setMethod("show", "CoordinateTransform", function(object) {
    cat("CoordinateTransform\n")
    cat("  type:", slot(object, "type"), "\n")
    cat("  ", slot(object, "input_cs"), " -> ",
        slot(object, "output_cs"), "\n", sep = "")
    if (slot(object, "type") == "affine") {
        cat("  affine:\n")
        print(slot(object, "affine"))
    }
})
