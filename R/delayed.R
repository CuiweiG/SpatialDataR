# R/delayed.R
# DelayedArray integration for lazy Zarr array access
# Enables out-of-memory image/label processing

#' @include AllClasses.R
#' @importFrom methods is
NULL

#' Load a Zarr array as DelayedArray
#'
#' Reads a Zarr array directory as a \code{DelayedArray}
#' for out-of-memory access. Falls back to in-memory
#' reading if \pkg{DelayedArray} or a Zarr backend is
#' not available.
#'
#' This is the recommended way to access large image and
#' label arrays, as it avoids loading entire raster data
#' into R memory.
#'
#' @param zarr_path Path to a Zarr array directory
#'   (containing \code{.zarray} metadata).
#' @return A \code{DelayedArray} if backends available,
#'   otherwise a regular array via \code{readZarrArray}.
#'
#' @export
#' @examples
#' store <- system.file("extdata", "xenium_mini.zarr",
#'     package = "SpatialDataR")
#' img_path <- file.path(store, "images",
#'     "morphology", "scale0")
#' \donttest{
#' arr <- readZarrDelayed(img_path)
#' }
readZarrDelayed <- function(zarr_path) {
    zarr_path <- normalizePath(zarr_path,
        mustWork = TRUE)

    ## Try Rarr's DelayedArray backend
    if (requireNamespace("Rarr", quietly = TRUE) &&
        requireNamespace("DelayedArray",
            quietly = TRUE)) {
        tryCatch({
            seed <- Rarr::ZarrArray(zarr_path)
            DelayedArray::DelayedArray(seed)
        }, error = function(e) {
            message(
                "DelayedArray failed, ",
                "falling back to in-memory: ",
                conditionMessage(e))
            readZarrArray(zarr_path)
        })
    } else {
        message(
            "Install Rarr + DelayedArray for ",
            "lazy access. Loading into memory.")
        readZarrArray(zarr_path)
    }
}

#' Load element as DelayedArray from SpatialData
#'
#' Convenience wrapper to load an image or label
#' element from a \code{SpatialData} object as a
#' \code{DelayedArray}.
#'
#' @param sd A \code{\linkS4class{SpatialData}} object.
#' @param element_name Character. Name of the image or
#'   label element.
#' @param scale Character. Multiscale level to load
#'   (default: \code{"scale0"}).
#'
#' @return A \code{DelayedArray} or regular array.
#'
#' @export
#' @examples
#' store <- system.file("extdata", "xenium_mini.zarr",
#'     package = "SpatialDataR")
#' sd <- readSpatialData(store)
#' \donttest{
#' img <- loadElement(sd, "morphology")
#' }
loadElement <- function(sd, element_name,
    scale = "scale0") {
    ## Check images
    elem <- NULL
    if (element_name %in% names(images(sd))) {
        elem <- images(sd)[[element_name]]
    } else if (element_name %in%
        names(spatialLabels(sd))) {
        elem <- spatialLabels(sd)[[element_name]]
    }

    if (is.null(elem)) {
        stop("Element '", element_name,
            "' not found in images or labels",
            call. = FALSE)
    }

    if (!is.list(elem) || !"path" %in% names(elem)) {
        stop("Element is already loaded (not a ",
            "path reference)", call. = FALSE)
    }

    arr_path <- file.path(elem$path, scale)
    if (!dir.exists(arr_path)) {
        ## Try without scale subdirectory
        arr_path <- elem$path
    }

    readZarrDelayed(arr_path)
}
