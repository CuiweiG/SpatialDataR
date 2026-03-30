# R/delayed.R
# DelayedArray integration for lazy Zarr array access
# Enables out-of-memory image/label processing

#' @include AllClasses.R
#' @importFrom methods is
NULL

#' Load a Zarr array as DelayedArray
#'
#' Reads a Zarr array directory as a \code{DelayedArray}
#' for out-of-memory access. For Zarr v2 stores, uses
#' \pkg{Rarr}'s \code{ZarrArray} seed for true lazy
#' chunk-level access. For Zarr v3 stores (or when Rarr
#' is not available), wraps the in-memory result in a
#' \code{DelayedArray} to provide a consistent interface.
#' Falls back to a regular array when neither
#' \code{DelayedArray} nor a Zarr backend is available.
#'
#' @param zarr_path Path to a Zarr array directory
#'   (containing \code{.zarray} or \code{zarr.json}
#'   metadata).
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

    has_rarr <- requireNamespace("Rarr", quietly = TRUE)
    has_da <- requireNamespace("DelayedArray",
        quietly = TRUE)
    is_v2 <- file.exists(
        file.path(zarr_path, ".zarray"))

    ## Zarr v2 + Rarr: true lazy access via ZarrArray seed
    if (is_v2 && has_rarr && has_da) {
        return(tryCatch({
            seed <- Rarr::ZarrArray(zarr_path)
            DelayedArray::DelayedArray(seed)
        }, error = function(e) {
            message(
                "Rarr DelayedArray failed, ",
                "falling back to in-memory: ",
                conditionMessage(e))
            readZarrArray(zarr_path)
        }))
    }

    ## Zarr v3 or no Rarr: read into memory, wrap in
    ## DelayedArray for consistent interface
    arr <- readZarrArray(zarr_path)

    if (has_da) {
        tryCatch(
            DelayedArray::DelayedArray(arr),
            error = function(e) arr
        )
    } else {
        message(
            "Install DelayedArray for lazy access. ",
            "Returning in-memory array.")
        arr
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
