# R/image-ops.R
# Image operations on Zarr arrays

#' @include AllClasses.R
NULL

#' Crop a Zarr image array by bounding box
#'
#' Loads a region of a Zarr image array specified by
#' pixel coordinates. Useful for extracting tissue
#' regions of interest without loading the full image.
#'
#' @param zarr_path Path to a Zarr array directory.
#' @param xmin,xmax Integer. Column (x) pixel range.
#' @param ymin,ymax Integer. Row (y) pixel range.
#'
#' @return A numeric array (cropped region).
#'
#' @export
#' @examples
#' store <- system.file("extdata", "xenium_mini.zarr",
#'     package = "SpatialDataR")
#' img_path <- file.path(store, "images",
#'     "morphology", "scale0")
#' if (requireNamespace("Rarr", quietly = TRUE) ||
#'     requireNamespace("pizzarr", quietly = TRUE)) {
#'     crop <- cropImage(img_path,
#'         xmin = 1, xmax = 10,
#'         ymin = 1, ymax = 10)
#'     dim(crop)
#' }
cropImage <- function(zarr_path, xmin, xmax,
    ymin, ymax) {
    arr <- readZarrArray(zarr_path)
    ndim <- length(dim(arr))

    if (ndim == 2L) {
        arr[ymin:ymax, xmin:xmax]
    } else if (ndim == 3L) {
        ## Assume CYX or YXC layout
        arr[, ymin:ymax, xmin:xmax]
    } else {
        stop("Unsupported array dimensions: ",
            ndim, call. = FALSE)
    }
}
