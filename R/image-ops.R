# R/image-ops.R
# Image operations on Zarr arrays

#' @include AllClasses.R
NULL

#' Crop a Zarr image array by bounding box
#'
#' Extracts a region of a Zarr image array specified by
#' pixel coordinates. When \pkg{Rarr} is available,
#' performs a partial read of only the required chunks
#' for memory efficiency on large images. Otherwise
#' falls back to a full read followed by subsetting.
#'
#' @param zarr_path Path to a Zarr array directory.
#' @param xmin,xmax Integer. Column (x) pixel range
#'   (1-indexed, inclusive).
#' @param ymin,ymax Integer. Row (y) pixel range
#'   (1-indexed, inclusive).
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

    zarr_path <- normalizePath(zarr_path,
        mustWork = TRUE)

    ## Try partial read via Rarr (chunk-aware, memory
    ## efficient for large images)
    zarray <- file.path(zarr_path, ".zarray")
    if (file.exists(zarray) &&
        requireNamespace("Rarr", quietly = TRUE)) {
        meta <- jsonlite::fromJSON(zarray,
            simplifyVector = TRUE)
        ndim <- length(meta[["shape"]])

        if (ndim == 2L) {
            idx <- list(ymin:ymax, xmin:xmax)
        } else if (ndim == 3L) {
            ## Full channel range; crop spatial dims
            nc <- meta[["shape"]][1L]
            idx <- list(seq_len(nc), ymin:ymax,
                xmin:xmax)
        } else {
            stop("Unsupported array dimensions: ",
                ndim, call. = FALSE)
        }

        return(tryCatch(
            Rarr::read_zarr_array(zarr_path,
                index = idx),
            error = function(e) {
                ## Fall through to full-read path
                .cropFull(zarr_path, xmin, xmax,
                    ymin, ymax)
            }
        ))
    }

    ## Fallback: full read then subset
    .cropFull(zarr_path, xmin, xmax, ymin, ymax)
}

#' Crop by full read then subset (fallback)
#' @param zarr_path Path to Zarr array.
#' @param xmin,xmax,ymin,ymax Pixel ranges.
#' @return Cropped array.
#' @keywords internal
.cropFull <- function(zarr_path, xmin, xmax,
    ymin, ymax) {
    arr <- readZarrArray(zarr_path)
    ndim <- length(dim(arr))

    if (ndim == 2L) {
        arr[ymin:ymax, xmin:xmax]
    } else if (ndim == 3L) {
        arr[, ymin:ymax, xmin:xmax]
    } else {
        stop("Unsupported array dimensions: ",
            ndim, call. = FALSE)
    }
}
