# R/transforms.R
# Coordinate transformation system

#' @include AllClasses.R
#' @include AllGenerics.R
#' @importFrom S4Vectors DataFrame
#' @importFrom methods new
NULL

#' Create a CoordinateTransform
#'
#' Constructs an affine or identity coordinate transformation.
#'
#' @param type Character. \code{"identity"} or \code{"affine"}.
#' @param affine Numeric matrix (3x3 for 2D affine). Ignored for
#'   identity transforms.
#' @param input_cs Character. Input coordinate system name.
#' @param output_cs Character. Output coordinate system name.
#'
#' @return A \code{\link{CoordinateTransform}} object.
#'
#' @export
#' @examples
#' # Identity transform
#' ct <- CoordinateTransform("identity")
#'
#' # Scale + translate
#' mat <- matrix(c(0.5, 0, 10, 0, 0.5, 20, 0, 0, 1),
#'     nrow = 3, byrow = TRUE)
#' ct <- CoordinateTransform("affine", affine = mat,
#'     input_cs = "pixels", output_cs = "microns")
#' ct
CoordinateTransform <- function(
    type = c("identity", "affine"),
    affine = diag(3),
    input_cs = "global",
    output_cs = "global"
) {
    type <- match.arg(type)
    if (type == "identity") affine <- diag(3)
    new("CoordinateTransform",
        type = type, affine = affine,
        input_cs = input_cs, output_cs = output_cs)
}

#' Apply coordinate transformation to a DataFrame
#'
#' Transforms x,y coordinates in a \code{DataFrame} using an
#' affine matrix.
#'
#' @param x A \code{DataFrame} with columns \code{x} and \code{y}.
#' @param transform A \code{\link{CoordinateTransform}}.
#' @param ... Additional arguments (unused).
#'
#' @return A \code{DataFrame} with transformed coordinates.
#'
#' @export
#' @rdname transformCoords
#' @examples
#' pts <- S4Vectors::DataFrame(x = c(100, 200), y = c(50, 150))
#' # Scale by 0.5
#' mat <- matrix(c(0.5, 0, 0, 0, 0.5, 0, 0, 0, 1),
#'     nrow = 3, byrow = TRUE)
#' ct <- CoordinateTransform("affine", affine = mat)
#' transformCoords(pts, ct)
setMethod("transformCoords",
    signature("DataFrame", "CoordinateTransform"),
    function(x, transform, ...) {

    if (slot(transform, "type") == "identity") return(x)

    aff <- slot(transform, "affine")
    xy <- cbind(x$x, x$y, 1)
    transformed <- xy %*% t(aff)
    x$x <- transformed[, 1]
    x$y <- transformed[, 2]
    x
})

#' Apply coordinate transformation to a matrix
#'
#' Transforms an Nx2 numeric matrix of (x, y) coordinates.
#'
#' @param x A numeric matrix with 2 columns (x, y).
#' @param transform A \code{\link{CoordinateTransform}}.
#' @param ... Additional arguments (unused).
#'
#' @return A numeric matrix with transformed coordinates.
#'
#' @export
#' @rdname transformCoords
#' @examples
#' coords <- matrix(c(100, 200, 50, 150), ncol = 2)
#' colnames(coords) <- c("x", "y")
#' mat <- matrix(c(2, 0, 0, 0, 2, 0, 0, 0, 1),
#'     nrow = 3, byrow = TRUE)
#' ct <- CoordinateTransform("affine", affine = mat)
#' transformCoords(coords, ct)
setMethod("transformCoords",
    signature("matrix", "CoordinateTransform"),
    function(x, transform, ...) {

    if (slot(transform, "type") == "identity") return(x)

    if (ncol(x) != 2L) {
        stop(
            "Matrix must have exactly 2 columns (x, y)",
            call. = FALSE
        )
    }

    aff <- slot(transform, "affine")
    xy <- cbind(x, 1)
    transformed <- xy %*% t(aff)
    result <- transformed[, seq_len(2), drop = FALSE]
    colnames(result) <- colnames(x)
    result
})

#' Parse transformation from Zarr metadata
#'
#' Reads coordinate transformation definitions from
#' SpatialData element \code{.zattrs} metadata.
#'
#' @param metadata List. Parsed \code{.zattrs} content.
#' @return A \code{\link{CoordinateTransform}} or \code{NULL}.
#'
#' @keywords internal
#' @examples
#' meta <- list(coordinateTransformations = list(
#'     list(type = "scale", scale = c(0.2125, 0.2125))))
#' SpatialDataR:::.parseTransform(meta)
.parseTransform <- function(metadata) {
    transforms <- metadata[["coordinateTransformations"]]
    if (is.null(transforms)) return(NULL)

    ## Take first transform
    tr <- if (is.list(transforms) && length(transforms) > 0L) {
        transforms[[1L]]
    } else {
        return(NULL)
    }

    tr_type <- tr[["type"]]
    if (is.null(tr_type)) return(NULL)

    if (tr_type == "identity") {
        CoordinateTransform("identity")
    } else if (tr_type == "affine") {
        mat <- tr[["affine"]]
        if (!is.null(mat)) {
            mat <- matrix(
                unlist(mat), nrow = 3, byrow = TRUE
            )
            CoordinateTransform("affine", affine = mat)
        } else {
            NULL
        }
    } else if (tr_type == "scale") {
        scales <- unlist(tr[["scale"]])
        n <- length(scales)
        mat <- diag(n + 1L)
        for (i in seq_len(n)) mat[i, i] <- scales[i]
        CoordinateTransform("affine", affine = mat)
    } else {
        NULL
    }
}
