# R/transforms.R
# Coordinate transformation system
# Supports identity, scale, affine, composition, and 3D transforms
# following the SpatialData / OME-NGFF specification.

#' @include AllClasses.R
#' @include AllGenerics.R
#' @importFrom S4Vectors DataFrame
#' @importFrom methods new
NULL

#' Create a CoordinateTransform
#'
#' Constructs an affine or identity coordinate transformation.
#' Supports 2D (3x3 matrix) and 3D (4x4 matrix) transforms.
#'
#' @param type Character. \code{"identity"} or \code{"affine"}.
#' @param affine Numeric matrix. 3x3 for 2D or 4x4 for 3D affine.
#'   Ignored for identity transforms.
#' @param input_cs Character. Input coordinate system name.
#' @param output_cs Character. Output coordinate system name.
#'
#' @return A \code{\link{CoordinateTransform}} object.
#'
#' @seealso \code{\link{composeTransforms}} for chaining,
#'   \code{\link{invertTransform}} for inversion.
#'
#' @references
#' OME-NGFF coordinate transformations specification.
#' \url{https://ngff.openmicroscopy.org/latest/}
#'
#' @export
#' @examples
#' # Identity transform
#' ct <- CoordinateTransform("identity")
#'
#' # 2D scale + translate
#' mat <- matrix(c(0.5, 0, 10, 0, 0.5, 20, 0, 0, 1),
#'     nrow = 3, byrow = TRUE)
#' ct <- CoordinateTransform("affine", affine = mat,
#'     input_cs = "pixels", output_cs = "microns")
#' ct
#'
#' # 3D affine (4x4)
#' ct3d <- CoordinateTransform("affine", affine = diag(4) * 0.5)
CoordinateTransform <- function(
    type = c("identity", "affine"),
    affine = diag(3),
    input_cs = "global",
    output_cs = "global"
) {
    type <- match.arg(type)
    if (type == "identity") affine <- diag(nrow(affine))
    new("CoordinateTransform",
        type = type, affine = affine,
        input_cs = input_cs, output_cs = output_cs)
}

#' Compose two coordinate transforms
#'
#' Chains two affine transforms: first applies \code{first}, then
#' \code{second}. The resulting affine is \code{second \%*\% first}.
#' Coordinate system labels are inherited: \code{input_cs} from
#' \code{first} and \code{output_cs} from \code{second}.
#'
#' @param first A \code{\link{CoordinateTransform}} applied first.
#' @param second A \code{\link{CoordinateTransform}} applied second.
#' @return A new \code{\link{CoordinateTransform}}.
#'
#' @references
#' Marconato L et al. (2025). SpatialData: an open and universal
#' data framework for spatial omics. \emph{Nat Methods} 22:58-62.
#' \doi{10.1038/s41592-024-02212-x}
#'
#' @export
#' @examples
#' # Scale then translate
#' s <- CoordinateTransform("affine",
#'     affine = diag(c(0.5, 0.5, 1)),
#'     input_cs = "pixels", output_cs = "scaled")
#' t <- CoordinateTransform("affine",
#'     affine = matrix(c(1,0,10, 0,1,20, 0,0,1),
#'         nrow = 3, byrow = TRUE),
#'     input_cs = "scaled", output_cs = "microns")
#' combined <- composeTransforms(s, t)
#' combined
composeTransforms <- function(first, second) {
    a1 <- slot(first, "affine")
    a2 <- slot(second, "affine")

    ## Pad dimensions if mismatched (2D + 3D)
    if (nrow(a1) != nrow(a2)) {
        n <- max(nrow(a1), nrow(a2))
        a1 <- .padAffine(a1, n)
        a2 <- .padAffine(a2, n)
    }

    composed <- a2 %*% a1
    tp <- if (all(composed == diag(nrow(composed)))) {
        "identity"
    } else {
        "affine"
    }
    CoordinateTransform(tp,
        affine = composed,
        input_cs = slot(first, "input_cs"),
        output_cs = slot(second, "output_cs"))
}

#' Invert a coordinate transform
#'
#' Computes the inverse affine transformation. Useful for mapping
#' from physical back to pixel coordinates.
#'
#' @param transform A \code{\link{CoordinateTransform}}.
#' @return A new \code{\link{CoordinateTransform}} with inverted
#'   affine and swapped coordinate system labels.
#'
#' @export
#' @examples
#' mat <- diag(c(0.2125, 0.2125, 1))
#' ct <- CoordinateTransform("affine", affine = mat,
#'     input_cs = "pixels", output_cs = "microns")
#' inv <- invertTransform(ct)
#' inv  # microns -> pixels
invertTransform <- function(transform) {
    if (slot(transform, "type") == "identity") {
        return(CoordinateTransform("identity",
            input_cs = slot(transform, "output_cs"),
            output_cs = slot(transform, "input_cs")))
    }
    inv <- solve(slot(transform, "affine"))
    CoordinateTransform("affine",
        affine = inv,
        input_cs = slot(transform, "output_cs"),
        output_cs = slot(transform, "input_cs"))
}

#' Pad affine matrix to larger dimension
#' @param mat Square matrix.
#' @param n Target dimension.
#' @return Padded square matrix.
#' @keywords internal
.padAffine <- function(mat, n) {
    if (nrow(mat) == n) return(mat)
    out <- diag(n)
    r <- nrow(mat)
    out[seq_len(r), seq_len(r)] <- mat
    out
}

## --- Methods --------------------------------------------------------

#' @rdname transformCoords
#' @export
setMethod("transformCoords",
    signature("DataFrame", "CoordinateTransform"),
    function(x, transform, ...) {
    if (slot(transform, "type") == "identity") return(x)
    aff <- slot(transform, "affine")
    xy <- cbind(x$x, x$y, 1)
    transformed <- xy %*% t(aff[seq_len(3), seq_len(3)])
    x$x <- transformed[, 1]
    x$y <- transformed[, 2]
    x
})

#' @rdname transformCoords
#' @export
setMethod("transformCoords",
    signature("matrix", "CoordinateTransform"),
    function(x, transform, ...) {
    if (slot(transform, "type") == "identity") return(x)
    nc <- ncol(x)
    if (!nc %in% c(2L, 3L)) {
        stop("Matrix must have 2 (x,y) or 3 (x,y,z) columns",
            call. = FALSE)
    }
    aff <- slot(transform, "affine")
    dim_needed <- nc + 1L
    if (nrow(aff) < dim_needed) {
        aff <- .padAffine(aff, dim_needed)
    }
    aug <- cbind(x, 1)
    sub <- aff[seq_len(dim_needed), seq_len(dim_needed)]
    transformed <- aug %*% t(sub)
    result <- transformed[, seq_len(nc), drop = FALSE]
    colnames(result) <- colnames(x)
    result
})

## --- Parsing --------------------------------------------------------

#' Parse transformation from Zarr metadata
#'
#' Reads coordinate transformation definitions from
#' SpatialData element \code{.zattrs} metadata. Supports
#' identity, affine, scale, translation, and sequence types
#' from the OME-NGFF specification.
#'
#' @param metadata List. Parsed \code{.zattrs} content.
#' @return A \code{\link{CoordinateTransform}} or \code{NULL}.
#'
#' @keywords internal
.parseTransform <- function(metadata) {
    transforms <- metadata[["coordinateTransformations"]]
    if (is.null(transforms) || length(transforms) == 0L) {
        return(NULL)
    }

    ## Parse each transform in the list
    parsed <- lapply(transforms, .parseSingleTransform)
    parsed <- Filter(Negate(is.null), parsed)
    if (length(parsed) == 0L) return(NULL)

    ## Compose sequence
    if (length(parsed) == 1L) return(parsed[[1L]])
    Reduce(composeTransforms, parsed)
}

#' Parse a single transform entry
#' @param tr List. One transform definition.
#' @return A \code{CoordinateTransform} or \code{NULL}.
#' @keywords internal
.parseSingleTransform <- function(tr) {
    tr_type <- tr[["type"]]
    if (is.null(tr_type)) return(NULL)

    if (tr_type == "identity") {
        CoordinateTransform("identity")
    } else if (tr_type == "affine") {
        mat <- tr[["affine"]]
        if (is.null(mat)) return(NULL)
        vals <- unlist(mat)
        n <- length(vals)
        dim <- as.integer(round(sqrt(n)))
        mat <- matrix(vals, nrow = dim, byrow = TRUE)
        CoordinateTransform("affine", affine = mat)
    } else if (tr_type == "scale") {
        scales <- unlist(tr[["scale"]])
        n <- length(scales)
        mat <- diag(n + 1L)
        for (i in seq_len(n)) mat[i, i] <- scales[i]
        CoordinateTransform("affine", affine = mat)
    } else if (tr_type == "translation") {
        offsets <- unlist(tr[["translation"]])
        n <- length(offsets)
        mat <- diag(n + 1L)
        mat[seq_len(n), n + 1L] <- offsets
        CoordinateTransform("affine", affine = mat)
    } else if (tr_type == "sequence") {
        ## Recursive: parse inner transforms
        inner <- tr[["transformations"]]
        if (is.null(inner)) return(NULL)
        parsed <- lapply(inner, .parseSingleTransform)
        parsed <- Filter(Negate(is.null), parsed)
        if (length(parsed) == 0L) return(NULL)
        if (length(parsed) == 1L) return(parsed[[1L]])
        Reduce(composeTransforms, parsed)
    } else {
        NULL
    }
}
