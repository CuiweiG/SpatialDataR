# R/aggregate.R
# Region-based aggregation of spatial elements
# Mirrors Python spatialdata.aggregate() but in pure R

#' @include AllClasses.R
#' @include AllGenerics.R
#' @importFrom S4Vectors DataFrame
NULL

#' Aggregate point features by region
#'
#' Counts or summarizes point-level features (e.g., transcripts)
#' within spatial regions defined by a shape table. This is the
#' R equivalent of Python \code{spatialdata.aggregate()}.
#'
#' For transcript data, this typically creates a cell-by-gene
#' count matrix from molecule coordinates and cell assignments.
#'
#' @param points A \code{DataFrame} with columns \code{x},
#'   \code{y}, and a feature column (e.g., \code{gene}).
#' @param regions A \code{DataFrame} with region identifiers
#'   and spatial coordinates. Must share a key column with
#'   \code{points} (e.g., \code{cell_id}).
#' @param feature_col Character. Column in \code{points}
#'   containing the feature to aggregate (default:
#'   \code{"gene"}).
#' @param region_col Character. Column shared between
#'   \code{points} and \code{regions} for joining (default:
#'   \code{"cell_id"}).
#' @param fun Character. Aggregation function: \code{"count"}
#'   (default), \code{"sum"}, or \code{"mean"}.
#'
#' @return A \code{DataFrame} in wide format: rows are regions,
#'   columns are features, values are aggregated counts.
#'
#' @references
#' Marconato L et al. (2024). SpatialData: an open and
#' universal data framework for spatial omics. \emph{Nat
#' Methods} 21:2196-2209.
#' \doi{10.1038/s41592-024-02212-x}
#'
#' @export
#' @examples
#' library(S4Vectors)
#' pts <- DataFrame(
#'     x = runif(20), y = runif(20),
#'     gene = sample(c("A", "B", "C"), 20, TRUE),
#'     cell_id = sample(1:5, 20, TRUE)
#' )
#' regions <- DataFrame(
#'     cell_id = 1:5,
#'     x = runif(5), y = runif(5)
#' )
#' counts <- aggregatePoints(pts, regions)
#' counts
aggregatePoints <- function(
    points,
    regions,
    feature_col = "gene",
    region_col = "cell_id",
    fun = c("count", "sum", "mean")
) {
    fun <- match.arg(fun)

    if (!feature_col %in% colnames(points)) {
        stop("Column '", feature_col,
            "' not found in points",
            call. = FALSE)
    }
    if (!region_col %in% colnames(points)) {
        stop("Column '", region_col,
            "' not found in points",
            call. = FALSE)
    }

    pts_df <- as.data.frame(points)
    regions_df <- as.data.frame(regions)
    mat <- .buildCountMatrix(
        pts_df[[feature_col]],
        pts_df[[region_col]],
        regions_df[[region_col]])
    if (fun == "mean") {
        mat <- .normalizeCounts(mat,
            pts_df[[region_col]])
    }
    result <- DataFrame(as.data.frame(mat))
    result[[region_col]] <- sort(unique(
        regions_df[[region_col]]))
    result
}

#' Build feature-by-region count matrix
#' @param features Character vector of feature IDs.
#' @param region_ids Vector of region assignments.
#' @param valid_regions Vector of valid region IDs.
#' @return Integer matrix.
#' @keywords internal
.buildCountMatrix <- function(features, region_ids,
    valid_regions) {
    ur <- sort(unique(as.character(valid_regions)))
    ## Filter to valid regions
    valid_mask <- as.character(region_ids) %in% ur
    feat_valid <- features[valid_mask]
    rid_valid <- as.character(region_ids[valid_mask])
    ## Vectorised tabulation
    uf <- sort(unique(feat_valid))
    tbl <- table(
        factor(rid_valid, levels = ur),
        factor(feat_valid, levels = uf))
    mat <- unclass(tbl)
    storage.mode(mat) <- "integer"
    dimnames(mat) <- list(ur, uf)
    mat
}

#' Normalize count matrix by region sizes
#' @param mat Count matrix.
#' @param region_ids Region assignments.
#' @return Normalized matrix.
#' @keywords internal
.normalizeCounts <- function(mat, region_ids) {
    sizes <- table(region_ids)
    for (r in rownames(mat)) {
        n <- as.integer(sizes[r])
        if (!is.na(n) && n > 0L)
            mat[r, ] <- mat[r, ] / n
    }
    mat
}
