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
    features <- pts_df[[feature_col]]
    region_ids <- pts_df[[region_col]]
    unique_features <- sort(unique(features))
    unique_regions <- sort(unique(
        as.data.frame(regions)[[region_col]]))

    ## Build count matrix
    mat <- matrix(0L,
        nrow = length(unique_regions),
        ncol = length(unique_features))
    rownames(mat) <- as.character(unique_regions)
    colnames(mat) <- unique_features

    for (i in seq_along(features)) {
        rid <- as.character(region_ids[i])
        fid <- features[i]
        if (rid %in% rownames(mat)) {
            mat[rid, fid] <- mat[rid, fid] + 1L
        }
    }

    if (fun == "mean" && region_col %in%
        colnames(points)) {
        region_sizes <- table(region_ids)
        for (r in rownames(mat)) {
            n <- as.integer(region_sizes[r])
            if (!is.na(n) && n > 0L) {
                mat[r, ] <- mat[r, ] / n
            }
        }
    }

    result <- DataFrame(mat)
    result[[region_col]] <- unique_regions
    result
}
