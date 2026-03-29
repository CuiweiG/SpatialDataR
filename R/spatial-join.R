# R/spatial-join.R
# Spatial join: assign points to regions by proximity

#' @include AllClasses.R
#' @importFrom S4Vectors DataFrame
NULL

#' Assign points to nearest regions
#'
#' For each point, finds the nearest region center and
#' assigns the point to that region. This enables
#' molecule-to-cell assignment without cell segmentation
#' polygons — using centroid proximity instead.
#'
#' For true point-in-polygon assignment, use
#' \pkg{sf} externally and merge the result.
#'
#' @param points A \code{DataFrame} with \code{x}, \code{y}.
#' @param regions A \code{DataFrame} with \code{x}, \code{y},
#'   and a region identifier column.
#' @param region_col Character. Column in \code{regions}
#'   containing region IDs (default: \code{"cell_id"}).
#' @param max_dist Numeric. Maximum distance for assignment.
#'   Points beyond this distance are assigned \code{NA}.
#'   Default: \code{Inf} (no limit).
#'
#' @return A copy of \code{points} with an added column
#'   \code{assigned_region}.
#'
#' @export
#' @examples
#' library(S4Vectors)
#' pts <- DataFrame(x = c(1, 2, 5), y = c(1, 2, 5))
#' regions <- DataFrame(
#'     x = c(1.5, 5.0), y = c(1.5, 5.0),
#'     cell_id = c("A", "B"))
#' result <- assignToRegions(pts, regions)
#' result$assigned_region
assignToRegions <- function(points, regions,
    region_col = "cell_id",
    max_dist = Inf) {

    if (!all(c("x", "y") %in% colnames(points)))
        stop("points must have x, y columns",
            call. = FALSE)
    if (!all(c("x", "y") %in% colnames(regions)))
        stop("regions must have x, y columns",
            call. = FALSE)
    if (!region_col %in% colnames(regions))
        stop("region_col not found", call. = FALSE)

    px <- points$x; py <- points$y
    rx <- regions$x; ry <- regions$y
    rids <- as.character(regions[[region_col]])

    assigned <- character(length(px))
    for (i in seq_along(px)) {
        dists <- sqrt((px[i] - rx)^2 +
            (py[i] - ry)^2)
        idx <- which.min(dists)
        if (dists[idx] <= max_dist) {
            assigned[i] <- rids[idx]
        } else {
            assigned[i] <- NA_character_
        }
    }

    out <- points
    out$assigned_region <- assigned
    out
}
