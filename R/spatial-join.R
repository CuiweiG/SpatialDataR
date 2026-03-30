# R/spatial-join.R
# Spatial join: assign points to regions by proximity or
# point-in-polygon testing

#' @include AllClasses.R
#' @include geometry.R
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

# ---- Point-in-polygon spatial join ----

#' Ray-casting point-in-polygon test (vectorized)
#'
#' Tests whether points are inside a polygon using the
#' ray-casting algorithm. Vectorized over points for
#' performance.
#'
#' @param px Numeric vector of point x coordinates.
#' @param py Numeric vector of point y coordinates.
#' @param poly_x Numeric vector of polygon vertex x coords.
#' @param poly_y Numeric vector of polygon vertex y coords.
#' @return Logical vector: TRUE if point is inside polygon.
#' @keywords internal
.pointInPolygon <- function(px, py, poly_x, poly_y) {
    n_verts <- length(poly_x)
    n_pts <- length(px)
    inside <- logical(n_pts)

    for (j in seq_len(n_verts)) {
        j_next <- if (j == n_verts) 1L else j + 1L
        yi <- poly_y[j]
        yj <- poly_y[j_next]
        xi <- poly_x[j]
        xj <- poly_x[j_next]

        ## Which points have y between yi and yj
        cond <- (yi > py) != (yj > py)
        if (!any(cond)) next

        ## x-intersection for qualifying points
        x_int <- xi + (py[cond] - yi) / (yj - yi) *
            (xj - xi)
        crosses <- px[cond] < x_int
        inside[cond] <- inside[cond] != crosses
    }
    inside
}

#' Compute bounding box of polygon vertices
#' @param poly_x Numeric vector of x coordinates.
#' @param poly_y Numeric vector of y coordinates.
#' @return Named numeric: xmin, xmax, ymin, ymax.
#' @keywords internal
.polyBBox <- function(poly_x, poly_y) {
    c(xmin = min(poly_x), xmax = max(poly_x),
      ymin = min(poly_y), ymax = max(poly_y))
}

#' Assign points to regions by point-in-polygon test
#'
#' For each point, tests membership in each polygon region
#' using the ray-casting algorithm. Points are assigned to
#' the first matching region. Uses bounding-box pre-filtering
#' for performance.
#'
#' @param points A \code{DataFrame} or \code{data.frame}
#'   with \code{x} and \code{y} columns.
#' @param regions A \code{DataFrame} or \code{data.frame}
#'   with a \code{geometry} column containing WKB-encoded
#'   polygons (as raw vectors).
#' @param region_names Optional character vector naming each
#'   region. If \code{NULL}, regions are numbered
#'   \code{"region_1"}, \code{"region_2"}, etc.
#'
#' @return Character vector of region assignments, with
#'   \code{NA} for points not inside any region. Length
#'   equals \code{nrow(points)}.
#'
#' @details
#' The algorithm:
#' \enumerate{
#'   \item Extracts polygon coordinates from WKB geometry.
#'   \item Pre-computes bounding boxes for each region.
#'   \item For each region, filters points by bounding box,
#'     then applies ray-casting only to candidates.
#'   \item Assigns each point to the first matching region.
#' }
#'
#' For the MERFISH anatomical dataset (3.7M points x 6
#' polygons), typical runtime is under 60 seconds.
#'
#' @export
#' @examples
#' library(S4Vectors)
#' ## Create simple polygon WKB (square 0-10)
#' mk_poly_wkb <- function(coords) {
#'     n_pts <- nrow(coords)
#'     wkb <- raw(0)
#'     wkb <- c(wkb, as.raw(0x01))
#'     wkb <- c(wkb, writeBin(3L, raw(), size = 4,
#'         endian = "little"))
#'     wkb <- c(wkb, writeBin(1L, raw(), size = 4,
#'         endian = "little"))
#'     wkb <- c(wkb, writeBin(n_pts, raw(), size = 4,
#'         endian = "little"))
#'     for (i in seq_len(n_pts)) {
#'         wkb <- c(wkb,
#'             writeBin(coords[i, 1], raw(), size = 8,
#'                 endian = "little"),
#'             writeBin(coords[i, 2], raw(), size = 8,
#'                 endian = "little"))
#'     }
#'     wkb
#' }
#' sq <- matrix(c(0,0, 10,0, 10,10, 0,10, 0,0),
#'     ncol = 2, byrow = TRUE)
#' pts <- DataFrame(x = c(5, 15), y = c(5, 5))
#' regions <- data.frame(
#'     geometry = I(list(mk_poly_wkb(sq))))
#' spatialJoin(pts, regions, region_names = "square")
spatialJoin <- function(points, regions,
    region_names = NULL) {

    ## Extract coordinates
    if (methods::is(points, "DataFrame")) {
        px <- as.numeric(points$x)
        py <- as.numeric(points$y)
    } else if (is.data.frame(points)) {
        px <- as.numeric(points$x)
        py <- as.numeric(points$y)
    } else {
        stop("points must be a DataFrame or data.frame ",
            "with x, y columns", call. = FALSE)
    }

    if (!all(c("x", "y") %in% colnames(points)))
        stop("points must have x, y columns",
            call. = FALSE)

    ## Get geometry column
    geom_col <- NULL
    if ("geometry" %in% colnames(regions)) {
        geom_col <- regions$geometry
    } else {
        stop("regions must have a 'geometry' column ",
            "with WKB data", call. = FALSE)
    }

    n_regions <- length(geom_col)
    n_pts <- length(px)

    if (is.null(region_names)) {
        region_names <- paste0("region_", seq_len(n_regions))
    }
    if (length(region_names) != n_regions)
        stop("region_names length must match number of ",
            "regions (", n_regions, ")", call. = FALSE)

    ## Parse all polygons and compute bounding boxes
    polys <- vector("list", n_regions)
    bboxes <- matrix(0, nrow = n_regions, ncol = 4L)
    colnames(bboxes) <- c("xmin", "xmax", "ymin", "ymax")

    for (r in seq_len(n_regions)) {
        wkb <- geom_col[[r]]
        parsed <- parseGeometry(wkb)
        geom_type <- geometryType(wkb)

        if (geom_type == "Polygon") {
            ## Use outer ring
            ring <- parsed[[1L]]
            polys[[r]] <- list(
                x = ring[, 1L], y = ring[, 2L])
        } else if (geom_type == "MultiPolygon") {
            ## Use first polygon's outer ring
            ring <- parsed[[1L]][[1L]]
            polys[[r]] <- list(
                x = ring[, 1L], y = ring[, 2L])
        } else {
            stop("Region ", r,
                " is not a Polygon or MultiPolygon",
                call. = FALSE)
        }
        bb <- .polyBBox(polys[[r]]$x, polys[[r]]$y)
        bboxes[r, ] <- bb
    }

    ## Assign points to regions
    assigned <- rep(NA_character_, n_pts)

    for (r in seq_len(n_regions)) {
        ## Bounding box filter
        candidates <- which(
            is.na(assigned) &
            px >= bboxes[r, "xmin"] &
            px <= bboxes[r, "xmax"] &
            py >= bboxes[r, "ymin"] &
            py <= bboxes[r, "ymax"])

        if (length(candidates) == 0L) next

        ## Ray-casting on candidates
        inside <- .pointInPolygon(
            px[candidates], py[candidates],
            polys[[r]]$x, polys[[r]]$y)

        assigned[candidates[inside]] <- region_names[r]
    }

    assigned
}
