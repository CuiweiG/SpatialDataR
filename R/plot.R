# R/plot.R
# Multi-layer spatial data visualization

#' @include AllClasses.R
#' @include AllGenerics.R
#' @include geometry.R
NULL

## Suppress R CMD check NOTEs for ggplot2 NSE
utils::globalVariables(c("x", "y", "group", ".data"))

#' Plot SpatialData elements with optional image overlay
#'
#' Creates a \pkg{ggplot2} visualization of SpatialData elements.
#' Supports image backgrounds, point overlays (transcripts or
#' point-type shapes), and polygon shapes.
#'
#' @param sdata A \code{\linkS4class{SpatialData}} object.
#' @param points Character or NULL. Name of a points element
#'   to plot (from \code{spatialPoints(sdata)}).
#' @param shapes Character or NULL. Name of a shapes element
#'   to plot (from \code{shapes(sdata)}).
#' @param image Character or NULL. Name of an image element
#'   to use as background (from \code{images(sdata)}).
#' @param color_by Character or NULL. Column name in the
#'   points data to map to color.
#' @param point_size Numeric. Size for \code{geom_point}
#'   (default 0.3).
#' @param point_alpha Numeric. Alpha for points (default 0.5).
#' @param shape_fill Character. Fill color for polygon shapes
#'   (default \code{"transparent"}).
#' @param shape_color Character. Border color for polygon
#'   shapes (default \code{"red"}).
#' @param max_points Integer. Maximum number of points to
#'   plot (subsampled if exceeded). Default 500000.
#' @param ... Additional parameters (currently unused).
#'
#' @return A \code{ggplot2} object.
#'
#' @export
#' @examples
#' ## plotSpatialData(sd, points = "transcripts",
#' ##     image = "morphology_focus")
plotSpatialData <- function(sdata,
    points = NULL, shapes = NULL,
    image = NULL, color_by = NULL,
    point_size = 0.3, point_alpha = 0.5,
    shape_fill = "transparent",
    shape_color = "red",
    max_points = 500000L,
    ...) {

    if (!requireNamespace("ggplot2", quietly = TRUE))
        stop("Install 'ggplot2' for plotting:\n",
            "  install.packages('ggplot2')",
            call. = FALSE)

    p <- ggplot2::ggplot()

    ## Image background layer
    if (!is.null(image)) {
        p <- .addImageLayer(p, sdata, image)
    }

    ## Points layer
    if (!is.null(points)) {
        p <- .addPointsLayer(p, sdata, points,
            color_by = color_by,
            point_size = point_size,
            point_alpha = point_alpha,
            max_points = max_points)
    }

    ## Shapes layer
    if (!is.null(shapes)) {
        p <- .addShapesLayer(p, sdata, shapes,
            fill = shape_fill,
            color = shape_color)
    }

    p <- p + ggplot2::coord_fixed() +
        ggplot2::theme_minimal() +
        ggplot2::labs(x = "x", y = "y")

    p
}

#' Add image background layer
#' @param p ggplot object.
#' @param sdata SpatialData object.
#' @param image_name Character. Name of image element.
#' @return Updated ggplot object.
#' @keywords internal
.addImageLayer <- function(p, sdata, image_name) {
    img_ref <- images(sdata)[[image_name]]
    if (is.null(img_ref))
        stop("Image '", image_name, "' not found",
            call. = FALSE)

    img_path <- img_ref$path

    ## Find the lowest-resolution pyramid level
    subdirs <- list.dirs(img_path,
        recursive = FALSE, full.names = FALSE)
    ## Filter to numeric scale levels
    levels <- suppressWarnings(as.integer(subdirs))
    levels <- levels[!is.na(levels)]
    if (length(levels) == 0L) {
        ## Single-scale image — read directly
        arr_path <- img_path
    } else {
        ## Use highest number = lowest resolution
        arr_path <- file.path(img_path,
            as.character(max(levels)))
    }

    arr <- readZarrArray(arr_path)

    ## arr may be [C, H, W] or [H, W] or [H, W, C]
    dims <- dim(arr)
    if (length(dims) == 3L) {
        if (dims[1L] <= 4L) {
            ## [C, H, W] — take first channel
            img_mat <- arr[1L, , ]
        } else if (dims[3L] <= 4L) {
            ## [H, W, C] — take first channel
            img_mat <- arr[, , 1L]
        } else {
            img_mat <- arr[1L, , ]
        }
    } else {
        img_mat <- arr
    }

    ## Normalize to [0, 1]
    img_range <- range(img_mat, na.rm = TRUE)
    if (img_range[2L] > img_range[1L]) {
        img_mat <- (img_mat - img_range[1L]) /
            (img_range[2L] - img_range[1L])
    }

    nrow_img <- nrow(img_mat)
    ncol_img <- ncol(img_mat)

    ## Create RGB raster — gray scale
    raster_mat <- matrix(
        grDevices::rgb(
            as.vector(img_mat),
            as.vector(img_mat),
            as.vector(img_mat)),
        nrow = nrow_img, ncol = ncol_img)

    ## Compute scale factor from pyramid level
    ## Image coordinates: row = y, col = x
    ## For OME-NGFF pyramids, scale factor = 2^level
    level_idx <- if (length(levels) > 0L) max(levels) else 0L
    scale <- 2^level_idx

    xmin <- 0
    xmax <- ncol_img * scale
    ymin <- 0
    ymax <- nrow_img * scale

    p + ggplot2::annotation_raster(
        raster_mat,
        xmin = xmin, xmax = xmax,
        ymin = -ymax, ymax = -ymin) +
        ggplot2::scale_y_continuous(
            trans = "reverse")
}

#' Add points layer
#' @param p ggplot object.
#' @param sdata SpatialData object.
#' @param points_name Character. Name of points element.
#' @param color_by Column name for coloring.
#' @param point_size Point size.
#' @param point_alpha Point alpha.
#' @param max_points Maximum points to plot.
#' @return Updated ggplot object.
#' @keywords internal
.addPointsLayer <- function(p, sdata, points_name,
    color_by = NULL,
    point_size = 0.3,
    point_alpha = 0.5,
    max_points = 500000L) {

    pts_data <- spatialPoints(sdata)[[points_name]]
    if (is.null(pts_data))
        stop("Points '", points_name, "' not found",
            call. = FALSE)

    ## Convert to data.frame for ggplot
    if (methods::is(pts_data, "DataFrame")) {
        df <- as.data.frame(pts_data)
    } else {
        stop("Points element must be a DataFrame",
            call. = FALSE)
    }

    if (!all(c("x", "y") %in% names(df)))
        stop("Points must have x, y columns",
            call. = FALSE)

    ## Subsample if too many
    if (nrow(df) > max_points) {
        idx <- sample.int(nrow(df), max_points)
        df <- df[idx, , drop = FALSE]
    }

    if (!is.null(color_by) && color_by %in% names(df)) {
        p + ggplot2::geom_point(
            data = df,
            ggplot2::aes(x = x, y = y,
                color = .data[[color_by]]),
            size = point_size,
            alpha = point_alpha)
    } else {
        p + ggplot2::geom_point(
            data = df,
            ggplot2::aes(x = x, y = y),
            size = point_size,
            alpha = point_alpha,
            color = "steelblue")
    }
}

#' Add shapes layer
#' @param p ggplot object.
#' @param sdata SpatialData object.
#' @param shapes_name Character. Name of shapes element.
#' @param fill Fill color for polygons.
#' @param color Border color for polygons.
#' @return Updated ggplot object.
#' @keywords internal
.addShapesLayer <- function(p, sdata, shapes_name,
    fill = "transparent",
    color = "red") {

    shape_elem <- shapes(sdata)[[shapes_name]]
    if (is.null(shape_elem))
        stop("Shapes '", shapes_name, "' not found",
            call. = FALSE)

    df <- loadShapeGeometry(shape_elem)

    ## Determine geometry type from first element
    geom_type <- geometryType(df$geometry[[1L]])

    if (geom_type == "Point") {
        ## Point shapes — plot as circles
        centroids <- geometryCentroids(df$geometry)
        centroids$radius <- if ("radius" %in% names(df)) {
            df$radius
        } else {
            rep(1, nrow(centroids))
        }
        p + ggplot2::geom_point(
            data = centroids,
            ggplot2::aes(x = x, y = y),
            size = 0.5, alpha = 0.3,
            color = color)
    } else {
        ## Polygon / MultiPolygon
        all_coords <- list()
        for (i in seq_along(df$geometry)) {
            parsed <- parseGeometry(df$geometry[[i]])
            if (geom_type == "Polygon") {
                ## parsed = list of rings; use outer ring
                coords <- parsed[[1L]]
                all_coords[[i]] <- data.frame(
                    x = coords[, 1L],
                    y = coords[, 2L],
                    group = i)
            } else {
                ## MultiPolygon: list of polygon ring lists
                for (pi in seq_along(parsed)) {
                    coords <- parsed[[pi]][[1L]]
                    all_coords[[length(all_coords) + 1L]] <-
                        data.frame(
                            x = coords[, 1L],
                            y = coords[, 2L],
                            group = paste0(i, ".", pi))
                }
            }
        }
        poly_df <- do.call(rbind, all_coords)
        p + ggplot2::geom_polygon(
            data = poly_df,
            ggplot2::aes(x = x, y = y, group = group),
            fill = fill, color = color,
            linewidth = 0.5)
    }
}
