# R/geometry.R
# WKB geometry parsing for GeoParquet-based shapes

#' @include AllClasses.R
NULL

# ---- Internal WKB parsing helpers ----

#' Read uint32 from raw bytes
#' @param raw_bytes Raw vector.
#' @param offset 1-based position to start reading.
#' @param endian "little" or "big".
#' @return Integer value.
#' @keywords internal
.readUint32 <- function(raw_bytes, offset, endian = "little") {
    readBin(raw_bytes[offset:(offset + 3L)],
        what = "integer", n = 1L, size = 4L,
        signed = TRUE, endian = endian)
}

#' Read float64 (double) from raw bytes
#' @param raw_bytes Raw vector.
#' @param offset 1-based position to start reading.
#' @param endian "little" or "big".
#' @return Numeric value.
#' @keywords internal
.readFloat64 <- function(raw_bytes, offset,
    endian = "little") {
    readBin(raw_bytes[offset:(offset + 7L)],
        what = "double", n = 1L, size = 8L,
        endian = endian)
}

#' Determine endianness from WKB byte
#' @param byte Single raw byte.
#' @return Character: "little" or "big".
#' @keywords internal
.wkbEndian <- function(byte) {
    if (as.integer(byte) == 1L) "little" else "big"
}

#' WKB type ID to type name
#' @param type_id Integer WKB type code.
#' @return Character type name.
#' @keywords internal
.wkbTypeName <- function(type_id) {
    switch(as.character(type_id),
        "1" = "Point",
        "2" = "LineString",
        "3" = "Polygon",
        "4" = "MultiPoint",
        "5" = "MultiLineString",
        "6" = "MultiPolygon",
        "7" = "GeometryCollection",
        paste0("Unknown(", type_id, ")"))
}

#' Parse a single WKB Point
#' @param raw_bytes Raw vector.
#' @param offset 1-based start of the point data (after type).
#' @param endian Endianness.
#' @return data.frame with x, y.
#' @keywords internal
.parseWKBPoint <- function(raw_bytes, offset,
    endian = "little") {
    x <- .readFloat64(raw_bytes, offset, endian)
    y <- .readFloat64(raw_bytes, offset + 8L, endian)
    data.frame(x = x, y = y)
}

#' Parse a single WKB Polygon
#' @param raw_bytes Raw vector.
#' @param offset 1-based start of polygon data (after type).
#' @param endian Endianness.
#' @return List of coordinate matrices (one per ring).
#' @keywords internal
.parseWKBPolygon <- function(raw_bytes, offset,
    endian = "little") {
    n_rings <- .readUint32(raw_bytes, offset, endian)
    pos <- offset + 4L
    rings <- vector("list", n_rings)
    for (i in seq_len(n_rings)) {
        n_pts <- .readUint32(raw_bytes, pos, endian)
        pos <- pos + 4L
        coords <- matrix(0.0, nrow = n_pts, ncol = 2L)
        colnames(coords) <- c("x", "y")
        for (j in seq_len(n_pts)) {
            coords[j, 1L] <- .readFloat64(
                raw_bytes, pos, endian)
            coords[j, 2L] <- .readFloat64(
                raw_bytes, pos + 8L, endian)
            pos <- pos + 16L
        }
        rings[[i]] <- coords
    }
    rings
}

#' Parse a single WKB MultiPolygon
#' @param raw_bytes Raw vector.
#' @param offset 1-based start (after type).
#' @param endian Endianness.
#' @return List of lists of coordinate matrices.
#' @keywords internal
.parseWKBMultiPolygon <- function(raw_bytes, offset,
    endian = "little") {
    n_polys <- .readUint32(raw_bytes, offset, endian)
    pos <- offset + 4L
    polygons <- vector("list", n_polys)
    for (i in seq_len(n_polys)) {
        ## Each polygon has its own WKB header
        poly_endian <- .wkbEndian(raw_bytes[pos])
        pos <- pos + 1L
        ## type_id (should be 3 = Polygon)
        pos <- pos + 4L
        rings <- .parseWKBPolygon(
            raw_bytes, pos, poly_endian)
        polygons[[i]] <- rings
        ## Advance pos past this polygon's data
        n_rings <- .readUint32(
            raw_bytes, pos, poly_endian)
        pos <- pos + 4L
        for (r in seq_len(n_rings)) {
            n_pts <- .readUint32(
                raw_bytes, pos, poly_endian)
            pos <- pos + 4L + n_pts * 16L
        }
    }
    polygons
}

#' Parse a single WKB geometry
#' @param wkb Raw vector.
#' @return Parsed geometry (data.frame for Point, list for
#'   Polygon/MultiPolygon).
#' @keywords internal
.parseOneWKB <- function(wkb) {
    endian <- .wkbEndian(wkb[1L])
    type_id <- .readUint32(wkb, 2L, endian)
    switch(as.character(type_id),
        "1" = .parseWKBPoint(wkb, 6L, endian),
        "3" = .parseWKBPolygon(wkb, 6L, endian),
        "6" = .parseWKBMultiPolygon(wkb, 6L, endian),
        stop("Unsupported WKB type: ", type_id,
            call. = FALSE))
}


# ---- Exported functions ----

#' Parse WKB geometry into coordinates
#'
#' Reads Well-Known Binary (WKB) geometry data and returns
#' coordinate structures. Supports Point, Polygon, and
#' MultiPolygon types as used in SpatialData GeoParquet
#' shapes.
#'
#' @param wkb Raw vector (single geometry) or list of raw
#'   vectors (multiple geometries).
#' @return For a single Point: \code{data.frame} with x, y.
#'   For a single Polygon: list of coordinate matrices (one
#'   per ring), each with x and y columns.
#'   For a single MultiPolygon: list of polygon ring lists.
#'   For a list input: list of parsed geometries.
#'
#' @export
#' @examples
#' # Create a WKB Point (little-endian, type=1, x=1.0, y=2.0)
#' wkb_pt <- c(
#'     as.raw(0x01),
#'     writeBin(1L, raw(), size = 4, endian = "little"),
#'     writeBin(1.0, raw(), size = 8, endian = "little"),
#'     writeBin(2.0, raw(), size = 8, endian = "little"))
#' parseGeometry(wkb_pt)
parseGeometry <- function(wkb) {
    if (is.raw(wkb)) {
        return(.parseOneWKB(wkb))
    }
    if (is.list(wkb)) {
        return(lapply(wkb, .parseOneWKB))
    }
    stop("wkb must be a raw vector or list of raw vectors",
        call. = FALSE)
}

#' Extract centroids from WKB geometries
#'
#' For Point geometries, returns the coordinates directly.
#' For Polygon geometries, computes the centroid as the
#' mean of the outer ring vertices. For MultiPolygon,
#' computes the mean of all vertices across all polygons.
#'
#' @param wkb List of raw vectors (WKB geometries).
#' @return A \code{data.frame} with \code{x} and \code{y}
#'   columns.
#'
#' @export
#' @examples
#' # Two WKB Points
#' mk_pt <- function(x, y) c(
#'     as.raw(0x01),
#'     writeBin(1L, raw(), size = 4, endian = "little"),
#'     writeBin(x, raw(), size = 8, endian = "little"),
#'     writeBin(y, raw(), size = 8, endian = "little"))
#' wkbs <- list(mk_pt(1.0, 2.0), mk_pt(3.0, 4.0))
#' geometryCentroids(wkbs)
geometryCentroids <- function(wkb) {
    if (!is.list(wkb))
        stop("wkb must be a list of raw vectors",
            call. = FALSE)

    n <- length(wkb)
    xs <- numeric(n)
    ys <- numeric(n)

    for (i in seq_len(n)) {
        g <- wkb[[i]]
        endian <- .wkbEndian(g[1L])
        type_id <- .readUint32(g, 2L, endian)

        if (type_id == 1L) {
            ## Point
            xs[i] <- .readFloat64(g, 6L, endian)
            ys[i] <- .readFloat64(g, 14L, endian)
        } else if (type_id == 3L) {
            ## Polygon — centroid of outer ring
            rings <- .parseWKBPolygon(g, 6L, endian)
            outer <- rings[[1L]]
            ## Exclude closing vertex if it duplicates first
            nr <- nrow(outer)
            if (nr > 1L &&
                outer[1L, 1L] == outer[nr, 1L] &&
                outer[1L, 2L] == outer[nr, 2L]) {
                outer <- outer[-nr, , drop = FALSE]
            }
            xs[i] <- mean(outer[, 1L])
            ys[i] <- mean(outer[, 2L])
        } else if (type_id == 6L) {
            ## MultiPolygon — mean of all vertices
            polys <- .parseWKBMultiPolygon(g, 6L, endian)
            all_x <- numeric(0L)
            all_y <- numeric(0L)
            for (poly in polys) {
                for (ring in poly) {
                    all_x <- c(all_x, ring[, 1L])
                    all_y <- c(all_y, ring[, 2L])
                }
            }
            xs[i] <- mean(all_x)
            ys[i] <- mean(all_y)
        } else {
            xs[i] <- NA_real_
            ys[i] <- NA_real_
        }
    }
    data.frame(x = xs, y = ys)
}

#' Get geometry type from WKB
#'
#' Reads the type code from a WKB binary geometry.
#'
#' @param wkb Raw vector containing WKB geometry.
#' @return Character string: \code{"Point"}, \code{"Polygon"},
#'   \code{"MultiPolygon"}, etc.
#'
#' @export
#' @examples
#' wkb_pt <- c(
#'     as.raw(0x01),
#'     writeBin(1L, raw(), size = 4, endian = "little"),
#'     writeBin(0.0, raw(), size = 8, endian = "little"),
#'     writeBin(0.0, raw(), size = 8, endian = "little"))
#' geometryType(wkb_pt)
geometryType <- function(wkb) {
    if (!is.raw(wkb))
        stop("wkb must be a raw vector", call. = FALSE)
    if (length(wkb) < 5L)
        stop("wkb too short to contain a geometry",
            call. = FALSE)
    endian <- .wkbEndian(wkb[1L])
    type_id <- .readUint32(wkb, 2L, endian)
    .wkbTypeName(type_id)
}

#' Load shapes from a SpatialData shape element
#'
#' Reads a GeoParquet shapes file and returns a
#' \code{data.frame} with a \code{geometry} column of raw
#' WKB vectors plus any additional columns (radius,
#' cell_id, etc.).
#'
#' @param shape_element Either a \code{DataFrame} (already
#'   loaded with geometry column) or a list with a
#'   \code{path} field pointing to a shapes directory
#'   containing \code{shapes.parquet}.
#' @return A \code{data.frame} with \code{geometry} (list
#'   of raw vectors) and other columns.
#'
#' @export
#' @examples
#' ## loadShapeGeometry(shapes(sd)[["cell_circles"]])
loadShapeGeometry <- function(shape_element) {
    if (is.data.frame(shape_element) ||
        methods::is(shape_element, "DataFrame")) {
        ## Already loaded as DataFrame
        df <- as.data.frame(shape_element)
        if (!"geometry" %in% names(df))
            stop("No 'geometry' column found",
                call. = FALSE)
        return(df)
    }
    if (is.list(shape_element) &&
        !is.null(shape_element$path)) {
        ## Path reference — load GeoParquet
        gpq <- list.files(shape_element$path,
            "^shapes[.]parquet$", full.names = TRUE)
        if (length(gpq) == 0L)
            stop("No shapes.parquet found in: ",
                shape_element$path, call. = FALSE)
        if (!requireNamespace("arrow", quietly = TRUE))
            stop("Install 'arrow' to read GeoParquet",
                call. = FALSE)
        df <- as.data.frame(
            arrow::read_parquet(gpq[1L]))
        return(df)
    }
    stop("Unrecognized shape_element format",
        call. = FALSE)
}
