# tests/testthat/test-geometry.R

# Helper: create WKB Point
.mk_wkb_point <- function(x, y) {
    c(as.raw(0x01),
      writeBin(1L, raw(), size = 4, endian = "little"),
      writeBin(x, raw(), size = 8, endian = "little"),
      writeBin(y, raw(), size = 8, endian = "little"))
}

# Helper: create WKB Polygon from coordinate matrix
# coords: matrix with 2 columns (x, y), first and last
#   row should be identical (closed ring)
.mk_wkb_polygon <- function(coords) {
    n_pts <- as.integer(nrow(coords))
    wkb <- c(
        as.raw(0x01),
        writeBin(3L, raw(), size = 4, endian = "little"),
        writeBin(1L, raw(), size = 4, endian = "little"),
        writeBin(n_pts, raw(), size = 4, endian = "little"))
    for (i in seq_len(n_pts)) {
        wkb <- c(wkb,
            writeBin(coords[i, 1], raw(), size = 8,
                endian = "little"),
            writeBin(coords[i, 2], raw(), size = 8,
                endian = "little"))
    }
    wkb
}

test_that("geometryType identifies Point", {
    wkb <- .mk_wkb_point(1.0, 2.0)
    expect_equal(geometryType(wkb), "Point")
})

test_that("geometryType identifies Polygon", {
    coords <- matrix(c(0, 0, 10, 0, 10, 10, 0, 10, 0, 0),
        ncol = 2, byrow = TRUE)
    wkb <- .mk_wkb_polygon(coords)
    expect_equal(geometryType(wkb), "Polygon")
})

test_that("geometryType errors on short input", {
    expect_error(geometryType(raw(3)), "too short")
})

test_that("geometryType errors on non-raw", {
    expect_error(geometryType("abc"), "raw vector")
})

test_that("parseGeometry Point returns data.frame", {
    wkb <- .mk_wkb_point(3.14, 2.72)
    result <- parseGeometry(wkb)
    expect_s3_class(result, "data.frame")
    expect_equal(result$x, 3.14)
    expect_equal(result$y, 2.72)
})

test_that("parseGeometry Polygon returns ring list", {
    coords <- matrix(c(0, 0, 10, 0, 10, 10, 0, 10, 0, 0),
        ncol = 2, byrow = TRUE)
    wkb <- .mk_wkb_polygon(coords)
    result <- parseGeometry(wkb)
    expect_type(result, "list")
    expect_equal(length(result), 1)
    ring <- result[[1]]
    expect_equal(nrow(ring), 5)
    expect_equal(colnames(ring), c("x", "y"))
    expect_equal(unname(ring[1, "x"]), 0)
    expect_equal(unname(ring[3, "x"]), 10)
})

test_that("parseGeometry list input returns list", {
    wkbs <- list(
        .mk_wkb_point(1.0, 2.0),
        .mk_wkb_point(3.0, 4.0))
    result <- parseGeometry(wkbs)
    expect_type(result, "list")
    expect_equal(length(result), 2)
    expect_equal(result[[1]]$x, 1.0)
    expect_equal(result[[2]]$y, 4.0)
})

test_that("parseGeometry errors on bad input", {
    expect_error(parseGeometry("not raw"), "raw vector")
})

test_that("geometryCentroids Point", {
    wkbs <- list(
        .mk_wkb_point(1.0, 2.0),
        .mk_wkb_point(3.0, 4.0),
        .mk_wkb_point(5.0, 6.0))
    result <- geometryCentroids(wkbs)
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 3)
    expect_equal(result$x, c(1.0, 3.0, 5.0))
    expect_equal(result$y, c(2.0, 4.0, 6.0))
})

test_that("geometryCentroids Polygon", {
    # Square 0-10 (centroid should be at 5, 5)
    coords <- matrix(c(0, 0, 10, 0, 10, 10, 0, 10, 0, 0),
        ncol = 2, byrow = TRUE)
    wkbs <- list(.mk_wkb_polygon(coords))
    result <- geometryCentroids(wkbs)
    expect_equal(nrow(result), 1)
    expect_equal(result$x, 5.0)
    expect_equal(result$y, 5.0)
})

test_that("geometryCentroids errors on non-list", {
    expect_error(geometryCentroids(raw(21)),
        "list of raw")
})

test_that("loadShapeGeometry path reference", {
    skip_if_not(
        dir.exists("C:/Users/win10/xenium_breast/data.zarr"),
        "Xenium data not available")
    skip_if_not(
        requireNamespace("arrow", quietly = TRUE),
        "arrow not installed")
    ref <- list(
        path = "C:/Users/win10/xenium_breast/data.zarr/shapes/cell_circles")
    df <- loadShapeGeometry(ref)
    expect_true(is.data.frame(df))
    expect_true("geometry" %in% names(df))
    expect_equal(nrow(df), 167780)
})

test_that("parseGeometry with real Xenium data", {
    skip_if_not(
        dir.exists("C:/Users/win10/xenium_breast/data.zarr"),
        "Xenium data not available")
    skip_if_not(
        requireNamespace("arrow", quietly = TRUE),
        "arrow not installed")
    ref <- list(
        path = "C:/Users/win10/xenium_breast/data.zarr/shapes/cell_circles")
    df <- loadShapeGeometry(ref)
    coords <- geometryCentroids(df$geometry)
    expect_equal(nrow(coords), 167780)
    expect_true(all(is.finite(coords$x)))
    expect_true(all(is.finite(coords$y)))
})

test_that("parseGeometry Polygon with real Xenium data", {
    skip_if_not(
        dir.exists("C:/Users/win10/xenium_breast/data.zarr"),
        "Xenium data not available")
    skip_if_not(
        requireNamespace("arrow", quietly = TRUE),
        "arrow not installed")
    ref <- list(
        path = "C:/Users/win10/xenium_breast/data.zarr/shapes/cell_boundaries")
    df <- loadShapeGeometry(ref)
    expect_equal(geometryType(df$geometry[[1]]), "Polygon")
    parsed <- parseGeometry(df$geometry[[1]])
    expect_type(parsed, "list")
    expect_true(nrow(parsed[[1]]) > 2)
})

test_that("parseGeometry MERFISH anatomical", {
    skip_if_not(
        dir.exists("C:/Users/win10/merfish_scverse/data.zarr"),
        "MERFISH data not available")
    skip_if_not(
        requireNamespace("arrow", quietly = TRUE),
        "arrow not installed")
    ref <- list(
        path = "C:/Users/win10/merfish_scverse/data.zarr/shapes/anatomical")
    df <- loadShapeGeometry(ref)
    expect_equal(nrow(df), 6)
    expect_equal(geometryType(df$geometry[[1]]), "Polygon")
})
