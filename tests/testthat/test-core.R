library(S4Vectors)

test_that("SpatialData class can be constructed", {
    sd <- new("SpatialData")
    expect_s4_class(sd, "SpatialData")
    expect_equal(length(images(sd)), 0L)
    expect_equal(length(spatialLabels(sd)), 0L)
    expect_equal(length(spatialPoints(sd)), 0L)
    expect_equal(length(shapes(sd)), 0L)
    expect_equal(length(tables(sd)), 0L)
})

test_that("readSpatialData reads mock Zarr store", {
    tmp <- tempfile()
    dir.create(file.path(tmp, "images", "he"), recursive = TRUE)
    dir.create(file.path(tmp, "points", "tx"), recursive = TRUE)
    writeLines('{"spatialdata_attrs": {"version": "0.1"}}',
               file.path(tmp, ".zattrs"))
    writeLines('{}', file.path(tmp, "images", "he", ".zattrs"))
    writeLines('{}', file.path(tmp, "points", "tx", ".zattrs"))

    sd <- readSpatialData(tmp)
    expect_s4_class(sd, "SpatialData")
    expect_equal(length(images(sd)), 1L)
    expect_equal(names(images(sd)), "he")
    expect_equal(length(spatialPoints(sd)), 1L)

    unlink(tmp, recursive = TRUE)
})

test_that("readSpatialData handles empty store", {
    tmp <- tempfile()
    dir.create(tmp)
    writeLines('{}', file.path(tmp, ".zattrs"))

    sd <- readSpatialData(tmp)
    expect_s4_class(sd, "SpatialData")
    expect_equal(length(images(sd)), 0L)

    unlink(tmp, recursive = TRUE)
})

test_that("CoordinateTransform identity works", {
    ct <- CoordinateTransform("identity")
    expect_s4_class(ct, "CoordinateTransform")
    expect_equal(slot(ct, "type"), "identity")

    pts <- DataFrame(x = c(100, 200), y = c(50, 150))
    result <- transformCoords(pts, ct)
    expect_equal(result$x, c(100, 200))
    expect_equal(result$y, c(50, 150))
})

test_that("CoordinateTransform affine scales", {
    mat <- matrix(c(0.5, 0, 0, 0, 0.5, 0, 0, 0, 1),
                  nrow = 3, byrow = TRUE)
    ct <- CoordinateTransform("affine", affine = mat)
    pts <- DataFrame(x = c(100, 200), y = c(50, 150))
    result <- transformCoords(pts, ct)
    expect_equal(result$x, c(50, 100))
    expect_equal(result$y, c(25, 75))
})

test_that("CoordinateTransform affine translates", {
    mat <- matrix(c(1, 0, 10, 0, 1, 20, 0, 0, 1),
                  nrow = 3, byrow = TRUE)
    ct <- CoordinateTransform("affine", affine = mat)
    pts <- DataFrame(x = c(0, 5), y = c(0, 5))
    result <- transformCoords(pts, ct)
    expect_equal(result$x, c(10, 15))
    expect_equal(result$y, c(20, 25))
})

test_that("show methods work", {
    sd <- new("SpatialData")
    expect_output(show(sd), "SpatialData object")

    ct <- CoordinateTransform("identity")
    expect_output(show(ct), "CoordinateTransform")
})
