library(S4Vectors)

test_that("SpatialData class can be constructed", {
    sd <- new("SpatialData")
    expect_s4_class(sd, "SpatialData")
    expect_equal(length(images(sd)), 0L)
    expect_equal(length(spatialLabels(sd)), 0L)
    expect_equal(length(spatialPoints(sd)), 0L)
    expect_equal(length(shapes(sd)), 0L)
    expect_equal(length(tables(sd)), 0L)
    expect_equal(length(coordinateSystems(sd)), 0L)
})

test_that("SpatialData validity rejects bad path", {
    expect_error(
        validObject(new("SpatialData", path = character(0))),
        "path"
    )
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

test_that("readSpatialData loads CSV points eagerly", {
    tmp <- tempfile()
    pts_dir <- file.path(tmp, "points", "transcripts")
    dir.create(pts_dir, recursive = TRUE)
    writeLines('{}', file.path(tmp, ".zattrs"))
    writeLines('{}', file.path(pts_dir, ".zattrs"))
    write.csv(
        data.frame(x = c(1, 2, 3), y = c(4, 5, 6),
            gene = c("A", "B", "C")),
        file.path(pts_dir, "data.csv"),
        row.names = FALSE
    )

    sd <- readSpatialData(tmp)
    pts <- spatialPoints(sd)[["transcripts"]]
    expect_s4_class(pts, "DataFrame")
    expect_equal(nrow(pts), 3L)
    expect_true("x" %in% colnames(pts))
    expect_true("y" %in% colnames(pts))

    unlink(tmp, recursive = TRUE)
})

test_that("readSpatialData loads CSV shapes eagerly", {
    tmp <- tempfile()
    shp_dir <- file.path(tmp, "shapes", "cells")
    dir.create(shp_dir, recursive = TRUE)
    writeLines('{}', file.path(tmp, ".zattrs"))
    writeLines('{}', file.path(shp_dir, ".zattrs"))
    write.csv(
        data.frame(cell_id = 1:5, x = runif(5), y = runif(5),
            radius = 0.3),
        file.path(shp_dir, "circles.csv"),
        row.names = FALSE
    )

    sd <- readSpatialData(tmp)
    shp <- shapes(sd)[["cells"]]
    expect_s4_class(shp, "DataFrame")
    expect_equal(nrow(shp), 5L)

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

test_that("transformCoords works on matrix", {
    mat <- matrix(c(2, 0, 0, 0, 3, 0, 0, 0, 1),
        nrow = 3, byrow = TRUE)
    ct <- CoordinateTransform("affine", affine = mat)
    coords <- matrix(c(10, 20, 5, 15), ncol = 2)
    colnames(coords) <- c("x", "y")
    result <- transformCoords(coords, ct)
    expect_equal(result[, 1], c(20, 40))
    expect_equal(result[, 2], c(15, 45))
})

test_that("transformCoords matrix rejects wrong dimensions", {
    ct <- CoordinateTransform("affine")
    bad <- matrix(1:8, ncol = 4)
    expect_error(transformCoords(bad, ct), "2.*3.*columns")
})

test_that("show methods work without errors", {
    sd <- new("SpatialData")
    expect_output(show(sd), "SpatialData object")

    ct <- CoordinateTransform("identity")
    expect_output(show(ct), "CoordinateTransform")
    expect_output(show(ct), "identity")

    ct2 <- CoordinateTransform("affine",
        affine = diag(3) * 0.5,
        input_cs = "px", output_cs = "um")
    expect_output(show(ct2), "affine")
    expect_output(show(ct2), "px -> um")
})

test_that("show method works on populated SpatialData", {
    tmp <- tempfile()
    dir.create(file.path(tmp, "images", "he"), recursive = TRUE)
    dir.create(file.path(tmp, "labels", "seg"), recursive = TRUE)
    writeLines('{"spatialdata_attrs": {"version": "0.1"}}',
        file.path(tmp, ".zattrs"))
    writeLines('{}', file.path(tmp, "images", "he", ".zattrs"))
    writeLines('{}', file.path(tmp, "labels", "seg", ".zattrs"))

    sd <- readSpatialData(tmp)
    expect_output(show(sd), "images\\(1\\)")
    expect_output(show(sd), "spatialLabels\\(1\\)")

    unlink(tmp, recursive = TRUE)
})
