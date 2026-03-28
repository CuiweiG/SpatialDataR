library(S4Vectors)

test_that("end-to-end: mini store reads all elements", {
    store <- system.file("extdata", "xenium_mini.zarr",
                          package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd <- readSpatialData(store)
    expect_s4_class(sd, "SpatialData")
    expect_true(length(images(sd)) >= 1)
    expect_true(length(spatialLabels(sd)) >= 1)
    expect_true(length(spatialPoints(sd)) >= 1)
    expect_true(length(shapes(sd)) >= 1)
    expect_true(length(tables(sd)) >= 1)
})

test_that("end-to-end: coordinate systems parsed", {
    store <- system.file("extdata", "xenium_mini.zarr",
                          package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd <- readSpatialData(store)
    cs <- coordinateSystems(sd)
    expect_true(length(cs) >= 1)
    expect_true("global" %in% names(cs))
})

test_that("end-to-end: transform chain works", {
    ct <- CoordinateTransform("affine",
        affine = matrix(c(0.2125, 0, 0, 0, 0.2125, 0, 0, 0, 1),
                        nrow = 3, byrow = TRUE))
    pts <- DataFrame(x = c(0, 100), y = c(0, 100))
    result <- transformCoords(pts, ct)
    expect_equal(result$x, c(0, 21.25), tolerance = 1e-6)
    expect_equal(result$y, c(0, 21.25), tolerance = 1e-6)
})

test_that("error: nonexistent path", {
    expect_error(readSpatialData("/nonexistent/zarr"))
})

test_that("error: readZarrArray on nonexistent", {
    expect_error(readZarrArray("/nonexistent/array"))
})

test_that("error: readParquetPoints without arrow", {
    skip_if(requireNamespace("arrow", quietly = TRUE))
    expect_error(readParquetPoints(tempdir()), "Install")
})

test_that("graceful: empty store", {
    tmp <- tempfile()
    dir.create(tmp)
    writeLines("{}", file.path(tmp, ".zattrs"))
    sd <- readSpatialData(tmp)
    expect_equal(length(images(sd)), 0L)
    unlink(tmp, recursive = TRUE)
})
