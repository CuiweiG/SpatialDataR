test_that("writeSpatialData roundtrip", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd <- readSpatialData(store)
    tmp <- tempfile(fileext = ".zarr")
    on.exit(unlink(tmp, recursive = TRUE))

    writeSpatialData(sd, tmp)
    expect_true(dir.exists(tmp))
    expect_true(file.exists(
        file.path(tmp, ".zattrs")))

    ## Read back
    sd2 <- readSpatialData(tmp)
    expect_s4_class(sd2, "SpatialData")
    expect_equal(length(sd2), length(sd))

    ## Points preserved
    pts1 <- spatialPoints(sd)[["transcripts"]]
    pts2 <- spatialPoints(sd2)[["transcripts"]]
    expect_equal(nrow(pts1), nrow(pts2))
})

test_that("writeSpatialData errors on existing", {
    tmp <- tempfile(fileext = ".zarr")
    dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE))

    sd <- new("SpatialData")
    expect_error(writeSpatialData(sd, tmp),
        "exists")
})

test_that("writeSpatialData overwrite works", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd <- readSpatialData(store)
    tmp <- tempfile(fileext = ".zarr")
    on.exit(unlink(tmp, recursive = TRUE))

    writeSpatialData(sd, tmp)
    writeSpatialData(sd, tmp, overwrite = TRUE)
    sd2 <- readSpatialData(tmp)
    expect_equal(length(sd2), length(sd))
})

test_that("writeSpatialData after bboxQuery", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd <- readSpatialData(store)
    sub <- bboxQuery(sd, xmin = 0, xmax = 2,
        ymin = 0, ymax = 2)
    tmp <- tempfile(fileext = ".zarr")
    on.exit(unlink(tmp, recursive = TRUE))

    writeSpatialData(sub, tmp)
    sd2 <- readSpatialData(tmp)

    n_orig <- nrow(
        spatialPoints(sd)[["transcripts"]])
    n_sub <- nrow(
        spatialPoints(sd2)[["transcripts"]])
    expect_true(n_sub < n_orig)
    expect_true(n_sub > 0L)
})

test_that("writeSpatialData validates result", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd <- readSpatialData(store)
    tmp <- tempfile(fileext = ".zarr")
    on.exit(unlink(tmp, recursive = TRUE))

    writeSpatialData(sd, tmp)
    v <- validateSpatialData(tmp)
    expect_true(v$valid)
})
