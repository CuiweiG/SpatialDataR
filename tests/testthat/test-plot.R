# tests/testthat/test-plot.R

test_that("plotSpatialData returns ggplot", {
    skip_if_not(
        requireNamespace("ggplot2", quietly = TRUE),
        "ggplot2 not installed")
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    skip_if(store == "", "Test data not available")
    sd <- readSpatialData(store)
    p <- plotSpatialData(sd)
    expect_s3_class(p, "ggplot")
})

test_that("plotSpatialData with points", {
    skip_if_not(
        requireNamespace("ggplot2", quietly = TRUE),
        "ggplot2 not installed")
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    skip_if(store == "", "Test data not available")
    sd <- readSpatialData(store)
    ## Only test if points are loaded as DataFrame
    pts <- spatialPoints(sd)
    skip_if(length(pts) == 0, "No points in test data")
    pt_name <- names(pts)[1]
    skip_if(!is(pts[[pt_name]], "DataFrame"),
        "Points not loaded as DataFrame")
    p <- plotSpatialData(sd, points = pt_name)
    expect_s3_class(p, "ggplot")
})

test_that("plotSpatialData with image background", {
    skip_if_not(
        requireNamespace("ggplot2", quietly = TRUE),
        "ggplot2 not installed")
    skip_if_not(
        dir.exists("C:/Users/win10/xenium_breast/data.zarr"),
        "Xenium data not available")
    sd <- readSpatialData(
        "C:/Users/win10/xenium_breast/data.zarr",
        elements = "images")
    p <- plotSpatialData(sd,
        image = "morphology_focus")
    expect_s3_class(p, "ggplot")
})

test_that("plotSpatialData with shapes", {
    skip_if_not(
        requireNamespace("ggplot2", quietly = TRUE),
        "ggplot2 not installed")
    skip_if_not(
        requireNamespace("arrow", quietly = TRUE),
        "arrow not installed")
    skip_if_not(
        dir.exists("C:/Users/win10/merfish_scverse/data.zarr"),
        "MERFISH data not available")
    sd <- readSpatialData(
        "C:/Users/win10/merfish_scverse/data.zarr",
        elements = "shapes")
    p <- plotSpatialData(sd,
        shapes = "anatomical")
    expect_s3_class(p, "ggplot")
})

test_that("plotSpatialData errors on missing element", {
    skip_if_not(
        requireNamespace("ggplot2", quietly = TRUE),
        "ggplot2 not installed")
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    skip_if(store == "", "Test data not available")
    sd <- readSpatialData(store)
    expect_error(
        plotSpatialData(sd, points = "nonexistent"),
        "not found")
})
