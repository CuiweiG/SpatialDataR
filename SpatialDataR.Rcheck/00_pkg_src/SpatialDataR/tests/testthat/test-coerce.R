library(S4Vectors)

test_that("elementSummary returns correct structure", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd <- readSpatialData(store)
    summary_df <- elementSummary(sd)

    expect_true(is.data.frame(summary_df))
    expect_true(nrow(summary_df) >= 5L)
    expect_true(all(c("type", "name", "loaded", "nrow") %in%
        colnames(summary_df)))
    expect_true("image" %in% summary_df$type)
    expect_true("points" %in% summary_df$type)
})

test_that("elementSummary on empty SpatialData", {
    sd <- new("SpatialData")
    summary_df <- elementSummary(sd)
    expect_equal(nrow(summary_df), 0L)
})

test_that("elementTransform extracts from image ref", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd <- readSpatialData(store)
    ct <- elementTransform(images(sd)[["morphology"]])
    expect_s4_class(ct, "CoordinateTransform")
    expect_equal(slot(ct, "type"), "affine")
})

test_that("elementTransform returns NULL for no transform", {
    elem <- list(path = "/tmp", metadata = list())
    expect_null(elementTransform(elem))
})

test_that("elementTransform on loaded DataFrame", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd <- readSpatialData(store)
    pts <- spatialPoints(sd)[["transcripts"]]
    ct <- elementTransform(pts)
    ## Points should have identity transform
    if (!is.null(ct)) {
        expect_s4_class(ct, "CoordinateTransform")
    }
})

test_that("coordinateSystemElements works", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd <- readSpatialData(store)
    cs_elems <- coordinateSystemElements(sd)
    expect_true(is.list(cs_elems))
})

test_that("[ operator subsets by name", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd <- readSpatialData(store)
    sub <- sd["transcripts"]
    expect_s4_class(sub, "SpatialData")
    expect_equal(length(spatialPoints(sub)), 1L)
    expect_equal(length(images(sub)), 0L)
})

test_that("[ operator with non-matching name returns empty", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd <- readSpatialData(store)
    sub <- sd["nonexistent"]
    expect_equal(length(images(sub)), 0L)
    expect_equal(length(spatialPoints(sub)), 0L)
})
