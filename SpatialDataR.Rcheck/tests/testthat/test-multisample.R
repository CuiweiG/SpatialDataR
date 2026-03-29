test_that("combineSpatialData merges two stores", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd1 <- readSpatialData(store)
    sd2 <- readSpatialData(store)
    combined <- combineSpatialData(sd1, sd2,
        sample_ids = c("tumor", "normal"))

    expect_s4_class(combined, "SpatialData")
    expect_equal(length(combined),
        length(sd1) + length(sd2))

    nms <- names(combined)
    expect_true(any(grepl("^tumor\\.", nms)))
    expect_true(any(grepl("^normal\\.", nms)))
})

test_that("combineSpatialData auto-names samples", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd1 <- readSpatialData(store)
    combined <- combineSpatialData(sd1, sd1)
    nms <- names(combined)
    expect_true(any(grepl("^sample1\\.", nms)))
    expect_true(any(grepl("^sample2\\.", nms)))
})

test_that("combineSpatialData errors on empty", {
    expect_error(combineSpatialData(),
        "At least one")
})

test_that("combineSpatialData single object passthrough", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd <- readSpatialData(store)
    result <- combineSpatialData(sd)
    expect_equal(length(result), length(sd))
})

test_that("filterSample extracts correctly", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd1 <- readSpatialData(store)
    sd2 <- readSpatialData(store)
    combined <- combineSpatialData(sd1, sd2,
        sample_ids = c("A", "B"))

    sdA <- filterSample(combined, "A")
    expect_s4_class(sdA, "SpatialData")
    expect_equal(length(sdA), length(sd1))

    ## Names should not have prefix
    nms <- names(sdA)
    expect_false(any(grepl("^A\\.", nms)))
})

test_that("filterSample non-matching returns empty", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd <- readSpatialData(store)
    combined <- combineSpatialData(sd,
        sample_ids = c("X"))
    empty <- filterSample(combined, "Z")
    expect_equal(length(empty), 0L)
})

test_that("coordinate systems are prefixed", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd1 <- readSpatialData(store)
    combined <- combineSpatialData(sd1, sd1,
        sample_ids = c("s1", "s2"))
    cs <- coordinateSystems(combined)
    expect_true(any(grepl("^s1\\.", names(cs))))
    expect_true(any(grepl("^s2\\.", names(cs))))
})
