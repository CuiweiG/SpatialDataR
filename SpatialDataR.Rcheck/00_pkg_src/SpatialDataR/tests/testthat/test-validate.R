test_that("validateSpatialData passes on mini store", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    result <- validateSpatialData(store)
    expect_true(result$valid)
    expect_equal(length(result$errors), 0L)
    expect_true(nrow(result$elements) >= 5L)
})

test_that("validateSpatialData detects missing .zattrs", {
    tmp <- tempfile()
    dir.create(tmp)
    result <- validateSpatialData(tmp)
    expect_false(result$valid)
    expect_true(any(grepl("Missing.*zattrs", result$errors)))
    unlink(tmp, recursive = TRUE)
})

test_that("validateSpatialData reports empty store", {
    tmp <- tempfile()
    dir.create(tmp)
    writeLines("{}", file.path(tmp, ".zattrs"))
    result <- validateSpatialData(tmp)
    expect_true(result$valid)
    expect_true(any(grepl("No elements",
        result$warnings)))
    unlink(tmp, recursive = TRUE)
})

test_that("validateSpatialData strict mode", {
    tmp <- tempfile()
    dir.create(tmp)
    writeLines("{}", file.path(tmp, ".zattrs"))
    result <- validateSpatialData(tmp, strict = TRUE)
    expect_false(result$valid)
    unlink(tmp, recursive = TRUE)
})

test_that("validateSpatialData checks element structure", {
    tmp <- tempfile()
    dir.create(file.path(tmp, "images", "he"),
        recursive = TRUE)
    writeLines("{}", file.path(tmp, ".zattrs"))
    ## No .zattrs in element
    result <- validateSpatialData(tmp)
    expect_true(any(grepl("missing .zattrs",
        result$warnings)))
    unlink(tmp, recursive = TRUE)
})

test_that("validate element reports on mini store", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    result <- validateSpatialData(store)
    elems <- result$elements
    ## All elements should have zattrs
    expect_true(all(elems$has_zattrs))
    ## Images should have data
    img_rows <- elems[elems$type == "images", ]
    expect_true(all(img_rows$has_data))
})
