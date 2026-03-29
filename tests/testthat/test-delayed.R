test_that("readZarrDelayed falls back gracefully", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    img_path <- file.path(store, "images",
        "morphology", "scale0")
    result <- suppressMessages(tryCatch(
        readZarrDelayed(img_path),
        error = function(e) NULL))
    if (requireNamespace("Rarr", quietly = TRUE)) {
        expect_true(!is.null(result))
    }
})

test_that("loadElement works on SpatialData", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd <- readSpatialData(store)
    result <- suppressMessages(tryCatch(
        loadElement(sd, "morphology"),
        error = function(e) NULL))
    if (requireNamespace("Rarr", quietly = TRUE)) {
        expect_true(!is.null(result))
    }
})

test_that("loadElement errors on missing", {
    sd <- new("SpatialData")
    expect_error(loadElement(sd, "nonexistent"),
        "not found")
})
