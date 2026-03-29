test_that("cropImage works on mini store", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")
    skip_if_not(
        requireNamespace("Rarr", quietly = TRUE) ||
        requireNamespace("pizzarr", quietly = TRUE))

    img_path <- file.path(store, "images",
        "morphology", "scale0")
    crop <- cropImage(img_path,
        xmin = 1, xmax = 10,
        ymin = 1, ymax = 10)
    expect_true(is.array(crop) || is.matrix(crop))
})
