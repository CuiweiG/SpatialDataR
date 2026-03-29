library(S4Vectors)

test_that("end-to-end: mini store reads all elements", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd <- readSpatialData(store)
    expect_s4_class(sd, "SpatialData")
    expect_true(length(images(sd)) >= 1L)
    expect_true(length(spatialLabels(sd)) >= 1L)
    expect_true(length(spatialPoints(sd)) >= 1L)
    expect_true(length(shapes(sd)) >= 1L)
    expect_true(length(tables(sd)) >= 1L)
})

test_that("end-to-end: points loaded as DataFrame", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd <- readSpatialData(store)
    pts <- spatialPoints(sd)[["transcripts"]]
    expect_s4_class(pts, "DataFrame")
    expect_true("x" %in% colnames(pts))
    expect_true("y" %in% colnames(pts))
    expect_true("gene" %in% colnames(pts))
    expect_true(nrow(pts) > 0L)
})

test_that("end-to-end: shapes loaded as DataFrame", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd <- readSpatialData(store)
    shp <- shapes(sd)[["cell_boundaries"]]
    expect_s4_class(shp, "DataFrame")
    expect_true("x" %in% colnames(shp))
    expect_true("radius" %in% colnames(shp))
})

test_that("end-to-end: tables loaded with obs/var", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd <- readSpatialData(store)
    tbl <- tables(sd)[["table"]]
    expect_true(is.list(tbl))
    expect_true("obs" %in% names(tbl))
    expect_true("var" %in% names(tbl))
    expect_s4_class(tbl$obs, "DataFrame")
    expect_true(nrow(tbl$obs) > 0L)
})

test_that("end-to-end: readCSVElement on points", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    csv_path <- file.path(store, "points", "transcripts",
        "transcripts.csv")
    df <- readCSVElement(csv_path)
    expect_s4_class(df, "DataFrame")
    expect_true(nrow(df) == 500L)
})

test_that("end-to-end: readSpatialTable on table", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    tbl <- readSpatialTable(
        file.path(store, "tables", "table"))
    expect_true(is.list(tbl))
    expect_s4_class(tbl$obs, "DataFrame")
    expect_s4_class(tbl$var, "DataFrame")
})

test_that("end-to-end: coordinate systems parsed", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd <- readSpatialData(store)
    cs <- coordinateSystems(sd)
    expect_true(length(cs) >= 1L)
    expect_true("global" %in% names(cs))
})

test_that("end-to-end: transform chain works", {
    ct <- CoordinateTransform("affine",
        affine = matrix(
            c(0.2125, 0, 0, 0, 0.2125, 0, 0, 0, 1),
            nrow = 3, byrow = TRUE))
    pts <- DataFrame(x = c(0, 100), y = c(0, 100))
    result <- transformCoords(pts, ct)
    expect_equal(result$x, c(0, 21.25), tolerance = 1e-6)
    expect_equal(result$y, c(0, 21.25), tolerance = 1e-6)
})

test_that("end-to-end: show works on loaded store", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd <- readSpatialData(store)
    expect_output(show(sd), "SpatialData object")
    expect_output(show(sd), "images\\(1\\)")
    expect_output(show(sd), "spatialLabels\\(1\\)")
    expect_output(show(sd), "spatialPoints\\(1\\)")
    expect_output(show(sd), "shapes\\(1\\)")
    expect_output(show(sd), "tables\\(1\\)")
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
