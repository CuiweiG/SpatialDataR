test_that("readSpatialData discovers all element types", {
    tmp <- tempfile()
    for (d in c("images/img1", "labels/seg1", "points/tx1",
                "shapes/poly1", "tables/tbl1")) {
        dir.create(file.path(tmp, d), recursive = TRUE)
        writeLines("{}", file.path(tmp, d, ".zattrs"))
    }
    writeLines('{"spatialdata_attrs": {"version": "0.1"}}',
               file.path(tmp, ".zattrs"))

    sd <- readSpatialData(tmp)
    expect_equal(length(images(sd)), 1L)
    expect_equal(length(spatialLabels(sd)), 1L)
    expect_equal(length(spatialPoints(sd)), 1L)
    expect_equal(length(shapes(sd)), 1L)
    expect_equal(length(tables(sd)), 1L)

    unlink(tmp, recursive = TRUE)
})

test_that("readSpatialData selective elements", {
    tmp <- tempfile()
    dir.create(file.path(tmp, "images", "img1"), recursive = TRUE)
    dir.create(file.path(tmp, "points", "tx1"), recursive = TRUE)
    writeLines("{}", file.path(tmp, "images", "img1", ".zattrs"))
    writeLines("{}", file.path(tmp, "points", "tx1", ".zattrs"))
    writeLines("{}", file.path(tmp, ".zattrs"))

    sd <- readSpatialData(tmp, elements = "images")
    expect_equal(length(images(sd)), 1L)
    expect_equal(length(spatialPoints(sd)), 0L)

    unlink(tmp, recursive = TRUE)
})

test_that("readSpatialData stores path", {
    tmp <- tempfile()
    dir.create(tmp)
    writeLines("{}", file.path(tmp, ".zattrs"))

    sd <- readSpatialData(tmp)
    expect_true(nchar(slot(sd, "path")) > 0)

    unlink(tmp, recursive = TRUE)
})

test_that("readSpatialData parses metadata", {
    tmp <- tempfile()
    dir.create(tmp)
    writeLines('{"spatialdata_attrs": {"version": "0.2"}}',
               file.path(tmp, ".zattrs"))

    sd <- readSpatialData(tmp)
    meta <- slot(sd, "metadata")
    expect_true("spatialdata_attrs" %in% names(meta))

    unlink(tmp, recursive = TRUE)
})

test_that("readSpatialData handles multiple elements per type", {
    tmp <- tempfile()
    for (d in c("images/he", "images/dapi", "images/morphology")) {
        dir.create(file.path(tmp, d), recursive = TRUE)
        writeLines("{}", file.path(tmp, d, ".zattrs"))
    }
    writeLines("{}", file.path(tmp, ".zattrs"))

    sd <- readSpatialData(tmp)
    expect_equal(length(images(sd)), 3L)
    expect_true(all(c("he", "dapi", "morphology") %in%
                        names(images(sd))))

    unlink(tmp, recursive = TRUE)
})

test_that("readSpatialData errors on non-existent path", {
    expect_error(readSpatialData("/nonexistent/path"))
})

test_that(".readZarrElements returns empty for missing dir", {
    result <- SpatialDataR:::.readZarrElements(
        "/nonexistent", "image")
    expect_equal(length(result), 0L)
})

test_that("readZarrArray errors without backend", {
    skip_if(requireNamespace("Rarr", quietly = TRUE) ||
            requireNamespace("pizzarr", quietly = TRUE))
    tmp <- tempfile()
    dir.create(tmp)
    expect_error(readZarrArray(tmp), "Install")
    unlink(tmp, recursive = TRUE)
})

test_that("readParquetPoints errors without arrow", {
    skip_if(requireNamespace("arrow", quietly = TRUE))
    expect_error(readParquetPoints(tempdir()), "Install")
})
