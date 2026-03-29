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
    expect_true(nchar(slot(sd, "path")) > 0L)

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
    for (d in c("images/he", "images/dapi", "images/morph")) {
        dir.create(file.path(tmp, d), recursive = TRUE)
        writeLines("{}", file.path(tmp, d, ".zattrs"))
    }
    writeLines("{}", file.path(tmp, ".zattrs"))

    sd <- readSpatialData(tmp)
    expect_equal(length(images(sd)), 3L)
    expect_true(all(c("he", "dapi", "morph") %in%
        names(images(sd))))

    unlink(tmp, recursive = TRUE)
})

test_that("readSpatialData errors on non-existent path", {
    expect_error(readSpatialData("/nonexistent/path"))
})

test_that(".readElementRefs returns empty for missing dir", {
    result <- SpatialDataR:::.readElementRefs(
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

test_that("readCSVElement reads single file", {
    tmp <- tempfile(fileext = ".csv")
    write.csv(data.frame(x = 1:3, y = 4:6), tmp,
        row.names = FALSE)
    df <- readCSVElement(tmp)
    expect_s4_class(df, "DataFrame")
    expect_equal(nrow(df), 3L)
    unlink(tmp)
})

test_that("readCSVElement reads directory", {
    tmp <- tempfile()
    dir.create(tmp)
    write.csv(data.frame(a = 1:2), file.path(tmp, "d.csv"),
        row.names = FALSE)
    df <- readCSVElement(tmp)
    expect_equal(nrow(df), 2L)
    unlink(tmp, recursive = TRUE)
})

test_that("readCSVElement errors on empty dir", {
    tmp <- tempfile()
    dir.create(tmp)
    expect_error(readCSVElement(tmp), "No .csv")
    unlink(tmp, recursive = TRUE)
})

test_that("readCSVElement errors on nonexistent", {
    expect_error(readCSVElement("/nonexistent"), "does not exist")
})

test_that("readSpatialTable reads CSV obs/var", {
    tmp <- tempfile()
    obs_dir <- file.path(tmp, "obs")
    var_dir <- file.path(tmp, "var")
    dir.create(obs_dir, recursive = TRUE)
    dir.create(var_dir, recursive = TRUE)
    writeLines('{}', file.path(tmp, ".zattrs"))
    writeLines('{"column-order": ["id"]}',
        file.path(obs_dir, ".zattrs"))
    writeLines('{"_index": "gene"}',
        file.path(var_dir, ".zattrs"))
    write.csv(data.frame(id = 1:5, type = rep("A", 5)),
        file.path(obs_dir, "obs.csv"), row.names = FALSE)
    write.csv(data.frame(gene = c("G1", "G2")),
        file.path(var_dir, "var.csv"), row.names = FALSE)

    result <- readSpatialTable(tmp)
    expect_true(is.list(result))
    expect_s4_class(result$obs, "DataFrame")
    expect_equal(nrow(result$obs), 5L)
    expect_s4_class(result$var, "DataFrame")

    unlink(tmp, recursive = TRUE)
})
