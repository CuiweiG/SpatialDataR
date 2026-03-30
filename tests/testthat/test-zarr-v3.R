# Tests for Zarr v3 reading support
# Requires real Xenium data at C:/Users/win10/xenium_breast/data.zarr

xenium_path <- "C:/Users/win10/xenium_breast/data.zarr"

test_that("readZarrV3Array reads int32 column", {
    skip_if_not(dir.exists(xenium_path),
        "Xenium test data not available")
    skip_if_not(requireNamespace("arrow", quietly = TRUE),
        "arrow package required")

    col_path <- file.path(xenium_path, "tables",
        "table", "obs", "cell_id")
    vals <- readZarrArray(col_path)
    expect_is(vals, "integer")
    expect_equal(length(vals), 167780L)
    expect_equal(vals[1], 1L)
    expect_equal(vals[167780], 167780L)
})

test_that("readZarrV3Array reads float64 column", {
    skip_if_not(dir.exists(xenium_path),
        "Xenium test data not available")
    skip_if_not(requireNamespace("arrow", quietly = TRUE),
        "arrow package required")

    col_path <- file.path(xenium_path, "tables",
        "table", "obs", "cell_area")
    vals <- readZarrArray(col_path)
    expect_is(vals, "numeric")
    expect_equal(length(vals), 167780L)
    expect_true(all(vals >= 0))
})

test_that("readZarrV3Array reads string array (vlen-utf8)", {
    skip_if_not(dir.exists(xenium_path),
        "Xenium test data not available")
    skip_if_not(requireNamespace("arrow", quietly = TRUE),
        "arrow package required")

    col_path <- file.path(xenium_path, "tables",
        "table", "var", "_index")
    vals <- readZarrArray(col_path)
    expect_is(vals, "character")
    expect_equal(length(vals), 313L)
    expect_true("EPCAM" %in% vals || nchar(vals[1]) > 0)
})

test_that("readZarrV3Array reads multi-dim array", {
    skip_if_not(dir.exists(xenium_path),
        "Xenium test data not available")
    skip_if_not(requireNamespace("arrow", quietly = TRUE),
        "arrow package required")

    sp_path <- file.path(xenium_path, "tables",
        "table", "obsm", "spatial")
    vals <- readZarrArray(sp_path)
    expect_equal(dim(vals), c(167780L, 2L))
    expect_is(vals, "matrix")
})

test_that(".discoverZarrColumns works for v3 obs", {
    skip_if_not(dir.exists(xenium_path),
        "Xenium test data not available")

    obs_path <- file.path(xenium_path, "tables",
        "table", "obs")
    cols <- SpatialDataR:::.discoverZarrColumns(obs_path)
    expect_true("cell_id" %in% cols)
    expect_true("cell_area" %in% cols)
    expect_true("region" %in% cols)
})

test_that(".readZarrV3Categorical works", {
    skip_if_not(dir.exists(xenium_path),
        "Xenium test data not available")
    skip_if_not(requireNamespace("arrow", quietly = TRUE),
        "arrow package required")

    region_path <- file.path(xenium_path, "tables",
        "table", "obs", "region")
    vals <- SpatialDataR:::.readZarrV3Categorical(
        region_path)
    expect_is(vals, "character")
    expect_equal(length(vals), 167780L)
    expect_equal(unique(vals), "cell_circles")
})

test_that(".readAnnDataX reads CSR sparse matrix", {
    skip_if_not(dir.exists(xenium_path),
        "Xenium test data not available")
    skip_if_not(requireNamespace("arrow", quietly = TRUE),
        "arrow package required")
    skip_if_not(requireNamespace("Matrix", quietly = TRUE),
        "Matrix package required")

    x_path <- file.path(xenium_path, "tables",
        "table", "X")
    mat <- SpatialDataR:::.readAnnDataX(x_path)
    expect_true(inherits(mat, "sparseMatrix"))
    expect_equal(dim(mat), c(167780L, 313L))
    expect_equal(Matrix::nnzero(mat), 10604415L)
})

test_that("readSpatialData works on full v3 store", {
    skip_if_not(dir.exists(xenium_path),
        "Xenium test data not available")
    skip_if_not(requireNamespace("arrow", quietly = TRUE),
        "arrow package required")

    sd <- readSpatialData(xenium_path)
    expect_s4_class(sd, "SpatialData")

    # Check tables
    tbl <- tables(sd)
    expect_true("table" %in% names(tbl))
    expect_equal(nrow(tbl[["table"]]$obs), 167780L)
    expect_equal(nrow(tbl[["table"]]$var), 313L)

    # Check images
    img <- images(sd)
    expect_equal(length(img), 2L)

    # Check points
    pts <- spatialPoints(sd)
    expect_true("transcripts" %in% names(pts))

    # Check shapes
    sh <- shapes(sd)
    expect_true(length(sh) >= 2L)
})
