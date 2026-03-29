## Tests for Python-native SpatialData stores
## Set SPATIALDATAR_PYTHON_STORE to the path of a Python-created
## SpatialData .zarr store to enable these tests.
## Created by: inst/scripts/create_python_native_store.py

test_that("read Python-native store (Zarr v3)", {
    store <- Sys.getenv("SPATIALDATAR_PYTHON_STORE", unset = "")
    if (store == "" || !dir.exists(store))
        skip("Python native store not available")

    sd <- readSpatialData(store)
    expect_s4_class(sd, "SpatialData")
    expect_true(length(sd) >= 2L)
})

test_that("validate Python-native store", {
    store <- Sys.getenv("SPATIALDATAR_PYTHON_STORE", unset = "")
    if (store == "" || !dir.exists(store)) skip("Python native store not available")

    v <- validateSpatialData(store)
    expect_true(v$valid)
})

test_that("read Parquet points from native store", {
    store <- Sys.getenv("SPATIALDATAR_PYTHON_STORE", unset = "")
    if (store == "" || !dir.exists(store)) skip("Python native store not available")
    skip_if_not_installed("arrow")

    sd <- readSpatialData(store)
    pts <- spatialPoints(sd)[["transcripts"]]
    expect_s4_class(pts, "DataFrame")
    expect_equal(nrow(pts), 100000L)
    expect_true("gene" %in% colnames(pts))
})

test_that("read GeoParquet shapes from native store", {
    store <- Sys.getenv("SPATIALDATAR_PYTHON_STORE", unset = "")
    if (store == "" || !dir.exists(store)) skip("Python native store not available")
    skip_if_not_installed("arrow")

    sd <- readSpatialData(store)
    shp <- shapes(sd)[["cell_boundaries"]]
    expect_s4_class(shp, "DataFrame")
    expect_equal(nrow(shp), 150L)
})

test_that("bboxQuery on Python-native store", {
    store <- Sys.getenv("SPATIALDATAR_PYTHON_STORE", unset = "")
    if (store == "" || !dir.exists(store)) skip("Python native store not available")
    skip_if_not_installed("arrow")

    sd <- readSpatialData(store)
    pts <- spatialPoints(sd)[["transcripts"]]
    cx <- median(pts$x); cy <- median(pts$y)
    sub <- bboxQuery(pts,
        xmin = cx - 200, xmax = cx + 200,
        ymin = cy - 200, ymax = cy + 200)
    expect_true(nrow(sub) > 0L)
    expect_true(nrow(sub) < nrow(pts))
})
