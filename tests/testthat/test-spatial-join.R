library(S4Vectors)

test_that("assignToRegions basic", {
    pts <- DataFrame(x = c(1, 2, 8), y = c(1, 2, 8))
    regions <- DataFrame(
        x = c(1.5, 8.0), y = c(1.5, 8.0),
        cell_id = c("A", "B"))
    result <- assignToRegions(pts, regions)
    expect_equal(result$assigned_region,
        c("A", "A", "B"))
})

test_that("assignToRegions max_dist", {
    pts <- DataFrame(x = c(0, 100), y = c(0, 100))
    regions <- DataFrame(
        x = 0, y = 0, cell_id = "A")
    result <- assignToRegions(pts, regions,
        max_dist = 10)
    expect_equal(result$assigned_region[1], "A")
    expect_true(is.na(result$assigned_region[2]))
})

test_that("assignToRegions errors", {
    expect_error(
        assignToRegions(
            DataFrame(a = 1), DataFrame(x = 1, y = 1)),
        "x, y")
})

# ---- spatialJoin tests ----

# Helper: create WKB Polygon from coordinate matrix
.mk_poly_wkb <- function(coords) {
    n_pts <- as.integer(nrow(coords))
    wkb <- c(
        as.raw(0x01),
        writeBin(3L, raw(), size = 4, endian = "little"),
        writeBin(1L, raw(), size = 4, endian = "little"),
        writeBin(n_pts, raw(), size = 4, endian = "little"))
    for (i in seq_len(n_pts)) {
        wkb <- c(wkb,
            writeBin(coords[i, 1], raw(), size = 8,
                endian = "little"),
            writeBin(coords[i, 2], raw(), size = 8,
                endian = "little"))
    }
    wkb
}

test_that("spatialJoin basic point-in-polygon", {
    sq <- matrix(c(0, 0, 10, 0, 10, 10, 0, 10, 0, 0),
        ncol = 2, byrow = TRUE)
    pts <- DataFrame(x = c(5, 15, 5), y = c(5, 5, 15))
    regions <- data.frame(
        geometry = I(list(.mk_poly_wkb(sq))))
    result <- spatialJoin(pts, regions,
        region_names = "inside")
    expect_equal(result[1], "inside")
    expect_true(is.na(result[2]))
    expect_true(is.na(result[3]))
})

test_that("spatialJoin multiple regions", {
    sq1 <- matrix(c(0, 0, 10, 0, 10, 10, 0, 10, 0, 0),
        ncol = 2, byrow = TRUE)
    sq2 <- matrix(c(20, 20, 30, 20, 30, 30, 20, 30, 20, 20),
        ncol = 2, byrow = TRUE)
    pts <- DataFrame(x = c(5, 25, 50), y = c(5, 25, 50))
    regions <- data.frame(
        geometry = I(list(
            .mk_poly_wkb(sq1),
            .mk_poly_wkb(sq2))))
    result <- spatialJoin(pts, regions,
        region_names = c("A", "B"))
    expect_equal(result, c("A", "B", NA_character_))
})

test_that("spatialJoin auto-names regions", {
    sq <- matrix(c(0, 0, 10, 0, 10, 10, 0, 10, 0, 0),
        ncol = 2, byrow = TRUE)
    pts <- DataFrame(x = 5, y = 5)
    regions <- data.frame(
        geometry = I(list(.mk_poly_wkb(sq))))
    result <- spatialJoin(pts, regions)
    expect_equal(result, "region_1")
})

test_that("spatialJoin errors on missing columns", {
    expect_error(
        spatialJoin(
            DataFrame(a = 1, b = 2),
            data.frame(geometry = I(list(raw(21))))),
        "x, y")
})

test_that("spatialJoin errors on missing geometry", {
    expect_error(
        spatialJoin(
            DataFrame(x = 1, y = 2),
            data.frame(a = 1)),
        "geometry")
})

test_that("spatialJoin with MERFISH data", {
    skip_if_not(
        dir.exists("C:/Users/win10/merfish_scverse/data.zarr"),
        "MERFISH data not available")
    skip_if_not(
        requireNamespace("arrow", quietly = TRUE),
        "arrow not installed")

    sd2 <- readSpatialData(
        "C:/Users/win10/merfish_scverse/data.zarr")
    pts <- spatialPoints(sd2)[["single_molecule"]]
    anat <- as.data.frame(arrow::read_parquet(
        "C:/Users/win10/merfish_scverse/data.zarr/shapes/anatomical/shapes.parquet"))

    t1 <- proc.time()
    result <- spatialJoin(pts, anat,
        region_names = c("VISp_I", "VISp_II/III",
            "VISp_IV", "VISp_V", "VISp_VI", "VISp_wm"))
    elapsed <- (proc.time() - t1)[3]

    expect_equal(length(result), nrow(pts))
    expect_true(elapsed < 60)
    ## Most points should be assigned
    assigned <- sum(!is.na(result))
    expect_true(assigned > 1000000)
    ## Check that all expected region names appear
    tbl <- table(result)
    expect_true(all(
        c("VISp_I", "VISp_IV") %in% names(tbl)))
})
