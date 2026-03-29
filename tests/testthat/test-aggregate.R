library(S4Vectors)

test_that("aggregatePoints counts correctly", {
    pts <- DataFrame(
        x = 1:6, y = 1:6,
        gene = c("A", "B", "A", "B", "A", "A"),
        cell_id = c(1, 1, 1, 2, 2, 2)
    )
    regions <- DataFrame(cell_id = 1:2, x = 0, y = 0)
    result <- aggregatePoints(pts, regions)

    expect_s4_class(result, "DataFrame")
    expect_true("A" %in% colnames(result))
    expect_true("B" %in% colnames(result))
    ## Cell 1: 2 A, 1 B
    r1 <- as.data.frame(result)
    row1 <- r1[r1$cell_id == 1, ]
    expect_equal(row1$A, 2L)
    expect_equal(row1$B, 1L)
    ## Cell 2: 2 A, 1 B
    row2 <- r1[r1$cell_id == 2, ]
    expect_equal(row2$A, 2L)
    expect_equal(row2$B, 1L)
})

test_that("aggregatePoints errors on missing column", {
    pts <- DataFrame(x = 1, y = 1, gene = "A")
    regions <- DataFrame(cell_id = 1)
    expect_error(
        aggregatePoints(pts, regions),
        "cell_id.*not found")
})

test_that("aggregatePoints with custom columns", {
    pts <- DataFrame(
        x = 1:4, y = 1:4,
        transcript = c("X", "Y", "X", "Y"),
        region = c("a", "a", "b", "b")
    )
    regions <- DataFrame(region = c("a", "b"))
    result <- aggregatePoints(pts, regions,
        feature_col = "transcript",
        region_col = "region")
    expect_equal(nrow(result), 2L)
    expect_true("X" %in% colnames(result))
})

test_that("aggregatePoints on mini store data", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd <- readSpatialData(store)
    pts <- spatialPoints(sd)[["transcripts"]]
    shp <- shapes(sd)[["cell_boundaries"]]

    if (!"cell_id" %in% colnames(pts)) skip("no cell_id")
    if (!"cell_id" %in% colnames(shp)) skip("no cell_id")

    counts <- aggregatePoints(pts, shp)
    expect_s4_class(counts, "DataFrame")
    expect_true(ncol(counts) > 1L)
})
