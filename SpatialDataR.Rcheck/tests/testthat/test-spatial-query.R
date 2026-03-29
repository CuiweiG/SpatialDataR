library(S4Vectors)

test_that("bboxQuery on DataFrame filters correctly", {
    pts <- DataFrame(
        x = c(1, 2, 3, 4, 5),
        y = c(5, 4, 3, 2, 1),
        gene = c("A", "B", "C", "D", "E")
    )
    sub <- bboxQuery(pts, xmin = 2, xmax = 4, ymin = 2, ymax = 4)
    expect_s4_class(sub, "DataFrame")
    expect_equal(nrow(sub), 3L)
    expect_true(all(sub$x >= 2 & sub$x <= 4))
    expect_true(all(sub$y >= 2 & sub$y <= 4))
})

test_that("bboxQuery returns 0 rows for empty region", {
    pts <- DataFrame(x = c(1, 2), y = c(1, 2))
    sub <- bboxQuery(pts, xmin = 10, xmax = 20, ymin = 10, ymax = 20)
    expect_equal(nrow(sub), 0L)
})

test_that("bboxQuery errors without x/y columns", {
    df <- DataFrame(a = 1:3, b = 4:6)
    expect_error(bboxQuery(df, 0, 1, 0, 1), "x.*y")
})

test_that("bboxQuery on SpatialData works", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd <- readSpatialData(store)
    full_n <- nrow(spatialPoints(sd)[["transcripts"]])

    sub_sd <- bboxQuery(sd, xmin = 0, xmax = 1,
        ymin = 0, ymax = 1)
    expect_s4_class(sub_sd, "SpatialData")

    sub_n <- nrow(spatialPoints(sub_sd)[["transcripts"]])
    expect_true(sub_n < full_n)
    expect_true(sub_n > 0L)
})

test_that("bboxQuery preserves images as references", {
    store <- system.file("extdata", "xenium_mini.zarr",
        package = "SpatialDataR")
    if (store == "") skip("mini store not found")

    sd <- readSpatialData(store)
    sub_sd <- bboxQuery(sd, xmin = 0, xmax = 1,
        ymin = 0, ymax = 1)

    img <- images(sub_sd)[["morphology"]]
    expect_true("bbox" %in% names(img))
    expect_equal(img$bbox$xmin, 0)
})
