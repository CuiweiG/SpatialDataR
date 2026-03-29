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
