library(S4Vectors)

test_that("CoordinateTransform constructor validates type", {
    expect_s4_class(CoordinateTransform("identity"),
                    "CoordinateTransform")
    expect_s4_class(CoordinateTransform("affine"),
                    "CoordinateTransform")
    expect_error(CoordinateTransform("invalid"))
})

test_that("identity transform returns unchanged coordinates", {
    ct <- CoordinateTransform("identity")
    pts <- DataFrame(x = c(1.5, 2.5, 3.5), y = c(4.0, 5.0, 6.0))
    result <- transformCoords(pts, ct)
    expect_equal(result$x, pts$x)
    expect_equal(result$y, pts$y)
})

test_that("affine scale works correctly", {
    mat <- diag(3)
    mat[1, 1] <- 2.0
    mat[2, 2] <- 3.0
    ct <- CoordinateTransform("affine", affine = mat)
    pts <- DataFrame(x = c(10, 20), y = c(5, 15))
    result <- transformCoords(pts, ct)
    expect_equal(result$x, c(20, 40))
    expect_equal(result$y, c(15, 45))
})

test_that("affine translation works correctly", {
    mat <- diag(3)
    mat[1, 3] <- 100
    mat[2, 3] <- 200
    ct <- CoordinateTransform("affine", affine = mat)
    pts <- DataFrame(x = c(0, 10), y = c(0, 20))
    result <- transformCoords(pts, ct)
    expect_equal(result$x, c(100, 110))
    expect_equal(result$y, c(200, 220))
})

test_that("affine combined scale+translate works", {
    mat <- matrix(c(0.5, 0, 10, 0, 0.5, 20, 0, 0, 1),
                  nrow = 3, byrow = TRUE)
    ct <- CoordinateTransform("affine", affine = mat,
        input_cs = "pixels", output_cs = "microns")
    pts <- DataFrame(x = c(100, 200), y = c(50, 150))
    result <- transformCoords(pts, ct)
    expect_equal(result$x, c(60, 110))
    expect_equal(result$y, c(45, 95))
    expect_equal(slot(ct, "input_cs"), "pixels")
    expect_equal(slot(ct, "output_cs"), "microns")
})

test_that("show method for CoordinateTransform works", {
    ct <- CoordinateTransform("affine",
        affine = diag(3) * 0.5,
        input_cs = "px", output_cs = "um")
    expect_output(show(ct), "affine")
    expect_output(show(ct), "px -> um")
})

test_that(".parseTransform handles identity", {
    meta <- list(coordinateTransformations = list(
        list(type = "identity")))
    ct <- SpatialDataR:::.parseTransform(meta)
    expect_s4_class(ct, "CoordinateTransform")
    expect_equal(slot(ct, "type"), "identity")
})

test_that(".parseTransform handles scale", {
    meta <- list(coordinateTransformations = list(
        list(type = "scale", scale = c(0.2125, 0.2125))))
    ct <- SpatialDataR:::.parseTransform(meta)
    expect_s4_class(ct, "CoordinateTransform")
    expect_equal(slot(ct, "type"), "affine")
    expect_equal(slot(ct, "affine")[1, 1], 0.2125)
})

test_that(".parseTransform returns NULL for empty", {
    expect_null(SpatialDataR:::.parseTransform(list()))
    expect_null(SpatialDataR:::.parseTransform(
        list(coordinateTransformations = list())))
})
