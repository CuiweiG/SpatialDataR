library(S4Vectors)

test_that("CoordinateTransform constructor validates type", {
    expect_s4_class(CoordinateTransform("identity"),
        "CoordinateTransform")
    expect_s4_class(CoordinateTransform("affine"),
        "CoordinateTransform")
    expect_error(CoordinateTransform("invalid"))
})

test_that("CoordinateTransform validity rejects bad type", {
    expect_error(
        new("CoordinateTransform", type = "bogus"),
        "identity.*affine"
    )
})

test_that("CoordinateTransform validity rejects non-square", {
    expect_error(
        new("CoordinateTransform",
            type = "affine",
            affine = matrix(1:6, nrow = 2)),
        "square"
    )
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

test_that("matrix method: identity returns unchanged", {
    ct <- CoordinateTransform("identity")
    coords <- matrix(c(1, 2, 3, 4), ncol = 2)
    result <- transformCoords(coords, ct)
    expect_equal(result, coords)
})

test_that("matrix method: scale works", {
    mat <- diag(3)
    mat[1, 1] <- 2.0
    mat[2, 2] <- 3.0
    ct <- CoordinateTransform("affine", affine = mat)
    coords <- matrix(c(10, 20, 5, 15), ncol = 2)
    colnames(coords) <- c("x", "y")
    result <- transformCoords(coords, ct)
    expect_equal(result[, 1], c(20, 40))
    expect_equal(result[, 2], c(15, 45))
    expect_equal(colnames(result), c("x", "y"))
})

test_that("matrix method: 3D coordinates work", {
    mat <- diag(4) * 2
    mat[4, 4] <- 1
    ct <- CoordinateTransform("affine", affine = mat)
    coords <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 3)
    colnames(coords) <- c("x", "y", "z")
    result <- transformCoords(coords, ct)
    expect_equal(result[, 1], c(2, 4))
    expect_equal(result[, 2], c(6, 8))
    expect_equal(result[, 3], c(10, 12))
})

test_that("matrix method: rejects 4-column matrix", {
    ct <- CoordinateTransform("affine")
    bad <- matrix(1:8, ncol = 4)
    expect_error(transformCoords(bad, ct), "2.*3.*columns")
})

test_that("3D identity transform with 4x4 matrix", {
    ct <- CoordinateTransform("identity", affine = diag(4))
    expect_equal(nrow(slot(ct, "affine")), 4L)
})

## composeTransforms
test_that("composeTransforms: scale then translate", {
    s <- CoordinateTransform("affine",
        affine = diag(c(2, 2, 1)),
        input_cs = "pixels", output_cs = "scaled")
    t <- CoordinateTransform("affine",
        affine = matrix(c(1,0,10, 0,1,20, 0,0,1),
            nrow = 3, byrow = TRUE),
        input_cs = "scaled", output_cs = "microns")
    combined <- composeTransforms(s, t)

    expect_equal(slot(combined, "input_cs"), "pixels")
    expect_equal(slot(combined, "output_cs"), "microns")

    pts <- DataFrame(x = c(5), y = c(10))
    result <- transformCoords(pts, combined)
    expect_equal(result$x, 20)  # 5*2 + 10
    expect_equal(result$y, 40)  # 10*2 + 20
})

test_that("composeTransforms: identity composed is same", {
    s <- CoordinateTransform("affine",
        affine = diag(c(0.5, 0.5, 1)))
    id <- CoordinateTransform("identity")
    c1 <- composeTransforms(s, id)
    c2 <- composeTransforms(id, s)
    expect_equal(slot(c1, "affine"), slot(s, "affine"))
    expect_equal(slot(c2, "affine"), slot(s, "affine"))
})

test_that("composeTransforms: mixed 2D+3D pads", {
    s2 <- CoordinateTransform("affine",
        affine = diag(c(2, 2, 1)))
    s3 <- CoordinateTransform("affine",
        affine = diag(c(3, 3, 3, 1)))
    combined <- composeTransforms(s2, s3)
    expect_equal(nrow(slot(combined, "affine")), 4L)
})

## invertTransform
test_that("invertTransform: scale", {
    ct <- CoordinateTransform("affine",
        affine = diag(c(0.5, 0.5, 1)),
        input_cs = "pixels", output_cs = "microns")
    inv <- invertTransform(ct)

    expect_equal(slot(inv, "input_cs"), "microns")
    expect_equal(slot(inv, "output_cs"), "pixels")

    pts <- DataFrame(x = c(50), y = c(25))
    result <- transformCoords(pts, inv)
    expect_equal(result$x, 100)
    expect_equal(result$y, 50)
})

test_that("invertTransform: roundtrip", {
    mat <- matrix(c(0.5, 0, 10, 0, 0.5, 20, 0, 0, 1),
        nrow = 3, byrow = TRUE)
    ct <- CoordinateTransform("affine", affine = mat)
    inv <- invertTransform(ct)

    pts <- DataFrame(x = c(100), y = c(200))
    fwd <- transformCoords(pts, ct)
    back <- transformCoords(fwd, inv)
    expect_equal(back$x, 100, tolerance = 1e-10)
    expect_equal(back$y, 200, tolerance = 1e-10)
})

test_that("invertTransform: identity", {
    ct <- CoordinateTransform("identity",
        input_cs = "a", output_cs = "b")
    inv <- invertTransform(ct)
    expect_equal(slot(inv, "type"), "identity")
    expect_equal(slot(inv, "input_cs"), "b")
    expect_equal(slot(inv, "output_cs"), "a")
})

## .parseTransform
test_that(".parseTransform handles identity", {
    meta <- list(coordinateTransformations = list(
        list(type = "identity")))
    ct <- SpatialDataR:::.parseTransform(meta)
    expect_s4_class(ct, "CoordinateTransform")
    expect_equal(slot(ct, "type"), "identity")
})

test_that(".parseTransform handles affine", {
    aff <- list(c(1, 0, 10), c(0, 1, 20), c(0, 0, 1))
    meta <- list(coordinateTransformations = list(
        list(type = "affine", affine = aff)))
    ct <- SpatialDataR:::.parseTransform(meta)
    expect_s4_class(ct, "CoordinateTransform")
    expect_equal(slot(ct, "affine")[1, 3], 10)
})

test_that(".parseTransform handles scale", {
    meta <- list(coordinateTransformations = list(
        list(type = "scale", scale = c(0.2125, 0.2125))))
    ct <- SpatialDataR:::.parseTransform(meta)
    expect_s4_class(ct, "CoordinateTransform")
    expect_equal(slot(ct, "affine")[1, 1], 0.2125)
})

test_that(".parseTransform handles translation", {
    meta <- list(coordinateTransformations = list(
        list(type = "translation",
            translation = c(100, 200))))
    ct <- SpatialDataR:::.parseTransform(meta)
    expect_s4_class(ct, "CoordinateTransform")
    expect_equal(slot(ct, "affine")[1, 3], 100)
    expect_equal(slot(ct, "affine")[2, 3], 200)
})

test_that(".parseTransform handles sequence", {
    meta <- list(coordinateTransformations = list(
        list(type = "sequence", transformations = list(
            list(type = "scale", scale = c(2, 2)),
            list(type = "translation",
                translation = c(10, 20))
        ))))
    ct <- SpatialDataR:::.parseTransform(meta)
    expect_s4_class(ct, "CoordinateTransform")
    pts <- DataFrame(x = c(5), y = c(10))
    result <- transformCoords(pts, ct)
    expect_equal(result$x, 20)  # 5*2 + 10
    expect_equal(result$y, 40)  # 10*2 + 20
})

test_that(".parseTransform composes multiple entries", {
    meta <- list(coordinateTransformations = list(
        list(type = "scale", scale = c(0.5, 0.5)),
        list(type = "translation", translation = c(10, 20))
    ))
    ct <- SpatialDataR:::.parseTransform(meta)
    expect_s4_class(ct, "CoordinateTransform")
})

test_that(".parseTransform returns NULL for empty", {
    expect_null(SpatialDataR:::.parseTransform(list()))
    expect_null(SpatialDataR:::.parseTransform(
        list(coordinateTransformations = list())))
})

test_that(".parseTransform returns NULL for unknown type", {
    meta <- list(coordinateTransformations = list(
        list(type = "rotation", angle = 45)))
    expect_null(SpatialDataR:::.parseTransform(meta))
})

test_that("show method for CoordinateTransform works", {
    ct <- CoordinateTransform("affine",
        affine = diag(3) * 0.5,
        input_cs = "px", output_cs = "um")
    expect_output(show(ct), "2D")
    expect_output(show(ct), "px -> um")

    ct3 <- CoordinateTransform("affine", affine = diag(4))
    expect_output(show(ct3), "3D")
})
