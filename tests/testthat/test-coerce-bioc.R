# tests/testthat/test-coerce-bioc.R
# Tests for Bioconductor bridge and spatial statistics

library(S4Vectors)

# ---- Bioconductor bridge tests (unit) ----

test_that(".extractIndex finds X_index column", {
    df <- DataFrame(X_index = c("a", "b"), val = 1:2)
    expect_equal(.extractIndex(df), c("a", "b"))
})

test_that(".extractIndex returns NULL when no index", {
    df <- DataFrame(val = 1:2)
    expect_null(.extractIndex(df))
})

test_that(".dropIndexCol removes X_index", {
    df <- DataFrame(X_index = c("a", "b"), val = 1:2)
    result <- .dropIndexCol(df)
    expect_false("X_index" %in% colnames(result))
    expect_true("val" %in% colnames(result))
})

test_that("toSingleCellExperiment errors without package", {
    skip_if(requireNamespace("SingleCellExperiment",
        quietly = TRUE))
    sd <- new("SpatialData")
    expect_error(toSingleCellExperiment(sd),
        "Install.*SingleCellExperiment")
})

test_that("toSingleCellExperiment errors with no tables", {
    skip_if_not_installed("SingleCellExperiment")
    skip_if_not_installed("SummarizedExperiment")
    sd <- new("SpatialData")
    expect_error(toSingleCellExperiment(sd),
        "No tables found")
})

test_that("toSpatialExperiment errors without package", {
    skip_if(requireNamespace("SpatialExperiment",
        quietly = TRUE))
    sd <- new("SpatialData")
    expect_error(toSpatialExperiment(sd),
        "Install.*SpatialExperiment")
})

# ---- Xenium integration tests ----

test_that("toSingleCellExperiment on Xenium data", {
    skip_if_not_installed("SingleCellExperiment")
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("Matrix")
    skip_if_not_installed("arrow")
    xenium <- "C:/Users/win10/xenium_breast/data.zarr"
    skip_if(!dir.exists(xenium),
        "Xenium breast data not available")

    sd <- readSpatialData(xenium)
    sce <- toSingleCellExperiment(sd)

    expect_s4_class(sce, "SingleCellExperiment")
    expect_true(ncol(sce) > 100000)  # 167780 cells
    expect_true(nrow(sce) > 300)     # 313 genes
    expect_true("EPCAM" %in% rownames(sce))

    cd <- SummarizedExperiment::colData(sce)
    expect_true("cell_id" %in% colnames(cd))
    expect_true("x_centroid" %in% colnames(cd))
    expect_true("y_centroid" %in% colnames(cd))
})

test_that("toSpatialExperiment on Xenium data", {
    skip_if_not_installed("SpatialExperiment")
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("Matrix")
    skip_if_not_installed("arrow")
    xenium <- "C:/Users/win10/xenium_breast/data.zarr"
    skip_if(!dir.exists(xenium),
        "Xenium breast data not available")

    sd <- readSpatialData(xenium)
    spe <- toSpatialExperiment(sd)

    expect_s4_class(spe, "SpatialExperiment")
    expect_true(ncol(spe) > 100000)
    expect_true(nrow(spe) > 300)

    sc <- SpatialExperiment::spatialCoords(spe)
    expect_equal(nrow(sc), ncol(spe))
    expect_equal(ncol(sc), 2L)
    expect_true(all(c("x", "y") %in% colnames(sc)))
})

# ---- Spatial statistics tests (unit) ----

test_that("spatialAutocorrelation detects spatial signal", {
    skip_if_not_installed("FNN")
    set.seed(42)
    n <- 200
    coords <- data.frame(x = runif(n), y = runif(n))
    expr <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
    colnames(expr) <- c("spatial", "random1", "random2")
    ## Add strong spatial signal
    expr[, 1] <- coords$x * 5 + rnorm(n, sd = 0.3)

    result <- spatialAutocorrelation(expr, coords)
    expect_equal(nrow(result), 3L)
    expect_true("gene" %in% names(result))
    expect_true("observed" %in% names(result))
    expect_true("p.value" %in% names(result))
    expect_true("adjusted.p" %in% names(result))

    ## Spatial gene should be significant
    sp_p <- result$p.value[result$gene == "spatial"]
    expect_true(sp_p < 0.01)

    ## Spatial gene I should be higher than random genes
    sp_I <- result$observed[result$gene == "spatial"]
    rand_I <- mean(result$observed[result$gene != "spatial"])
    expect_true(sp_I > rand_I)
})

test_that("spatialAutocorrelation works with gene subset", {
    skip_if_not_installed("FNN")
    set.seed(42)
    n <- 100
    coords <- data.frame(x = runif(n), y = runif(n))
    expr <- matrix(rnorm(n * 5), nrow = n, ncol = 5)
    colnames(expr) <- paste0("g", 1:5)

    result <- spatialAutocorrelation(expr, coords,
        genes = c("g1", "g3"))
    expect_equal(nrow(result), 2L)
    expect_equal(result$gene, c("g1", "g3"))
})

test_that("spatialAutocorrelation errors on mismatched dims", {
    skip_if_not_installed("FNN")
    expr <- matrix(rnorm(20), nrow = 10)
    coords <- data.frame(x = runif(5), y = runif(5))
    expect_error(spatialAutocorrelation(expr, coords),
        "expr_mat rows")
})

test_that("spatialDiffExpression detects DE genes", {
    set.seed(42)
    n <- 200
    expr <- matrix(rpois(n * 5, lambda = 5),
        nrow = n, ncol = 5)
    colnames(expr) <- paste0("gene", 1:5)
    ## Make gene1 differentially expressed
    expr[1:100, 1] <- rpois(100, lambda = 20)

    result <- spatialDiffExpression(expr,
        group1 = 1:100, group2 = 101:200)

    expect_equal(nrow(result), 5L)
    expect_true("log2FC" %in% names(result))
    expect_true("adjusted.p" %in% names(result))

    ## gene1 should be significant
    de_p <- result$adjusted.p[result$gene == "gene1"]
    expect_true(de_p < 0.01)

    ## gene1 should have positive log2FC
    de_lfc <- result$log2FC[result$gene == "gene1"]
    expect_true(de_lfc > 1)
})

test_that("spatialDiffExpression works with t.test", {
    set.seed(42)
    n <- 100
    expr <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
    colnames(expr) <- paste0("g", 1:3)

    result <- spatialDiffExpression(expr,
        group1 = 1:50, group2 = 51:100, method = "t.test")
    expect_equal(nrow(result), 3L)
})

test_that("toPointPattern creates ppp object", {
    skip_if_not_installed("spatstat.geom")
    pts <- data.frame(
        x = runif(100), y = runif(100),
        gene = sample(c("A", "B"), 100, TRUE)
    )
    pp <- toPointPattern(pts, marks = "gene")
    expect_true(inherits(pp, "ppp"))
    expect_equal(pp$n, 100L)
})

test_that("toPointPattern works without marks", {
    skip_if_not_installed("spatstat.geom")
    pts <- data.frame(x = runif(50), y = runif(50))
    pp <- toPointPattern(pts)
    expect_true(inherits(pp, "ppp"))
    expect_equal(pp$n, 50L)
})

test_that("toPointPattern errors without x/y", {
    skip_if_not_installed("spatstat.geom")
    pts <- data.frame(a = 1:10, b = 1:10)
    expect_error(toPointPattern(pts), "x.*y")
})

# ---- Xenium Moran's I integration test ----

test_that("Moran's I on Xenium subsampled data", {
    skip_if_not_installed("SingleCellExperiment")
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("Matrix")
    skip_if_not_installed("FNN")
    skip_if_not_installed("arrow")
    xenium <- "C:/Users/win10/xenium_breast/data.zarr"
    skip_if(!dir.exists(xenium),
        "Xenium breast data not available")

    sd <- readSpatialData(xenium)
    sce <- toSingleCellExperiment(sd)

    cd <- SummarizedExperiment::colData(sce)
    coords <- data.frame(
        x = cd[["x_centroid"]],
        y = cd[["y_centroid"]]
    )

    set.seed(42)
    idx <- sample(ncol(sce), 2000)
    sub_expr <- t(as.matrix(
        SummarizedExperiment::assay(sce[, idx], "counts")
    ))
    sub_coords <- coords[idx, ]

    morans <- spatialAutocorrelation(sub_expr, sub_coords)

    ## Known spatially variable genes
    expect_true(
        morans$p.value[morans$gene == "EPCAM"] < 0.05)
    expect_true(
        morans$p.value[morans$gene == "LUM"] < 0.05)
    expect_true(
        morans$p.value[morans$gene == "PTPRC"] < 0.05)
})

# ---- spatstat integration on Xenium ----

test_that("toPointPattern on Xenium transcripts", {
    skip_if_not_installed("spatstat.geom")
    skip_if_not_installed("arrow")
    xenium <- "C:/Users/win10/xenium_breast/data.zarr"
    skip_if(!dir.exists(xenium),
        "Xenium breast data not available")

    sd <- readSpatialData(xenium)
    pts <- spatialPoints(sd)[["transcripts"]]
    set.seed(42)
    sub_pts <- pts[sample(nrow(pts), 50000), ]
    pp <- toPointPattern(sub_pts, marks = "feature_name")

    expect_true(inherits(pp, "ppp"))
    expect_equal(pp$n, 50000L)
})
