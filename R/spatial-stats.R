# R/spatial-stats.R
# Spatial statistics: Moran's I, differential expression, spatstat bridge

#' @include AllClasses.R
#' @importFrom methods is
#' @importFrom stats pnorm p.adjust wilcox.test t.test
NULL

#' Compute Moran's I spatial autocorrelation
#'
#' Tests whether gene expression is spatially clustered,
#' dispersed, or random. Uses cell spatial coordinates and
#' a k-nearest-neighbor weight matrix.
#'
#' No Python equivalent — leverages R's statistical computing
#' strength for spatial analysis of gene expression.
#'
#' @param expr_mat Numeric matrix (cells x genes) or
#'   \code{dgCMatrix}. Rows are cells, columns are genes.
#' @param coords A \code{data.frame} with \code{x} and
#'   \code{y} columns (cell coordinates). Must have the same
#'   number of rows as \code{expr_mat}.
#' @param genes Character vector of gene names to test.
#'   Default: all columns of \code{expr_mat}.
#' @param k Integer. Number of nearest neighbors for
#'   spatial weight matrix. Default: 10.
#' @param alternative Character. \code{"greater"} (clustered),
#'   \code{"less"} (dispersed), or \code{"two.sided"}.
#'
#' @return A \code{data.frame} with columns:
#'   \describe{
#'     \item{gene}{Gene name}
#'     \item{observed}{Observed Moran's I statistic}
#'     \item{expected}{Expected I under null hypothesis}
#'     \item{sd}{Standard deviation under null}
#'     \item{p.value}{Raw p-value}
#'     \item{adjusted.p}{BH-adjusted p-value}
#'   }
#'
#' @details
#' Moran's I is computed as:
#' \deqn{I = \frac{n}{S_0} \cdot
#'   \frac{\sum_{i,j} w_{ij} z_i z_j}{\sum_i z_i^2}}
#' where \eqn{z = x - \bar{x}}, \eqn{w_{ij}} are spatial
#' weights from k-NN, and \eqn{S_0 = \sum w_{ij}}.
#'
#' The implementation is vectorized across genes for
#' efficiency. For large datasets (>5000 cells), consider
#' subsampling.
#'
#' @references
#' Moran PAP (1950). Notes on continuous stochastic phenomena.
#' \emph{Biometrika} 37:17-23.
#'
#' @export
#' @examples
#' set.seed(42)
#' n <- 200
#' coords <- data.frame(x = runif(n), y = runif(n))
#' expr <- matrix(rnorm(n * 5), nrow = n, ncol = 5)
#' colnames(expr) <- paste0("gene", 1:5)
#' ## Add spatial signal to gene1
#' expr[, 1] <- coords$x * 3 + rnorm(n, sd = 0.5)
#' result <- spatialAutocorrelation(expr, coords)
#' result
spatialAutocorrelation <- function(
    expr_mat,
    coords,
    genes = NULL,
    k = 10L,
    alternative = c("greater", "less", "two.sided")
) {
    alternative <- match.arg(alternative)
    k <- as.integer(k)

    ## Validate inputs
    if (!is.data.frame(coords))
        coords <- as.data.frame(coords)
    if (!all(c("x", "y") %in% names(coords)))
        stop("coords must have 'x' and 'y' columns",
            call. = FALSE)

    ## Convert sparse matrix to dense if needed
    if (methods::is(expr_mat, "sparseMatrix")) {
        expr_mat <- as.matrix(expr_mat)
    }
    if (!is.matrix(expr_mat))
        expr_mat <- as.matrix(expr_mat)

    n <- nrow(expr_mat)
    if (n != nrow(coords))
        stop("expr_mat rows (", n,
            ") != coords rows (", nrow(coords), ")",
            call. = FALSE)

    ## Subset genes if specified
    if (!is.null(genes)) {
        genes <- intersect(genes, colnames(expr_mat))
        if (length(genes) == 0L)
            stop("None of the specified genes found in ",
                "expr_mat", call. = FALSE)
        expr_mat <- expr_mat[, genes, drop = FALSE]
    }

    if (is.null(colnames(expr_mat)))
        colnames(expr_mat) <- paste0("V", seq_len(ncol(expr_mat)))

    ## Build k-NN spatial weight matrix
    knn <- .buildKNN(coords, k)
    nn_idx <- knn$nn.index  # n x k matrix of neighbor indices

    ## Compute S0 (sum of all weights; binary weights = n * k)
    S0 <- as.numeric(n) * as.numeric(k)

    ## Expected Moran's I
    EI <- -1.0 / (n - 1L)

    ## Vectorized Moran's I computation
    n_genes <- ncol(expr_mat)
    gene_names <- colnames(expr_mat)
    observed <- numeric(n_genes)
    sds <- numeric(n_genes)

    for (g in seq_len(n_genes)) {
        x <- expr_mat[, g]
        z <- x - mean(x)
        ss <- sum(z * z)
        if (ss == 0) {
            observed[g] <- NA_real_
            sds[g] <- NA_real_
            next
        }

        ## Compute numerator: sum_ij w_ij * z_i * z_j
        ## Using the k-NN index: for each cell i, sum z[j]
        ## for its k neighbors
        z_neighbors <- matrix(z[nn_idx], nrow = n, ncol = k)
        wz <- rowSums(z_neighbors)
        numerator <- sum(z * wz)

        observed[g] <- (n / S0) * (numerator / ss)

        ## Standard deviation under randomization assumption
        sds[g] <- .moranVariance(n, k, S0)
    }

    ## P-values
    z_scores <- (observed - EI) / sds
    p_values <- switch(alternative,
        greater = pnorm(z_scores, lower.tail = FALSE),
        less = pnorm(z_scores, lower.tail = TRUE),
        two.sided = 2 * pnorm(abs(z_scores),
            lower.tail = FALSE)
    )

    data.frame(
        gene = gene_names,
        observed = observed,
        expected = EI,
        sd = sds,
        p.value = p_values,
        adjusted.p = p.adjust(p_values, method = "BH"),
        stringsAsFactors = FALSE
    )
}

#' Build k-nearest-neighbor index
#' @param coords data.frame with x, y.
#' @param k Number of neighbors.
#' @return List with nn.index and nn.dist matrices.
#' @keywords internal
.buildKNN <- function(coords, k) {
    coord_mat <- as.matrix(coords[, c("x", "y")])
    if (requireNamespace("FNN", quietly = TRUE)) {
        FNN::get.knn(coord_mat, k = k)
    } else {
        ## Fallback: use dist() — slower but no dependency
        message("Install 'FNN' for faster k-NN computation")
        d <- as.matrix(dist(coord_mat))
        n <- nrow(d)
        nn_idx <- matrix(0L, nrow = n, ncol = k)
        nn_dist <- matrix(0.0, nrow = n, ncol = k)
        for (i in seq_len(n)) {
            ord <- order(d[i, ])
            ## Skip self (first entry)
            idx <- ord[2:(k + 1L)]
            nn_idx[i, ] <- idx
            nn_dist[i, ] <- d[i, idx]
        }
        list(nn.index = nn_idx, nn.dist = nn_dist)
    }
}

#' Compute asymptotic standard deviation of Moran's I
#' @param n Number of observations.
#' @param k Number of neighbors.
#' @param S0 Sum of weights.
#' @return Numeric standard deviation.
#' @keywords internal
.moranVariance <- function(n, k, S0) {
    ## Under normality assumption with row-standardized weights:
    ## Var(I) = [n * ((n^2 - 3n + 3) * S1 - n * S2 + 3 * S0^2)
    ##          - b2 * ((n^2 - n) * S1 - 2n * S2 + 6 * S0^2)]
    ##          / [(n - 1)(n - 2)(n - 3) * S0^2]  - EI^2
    ## For binary k-NN weights:
    ##   S0 = n * k
    ##   S1 = 2 * n * k (each pair counted twice for symmetric)
    ##   But k-NN is NOT symmetric. Use approximation.
    ## Simplified asymptotic: sd ≈ 1/sqrt(n) for large n
    ## More precise: use the permutation-based variance
    EI <- -1.0 / (n - 1)

    ## For k-NN binary weights (asymmetric):
    ## Use the randomization variance formula
    ## Under randomization: Var(I) depends on S1, S2, S0
    ## S0 = n*k (total weight)
    ## S1 = sum_ij (w_ij + w_ji)^2
    ##    For asymmetric k-NN: most pairs have w_ij + w_ji = 1,
    ##    some mutual neighbors have sum = 2
    ##    Approximate S1 ≈ 2 * n * k (conservative)
    S1 <- 2.0 * S0
    ## S2 = sum_i (sum_j w_ij + sum_j w_ji)^2
    ##    Each row sum is k (out-degree). In-degree varies.
    ##    Approximate: S2 ≈ n * (2*k)^2 = 4*n*k^2
    S2 <- 4.0 * n * k^2

    nf <- as.numeric(n)
    kf <- as.numeric(k)
    S0f <- as.numeric(S0)

    ## Randomization variance (assumes no normality)
    A <- nf * ((nf * nf - 3 * nf + 3) * S1 -
        nf * S2 + 3 * S0f * S0f)
    ## Assume b2 (kurtosis) = 3 for normal data
    B <- 3.0 * ((nf * nf - nf) * S1 -
        2 * nf * S2 + 6 * S0f * S0f)
    C <- (nf - 1) * (nf - 2) * (nf - 3) * S0f * S0f

    var_I <- (A - B) / C - EI * EI
    if (var_I <= 0) var_I <- 1.0 / nf  # fallback
    sqrt(var_I)
}

#' Spatial differential expression between regions
#'
#' Identifies genes whose expression differs significantly
#' between two groups of cells. Computes log2 fold changes,
#' p-values, and detection rates.
#'
#' @param expr_mat Numeric matrix (cells x genes) or
#'   \code{dgCMatrix}. Rows are cells, columns are genes.
#' @param group1 Integer or logical vector indexing cells
#'   in group 1.
#' @param group2 Integer or logical vector indexing cells
#'   in group 2.
#' @param method Character. Test method: \code{"wilcox"}
#'   (default) or \code{"t.test"}.
#'
#' @return A \code{data.frame} with columns:
#'   \describe{
#'     \item{gene}{Gene name}
#'     \item{log2FC}{Log2 fold change (group1 / group2)}
#'     \item{p.value}{Raw p-value}
#'     \item{adjusted.p}{BH-adjusted p-value}
#'     \item{mean1}{Mean expression in group 1}
#'     \item{mean2}{Mean expression in group 2}
#'     \item{pct1}{Detection rate in group 1}
#'     \item{pct2}{Detection rate in group 2}
#'   }
#'
#' @export
#' @examples
#' set.seed(42)
#' n <- 100
#' expr <- matrix(rpois(n * 10, lambda = 5),
#'     nrow = n, ncol = 10)
#' colnames(expr) <- paste0("gene", 1:10)
#' ## Make gene1 differentially expressed
#' expr[1:50, 1] <- rpois(50, lambda = 20)
#' result <- spatialDiffExpression(
#'     expr, group1 = 1:50, group2 = 51:100)
#' head(result)
spatialDiffExpression <- function(
    expr_mat,
    group1,
    group2,
    method = c("wilcox", "t.test")
) {
    method <- match.arg(method)

    ## Convert sparse matrix
    if (methods::is(expr_mat, "sparseMatrix"))
        expr_mat <- as.matrix(expr_mat)
    if (!is.matrix(expr_mat))
        expr_mat <- as.matrix(expr_mat)

    if (is.null(colnames(expr_mat)))
        colnames(expr_mat) <- paste0("V", seq_len(ncol(expr_mat)))

    mat1 <- expr_mat[group1, , drop = FALSE]
    mat2 <- expr_mat[group2, , drop = FALSE]

    n_genes <- ncol(expr_mat)
    gene_names <- colnames(expr_mat)

    mean1 <- colMeans(mat1)
    mean2 <- colMeans(mat2)
    pct1 <- colMeans(mat1 > 0)
    pct2 <- colMeans(mat2 > 0)

    ## Pseudocount for log2FC
    pseudo <- 1
    log2fc <- log2((mean1 + pseudo) / (mean2 + pseudo))

    ## P-values
    p_values <- numeric(n_genes)
    test_fn <- if (method == "wilcox") {
        function(x, y) {
            suppressWarnings(
                wilcox.test(x, y)$p.value
            )
        }
    } else {
        function(x, y) {
            suppressWarnings(
                t.test(x, y)$p.value
            )
        }
    }

    for (g in seq_len(n_genes)) {
        p_values[g] <- tryCatch(
            test_fn(mat1[, g], mat2[, g]),
            error = function(e) NA_real_
        )
    }

    data.frame(
        gene = gene_names,
        log2FC = log2fc,
        p.value = p_values,
        adjusted.p = p.adjust(p_values, method = "BH"),
        mean1 = mean1,
        mean2 = mean2,
        pct1 = pct1,
        pct2 = pct2,
        stringsAsFactors = FALSE
    )
}

#' Convert points to spatstat point pattern
#'
#' Creates a \code{spatstat.geom::ppp} object from spatial
#' coordinates, enabling the full spatstat toolkit: Ripley's K,
#' pair correlation, intensity estimation, etc.
#'
#' @param points A \code{data.frame} or \code{DataFrame} with
#'   \code{x} and \code{y} columns.
#' @param marks Character string naming a column to use as
#'   marks (e.g., gene names), or \code{NULL} for unmarked.
#' @param window An \code{owin} object, or \code{NULL} to use
#'   the bounding box of the points.
#'
#' @return A \code{spatstat.geom::ppp} object.
#'
#' @export
#' @examples
#' ## Requires spatstat.geom
#' pts <- data.frame(
#'     x = runif(100), y = runif(100),
#'     gene = sample(c("A", "B"), 100, TRUE))
#' if (requireNamespace("spatstat.geom", quietly = TRUE)) {
#'     pp <- toPointPattern(pts, marks = "gene")
#'     pp
#' }
toPointPattern <- function(points, marks = NULL,
    window = NULL) {
    if (!requireNamespace("spatstat.geom",
            quietly = TRUE)) {
        stop(
            "Install 'spatstat.geom' for point pattern ",
            "analysis:\n",
            "  install.packages('spatstat.geom')",
            call. = FALSE
        )
    }

    df <- as.data.frame(points)
    if (!all(c("x", "y") %in% names(df)))
        stop("points must have 'x' and 'y' columns",
            call. = FALSE)

    x <- df$x
    y <- df$y

    ## Build window
    if (is.null(window)) {
        xr <- range(x, na.rm = TRUE)
        yr <- range(y, na.rm = TRUE)
        ## Add small buffer to avoid points on boundary
        dx <- diff(xr) * 0.01
        dy <- diff(yr) * 0.01
        window <- spatstat.geom::owin(
            xrange = xr + c(-dx, dx),
            yrange = yr + c(-dy, dy)
        )
    }

    ## Build marks
    mk <- NULL
    if (!is.null(marks) && is.character(marks) &&
        length(marks) == 1L) {
        if (marks %in% names(df)) {
            mk <- factor(df[[marks]])
        }
    }

    spatstat.geom::ppp(
        x = x, y = y,
        window = window,
        marks = mk
    )
}
