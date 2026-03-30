# R/coerce-bioc.R
# Bioconductor bridge: SpatialData â†?SingleCellExperiment / SpatialExperiment

#' @include AllClasses.R
#' @include AllGenerics.R
#' @importFrom methods is slot new
#' @importFrom S4Vectors DataFrame SimpleList
NULL

#' Convert SpatialData to SingleCellExperiment
#'
#' One-step conversion from a SpatialData object to a
#' \code{SingleCellExperiment}, ready for downstream
#' Bioconductor analysis (scran, scater, etc.).
#'
#' If the table contains a stored X matrix (e.g., CSR sparse
#' from Zarr v3), it is used directly. Otherwise, if
#' \code{points_name} and \code{shapes_name} are provided,
#' \code{\link{aggregatePoints}} is called to build a count
#' matrix from transcript coordinates.
#'
#' @param sdata A \code{\linkS4class{SpatialData}} object.
#' @param table_name Character. Name of the table element to
#'   use. Default: first available table.
#' @param points_name Character. Name of the points element
#'   for building a count matrix via
#'   \code{\link{aggregatePoints}} when the table lacks an X
#'   matrix. Default: \code{NULL} (auto-detect).
#' @param shapes_name Character. Name of the shapes element
#'   for cell regions. Default: \code{NULL} (auto-detect).
#' @param feature_col Character. Column in points containing
#'   gene names. Default: \code{"feature_name"}.
#' @param region_col Character. Column shared between points
#'   and shapes for joining. Default: \code{"cell_id"}.
#'
#' @return A \code{SingleCellExperiment} object with a
#'   \code{counts} assay, \code{colData} from obs, and
#'   \code{rowData} from var.
#'
#' @references
#' Marconato L et al. (2025). SpatialData: an open and
#' universal data framework for spatial omics. \emph{Nat
#' Methods} 22:58-62.
#' \doi{10.1038/s41592-024-02212-x}
#'
#' @export
#' @examples
#' ## Requires a SpatialData Zarr store and SingleCellExperiment
#' ## sd <- readSpatialData("path/to/data.zarr")
#' ## sce <- toSingleCellExperiment(sd)
toSingleCellExperiment <- function(
    sdata,
    table_name = NULL,
    points_name = NULL,
    shapes_name = NULL,
    feature_col = "feature_name",
    region_col = "cell_id"
) {
    if (!requireNamespace("SingleCellExperiment",
            quietly = TRUE)) {
        stop(
            "Install 'SingleCellExperiment' for this ",
            "conversion:\n",
            "  BiocManager::install('SingleCellExperiment')",
            call. = FALSE
        )
    }
    if (!requireNamespace("SummarizedExperiment",
            quietly = TRUE)) {
        stop(
            "Install 'SummarizedExperiment':\n",
            "  BiocManager::install('SummarizedExperiment')",
            call. = FALSE
        )
    }

    ## Resolve table
    tbls <- tables(sdata)
    if (length(tbls) == 0L)
        stop("No tables found in SpatialData object",
            call. = FALSE)
    if (is.null(table_name)) table_name <- names(tbls)[1L]
    tbl <- tbls[[table_name]]
    if (is.null(tbl))
        stop("Table '", table_name, "' not found",
            call. = FALSE)

    ## Extract obs and var
    obs <- var_df <- DataFrame()
    x_mat <- NULL
    coords_mat <- NULL

    if (is.list(tbl) && "obs" %in% names(tbl)) {
        obs <- tbl[["obs"]]
        var_df <- tbl[["var"]]
    }

    ## Try to read X matrix from Zarr store
    table_path <- file.path(slot(sdata, "path"),
        "tables", table_name)
    x_path <- file.path(table_path, "X")
    if (dir.exists(x_path)) {
        x_mat <- tryCatch(
            .readAnnDataX(x_path),
            error = function(e) NULL
        )
    }

    ## Try to read spatial coordinates from obsm/spatial
    obsm_path <- file.path(table_path, "obsm", "spatial")
    if (dir.exists(obsm_path)) {
        coords_mat <- tryCatch(
            .readAnnDataX(obsm_path),
            error = function(e) NULL
        )
    }

    ## Fallback: build count matrix from points + shapes
    if (is.null(x_mat)) {
        x_mat <- .buildCountMatrixFromPoints(
            sdata, points_name, shapes_name,
            feature_col, region_col
        )
    }

    if (is.null(x_mat))
        stop("Could not obtain expression matrix: ",
            "no X matrix in table and no points/shapes ",
            "available for aggregation", call. = FALSE)

    ## Ensure matrix is genes x cells (features x observations)
    ## AnnData stores cells x genes, so transpose
    if (nrow(x_mat) == nrow(obs) && ncol(x_mat) != nrow(obs)) {
        x_mat <- Matrix::t(x_mat)
    }

    ## Set row/col names
    if (nrow(var_df) > 0L && nrow(var_df) == nrow(x_mat)) {
        ## Use _index column for gene names if present
        gene_names <- .extractIndex(var_df)
        if (!is.null(gene_names)) {
            rownames(x_mat) <- gene_names
        }
    }
    if (nrow(obs) > 0L && ncol(x_mat) == nrow(obs)) {
        cell_names <- .extractIndex(obs)
        if (!is.null(cell_names)) {
            colnames(x_mat) <- cell_names
        }
    }

    ## Clean up obs: drop _index column from colData
    obs <- .dropIndexCol(obs)
    var_df <- .dropIndexCol(var_df)

    ## Build SingleCellExperiment
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = x_mat),
        colData = if (nrow(obs) > 0L) obs else NULL,
        rowData = if (nrow(var_df) > 0L) var_df else NULL
    )

    ## Attach spatial coords as internal colData if available
    if (!is.null(coords_mat) &&
        nrow(coords_mat) == ncol(sce)) {
        colnames(coords_mat) <- c("x_centroid", "y_centroid")
        cd <- SummarizedExperiment::colData(sce)
        cd$x_centroid <- coords_mat[, 1L]
        cd$y_centroid <- coords_mat[, 2L]
        SummarizedExperiment::colData(sce) <- cd
    }

    sce
}

#' Convert SpatialData to SpatialExperiment
#'
#' Like \code{\link{toSingleCellExperiment}} but returns a
#' \code{SpatialExperiment} with spatial coordinates attached.
#' Coordinates are obtained from \code{obsm/spatial} in the
#' table, or computed from cell shape centroids.
#'
#' @inheritParams toSingleCellExperiment
#' @param ... Additional arguments passed to
#'   \code{\link{toSingleCellExperiment}}.
#'
#' @return A \code{SpatialExperiment} object.
#'
#' @references
#' Marconato L et al. (2025). SpatialData: an open and
#' universal data framework for spatial omics. \emph{Nat
#' Methods} 22:58-62.
#' \doi{10.1038/s41592-024-02212-x}
#'
#' @export
#' @examples
#' ## Requires a SpatialData Zarr store and SpatialExperiment
#' ## sd <- readSpatialData("path/to/data.zarr")
#' ## spe <- toSpatialExperiment(sd)
toSpatialExperiment <- function(
    sdata,
    table_name = NULL,
    points_name = NULL,
    shapes_name = NULL,
    feature_col = "feature_name",
    region_col = "cell_id",
    ...
) {
    if (!requireNamespace("SpatialExperiment",
            quietly = TRUE)) {
        stop(
            "Install 'SpatialExperiment' for this ",
            "conversion:\n",
            "  BiocManager::install('SpatialExperiment')",
            call. = FALSE
        )
    }
    if (!requireNamespace("SummarizedExperiment",
            quietly = TRUE)) {
        stop(
            "Install 'SummarizedExperiment':\n",
            "  BiocManager::install('SummarizedExperiment')",
            call. = FALSE
        )
    }

    ## Resolve table
    tbls <- tables(sdata)
    if (length(tbls) == 0L)
        stop("No tables found in SpatialData object",
            call. = FALSE)
    if (is.null(table_name)) table_name <- names(tbls)[1L]
    tbl <- tbls[[table_name]]
    if (is.null(tbl))
        stop("Table '", table_name, "' not found",
            call. = FALSE)

    ## Extract obs and var
    obs <- var_df <- DataFrame()
    x_mat <- NULL
    coords_mat <- NULL

    if (is.list(tbl) && "obs" %in% names(tbl)) {
        obs <- tbl[["obs"]]
        var_df <- tbl[["var"]]
    }

    ## Try to read X matrix from Zarr store
    table_path <- file.path(slot(sdata, "path"),
        "tables", table_name)
    x_path <- file.path(table_path, "X")
    if (dir.exists(x_path)) {
        x_mat <- tryCatch(
            .readAnnDataX(x_path),
            error = function(e) NULL
        )
    }

    ## Try to read spatial coordinates from obsm/spatial
    obsm_path <- file.path(table_path, "obsm", "spatial")
    if (dir.exists(obsm_path)) {
        coords_mat <- tryCatch(
            .readAnnDataX(obsm_path),
            error = function(e) NULL
        )
    }

    ## Fallback: build from points + shapes
    if (is.null(x_mat)) {
        x_mat <- .buildCountMatrixFromPoints(
            sdata, points_name, shapes_name,
            feature_col, region_col
        )
    }

    if (is.null(x_mat))
        stop("Could not obtain expression matrix",
            call. = FALSE)

    ## Transpose: AnnData is cells x genes, we need genes x cells
    if (nrow(x_mat) == nrow(obs) && ncol(x_mat) != nrow(obs)) {
        x_mat <- Matrix::t(x_mat)
    }

    ## Set row/col names
    gene_names <- NULL
    cell_names <- NULL
    if (nrow(var_df) > 0L && nrow(var_df) == nrow(x_mat)) {
        gene_names <- .extractIndex(var_df)
        if (!is.null(gene_names))
            rownames(x_mat) <- gene_names
    }
    if (nrow(obs) > 0L && ncol(x_mat) == nrow(obs)) {
        cell_names <- .extractIndex(obs)
        if (!is.null(cell_names))
            colnames(x_mat) <- cell_names
    }

    ## Clean up colData
    obs <- .dropIndexCol(obs)
    var_df <- .dropIndexCol(var_df)

    ## If no coords from obsm, try computing from shapes
    if (is.null(coords_mat)) {
        coords_mat <- .getCoordsFromShapes(
            sdata, shapes_name, ncol(x_mat)
        )
    }

    ## Build spatial coordinates matrix
    sp_coords <- NULL
    if (!is.null(coords_mat) &&
        nrow(coords_mat) == ncol(x_mat)) {
        sp_coords <- as.matrix(coords_mat)
        colnames(sp_coords) <- c("x", "y")
    }

    ## Build SpatialExperiment
    spe_args <- list(
        assays = list(counts = x_mat),
        colData = if (nrow(obs) > 0L) obs else NULL,
        rowData = if (nrow(var_df) > 0L) var_df else NULL
    )
    if (!is.null(sp_coords)) {
        spe_args$spatialCoords <- sp_coords
    }

    do.call(SpatialExperiment::SpatialExperiment, spe_args)
}

# ---- Internal helpers ----

#' Extract index column values from a DataFrame
#' @param df A DataFrame that may have an X_index or _index column.
#' @return Character vector or NULL.
#' @keywords internal
.extractIndex <- function(df) {
    if ("X_index" %in% colnames(df)) {
        vals <- as.character(df[["X_index"]])
        if (length(vals) > 0L && !all(is.na(vals)))
            return(vals)
    }
    if ("_index" %in% colnames(df)) {
        vals <- as.character(df[["_index"]])
        if (length(vals) > 0L && !all(is.na(vals)))
            return(vals)
    }
    NULL
}

#' Remove index columns from a DataFrame for colData/rowData
#' @param df A DataFrame.
#' @return DataFrame with index columns removed.
#' @keywords internal
.dropIndexCol <- function(df) {
    drop <- intersect(colnames(df), c("X_index", "_index"))
    if (length(drop) > 0L && nrow(df) > 0L) {
        keep <- setdiff(colnames(df), drop)
        if (length(keep) > 0L) {
            df <- df[, keep, drop = FALSE]
        }
    }
    df
}

#' Build count matrix from points and shapes via aggregation
#' @param sdata SpatialData object.
#' @param points_name Name of points element (or NULL for auto).
#' @param shapes_name Name of shapes element (or NULL for auto).
#' @param feature_col Feature column name.
#' @param region_col Region column name.
#' @return A sparse or dense matrix, or NULL.
#' @keywords internal
.buildCountMatrixFromPoints <- function(
    sdata, points_name, shapes_name,
    feature_col, region_col
) {
    pts_list <- spatialPoints(sdata)
    shp_list <- shapes(sdata)

    if (length(pts_list) == 0L || length(shp_list) == 0L)
        return(NULL)

    ## Auto-detect names
    if (is.null(points_name)) points_name <- names(pts_list)[1L]
    if (is.null(shapes_name)) shapes_name <- names(shp_list)[1L]

    pts <- pts_list[[points_name]]
    shp <- shp_list[[shapes_name]]

    if (!is(pts, "DataFrame") || !is(shp, "DataFrame"))
        return(NULL)

    ## Check required columns
    if (!feature_col %in% colnames(pts)) return(NULL)
    if (!region_col %in% colnames(pts)) return(NULL)

    ## Use aggregatePoints
    counts_df <- tryCatch(
        aggregatePoints(pts, shp,
            feature_col = feature_col,
            region_col = region_col),
        error = function(e) NULL
    )
    if (is.null(counts_df)) return(NULL)

    ## Convert to matrix
    region_ids <- counts_df[[region_col]]
    mat_cols <- setdiff(colnames(counts_df), region_col)
    mat <- as.matrix(as.data.frame(counts_df[, mat_cols]))
    rownames(mat) <- as.character(region_ids)
    ## Return as cells x genes (will be transposed later)
    mat
}

#' Get spatial coordinates from shape centroids
#' @param sdata SpatialData object.
#' @param shapes_name Name of shapes element (or NULL for auto).
#' @param n_cells Expected number of cells.
#' @return A matrix with x, y columns or NULL.
#' @keywords internal
.getCoordsFromShapes <- function(sdata, shapes_name, n_cells) {
    shp_list <- shapes(sdata)
    if (length(shp_list) == 0L) return(NULL)

    if (is.null(shapes_name)) {
        ## Try to find shapes with matching row count
        for (nm in names(shp_list)) {
            s <- shp_list[[nm]]
            if (is(s, "DataFrame") && nrow(s) == n_cells &&
                "geometry" %in% colnames(s)) {
                shapes_name <- nm
                break
            }
        }
        if (is.null(shapes_name)) return(NULL)
    }

    shp <- shp_list[[shapes_name]]
    if (!is(shp, "DataFrame") ||
        !"geometry" %in% colnames(shp))
        return(NULL)

    geom <- shp[["geometry"]]
    if (!is.list(geom)) return(NULL)

    tryCatch(
        as.matrix(geometryCentroids(geom)),
        error = function(e) NULL
    )
}
