# R/read-elements.R
# Element-specific readers for Zarr arrays and Parquet tables

#' @include AllClasses.R
#' @importFrom S4Vectors DataFrame SimpleList
NULL

#' Read a Zarr array lazily
#'
#' Reads a Zarr array as a \code{DelayedArray} for out-of-memory
#' access. Uses \pkg{Rarr} or \pkg{pizzarr} as backend.
#'
#' @param zarr_path Path to a Zarr array directory.
#' @return A \code{DelayedArray} or numeric matrix.
#'
#' @export
#' @examples
#' # Requires a Zarr array on disk
#' # readZarrArray("path/to/images/morphology/scale0")
readZarrArray <- function(zarr_path) {
    zarr_path <- normalizePath(zarr_path, mustWork = TRUE)

    if (requireNamespace("Rarr", quietly = TRUE)) {
        Rarr::read_zarr_array(zarr_path)
    } else if (requireNamespace("pizzarr", quietly = TRUE)) {
        store <- pizzarr::DirectoryStore$new(zarr_path)
        arr <- pizzarr::zarr_open(store, mode = "r")
        arr$get_item(".")$as.array()
    } else {
        stop("Install 'Rarr' or 'pizzarr' to read Zarr arrays: ",
             "BiocManager::install('Rarr')")
    }
}

#' Read a Parquet-backed point table
#'
#' Reads spatial points stored as Parquet files (used by
#' SpatialData for transcript coordinates).
#'
#' @param parquet_path Path to a \code{.parquet} file or directory.
#' @return A \code{DataFrame} with x, y coordinates and metadata.
#'
#' @export
readParquetPoints <- function(parquet_path) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
        stop("Install 'arrow' to read Parquet files: ",
             "install.packages('arrow')")
    }

    pq_files <- list.files(parquet_path, pattern = "[.]parquet$",
                            full.names = TRUE, recursive = TRUE)
    if (length(pq_files) == 0L) {
        ## Maybe it's a single file
        if (file.exists(parquet_path) &&
            grepl("[.]parquet$", parquet_path)) {
            pq_files <- parquet_path
        } else {
            stop("No .parquet files found in: ", parquet_path)
        }
    }

    dfs <- lapply(pq_files, function(f) {
        as.data.frame(arrow::read_parquet(f))
    })
    combined <- do.call(rbind, dfs)
    DataFrame(combined)
}

#' Convert a SpatialData table to SpatialExperiment
#'
#' Reads an AnnData-formatted table from a Zarr store and
#' converts to a \code{SpatialExperiment} object.
#'
#' @param table_path Path to a table Zarr group (AnnData format).
#' @return A \code{SpatialExperiment} object, or \code{NULL} if
#'   conversion fails.
#'
#' @export
readSpatialTable <- function(table_path) {
    table_path <- normalizePath(table_path, mustWork = TRUE)

    ## Read obs (cell metadata)
    obs_path <- file.path(table_path, "obs")
    obs <- if (dir.exists(obs_path)) {
        .readZarrDataFrame(obs_path)
    } else {
        DataFrame(row.names = character())
    }

    ## Read X (expression matrix)
    x_path <- file.path(table_path, "X")
    x_mat <- if (dir.exists(x_path)) {
        tryCatch(readZarrArray(x_path),
                 error = function(e) NULL)
    } else NULL

    ## Read var (gene metadata)
    var_path <- file.path(table_path, "var")
    var_df <- if (dir.exists(var_path)) {
        .readZarrDataFrame(var_path)
    } else NULL

    ## Read spatial coordinates from obsm
    obsm_path <- file.path(table_path, "obsm", "spatial")
    coords <- if (dir.exists(obsm_path)) {
        tryCatch(readZarrArray(obsm_path),
                 error = function(e) NULL)
    } else NULL

    ## Build SpatialExperiment
    if (!is.null(x_mat) &&
        requireNamespace("SpatialExperiment", quietly = TRUE)) {
        se_args <- list(assays = list(counts = t(x_mat)),
                        colData = obs)
        if (!is.null(coords) && nrow(coords) == nrow(obs)) {
            se_args$spatialCoords <- coords
        }
        tryCatch(
            do.call(SpatialExperiment::SpatialExperiment, se_args),
            error = function(e) {
                message("SpatialExperiment creation failed: ",
                        conditionMessage(e))
                NULL
            }
        )
    } else {
        NULL
    }
}

#' @keywords internal
.readZarrDataFrame <- function(zarr_group_path) {
    ## Read column names from .zattrs
    zattrs <- file.path(zarr_group_path, ".zattrs")
    if (!file.exists(zattrs)) return(DataFrame())

    meta <- jsonlite::fromJSON(zattrs, simplifyVector = FALSE)
    col_order <- meta[["column-order"]]
    if (is.null(col_order)) {
        col_order <- meta[["_index"]]
    }
    if (is.null(col_order)) {
        ## Try to discover columns from subdirectories
        col_order <- list.dirs(zarr_group_path, recursive = FALSE,
                               full.names = FALSE)
        col_order <- col_order[col_order != ".zattrs"]
    }

    if (length(col_order) == 0L) return(DataFrame())

    cols <- lapply(col_order, function(col) {
        col_path <- file.path(zarr_group_path, col)
        if (dir.exists(col_path)) {
            tryCatch(as.vector(readZarrArray(col_path)),
                     error = function(e) NA)
        } else NA
    })
    names(cols) <- col_order
    cols <- cols[!vapply(cols, function(x) all(is.na(x)), logical(1))]
    if (length(cols) == 0L) return(DataFrame())
    DataFrame(cols)
}
