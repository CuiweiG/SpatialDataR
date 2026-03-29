# R/read-elements.R
# Element-specific readers for Zarr arrays, CSV, and Parquet

#' @include AllClasses.R
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom jsonlite fromJSON
NULL

#' Read a Zarr array into memory
#'
#' Reads a Zarr array using \pkg{Rarr} or \pkg{pizzarr} as backend.
#' The full array is loaded into memory as a numeric matrix or
#' array.
#'
#' @param zarr_path Path to a Zarr array directory (containing
#'   \code{.zarray} metadata and chunk files).
#' @return A numeric matrix or array.
#'
#' @export
#' @examples
#' store <- system.file("extdata", "xenium_mini.zarr",
#'     package = "SpatialDataR")
#' img_path <- file.path(store, "images", "morphology", "scale0")
#' if (requireNamespace("Rarr", quietly = TRUE) ||
#'     requireNamespace("pizzarr", quietly = TRUE)) {
#'     arr <- readZarrArray(img_path)
#'     dim(arr)
#' }
readZarrArray <- function(zarr_path) {
    zarr_path <- normalizePath(zarr_path, mustWork = TRUE)

    if (requireNamespace("Rarr", quietly = TRUE)) {
        Rarr::read_zarr_array(zarr_path)
    } else if (requireNamespace("pizzarr", quietly = TRUE)) {
        store <- pizzarr::DirectoryStore$new(zarr_path)
        arr <- pizzarr::zarr_open(store, mode = "r")
        arr$get_item(".")$as.array()
    } else {
        stop(
            "Install 'Rarr' or 'pizzarr' to read Zarr arrays:\n",
            "  BiocManager::install('Rarr')",
            call. = FALSE
        )
    }
}

#' Read a Parquet-backed point table
#'
#' Reads spatial points stored as Parquet files (used by
#' SpatialData for transcript coordinates).
#'
#' @param parquet_path Path to a \code{.parquet} file or a
#'   directory containing \code{.parquet} files.
#' @return A \code{DataFrame} with x, y coordinates and metadata.
#'
#' @export
#' @examples
#' ## Parquet reading requires the arrow package
#' ## readParquetPoints("path/to/transcripts/")
#' ## Create a small mock example instead
#' df <- S4Vectors::DataFrame(
#'     x = c(1.5, 2.3, 4.1),
#'     y = c(0.8, 3.2, 1.7),
#'     gene = c("EPCAM", "VIM", "KRT18")
#' )
#' df
readParquetPoints <- function(parquet_path) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
        stop(
            "Install 'arrow' to read Parquet files:\n",
            "  install.packages('arrow')",
            call. = FALSE
        )
    }

    pq_files <- list.files(
        parquet_path,
        pattern = "[.]parquet$",
        full.names = TRUE,
        recursive = TRUE
    )
    if (length(pq_files) == 0L) {
        ## Maybe a single file path
        if (file.exists(parquet_path) &&
            grepl("[.]parquet$", parquet_path)) {
            pq_files <- parquet_path
        } else {
            stop(
                "No .parquet files found in: ",
                parquet_path,
                call. = FALSE
            )
        }
    }

    dfs <- lapply(pq_files, function(f) {
        as.data.frame(arrow::read_parquet(f))
    })
    combined <- do.call(rbind, dfs)
    DataFrame(combined)
}

#' Read CSV-backed point or shape data
#'
#' Reads spatial data stored as CSV files within a SpatialData
#' element directory. This is a fallback for stores that use CSV
#' instead of Parquet for points or shapes.
#'
#' @param csv_dir Path to a directory containing \code{.csv} files,
#'   or a single \code{.csv} file path.
#' @return A \code{DataFrame}.
#'
#' @export
#' @examples
#' store <- system.file("extdata", "xenium_mini.zarr",
#'     package = "SpatialDataR")
#' pts_csv <- file.path(store, "points", "transcripts",
#'     "transcripts.csv")
#' df <- readCSVElement(pts_csv)
#' head(df)
readCSVElement <- function(csv_dir) {
    if (file.exists(csv_dir) && grepl("[.]csv$", csv_dir)) {
        csv_files <- csv_dir
    } else if (dir.exists(csv_dir)) {
        csv_files <- list.files(
            csv_dir,
            pattern = "[.]csv$",
            full.names = TRUE,
            recursive = TRUE
        )
    } else {
        stop(
            "Path does not exist: ", csv_dir,
            call. = FALSE
        )
    }

    if (length(csv_files) == 0L) {
        stop(
            "No .csv files found in: ", csv_dir,
            call. = FALSE
        )
    }

    dfs <- lapply(csv_files, function(f) {
        utils::read.csv(f, stringsAsFactors = FALSE)
    })
    combined <- do.call(rbind, dfs)
    DataFrame(combined)
}

#' Read a SpatialData annotation table
#'
#' Reads an AnnData-formatted table from a SpatialData Zarr
#' store. By default returns a list with \code{obs} and
#' \code{var} DataFrames. Optionally converts to
#' \code{SpatialExperiment} when expression data and the
#' package are available.
#'
#' @param table_path Path to a table Zarr group.
#' @param as Character. Return type: \code{"list"} (default,
#'   always works) or \code{"SpatialExperiment"} (requires
#'   the \pkg{SpatialExperiment} package and an X matrix).
#' @return A list with \code{obs} and \code{var} DataFrames
#'   (when \code{as = "list"}), or a \code{SpatialExperiment}
#'   (when \code{as = "SpatialExperiment"} and requirements
#'   are met).
#'
#' @export
#' @examples
#' store <- system.file("extdata", "xenium_mini.zarr",
#'     package = "SpatialDataR")
#' tbl <- readSpatialTable(
#'     file.path(store, "tables", "table"))
#' names(tbl)
#' head(tbl$obs)
readSpatialTable <- function(table_path,
    as = c("list", "SpatialExperiment")) {
    as <- match.arg(as)
    table_path <- normalizePath(table_path,
        mustWork = TRUE)
    obs <- .readAnnDataGroup(
        file.path(table_path, "obs"))
    var_df <- .readAnnDataGroup(
        file.path(table_path, "var"))

    if (as == "list") {
        return(list(obs = obs, var = var_df))
    }

    x_mat <- .readAnnDataX(
        file.path(table_path, "X"))
    coords <- .readAnnDataX(
        file.path(table_path, "obsm", "spatial"))

    if (!is.null(x_mat) && requireNamespace(
        "SpatialExperiment", quietly = TRUE)) {
        .buildSpatialExperiment(
            x_mat, obs, var_df, coords)
    } else {
        message("Returning list: X matrix or ",
            "SpatialExperiment not available.")
        list(obs = obs, var = var_df)
    }
}

#' Build SpatialExperiment from table components
#' @param x_mat Expression matrix.
#' @param obs Cell metadata DataFrame.
#' @param var_df Gene metadata DataFrame.
#' @param coords Spatial coordinates matrix or NULL.
#' @return A \code{SpatialExperiment} or list fallback.
#' @keywords internal
.buildSpatialExperiment <- function(x_mat, obs, var_df, coords) {
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
            list(obs = obs, var = var_df)
        }
    )
}

#' Read an AnnData obs/var group (Zarr arrays or CSV)
#' @param group_path Path to the obs or var Zarr group.
#' @return A \code{DataFrame}.
#' @keywords internal
.readAnnDataGroup <- function(group_path) {
    if (!dir.exists(group_path)) return(DataFrame())

    ## Strategy 1: CSV files present
    csv_files <- list.files(group_path,
        pattern = "[.]csv$", full.names = TRUE)
    if (length(csv_files) > 0L) {
        df <- utils::read.csv(csv_files[1L],
            stringsAsFactors = FALSE)
        return(DataFrame(df))
    }

    ## Strategy 2: Zarr-based columnar storage
    col_order <- .discoverZarrColumns(group_path)
    if (length(col_order) == 0L) return(DataFrame())
    .readZarrColumns(group_path, col_order)
}

#' Discover column order from Zarr group .zattrs
#' @param group_path Path to a Zarr group.
#' @return Character vector of column names.
#' @keywords internal
.discoverZarrColumns <- function(group_path) {
    zattrs <- file.path(group_path, ".zattrs")
    if (!file.exists(zattrs)) return(character())

    meta <- jsonlite::fromJSON(zattrs, simplifyVector = FALSE)
    col_order <- meta[["column-order"]]
    if (is.null(col_order)) {
        idx <- meta[["_index"]]
        if (!is.null(idx)) return(list(idx))
    }
    if (is.null(col_order)) {
        subdirs <- list.dirs(group_path,
            recursive = FALSE, full.names = FALSE)
        col_order <- subdirs[subdirs != ".zattrs"]
    }
    col_order
}

#' Read Zarr columnar arrays into a DataFrame
#' @param group_path Path to a Zarr group.
#' @param col_order Character vector of column names.
#' @return A \code{DataFrame}.
#' @keywords internal
.readZarrColumns <- function(group_path, col_order) {
    cols <- lapply(col_order, function(col) {
        col_path <- file.path(group_path, col)
        if (dir.exists(col_path)) {
            tryCatch(as.vector(readZarrArray(col_path)),
                error = function(e) NA)
        } else {
            NA
        }
    })
    names(cols) <- col_order
    valid <- vapply(cols,
        function(x) !all(is.na(x)), logical(1))
    cols <- cols[valid]
    if (length(cols) == 0L) return(DataFrame())
    DataFrame(cols)
}

#' Read AnnData X matrix (Zarr or CSV)
#' @param x_path Path to the X matrix directory or CSV file.
#' @return A numeric matrix or \code{NULL}.
#' @keywords internal
.readAnnDataX <- function(x_path) {
    if (!file.exists(x_path) && !dir.exists(x_path)) {
        return(NULL)
    }

    ## CSV
    if (file.exists(x_path) && grepl("[.]csv$", x_path)) {
        mat <- as.matrix(
            utils::read.csv(x_path, row.names = 1)
        )
        return(mat)
    }

    ## Directory — could be Zarr array or contain CSV
    if (dir.exists(x_path)) {
        csv <- list.files(x_path, pattern = "[.]csv$",
            full.names = TRUE)
        if (length(csv) > 0L) {
            return(as.matrix(
                utils::read.csv(csv[1L], row.names = 1)
            ))
        }
        ## Try Zarr
        tryCatch(
            readZarrArray(x_path),
            error = function(e) NULL
        )
    } else {
        NULL
    }
}
