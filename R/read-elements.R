# R/read-elements.R
# Element-specific readers for Zarr arrays, CSV, and Parquet

#' @include AllClasses.R
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom jsonlite fromJSON
#' @importFrom utils read.csv
NULL

# ---- Zarr v3 low-level helpers (internal) ----

#' Decompress raw bytes with zstd via arrow
#'
#' Uses the \pkg{arrow} package's Codec interface to decompress
#' zstd-compressed raw bytes.
#'
#' @param raw_bytes Raw vector of compressed data.
#' @param max_size Integer. Maximum decompressed size to read.
#' @return Raw vector of decompressed data.
#' @keywords internal
.zstdDecompress <- function(raw_bytes, max_size = 100000000L) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
        stop(
            "Install 'arrow' to read zstd-compressed ",
            "Zarr v3 arrays:\n",
            "  install.packages('arrow')",
            call. = FALSE
        )
    }
    codec <- arrow::Codec$create("zstd")
    buf <- arrow::buffer(raw_bytes)
    reader <- arrow::BufferReader$create(buf)
    stream <- arrow::CompressedInputStream$create(
        reader, codec
    )
    out <- stream$Read(max_size)
    as.raw(out)
}

#' Parse vlen-utf8 encoded raw bytes into character vector
#'
#' The vlen-utf8 codec stores: 4-byte little-endian count,
#' then for each string: 4-byte length + UTF-8 bytes.
#'
#' @param raw_bytes Raw vector of decompressed vlen-utf8 data.
#' @return Character vector.
#' @keywords internal
.parseVlenUtf8 <- function(raw_bytes) {
    n <- readBin(raw_bytes[1:4], "integer",
        n = 1L, size = 4L, endian = "little")
    if (n == 0L) return(character(0L))
    pos <- 5L
    strings <- character(n)
    for (i in seq_len(n)) {
        str_len <- readBin(
            raw_bytes[pos:(pos + 3L)], "integer",
            n = 1L, size = 4L, endian = "little"
        )
        pos <- pos + 4L
        if (str_len > 0L) {
            strings[i] <- rawToChar(
                raw_bytes[pos:(pos + str_len - 1L)]
            )
            pos <- pos + str_len
        }
    }
    strings
}

#' Get R type info for a Zarr v3 data_type string
#'
#' @param dtype Character. Zarr v3 data type string.
#' @return A list with \code{what}, \code{size}, \code{signed}.
#' @keywords internal
## Convert raw int64 bytes to R numeric
## R has no native int64; we read as pairs of uint32 and combine
## Accurate for values up to 2^53 (9e15), which covers all counts
.readInt64Raw <- function(raw_data, n, endian = "little") {
    ## Read as signed int32 pairs (R doesn't support unsigned int32)
    all_ints <- readBin(raw_data, what = "integer", n = 2L * n,
                        size = 4L, signed = TRUE, endian = endian)
    lo_vals <- all_ints[seq(1L, 2L * n, by = 2L)]
    hi_vals <- all_ints[seq(2L, 2L * n, by = 2L)]
    ## Convert signed int32 to unsigned via: if < 0, add 2^32
    lo_unsigned <- ifelse(lo_vals < 0L,
                          as.numeric(lo_vals) + 4294967296,
                          as.numeric(lo_vals))
    hi_unsigned <- ifelse(hi_vals < 0L,
                          as.numeric(hi_vals) + 4294967296,
                          as.numeric(hi_vals))
    ## Combine: value = hi * 2^32 + lo
    result <- hi_unsigned * 4294967296 + lo_unsigned
    result
}

.zarrV3TypeInfo <- function(dtype) {
    switch(dtype,
        int8    = list(what = "integer", size = 1L,
            signed = TRUE),
        int16   = list(what = "integer", size = 2L,
            signed = TRUE),
        int32   = list(what = "integer", size = 4L,
            signed = TRUE),
        int64   = list(what = "integer", size = 8L,
            signed = TRUE, is_int64 = TRUE),
        uint8   = list(what = "integer", size = 1L,
            signed = FALSE),
        uint16  = list(what = "integer", size = 2L,
            signed = FALSE),
        float32 = list(what = "double",  size = 4L,
            signed = TRUE),
        float64 = list(what = "double",  size = 8L,
            signed = TRUE),
        stop("Unsupported Zarr v3 data type: ", dtype,
            call. = FALSE)
    )
}

#' Read a Zarr v3 array from disk
#'
#' Reads a Zarr v3 array directory containing \code{zarr.json}
#' metadata and chunk files under \code{c/}. Supports numeric
#' data types with \code{bytes} + \code{zstd} codecs, and
#' string arrays with \code{vlen-utf8} + \code{zstd} codecs.
#'
#' @param zarr_path Path to a Zarr v3 array directory
#'   (containing \code{zarr.json}).
#' @return A numeric vector, matrix, or character vector
#'   depending on the array's data type and shape.
#'
#' @keywords internal
readZarrV3Array <- function(zarr_path) {
    zarr_json <- file.path(zarr_path, "zarr.json")
    meta <- jsonlite::fromJSON(zarr_json,
        simplifyVector = FALSE)

    shape <- as.integer(unlist(meta[["shape"]]))
    dtype <- meta[["data_type"]]
    fill_value <- meta[["fill_value"]]
    chunk_shape <- as.integer(unlist(
        meta[["chunk_grid"]][["configuration"]][["chunk_shape"]]
    ))
    codecs <- meta[["codecs"]]

    ## Determine codec pipeline
    codec_names <- vapply(codecs, `[[`, character(1L),
        "name")
    has_zstd <- "zstd" %in% codec_names
    is_vlen <- "vlen-utf8" %in% codec_names
    is_string <- dtype == "string" || is_vlen

    ## Determine endianness
    endian <- "little"
    bytes_codec <- codecs[codec_names == "bytes"]
    if (length(bytes_codec) > 0L) {
        cfg <- bytes_codec[[1L]][["configuration"]]
        if (!is.null(cfg) && !is.null(cfg[["endian"]])) {
            endian <- cfg[["endian"]]
        }
    }

    total_elements <- prod(shape)
    if (total_elements == 0L) {
        if (is_string) return(character(0L))
        return(numeric(0L))
    }

    ## Handle scalar arrays (shape = [])
    if (length(shape) == 0L || all(shape == 0L)) {
        chunk_file <- file.path(zarr_path, "c")
        if (file.exists(chunk_file) && !dir.exists(chunk_file)) {
            raw_data <- readBin(chunk_file, "raw",
                n = file.info(chunk_file)$size)
            if (has_zstd) {
                raw_data <- .zstdDecompress(raw_data,
                    max_size = 10000L)
            }
            if (is_string) {
                return(.parseVlenUtf8(raw_data))
            }
        }
        if (is_string) return("")
        return(fill_value)
    }

    ## Calculate number of chunks per dimension
    ndim <- length(shape)
    n_chunks <- ceiling(shape / chunk_shape)

    ## Generate all chunk indices
    if (ndim == 1L) {
        chunk_indices <- as.list(seq(0L,
            n_chunks[1L] - 1L))
    } else {
        idx_lists <- lapply(seq_len(ndim), function(d) {
            seq(0L, n_chunks[d] - 1L)
        })
        chunk_indices <- as.matrix(
            do.call(expand.grid, idx_lists)
        )
    }

    ## For string arrays, concatenate
    if (is_string) {
        all_strings <- character(0L)
        if (ndim == 1L) {
            for (ci in seq_along(chunk_indices)) {
                idx <- chunk_indices[[ci]]
                chunk_file <- file.path(
                    zarr_path, "c", as.character(idx)
                )
                if (!file.exists(chunk_file)) {
                    ## Fill with empty strings
                    n_fill <- min(chunk_shape[1L],
                        shape[1L] - idx * chunk_shape[1L])
                    all_strings <- c(all_strings,
                        rep("", n_fill))
                    next
                }
                raw_data <- readBin(chunk_file, "raw",
                    n = file.info(chunk_file)$size)
                if (has_zstd) {
                    raw_data <- .zstdDecompress(raw_data)
                }
                strings <- .parseVlenUtf8(raw_data)
                all_strings <- c(all_strings, strings)
            }
        }
        ## Trim to exact shape
        if (length(all_strings) > total_elements) {
            all_strings <- all_strings[
                seq_len(total_elements)]
        }
        return(all_strings)
    }

    ## Numeric arrays
    type_info <- .zarrV3TypeInfo(dtype)
    element_bytes <- type_info$size

    ## Allocate output
    ## int64 values are stored as R numeric (double)
    result_what <- if (isTRUE(type_info$is_int64)) "numeric"
                   else type_info$what
    if (ndim == 1L) {
        result <- vector(result_what,
            total_elements)
        if (!is.null(fill_value) &&
            !identical(fill_value, 0) &&
            !identical(fill_value, 0.0)) {
            result[] <- fill_value
        }

        for (ci in seq_along(chunk_indices)) {
            idx <- chunk_indices[[ci]]
            chunk_file <- file.path(
                zarr_path, "c", as.character(idx)
            )
            start <- idx * chunk_shape[1L] + 1L
            n_elements <- min(chunk_shape[1L],
                shape[1L] - idx * chunk_shape[1L])

            if (!file.exists(chunk_file)) next

            raw_data <- readBin(chunk_file, "raw",
                n = file.info(chunk_file)$size)
            if (has_zstd) {
                expected <- as.integer(
                    n_elements * element_bytes)
                raw_data <- .zstdDecompress(raw_data,
                    max_size = expected + 1024L)
            }
            if (isTRUE(type_info$is_int64)) {
                vals <- .readInt64Raw(raw_data, n_elements,
                    endian = endian)
            } else {
                vals <- readBin(raw_data,
                    what = type_info$what,
                    n = n_elements,
                    size = type_info$size,
                    signed = type_info$signed,
                    endian = endian
                )
            }
            result[start:(start + length(vals) - 1L)] <-
                vals
        }
        return(result)

    } else {
        ## Multi-dimensional array
        result <- array(
            if (result_what == "integer") 0L else 0.0,
            dim = shape
        )

        n_chunk_rows <- nrow(chunk_indices)
        for (ri in seq_len(n_chunk_rows)) {
            idx <- as.integer(chunk_indices[ri, ])
            chunk_file <- file.path(zarr_path, "c",
                paste(idx, collapse = .Platform$file.sep)
            )

            if (!file.exists(chunk_file)) next

            ## Calculate this chunk's actual size
            actual_size <- integer(ndim)
            for (d in seq_len(ndim)) {
                actual_size[d] <- min(
                    chunk_shape[d],
                    shape[d] - idx[d] * chunk_shape[d]
                )
            }
            n_elements <- prod(actual_size)

            raw_data <- readBin(chunk_file, "raw",
                n = file.info(chunk_file)$size)
            if (has_zstd) {
                expected <- as.integer(
                    n_elements * element_bytes)
                raw_data <- .zstdDecompress(raw_data,
                    max_size = expected + 1024L)
            }
            if (isTRUE(type_info$is_int64)) {
                vals <- .readInt64Raw(raw_data, n_elements,
                    endian = endian)
            } else {
                vals <- readBin(raw_data,
                    what = type_info$what,
                    n = n_elements,
                    size = type_info$size,
                    signed = type_info$signed,
                    endian = endian
                )
            }

            ## Place values into result array
            ## Zarr uses C-order (row-major)
            ranges <- lapply(seq_len(ndim), function(d) {
                start <- idx[d] * chunk_shape[d] + 1L
                end <- start + actual_size[d] - 1L
                start:end
            })

            ## Create index call
            chunk_arr <- array(vals, dim = actual_size)
            idx_call <- lapply(ranges, identity)
            result <- do.call(`[<-`,
                c(list(result), idx_call,
                    list(value = chunk_arr)))
        }
        return(result)
    }
}

#' Read a Zarr array into memory
#'
#' Reads a Zarr array from a directory. Automatically detects
#' Zarr v2 (\code{.zarray} metadata) or Zarr v3
#' (\code{zarr.json} metadata) format. For v2, uses \pkg{Rarr}
#' or \pkg{pizzarr}. For v3, uses built-in reader with
#' \pkg{arrow} for zstd decompression.
#'
#' @param zarr_path Path to a Zarr array directory (containing
#'   \code{.zarray} or \code{zarr.json} metadata and chunk
#'   files).
#' @return A numeric vector, matrix, array, or character vector.
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

    ## Zarr v3: zarr.json with node_type = "array"
    zarr_json <- file.path(zarr_path, "zarr.json")
    if (file.exists(zarr_json)) {
        meta <- jsonlite::fromJSON(zarr_json,
            simplifyVector = FALSE)
        if (identical(meta[["node_type"]], "array") &&
            identical(meta[["zarr_format"]], 3L)) {
            return(readZarrV3Array(zarr_path))
        }
    }

    ## Zarr v2: .zarray
    if (requireNamespace("Rarr", quietly = TRUE)) {
        Rarr::read_zarr_array(zarr_path)
    } else if (requireNamespace("pizzarr", quietly = TRUE)) {
        store <- pizzarr::DirectoryStore$new(zarr_path)
        arr <- pizzarr::zarr_open(store, mode = "r")
        arr$get_item(".")$as.array()
    } else {
        stop(
            "Install 'Rarr' or 'pizzarr' to read Zarr v2 ",
            "arrays:\n  BiocManager::install('Rarr')",
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

#' Discover column order from Zarr group metadata
#'
#' Reads column-order and index information from either
#' \code{.zattrs} (Zarr v2) or \code{zarr.json} (Zarr v3).
#'
#' @param group_path Path to a Zarr group.
#' @return Character vector of column names.
#' @keywords internal
.discoverZarrColumns <- function(group_path) {
    ## Try Zarr v2: .zattrs
    zattrs <- file.path(group_path, ".zattrs")
    if (file.exists(zattrs)) {
        meta <- jsonlite::fromJSON(zattrs,
            simplifyVector = FALSE)
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
        return(col_order)
    }

    ## Try Zarr v3: zarr.json with attributes
    zarr_json <- file.path(group_path, "zarr.json")
    if (file.exists(zarr_json)) {
        meta <- jsonlite::fromJSON(zarr_json,
            simplifyVector = FALSE)
        attrs <- meta[["attributes"]]
        if (!is.null(attrs)) {
            col_order <- attrs[["column-order"]]
            idx <- attrs[["_index"]]
            if (!is.null(col_order)) {
                ## Include _index column if present
                if (!is.null(idx) &&
                    !(idx %in% col_order)) {
                    col_order <- c(idx, col_order)
                }
                return(col_order)
            }
            if (!is.null(idx)) return(list(idx))
        }
    }

    ## Fallback: scan subdirectories
    subdirs <- list.dirs(group_path,
        recursive = FALSE, full.names = FALSE)
    ## Exclude metadata files and non-array dirs
    subdirs <- subdirs[!subdirs %in%
        c(".zattrs", "zarr.json", "c")]
    if (length(subdirs) > 0L) return(subdirs)

    character()
}

#' Read a Zarr v3 categorical column
#'
#' Categorical columns in AnnData/Zarr v3 are stored as a
#' group with \code{categories/} (string array) and
#' \code{codes/} (integer array) sub-arrays.
#'
#' @param col_path Path to the categorical column directory.
#' @return A character vector (factor resolved to strings).
#' @keywords internal
.readZarrV3Categorical <- function(col_path) {
    cats_path <- file.path(col_path, "categories")
    codes_path <- file.path(col_path, "codes")

    if (!dir.exists(cats_path) || !dir.exists(codes_path)) {
        return(NA)
    }

    categories <- readZarrArray(cats_path)
    codes <- readZarrArray(codes_path)

    ## Map codes to categories (0-indexed)
    categories[codes + 1L]
}

#' Read Zarr columnar arrays into a DataFrame
#' @param group_path Path to a Zarr group.
#' @param col_order Character vector of column names.
#' @return A \code{DataFrame}.
#' @keywords internal
.readZarrColumns <- function(group_path, col_order) {
    cols <- lapply(col_order, function(col) {
        col_path <- file.path(group_path, col)
        if (!dir.exists(col_path)) return(NA)

        ## Check if this is a categorical column (v3)
        zarr_json <- file.path(col_path, "zarr.json")
        if (file.exists(zarr_json)) {
            meta <- jsonlite::fromJSON(zarr_json,
                simplifyVector = FALSE)
            attrs <- meta[["attributes"]]
            enc_type <- if (!is.null(attrs)) {
                attrs[["encoding-type"]]
            } else {
                NULL
            }
            ## Categorical: group with categories + codes
            if (identical(enc_type, "categorical") ||
                identical(meta[["node_type"]], "group")) {
                cats_dir <- file.path(col_path,
                    "categories")
                if (dir.exists(cats_dir)) {
                    return(tryCatch(
                        .readZarrV3Categorical(col_path),
                        error = function(e) NA
                    ))
                }
            }
        }

        ## Regular array column
        tryCatch(as.vector(readZarrArray(col_path)),
            error = function(e) NA)
    })
    names(cols) <- col_order
    valid <- vapply(cols,
        function(x) !all(is.na(x)), logical(1))
    cols <- cols[valid]
    if (length(cols) == 0L) return(DataFrame())
    DataFrame(cols)
}

#' Read AnnData X matrix (Zarr, CSV, or sparse)
#' @param x_path Path to the X matrix directory or CSV file.
#' @return A numeric matrix, sparse matrix, or \code{NULL}.
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

    ## Directory — check for Zarr v3 sparse matrix
    if (dir.exists(x_path)) {
        ## Check for CSR/CSC sparse matrix (v3)
        zarr_json <- file.path(x_path, "zarr.json")
        if (file.exists(zarr_json)) {
            meta <- jsonlite::fromJSON(zarr_json,
                simplifyVector = FALSE)
            attrs <- meta[["attributes"]]
            enc_type <- if (!is.null(attrs)) {
                attrs[["encoding-type"]]
            } else {
                NULL
            }

            if (identical(enc_type, "csr_matrix") ||
                identical(enc_type, "csc_matrix")) {
                return(.readZarrV3Sparse(
                    x_path, attrs, enc_type))
            }
        }

        ## CSV files in directory
        csv <- list.files(x_path, pattern = "[.]csv$",
            full.names = TRUE)
        if (length(csv) > 0L) {
            return(as.matrix(
                utils::read.csv(csv[1L], row.names = 1)
            ))
        }

        ## Try Zarr dense array
        tryCatch(
            readZarrArray(x_path),
            error = function(e) NULL
        )
    } else {
        NULL
    }
}

#' Read a Zarr v3 sparse matrix (CSR or CSC)
#'
#' Reads data, indices, and indptr sub-arrays and constructs
#' a \code{Matrix::sparseMatrix}.
#'
#' @param x_path Path to the sparse matrix Zarr group.
#' @param attrs List. Parsed attributes from zarr.json.
#' @param enc_type Character. Either \code{"csr_matrix"} or
#'   \code{"csc_matrix"}.
#' @return A sparse matrix (\code{dgCMatrix} or
#'   \code{dgRMatrix}), or \code{NULL} on failure.
#' @keywords internal
.readZarrV3Sparse <- function(x_path, attrs, enc_type) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        message("Install 'Matrix' for sparse matrix support")
        return(NULL)
    }

    shape <- as.integer(unlist(attrs[["shape"]]))
    nrow_mat <- shape[1L]
    ncol_mat <- shape[2L]

    data_arr <- tryCatch(
        readZarrArray(file.path(x_path, "data")),
        error = function(e) NULL
    )
    indices_arr <- tryCatch(
        readZarrArray(file.path(x_path, "indices")),
        error = function(e) NULL
    )
    indptr_arr <- tryCatch(
        readZarrArray(file.path(x_path, "indptr")),
        error = function(e) NULL
    )

    if (is.null(data_arr) || is.null(indices_arr) ||
        is.null(indptr_arr)) {
        return(NULL)
    }

    if (identical(enc_type, "csr_matrix")) {
        ## CSR: indptr indexes rows, indices are column indices
        Matrix::sparseMatrix(
            i = rep(seq_len(nrow_mat),
                diff(as.integer(indptr_arr))),
            j = as.integer(indices_arr) + 1L,
            x = as.numeric(data_arr),
            dims = c(nrow_mat, ncol_mat)
        )
    } else {
        ## CSC: indptr indexes columns, indices are row indices
        Matrix::sparseMatrix(
            i = as.integer(indices_arr) + 1L,
            j = rep(seq_len(ncol_mat),
                diff(as.integer(indptr_arr))),
            x = as.numeric(data_arr),
            dims = c(nrow_mat, ncol_mat)
        )
    }
}
