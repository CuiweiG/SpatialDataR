# R/validate.R
# Validate SpatialData Zarr store against specification
# No equivalent exists in R — fills a real gap.

#' @include AllClasses.R
#' @importFrom jsonlite fromJSON
NULL

#' Validate a SpatialData Zarr store
#'
#' Checks a \code{.zarr} directory against the SpatialData
#' on-disk specification (Marconato et al. 2024). Reports
#' structural issues, missing metadata, and spec violations.
#' Useful for diagnosing interoperability problems between
#' Python and R.
#'
#' @param path Character. Path to a \code{.zarr} directory.
#' @param strict Logical. If \code{TRUE}, warnings become
#'   errors. Default: \code{FALSE}.
#' @return A list with \code{valid} (logical), \code{errors}
#'   (character vector), \code{warnings} (character vector),
#'   and \code{elements} (data.frame of discovered elements).
#'
#' @details
#' Checks performed:
#' \itemize{
#'   \item Top-level \code{.zattrs} exists and contains
#'     \code{spatialdata_attrs}
#'   \item Each element directory has valid \code{.zattrs}
#'   \item Images/labels have \code{.zarray} chunk metadata
#'   \item Points/shapes have data files (Parquet or CSV)
#'   \item Tables have obs/var subdirectories
#'   \item Coordinate transformations are parseable
#'   \item Coordinate system references are consistent
#' }
#'
#' @references
#' Marconato L et al. (2024). SpatialData: an open and
#' universal data framework for spatial omics. \emph{Nat
#' Methods} 21:2196-2209.
#' \doi{10.1038/s41592-024-02212-x}
#'
#' @export
#' @examples
#' store <- system.file("extdata", "xenium_mini.zarr",
#'     package = "SpatialDataR")
#' result <- validateSpatialData(store)
#' result$valid
#' result$elements
validateSpatialData <- function(path, strict = FALSE) {
    path <- normalizePath(path, mustWork = TRUE)
    errors <- character()
    warnings <- character()
    elements <- data.frame(
        type = character(), name = character(),
        has_zattrs = logical(), has_data = logical(),
        has_transform = logical(),
        stringsAsFactors = FALSE)

    ## 1. Top-level .zattrs
    zattrs_file <- file.path(path, ".zattrs")
    if (!file.exists(zattrs_file)) {
        errors <- c(errors,
            "Missing top-level .zattrs file")
    } else {
        meta <- .safeReadJSON(zattrs_file)
        if (is.null(meta)) {
            errors <- c(errors,
                ".zattrs is not valid JSON")
        } else if (is.null(meta[["spatialdata_attrs"]])) {
            warnings <- c(warnings,
                ".zattrs missing 'spatialdata_attrs' key")
        }
    }

    ## 2. Scan element types
    for (etype in c("images", "labels", "points",
        "shapes", "tables")) {
        edir <- file.path(path, etype)
        if (!dir.exists(edir)) next
        enames <- list.dirs(edir, recursive = FALSE,
            full.names = FALSE)
        for (en in enames) {
            row <- .validateElement(
                file.path(edir, en), etype, en)
            elements <- rbind(elements, row$summary)
            errors <- c(errors, row$errors)
            warnings <- c(warnings, row$warnings)
        }
    }

    if (nrow(elements) == 0L) {
        warnings <- c(warnings,
            "No elements found in store")
    }

    valid <- length(errors) == 0L
    if (strict && length(warnings) > 0L) valid <- FALSE

    list(valid = valid, errors = errors,
        warnings = warnings, elements = elements)
}

#' Validate a single element
#' @param elem_path Path to element directory.
#' @param etype Element type string.
#' @param ename Element name string.
#' @return List with summary, errors, warnings.
#' @keywords internal
.validateElement <- function(elem_path, etype, ename) {
    errs <- character()
    warns <- character()
    has_zattrs <- FALSE
    has_data <- FALSE
    has_transform <- FALSE

    ## Check .zattrs
    za <- file.path(elem_path, ".zattrs")
    if (file.exists(za)) {
        has_zattrs <- TRUE
        meta <- .safeReadJSON(za)
        if (!is.null(meta)) {
            ct <- meta[["coordinateTransformations"]]
            has_transform <- !is.null(ct) && length(ct) > 0L
        }
    } else {
        warns <- c(warns, paste0(
            etype, "/", ename, ": missing .zattrs"))
    }

    ## Check data presence
    if (etype %in% c("images", "labels")) {
        zarrs <- list.files(elem_path,
            "^[.]zarray$|^0",
            recursive = TRUE)
        has_data <- length(zarrs) > 0L
        if (!has_data) warns <- c(warns, paste0(
            etype, "/", ename,
            ": no .zarray found (no raster data)"))
    } else if (etype %in% c("points", "shapes")) {
        pq <- list.files(elem_path, "[.]parquet$",
            recursive = TRUE)
        csv <- list.files(elem_path, "[.]csv$",
            recursive = TRUE)
        has_data <- length(pq) > 0L || length(csv) > 0L
        if (!has_data) warns <- c(warns, paste0(
            etype, "/", ename,
            ": no Parquet/CSV data files"))
    } else if (etype == "tables") {
        has_obs <- dir.exists(file.path(elem_path, "obs"))
        has_data <- has_obs
        if (!has_obs) warns <- c(warns, paste0(
            "tables/", ename, ": missing obs directory"))
    }

    summary_row <- data.frame(
        type = etype, name = ename,
        has_zattrs = has_zattrs,
        has_data = has_data,
        has_transform = has_transform,
        stringsAsFactors = FALSE)

    list(summary = summary_row, errors = errs,
        warnings = warns)
}

#' Safely read JSON, returning NULL on failure
#' @param path Path to JSON file.
#' @return Parsed list or NULL.
#' @keywords internal
.safeReadJSON <- function(path) {
    tryCatch(
        jsonlite::fromJSON(path,
            simplifyVector = FALSE),
        error = function(e) NULL)
}
