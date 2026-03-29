# R/write-zarr.R
# Write SpatialData objects back to Zarr format
# Enables roundtrip: Python -> R -> Python

#' @include AllClasses.R
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom jsonlite toJSON fromJSON
#' @importFrom methods slot is
#' @importFrom utils write.csv
NULL

#' Write a SpatialData object to Zarr format
#'
#' Writes a \code{SpatialData} object to a SpatialData-formatted
#' Zarr directory. This enables roundtrip interoperability:
#' data read from Python can be modified in R and written back
#' for use in the Python spatialdata ecosystem.
#'
#' @param sd A \code{\linkS4class{SpatialData}} object.
#' @param path Character. Output path for the \code{.zarr}
#'   directory. Will be created if it does not exist.
#' @param overwrite Logical. If \code{TRUE}, overwrite existing
#'   store. Default: \code{FALSE}.
#'
#' @return Invisible \code{path}.
#'
#' @details
#' The writer creates a spec-compliant SpatialData Zarr store:
#' \itemize{
#'   \item Top-level \code{.zattrs} with \code{spatialdata_attrs}
#'     and coordinate system definitions
#'   \item Points/shapes written as CSV (compatible with
#'     SpatialData readers)
#'   \item Tables written as AnnData-format obs/var CSV
#'   \item Image/label path references are copied if accessible
#'   \item Element-level \code{.zattrs} with coordinate
#'     transformations preserved
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
#' sd <- readSpatialData(store)
#'
#' # Subset and write back
#' sub <- bboxQuery(sd, xmin = 0, xmax = 2,
#'     ymin = 0, ymax = 2)
#' tmp <- tempfile(fileext = ".zarr")
#' writeSpatialData(sub, tmp)
#'
#' # Verify roundtrip
#' sd2 <- readSpatialData(tmp)
#' sd2
#' unlink(tmp, recursive = TRUE)
writeSpatialData <- function(sd, path,
    overwrite = FALSE) {
    if (dir.exists(path) && !overwrite) {
        stop("Directory exists: ", path,
            ". Use overwrite=TRUE.", call. = FALSE)
    }
    if (dir.exists(path) && overwrite) {
        unlink(path, recursive = TRUE)
    }
    dir.create(path, recursive = TRUE)

    ## Write top-level .zattrs
    meta <- slot(sd, "metadata")
    cs <- coordinateSystems(sd)
    if (length(meta) == 0L) {
        meta <- list(spatialdata_attrs = list(
            version = "0.1"))
    }
    if (length(cs) > 0L) {
        meta$spatialdata_attrs$coordinate_systems <- cs
    }
    writeLines(
        jsonlite::toJSON(meta, auto_unbox = TRUE,
            pretty = TRUE),
        file.path(path, ".zattrs"))

    ## Write points
    .writeDataFrameElements(
        spatialPoints(sd), path, "points")

    ## Write shapes
    .writeDataFrameElements(
        shapes(sd), path, "shapes")

    ## Write tables
    .writeTableElements(tables(sd), path)

    ## Write image/label refs (copy if possible)
    .writeRefElements(images(sd), path, "images")
    .writeRefElements(spatialLabels(sd), path, "labels")

    invisible(path)
}

#' Write DataFrame elements (points/shapes)
#' @param sl SimpleList of DataFrames or refs.
#' @param base_path Zarr store root.
#' @param type_dir "points" or "shapes".
#' @return Invisible \code{NULL}.
#' @keywords internal
.writeDataFrameElements <- function(sl, base_path,
    type_dir) {
    if (length(sl) == 0L) return(invisible(NULL))
    for (nm in names(sl)) {
        elem <- sl[[nm]]
        elem_dir <- file.path(base_path, type_dir, nm)
        dir.create(elem_dir, recursive = TRUE)

        ## Write .zattrs
        meta <- NULL
        if (is(elem, "DataFrame")) {
            meta <- attr(elem, "spatialdata_metadata")
            utils::write.csv(
                as.data.frame(elem),
                file.path(elem_dir,
                    paste0(nm, ".csv")),
                row.names = FALSE)
        } else if (is.list(elem) &&
            "metadata" %in% names(elem)) {
            meta <- elem$metadata
        }
        if (is.null(meta)) meta <- list()
        if (length(meta) == 0L) {
            meta <- list(
                coordinateTransformations = list(
                    list(type = "identity")))
        }
        writeLines(
            jsonlite::toJSON(meta, auto_unbox = TRUE,
                pretty = TRUE),
            file.path(elem_dir, ".zattrs"))
    }
}

#' Write table elements
#' @param sl SimpleList of table objects.
#' @param base_path Zarr store root.
#' @return Invisible \code{NULL}.
#' @keywords internal
.writeTableElements <- function(sl, base_path) {
    if (length(sl) == 0L) return(invisible(NULL))
    for (nm in names(sl)) {
        tbl <- sl[[nm]]
        tbl_dir <- file.path(base_path, "tables", nm)

        if (is.list(tbl) && "obs" %in% names(tbl)) {
            obs_dir <- file.path(tbl_dir, "obs")
            var_dir <- file.path(tbl_dir, "var")
            dir.create(obs_dir, recursive = TRUE)
            dir.create(var_dir, recursive = TRUE)

            if (is(tbl$obs, "DataFrame") &&
                nrow(tbl$obs) > 0L) {
                utils::write.csv(
                    as.data.frame(tbl$obs),
                    file.path(obs_dir, "obs.csv"),
                    row.names = FALSE)
            }
            if (is(tbl$var, "DataFrame") &&
                nrow(tbl$var) > 0L) {
                utils::write.csv(
                    as.data.frame(tbl$var),
                    file.path(var_dir, "var.csv"),
                    row.names = FALSE)
            }
            writeLines(
                '{"encoding-type": "anndata"}',
                file.path(tbl_dir, ".zattrs"))
            writeLines('{}',
                file.path(obs_dir, ".zattrs"))
            writeLines('{}',
                file.path(var_dir, ".zattrs"))
        }
    }
}

#' Write image/label references
#' @param sl SimpleList of element refs.
#' @param base_path Zarr store root.
#' @param type_dir "images" or "labels".
#' @return Invisible \code{NULL}.
#' @keywords internal
.writeRefElements <- function(sl, base_path,
    type_dir) {
    if (length(sl) == 0L) return(invisible(NULL))
    for (nm in names(sl)) {
        elem <- sl[[nm]]
        if (!is.list(elem) ||
            !"path" %in% names(elem)) next
        src <- elem$path
        dst <- file.path(base_path, type_dir, nm)
        if (dir.exists(src)) {
            dir.create(dst, recursive = TRUE)
            file.copy(
                list.files(src, full.names = TRUE,
                    recursive = TRUE),
                dst, recursive = TRUE)
        }
        ## Write .zattrs
        meta <- elem$metadata
        if (is.null(meta)) meta <- list()
        writeLines(
            jsonlite::toJSON(meta, auto_unbox = TRUE,
                pretty = TRUE),
            file.path(dst, ".zattrs"))
    }
}
