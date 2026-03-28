# R/read-zarr.R
# Core Zarr reading �?the key innovation: no Python dependency

#' @include AllClasses.R
#' @include AllGenerics.R
#' @importFrom jsonlite fromJSON
#' @importFrom S4Vectors SimpleList DataFrame
#' @importFrom methods new
NULL

#' Read a SpatialData Zarr store into R
#'
#' Reads a SpatialData-formatted \code{.zarr} directory and returns
#' a \code{\linkS4class{SpatialData}} object with lazy references to all
#' elements. No data is loaded into memory until explicitly accessed.
#'
#' @param path Character. Path to a \code{.zarr} directory.
#' @param elements Character vector. Which elements to read:
#'   \code{"images"}, \code{"labels"}, \code{"points"},
#'   \code{"shapes"}, \code{"tables"}. Default: all present.
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{\linkS4class{SpatialData}} object.
#'
#' @details
#' The function discovers elements by reading the top-level
#' \code{.zattrs} metadata. Each element type is read lazily:
#' images and labels as paths (converted to \code{DelayedArray}
#' on access), points and shapes as \code{DataFrame} objects.
#'
#' Zarr array reading uses \pkg{Rarr} (if available) or
#' \pkg{pizzarr} as backends. For Parquet-backed points,
#' \pkg{arrow} is used.
#'
#' @references
#' Marconato L et al. (2024). SpatialData: an open and universal
#' data framework for spatial omics. \emph{Nat Methods} 21:2196-2209.
#' \doi{10.1038/s41592-024-02212-x}
#'
#' @export
#' @examples
#' # Create a minimal mock SpatialData store for demo
#' tmp <- tempfile()
#' dir.create(file.path(tmp, "images", "morphology"), recursive = TRUE)
#' writeLines('{"spatialdata_attrs": {"version": "0.1"}}',
#'     file.path(tmp, ".zattrs"))
#' writeLines('{}', file.path(tmp, "images", "morphology", ".zattrs"))
#' sd <- readSpatialData(tmp)
#' sd
#' unlink(tmp, recursive = TRUE)
readSpatialData <- function(path, elements = NULL, ...) {

    path <- normalizePath(path, mustWork = TRUE)

    ## Read top-level metadata
    zattrs_file <- file.path(path, ".zattrs")
    metadata <- if (file.exists(zattrs_file)) {
        fromJSON(zattrs_file, simplifyVector = FALSE)
    } else {
        list()
    }

    ## Discover elements
    all_dirs <- list.dirs(path, recursive = FALSE,
                          full.names = FALSE)
    element_types <- c("images", "labels", "points",
                       "shapes", "tables")
    present <- intersect(all_dirs, element_types)

    if (!is.null(elements)) {
        present <- intersect(present, elements)
    }

    ## Read each element type
    img_list <- if ("images" %in% present) {
        .readZarrElements(file.path(path, "images"), "image")
    } else SimpleList()

    lbl_list <- if ("labels" %in% present) {
        .readZarrElements(file.path(path, "labels"), "label")
    } else SimpleList()

    pts_list <- if ("points" %in% present) {
        .readZarrElements(file.path(path, "points"), "points")
    } else SimpleList()

    shp_list <- if ("shapes" %in% present) {
        .readZarrElements(file.path(path, "shapes"), "shapes")
    } else SimpleList()

    tbl_list <- if ("tables" %in% present) {
        .readZarrElements(file.path(path, "tables"), "table")
    } else SimpleList()

    ## Extract coordinate systems from metadata
    cs <- .extractCoordinateSystems(metadata)

    new("SpatialData",
        images             = img_list,
        labels             = lbl_list,
        points             = pts_list,
        shapes             = shp_list,
        tables             = tbl_list,
        coordinate_systems = cs,
        metadata           = metadata,
        path               = path)
}

#' @keywords internal
.readZarrElements <- function(dir_path, type) {
    if (!dir.exists(dir_path)) return(SimpleList())

    element_names <- list.dirs(dir_path, recursive = FALSE,
                               full.names = FALSE)
    if (length(element_names) == 0L) return(SimpleList())

    elements <- lapply(element_names, function(name) {
        elem_path <- file.path(dir_path, name)

        ## Read element metadata
        zattrs <- file.path(elem_path, ".zattrs")
        meta <- if (file.exists(zattrs)) {
            fromJSON(zattrs, simplifyVector = FALSE)
        } else list()

        ## Store as lazy reference
        list(path = elem_path,
             type = type,
             name = name,
             metadata = meta,
             loaded = FALSE)
    })
    names(elements) <- element_names
    SimpleList(elements)
}

#' @keywords internal
.extractCoordinateSystems <- function(metadata) {
    sd_attrs <- metadata[["spatialdata_attrs"]]
    if (is.null(sd_attrs)) return(list())

    cs <- sd_attrs[["coordinate_systems"]]
    if (is.null(cs)) return(list())

    ## Return as-is (named list of coordinate system definitions)
    ## Users can construct CoordinateTransform objects manually
    ## or use .parseTransform() on element-level metadata
    cs
}
