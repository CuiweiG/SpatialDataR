# R/read-zarr.R
# Core Zarr store reading — no Python dependency

#' @include AllClasses.R
#' @include AllGenerics.R
#' @importFrom jsonlite fromJSON
#' @importFrom S4Vectors SimpleList DataFrame
#' @importFrom methods new
NULL

#' Read a SpatialData Zarr store into R
#'
#' Reads a SpatialData-formatted \code{.zarr} directory and returns
#' a \code{\linkS4class{SpatialData}} object. Each element type is
#' discovered from the directory structure and \code{.zattrs}
#' metadata. Points and shapes stored as CSV or Parquet are
#' loaded as \code{DataFrame} objects. Images and labels are
#' stored as path references — use \code{\link{readZarrArray}} to
#' load them into memory.
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
#' The function discovers elements by scanning subdirectories
#' of the Zarr store. Images and labels are stored as lightweight
#' path references (a list with \code{path}, \code{name}, and
#' \code{metadata} fields). Use \code{\link{readZarrArray}} to
#' load the actual array data from the stored path.
#'
#' Points and shapes are eagerly loaded when CSV or Parquet files
#' are present. Tables are loaded via \code{\link{readSpatialTable}}
#' when the obs/var structure is detected.
#'
#' @references
#' Marconato L et al. (2024). SpatialData: an open and universal
#' data framework for spatial omics. \emph{Nat Methods} 21:2196-2209.
#' \doi{10.1038/s41592-024-02212-x}
#'
#' @export
#' @examples
#' store <- system.file("extdata", "xenium_mini.zarr",
#'     package = "SpatialDataR")
#' sd <- readSpatialData(store)
#' sd
#'
#' # Access element references
#' images(sd)
#' spatialPoints(sd)
#'
#' # Selective reading
#' sd2 <- readSpatialData(store, elements = c("images", "labels"))
#' sd2
readSpatialData <- function(path, elements = NULL, ...) {
    path <- normalizePath(path, mustWork = TRUE)
    metadata <- .readZattrs(path)

    present <- .discoverElements(path, elements)

    .get <- function(type, reader, ...) {
        if (type %in% present) reader(file.path(path, type), ...)
        else SimpleList()
    }

    new("SpatialData",
        images             = .get("images", .readElementRefs, "image"),
        labels             = .get("labels", .readElementRefs, "label"),
        points             = .get("points", .readPointsOrShapes, "points"),
        shapes             = .get("shapes", .readPointsOrShapes, "shapes"),
        tables             = .get("tables", .readTablesGroup),
        coordinate_systems = .extractCoordinateSystems(metadata),
        metadata           = metadata,
        path               = path)
}

#' Read .zattrs from a Zarr directory
#' @param path Path to directory containing .zattrs.
#' @return A list of parsed metadata.
#' @keywords internal
.readZattrs <- function(path) {
    zattrs_file <- file.path(path, ".zattrs")
    if (file.exists(zattrs_file)) {
        fromJSON(zattrs_file, simplifyVector = FALSE)
    } else {
        list()
    }
}

#' Discover element types in a Zarr store
#' @param path Path to the Zarr store root.
#' @param elements Optional character vector filter.
#' @return Character vector of present element type names.
#' @keywords internal
.discoverElements <- function(path, elements = NULL) {
    all_dirs <- list.dirs(path, recursive = FALSE,
        full.names = FALSE)
    present <- intersect(all_dirs,
        c("images", "labels", "points", "shapes", "tables"))
    if (!is.null(elements)) intersect(present, elements)
    else present
}

#' Read element path references (images/labels)
#' @param dir_path Path to the element type directory.
#' @param type Character label for element type.
#' @return A \code{SimpleList} of element descriptor lists.
#' @keywords internal
.readElementRefs <- function(dir_path, type) {
    if (!dir.exists(dir_path)) return(SimpleList())

    element_names <- list.dirs(
        dir_path, recursive = FALSE, full.names = FALSE
    )
    if (length(element_names) == 0L) return(SimpleList())

    elements <- lapply(element_names, function(name) {
        elem_path <- file.path(dir_path, name)
        zattrs <- file.path(elem_path, ".zattrs")
        meta <- if (file.exists(zattrs)) {
            fromJSON(zattrs, simplifyVector = FALSE)
        } else {
            list()
        }
        list(
            path = elem_path,
            type = type,
            name = name,
            metadata = meta
        )
    })
    names(elements) <- element_names
    SimpleList(elements)
}

#' Read points or shapes — CSV/Parquet eager load
#' @param dir_path Path to the points or shapes directory.
#' @param type Character label for element type.
#' @return A \code{SimpleList} of \code{DataFrame} or descriptor lists.
#' @keywords internal
.readPointsOrShapes <- function(dir_path, type) {
    if (!dir.exists(dir_path)) return(SimpleList())
    element_names <- list.dirs(
        dir_path, recursive = FALSE, full.names = FALSE)
    if (length(element_names) == 0L) return(SimpleList())

    elements <- lapply(element_names, function(name) {
        elem_path <- file.path(dir_path, name)
        meta <- .readZattrs(elem_path)
        data <- .loadTabularElement(elem_path)
        if (!is.null(data)) {
            attr(data, "spatialdata_metadata") <- meta
            data
        } else {
            list(path = elem_path, type = type,
                name = name, metadata = meta)
        }
    })
    names(elements) <- element_names
    SimpleList(elements)
}

#' Load tabular data from an element directory
#' @param elem_path Path to element directory.
#' @return A \code{DataFrame} or \code{NULL}.
#' @keywords internal
.loadTabularElement <- function(elem_path) {
    pq <- list.files(elem_path, "[.]parquet$",
        full.names = TRUE, recursive = TRUE)
    csv <- list.files(elem_path, "[.]csv$",
        full.names = TRUE, recursive = TRUE)

    data <- NULL
    if (length(pq) > 0L &&
        requireNamespace("arrow", quietly = TRUE)) {
        data <- tryCatch(readParquetPoints(elem_path),
            error = function(e) NULL)
    }
    if (is.null(data) && length(csv) > 0L) {
        data <- tryCatch(readCSVElement(csv[1L]),
            error = function(e) NULL)
    }
    data
}

#' Read tables group — attempt SpatialExperiment conversion
#' @param dir_path Path to the tables directory.
#' @return A \code{SimpleList} of table objects (list or SpatialExperiment).
#' @keywords internal
.readTablesGroup <- function(dir_path) {
    if (!dir.exists(dir_path)) return(SimpleList())

    table_names <- list.dirs(
        dir_path, recursive = FALSE, full.names = FALSE
    )
    if (length(table_names) == 0L) return(SimpleList())

    tables <- lapply(table_names, function(name) {
        tryCatch(
            readSpatialTable(file.path(dir_path, name)),
            error = function(e) {
                ## Fallback to reference
                list(
                    path = file.path(dir_path, name),
                    type = "table",
                    name = name,
                    error = conditionMessage(e)
                )
            }
        )
    })
    names(tables) <- table_names
    SimpleList(tables)
}

#' Extract coordinate systems from store metadata
#' @param metadata List. Parsed top-level \code{.zattrs}.
#' @return A named list of coordinate system definitions.
#' @keywords internal
.extractCoordinateSystems <- function(metadata) {
    sd_attrs <- metadata[["spatialdata_attrs"]]
    if (is.null(sd_attrs)) return(list())

    cs <- sd_attrs[["coordinate_systems"]]
    if (is.null(cs)) return(list())

    cs
}
