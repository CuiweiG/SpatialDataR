# R/read-zarr.R
# Core Zarr store reading â€?no Python dependency

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
#' stored as path references â€?use \code{\link{readZarrArray}} to
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
#' Marconato L et al. (2025). SpatialData: an open and universal
#' data framework for spatial omics. \emph{Nat Methods} 22:58-62.
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
    ## Zarr v2: .zattrs
    zattrs_v2 <- file.path(path, ".zattrs")
    if (file.exists(zattrs_v2)) {
        return(fromJSON(zattrs_v2,
            simplifyVector = FALSE))
    }
    ## Zarr v3: zarr.json
    zarr_json <- file.path(path, "zarr.json")
    if (file.exists(zarr_json)) {
        meta <- fromJSON(zarr_json,
            simplifyVector = FALSE)
        ## Extract attributes from Zarr v3 format
        attrs <- meta[["attributes"]]
        if (!is.null(attrs)) return(attrs)
        return(meta)
    }
    list()
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
        meta <- .readZattrs(elem_path)
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

#' Read points or shapes â€?CSV/Parquet eager load
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
        ## GeoParquet single file (shapes)
        if (is.null(data)) {
            gpq <- list.files(elem_path,
                "^shapes[.]parquet$",
                full.names = TRUE)
            if (length(gpq) > 0L &&
                requireNamespace("arrow",
                    quietly = TRUE)) {
                data <- tryCatch({
                    df <- as.data.frame(
                        arrow::read_parquet(gpq[1L]))
                    DataFrame(df)
                }, error = function(e) NULL)
            }
        }
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
    ## Check for Parquet files (incl. subdirectories)
    pq <- list.files(elem_path, "[.]parquet$",
        full.names = TRUE, recursive = TRUE)
    ## Also check for .parquet directories
    pq_dirs <- list.dirs(elem_path,
        recursive = FALSE, full.names = TRUE)
    pq_dirs <- pq_dirs[grepl("[.]parquet$", pq_dirs)]

    csv <- list.files(elem_path, "[.]csv$",
        full.names = TRUE, recursive = TRUE)

    data <- NULL
    ## Try Parquet files first
    if ((length(pq) > 0L || length(pq_dirs) > 0L) &&
        requireNamespace("arrow", quietly = TRUE)) {
        target <- if (length(pq) > 0L) {
            elem_path
        } else {
            pq_dirs[1L]
        }
        data <- tryCatch(
            readParquetPoints(target),
            error = function(e) NULL)
    }
    ## Fallback to CSV
    if (is.null(data) && length(csv) > 0L) {
        data <- tryCatch(readCSVElement(csv[1L]),
            error = function(e) NULL)
    }
    data
}

#' Read tables group â€?attempt SpatialExperiment conversion
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
    if (is.null(sd_attrs)) {
        ## Zarr v3: might be nested under attributes
        sd_attrs <- metadata[["attributes"]]
        if (!is.null(sd_attrs)) {
            sd_attrs <- sd_attrs[["spatialdata_attrs"]]
        }
    }
    if (is.null(sd_attrs)) return(list())
    cs <- sd_attrs[["coordinate_systems"]]
    if (is.null(cs)) return(list())
    cs
}
