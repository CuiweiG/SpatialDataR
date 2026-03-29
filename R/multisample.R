# R/multisample.R
# Multi-sample operations on SpatialData objects
# Key gap: no R package can combine SpatialData stores

#' @include AllClasses.R
#' @include AllGenerics.R
#' @importFrom S4Vectors SimpleList DataFrame
#' @importFrom methods new slot
NULL

#' Combine multiple SpatialData objects
#'
#' Merges elements from multiple \code{SpatialData} objects
#' into a single container, with automatic name disambiguation.
#' Useful for multi-sample or multi-region spatial analyses.
#'
#' Element names are prefixed with \code{sample_ids} to avoid
#' collisions. Coordinate systems from all samples are merged.
#'
#' @param ... \code{SpatialData} objects to combine.
#' @param sample_ids Character vector of sample identifiers.
#'   If \code{NULL}, defaults to \code{"sample1"},
#'   \code{"sample2"}, etc.
#'
#' @return A new \code{SpatialData} object containing all
#'   elements from all inputs.
#'
#' @export
#' @examples
#' store <- system.file("extdata", "xenium_mini.zarr",
#'     package = "SpatialDataR")
#' sd1 <- readSpatialData(store)
#' sd2 <- readSpatialData(store)
#' combined <- combineSpatialData(sd1, sd2,
#'     sample_ids = c("tumor", "normal"))
#' names(combined)
#' length(combined)
combineSpatialData <- function(..., sample_ids = NULL) {
    objects <- list(...)
    n <- length(objects)
    if (n == 0L) {
        stop("At least one SpatialData object required",
            call. = FALSE)
    }
    if (n == 1L) return(objects[[1L]])

    if (is.null(sample_ids)) {
        sample_ids <- paste0("sample", seq_len(n))
    }
    if (length(sample_ids) != n) {
        stop("Length of 'sample_ids' must match number ",
            "of SpatialData objects",
            call. = FALSE)
    }

    .merge <- function(accessor) {
        merged <- SimpleList()
        for (i in seq_len(n)) {
            sl <- accessor(objects[[i]])
            if (length(sl) == 0L) next
            names(sl) <- paste0(
                sample_ids[i], ".", names(sl))
            merged <- c(merged, sl)
        }
        merged
    }
    all_cs <- .mergeCoordSystems(
        objects, sample_ids, n)
    new("SpatialData",
        images = .merge(images),
        labels = .merge(spatialLabels),
        points = .merge(spatialPoints),
        shapes = .merge(shapes),
        tables = .merge(tables),
        coordinate_systems = all_cs,
        metadata = list(sample_ids = sample_ids),
        path = NA_character_)
}

#' Merge coordinate systems with sample prefixes
#' @param objects List of SpatialData.
#' @param sample_ids Character vector.
#' @param n Number of objects.
#' @return Named list.
#' @keywords internal
.mergeCoordSystems <- function(objects, sample_ids, n) {
    cs <- list()
    for (i in seq_len(n)) {
        for (nm in names(coordinateSystems(objects[[i]]))) {
            cs[[paste0(sample_ids[i], ".", nm)]] <-
                coordinateSystems(objects[[i]])[[nm]]
        }
    }
    cs
}

#' Filter SpatialData by sample
#'
#' Extracts elements belonging to a specific sample from a
#' combined \code{SpatialData} object (created by
#' \code{\link{combineSpatialData}}).
#'
#' @param x A \code{SpatialData} object.
#' @param sample_id Character. Sample identifier to extract.
#'
#' @return A \code{SpatialData} object with only elements
#'   from the specified sample.
#'
#' @export
#' @examples
#' store <- system.file("extdata", "xenium_mini.zarr",
#'     package = "SpatialDataR")
#' sd1 <- readSpatialData(store)
#' sd2 <- readSpatialData(store)
#' combined <- combineSpatialData(sd1, sd2,
#'     sample_ids = c("A", "B"))
#' sdA <- filterSample(combined, "A")
#' names(sdA)
filterSample <- function(x, sample_id) {
    prefix <- paste0(sample_id, ".")
    .filt <- function(sl) {
        keep <- grep(paste0("^", prefix), names(sl))
        if (length(keep) == 0L) return(SimpleList())
        out <- sl[keep]
        names(out) <- sub(paste0("^", prefix), "",
            names(out))
        out
    }
    cs <- coordinateSystems(x)
    cs_keep <- grep(paste0("^", prefix), names(cs))
    cs_sub <- cs[cs_keep]
    if (length(cs_sub) > 0L) {
        names(cs_sub) <- sub(paste0("^", prefix), "",
            names(cs_sub))
    }

    new("SpatialData",
        images = .filt(images(x)),
        labels = .filt(spatialLabels(x)),
        points = .filt(spatialPoints(x)),
        shapes = .filt(shapes(x)),
        tables = .filt(tables(x)),
        coordinate_systems = cs_sub,
        metadata = slot(x, "metadata"),
        path = slot(x, "path"))
}
