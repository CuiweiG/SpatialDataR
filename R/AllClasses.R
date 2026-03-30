# R/AllClasses.R
# S4 class definitions for SpatialDataR

#' @importFrom methods setClass setValidity new is validObject slot
#' @importFrom S4Vectors DataFrame SimpleList
NULL

#' SpatialData: Container for multi-modal spatial omics data
#'
#' A lightweight R-native representation of a SpatialData Zarr
#' store. Holds references to images, labels, points, shapes,
#' and annotation tables. Element data is loaded on explicit
#' request via element-level reader functions.
#'
#' @slot images \code{SimpleList} of image element descriptors.
#'   Each entry is a list with \code{path}, \code{type},
#'   \code{name}, and \code{metadata}.
#' @slot labels \code{SimpleList} of label/segmentation descriptors.
#' @slot points \code{SimpleList} of point element descriptors or
#'   \code{DataFrame} objects (when loaded via CSV/Parquet).
#' @slot shapes \code{SimpleList} of shape descriptors or
#'   \code{DataFrame} objects (when loaded via CSV/Parquet).
#' @slot tables \code{SimpleList} of annotation table descriptors
#'   or \code{SpatialExperiment} objects (when successfully loaded).
#' @slot coordinate_systems \code{list} of coordinate system
#'   definitions parsed from the top-level \code{.zattrs}.
#' @slot metadata \code{list} of additional Zarr store metadata.
#' @slot path Character. Path to the source Zarr store.
#'
#' @references
#' Marconato L et al. (2025). SpatialData: an open and universal
#' data framework for spatial omics. \emph{Nat Methods} 22:58-62.
#' \doi{10.1038/s41592-024-02212-x}
#'
#' @name SpatialData-class
#' @exportClass SpatialData
#' @examples
#' showClass("SpatialData")
.SpatialData <- setClass("SpatialData",
    slots = list(
        images             = "SimpleList",
        labels             = "SimpleList",
        points             = "SimpleList",
        shapes             = "SimpleList",
        tables             = "SimpleList",
        coordinate_systems = "list",
        metadata           = "list",
        path               = "character"
    ),
    prototype = list(
        images             = S4Vectors::SimpleList(),
        labels             = S4Vectors::SimpleList(),
        points             = S4Vectors::SimpleList(),
        shapes             = S4Vectors::SimpleList(),
        tables             = S4Vectors::SimpleList(),
        coordinate_systems = list(),
        metadata           = list(),
        path               = NA_character_
    )
)

setValidity("SpatialData", function(object) {
    msg <- character()

    ## path must be length 1
    if (length(slot(object, "path")) != 1L) {
        msg <- c(msg, "'path' must be a single character string")
    }

    ## coordinate_systems must be a named list (or empty)
    cs <- slot(object, "coordinate_systems")
    if (length(cs) > 0L && is.null(names(cs))) {
        msg <- c(msg,
            "'coordinate_systems' must be a named list")
    }

    if (length(msg) > 0L) msg else TRUE
})

#' CoordinateTransform: Spatial coordinate transformation
#'
#' Represents an affine transformation between coordinate systems.
#'
#' @slot type Character. \code{"identity"} or \code{"affine"}.
#' @slot affine Numeric matrix (3x3 for 2D).
#' @slot input_cs Character. Input coordinate system name.
#' @slot output_cs Character. Output coordinate system name.
#'
#' @references
#' Marconato L et al. (2025). SpatialData: an open and universal
#' data framework for spatial omics. \emph{Nat Methods} 22:58-62.
#' \doi{10.1038/s41592-024-02212-x}
#'
#' @name CoordinateTransform-class
#' @exportClass CoordinateTransform
#' @examples
#' ct <- CoordinateTransform("affine", affine = diag(3) * 0.5)
#' ct
.CoordinateTransform <- setClass("CoordinateTransform",
    slots = list(
        type      = "character",
        affine    = "matrix",
        input_cs  = "character",
        output_cs = "character"
    ),
    prototype = list(
        type      = "identity",
        affine    = diag(3),
        input_cs  = "global",
        output_cs = "global"
    )
)

setValidity("CoordinateTransform", function(object) {
    msg <- character()

    tp <- slot(object, "type")
    if (!tp %in% c("identity", "affine")) {
        msg <- c(msg,
            "'type' must be 'identity' or 'affine'")
    }

    aff <- slot(object, "affine")
    if (!is.numeric(aff) || nrow(aff) != ncol(aff)) {
        msg <- c(msg, "'affine' must be a square numeric matrix")
    }

    if (length(msg) > 0L) msg else TRUE
})
