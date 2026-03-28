# R/AllClasses.R
# S4 class definitions for SpatialDataR

#' @importFrom methods setClass setValidity new is validObject slot
#' @importFrom S4Vectors DataFrame SimpleList
NULL

#' SpatialData: Container for multi-modal spatial omics data
#'
#' A lightweight R-native representation of a SpatialData Zarr
#' store. Holds references to images, labels, points, shapes,
#' and annotation tables without loading full data into memory.
#'
#' @slot images \code{SimpleList} of image references.
#' @slot labels \code{SimpleList} of label/segmentation references.
#' @slot points \code{SimpleList} of point DataFrames.
#' @slot shapes \code{SimpleList} of shape DataFrames.
#' @slot tables \code{SimpleList} of annotation tables.
#' @slot coordinate_systems \code{list} of coordinate systems.
#' @slot metadata \code{list} of additional metadata.
#' @slot path Character. Path to the source Zarr store.
#'
#' @references
#' Marconato L et al. (2024). SpatialData: an open and universal
#' data framework for spatial omics. \emph{Nat Methods} 21:2196-2209.
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
#' Marconato L et al. (2024). SpatialData: an open and universal
#' data framework for spatial omics. \emph{Nat Methods} 21:2196-2209.
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
