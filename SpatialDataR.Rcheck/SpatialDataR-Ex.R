pkgname <- "SpatialDataR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "SpatialDataR-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('SpatialDataR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("CoordinateTransform-class")
### * CoordinateTransform-class

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CoordinateTransform-class
### Title: CoordinateTransform: Spatial coordinate transformation
### Aliases: CoordinateTransform-class .CoordinateTransform
###   show,CoordinateTransform-method

### ** Examples

ct <- CoordinateTransform("affine", affine = diag(3) * 0.5)
ct



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CoordinateTransform-class", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("CoordinateTransform")
### * CoordinateTransform

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CoordinateTransform
### Title: Create a CoordinateTransform
### Aliases: CoordinateTransform

### ** Examples

# Identity transform
ct <- CoordinateTransform("identity")

# 2D scale + translate
mat <- matrix(c(0.5, 0, 10, 0, 0.5, 20, 0, 0, 1),
    nrow = 3, byrow = TRUE)
ct <- CoordinateTransform("affine", affine = mat,
    input_cs = "pixels", output_cs = "microns")
ct

# 3D affine (4x4)
ct3d <- CoordinateTransform("affine", affine = diag(4) * 0.5)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CoordinateTransform", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("SpatialData-accessors")
### * SpatialData-accessors

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: images
### Title: Access images from a SpatialData object
### Aliases: images spatialLabels spatialPoints shapes tables
###   coordinateSystems images,SpatialData-method
###   spatialLabels,SpatialData-method spatialPoints,SpatialData-method
###   shapes,SpatialData-method tables,SpatialData-method
###   coordinateSystems,SpatialData-method

### ** Examples

sd <- new("SpatialData")
images(sd)
spatialLabels(sd)
spatialPoints(sd)
shapes(sd)
tables(sd)
coordinateSystems(sd)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("SpatialData-accessors", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("SpatialData-class")
### * SpatialData-class

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: SpatialData-class
### Title: SpatialData: Container for multi-modal spatial omics data
### Aliases: SpatialData-class .SpatialData [,SpatialData,character-method
###   length,SpatialData-method names,SpatialData-method
###   show,SpatialData-method

### ** Examples

showClass("SpatialData")
store <- system.file("extdata", "xenium_mini.zarr",
    package = "SpatialDataR")
sd <- readSpatialData(store)
sd_sub <- sd["transcripts"]
sd_sub
sd <- new("SpatialData")
length(sd)
sd <- new("SpatialData")
names(sd)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("SpatialData-class", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("SpatialDataR-package")
### * SpatialDataR-package

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: SpatialDataR-package
### Title: SpatialDataR: Native R Interface to the SpatialData Zarr Format
### Aliases: SpatialDataR SpatialDataR-package
### Keywords: package

### ** Examples

store <- system.file("extdata", "xenium_mini.zarr",
    package = "SpatialDataR")
sd <- readSpatialData(store)
sd

# Element summary
elementSummary(sd)

# Bounding box query
library(S4Vectors)
pts <- spatialPoints(sd)[["transcripts"]]
sub <- bboxQuery(pts, xmin = 0, xmax = 2, ymin = 0, ymax = 2)
nrow(sub)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("SpatialDataR-package", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("aggregatePoints")
### * aggregatePoints

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: aggregatePoints
### Title: Aggregate point features by region
### Aliases: aggregatePoints

### ** Examples

library(S4Vectors)
pts <- DataFrame(
    x = runif(20), y = runif(20),
    gene = sample(c("A", "B", "C"), 20, TRUE),
    cell_id = sample(1:5, 20, TRUE)
)
regions <- DataFrame(
    cell_id = 1:5,
    x = runif(5), y = runif(5)
)
counts <- aggregatePoints(pts, regions)
counts



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("aggregatePoints", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bboxQuery")
### * bboxQuery

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bboxQuery
### Title: Bounding box spatial query
### Aliases: bboxQuery bboxQuery,DataFrame-method
###   bboxQuery,SpatialData-method

### ** Examples

library(S4Vectors)
pts <- DataFrame(
    x = c(1, 2, 3, 4, 5),
    y = c(5, 4, 3, 2, 1),
    gene = c("A", "B", "C", "D", "E")
)
# Query points in [2,4] x [2,4]
sub <- bboxQuery(pts, xmin = 2, xmax = 4, ymin = 2, ymax = 4)
sub



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bboxQuery", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("combineSpatialData")
### * combineSpatialData

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: combineSpatialData
### Title: Combine multiple SpatialData objects
### Aliases: combineSpatialData

### ** Examples

store <- system.file("extdata", "xenium_mini.zarr",
    package = "SpatialDataR")
sd1 <- readSpatialData(store)
sd2 <- readSpatialData(store)
combined <- combineSpatialData(sd1, sd2,
    sample_ids = c("tumor", "normal"))
names(combined)
length(combined)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("combineSpatialData", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("composeTransforms")
### * composeTransforms

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: composeTransforms
### Title: Compose two coordinate transforms
### Aliases: composeTransforms

### ** Examples

# Scale then translate
s <- CoordinateTransform("affine",
    affine = diag(c(0.5, 0.5, 1)),
    input_cs = "pixels", output_cs = "scaled")
t <- CoordinateTransform("affine",
    affine = matrix(c(1,0,10, 0,1,20, 0,0,1),
        nrow = 3, byrow = TRUE),
    input_cs = "scaled", output_cs = "microns")
combined <- composeTransforms(s, t)
combined



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("composeTransforms", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("coordinateSystemElements")
### * coordinateSystemElements

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: coordinateSystemElements
### Title: List coordinate systems with their elements
### Aliases: coordinateSystemElements

### ** Examples

store <- system.file("extdata", "xenium_mini.zarr",
    package = "SpatialDataR")
sd <- readSpatialData(store)
coordinateSystemElements(sd)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("coordinateSystemElements", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dot-parseTransform")
### * dot-parseTransform

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: .parseTransform
### Title: Parse transformation from Zarr metadata
### Aliases: .parseTransform
### Keywords: internal

### ** Examples

meta <- list(coordinateTransformations = list(
    list(type = "scale", scale = c(0.2125, 0.2125))))
SpatialDataR:::.parseTransform(meta)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dot-parseTransform", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("elementSummary")
### * elementSummary

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: elementSummary
### Title: Coerce SpatialData to a summary data.frame
### Aliases: elementSummary

### ** Examples

store <- system.file("extdata", "xenium_mini.zarr",
    package = "SpatialDataR")
sd <- readSpatialData(store)
elementSummary(sd)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("elementSummary", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("elementTransform")
### * elementTransform

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: elementTransform
### Title: Extract element transform from metadata
### Aliases: elementTransform

### ** Examples

store <- system.file("extdata", "xenium_mini.zarr",
    package = "SpatialDataR")
sd <- readSpatialData(store)
ct <- elementTransform(images(sd)[["morphology"]])
ct



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("elementTransform", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("filterSample")
### * filterSample

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: filterSample
### Title: Filter SpatialData by sample
### Aliases: filterSample

### ** Examples

store <- system.file("extdata", "xenium_mini.zarr",
    package = "SpatialDataR")
sd1 <- readSpatialData(store)
sd2 <- readSpatialData(store)
combined <- combineSpatialData(sd1, sd2,
    sample_ids = c("A", "B"))
sdA <- filterSample(combined, "A")
names(sdA)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("filterSample", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("invertTransform")
### * invertTransform

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: invertTransform
### Title: Invert a coordinate transform
### Aliases: invertTransform

### ** Examples

mat <- diag(c(0.2125, 0.2125, 1))
ct <- CoordinateTransform("affine", affine = mat,
    input_cs = "pixels", output_cs = "microns")
inv <- invertTransform(ct)
inv  # microns -> pixels



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("invertTransform", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("readCSVElement")
### * readCSVElement

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: readCSVElement
### Title: Read CSV-backed point or shape data
### Aliases: readCSVElement

### ** Examples

store <- system.file("extdata", "xenium_mini.zarr",
    package = "SpatialDataR")
pts_csv <- file.path(store, "points", "transcripts",
    "transcripts.csv")
df <- readCSVElement(pts_csv)
head(df)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("readCSVElement", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("readParquetPoints")
### * readParquetPoints

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: readParquetPoints
### Title: Read a Parquet-backed point table
### Aliases: readParquetPoints

### ** Examples

## Parquet reading requires the arrow package
## readParquetPoints("path/to/transcripts/")
## Create a small mock example instead
df <- S4Vectors::DataFrame(
    x = c(1.5, 2.3, 4.1),
    y = c(0.8, 3.2, 1.7),
    gene = c("EPCAM", "VIM", "KRT18")
)
df



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("readParquetPoints", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("readSpatialData")
### * readSpatialData

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: readSpatialData
### Title: Read a SpatialData Zarr store into R
### Aliases: readSpatialData

### ** Examples

store <- system.file("extdata", "xenium_mini.zarr",
    package = "SpatialDataR")
sd <- readSpatialData(store)
sd

# Access element references
images(sd)
spatialPoints(sd)

# Selective reading
sd2 <- readSpatialData(store, elements = c("images", "labels"))
sd2



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("readSpatialData", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("readSpatialTable")
### * readSpatialTable

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: readSpatialTable
### Title: Convert a SpatialData table to SpatialExperiment
### Aliases: readSpatialTable

### ** Examples

store <- system.file("extdata", "xenium_mini.zarr",
    package = "SpatialDataR")
tbl <- readSpatialTable(file.path(store, "tables", "table"))
tbl



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("readSpatialTable", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("readZarrArray")
### * readZarrArray

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: readZarrArray
### Title: Read a Zarr array into memory
### Aliases: readZarrArray

### ** Examples

store <- system.file("extdata", "xenium_mini.zarr",
    package = "SpatialDataR")
img_path <- file.path(store, "images", "morphology", "scale0")
if (requireNamespace("Rarr", quietly = TRUE) ||
    requireNamespace("pizzarr", quietly = TRUE)) {
    arr <- readZarrArray(img_path)
    dim(arr)
}



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("readZarrArray", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("transformCoords")
### * transformCoords

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: transformCoords
### Title: Apply coordinate transformation
### Aliases: transformCoords
###   transformCoords,DataFrame,CoordinateTransform-method
###   transformCoords,matrix,CoordinateTransform-method

### ** Examples

library(S4Vectors)
pts <- DataFrame(x = c(100, 200), y = c(50, 150))
ct <- CoordinateTransform("affine",
    affine = diag(c(0.5, 0.5, 1)))
transformCoords(pts, ct)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("transformCoords", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("validateSpatialData")
### * validateSpatialData

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: validateSpatialData
### Title: Validate a SpatialData Zarr store
### Aliases: validateSpatialData

### ** Examples

store <- system.file("extdata", "xenium_mini.zarr",
    package = "SpatialDataR")
result <- validateSpatialData(store)
result$valid
result$elements



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("validateSpatialData", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
