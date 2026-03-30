.libPaths(c("C:/Users/win10/R/win-library/4.4", .libPaths()))
devtools::load_all("C:/Users/win10/SpatialDataR", quiet = TRUE)

sd <- readSpatialData("C:/Users/win10/merfish_scverse/data.zarr")

## Check anatomical shapes — these should be layer boundary polygons
anat <- shapes(sd)[["anatomical"]]
cat("Anatomical class:", class(anat)[1], "\n")
cat("Rows:", nrow(anat), "\n")
cat("Cols:", paste(names(anat), collapse = ", "), "\n\n")

## Check geometry column
adf <- as.data.frame(anat)
for (i in seq_len(nrow(adf))) {
    g <- adf$geometry[[i]]
    raw <- as.raw(g)
    cat("Row", i, ": geometry size =", length(raw), "bytes")
    ## WKB type: bytes 2-5 (uint32)
    wkb_type <- readBin(raw[2:5], "integer", size = 4, endian = "little")
    cat(", WKB type =", wkb_type)
    ## WKB types: 1=Point, 2=LineString, 3=Polygon, 4=MultiPoint,
    ##            5=MultiLineString, 6=MultiPolygon
    type_name <- switch(as.character(wkb_type),
        "1" = "Point", "2" = "LineString", "3" = "Polygon",
        "4" = "MultiPoint", "5" = "MultiLineString",
        "6" = "MultiPolygon", paste0("Unknown(", wkb_type, ")"))
    cat(" (", type_name, ")\n")
}

## Try to parse polygon coordinates for the first shape
cat("\n=== Parsing first polygon ===\n")
g1 <- as.raw(adf$geometry[[1]])
endian <- g1[1]  # 1 = little-endian
wkb_type <- readBin(g1[2:5], "integer", size = 4, endian = "little")
cat("Type:", wkb_type, "\n")

if (wkb_type == 3) {  # Polygon
    n_rings <- readBin(g1[6:9], "integer", size = 4, endian = "little")
    cat("Rings:", n_rings, "\n")
    offset <- 10
    for (r in seq_len(n_rings)) {
        n_pts <- readBin(g1[offset:(offset+3)], "integer", size = 4, endian = "little")
        cat("Ring", r, ":", n_pts, "points\n")
        offset <- offset + 4
        coords <- matrix(0, nrow = n_pts, ncol = 2)
        for (p in seq_len(n_pts)) {
            x <- readBin(g1[offset:(offset+7)], "double", size = 8, endian = "little")
            y <- readBin(g1[(offset+8):(offset+15)], "double", size = 8, endian = "little")
            coords[p, ] <- c(x, y)
            offset <- offset + 16
        }
        cat("  x range:", range(coords[, 1]), "\n")
        cat("  y range:", range(coords[, 2]), "\n")
    }
} else if (wkb_type == 6) {  # MultiPolygon
    n_polys <- readBin(g1[6:9], "integer", size = 4, endian = "little")
    cat("Polygons:", n_polys, "\n")
}

## Also check cells shapes for comparison
cells <- shapes(sd)[["cells"]]
cat("\n=== Cells ===\n")
cat("Rows:", nrow(cells), "\n")
cat("Cols:", paste(names(cells), collapse = ", "), "\n")
