#!/usr/bin/env Rscript
# ============================================================================
# Fig 5: Cross-platform — Xenium (10x) + MERFISH (Allen Institute)
# Panel b: use anatomical polygons from shapes/anatomical (proper approach)
# ============================================================================

.libPaths(c("C:/Users/win10/R/win-library/4.4",
            "C:/Users/win10/AppData/Local/R/win-library/4.4",
            .libPaths()))

suppressPackageStartupMessages({
    devtools::load_all("C:/Users/win10/SpatialDataR", quiet = TRUE)
    library(S4Vectors)
    library(ggplot2)
    library(patchwork)
    library(viridis)
    library(arrow)
})

setwd("C:/Users/win10/SpatialDataR")
od <- "man/figures"

th <- function(base = 9) {
    theme_classic(base_size = base, base_family = "sans") +
    theme(
        axis.line        = element_line(linewidth = 0.4, colour = "black"),
        axis.ticks       = element_line(linewidth = 0.3, colour = "black"),
        axis.title       = element_text(size = rel(1), colour = "black"),
        axis.text        = element_text(size = rel(0.9), colour = "grey20"),
        plot.title       = element_text(size = rel(1.1), face = "bold",
                                        hjust = 0, margin = margin(b = 1)),
        plot.subtitle    = element_text(size = rel(0.75), colour = "grey30",
                                        hjust = 0, margin = margin(b = 3)),
        plot.margin      = margin(4, 6, 4, 4),
        legend.title     = element_text(size = rel(0.85), face = "bold"),
        legend.text      = element_text(size = rel(0.75)),
        legend.key.size  = unit(0.3, "cm"),
        legend.background = element_blank()
    )
}

## ---- Parse WKB Polygon geometry -------------------------------------------
.parseWKBPolygon <- function(raw_bytes) {
    ## WKB: 1 byte endian + 4 bytes type + 4 bytes nrings
    ## Each ring: 4 bytes npoints + npoints * (8 bytes x + 8 bytes y)
    endian <- ifelse(raw_bytes[1] == 1, "little", "big")
    n_rings <- readBin(raw_bytes[6:9], "integer", size = 4, endian = endian)
    offset <- 10
    coords_list <- list()
    for (r in seq_len(n_rings)) {
        n_pts <- readBin(raw_bytes[offset:(offset + 3)],
                          "integer", size = 4, endian = endian)
        offset <- offset + 4
        coords <- matrix(0, nrow = n_pts, ncol = 2)
        for (p in seq_len(n_pts)) {
            coords[p, 1] <- readBin(raw_bytes[offset:(offset + 7)],
                                     "double", size = 8, endian = endian)
            coords[p, 2] <- readBin(raw_bytes[(offset + 8):(offset + 15)],
                                     "double", size = 8, endian = endian)
            offset <- offset + 16
        }
        coords_list[[r]] <- coords
    }
    coords_list[[1]]  # Return outer ring
}

## ======================================================================
## READ BOTH DATASETS
## ======================================================================
cat("=== Reading Xenium ===\n")
sd_xenium <- readSpatialData("C:/Users/win10/xenium_breast/data.zarr")
print(sd_xenium)

cat("\n=== Reading MERFISH ===\n")
sd_merfish <- readSpatialData("C:/Users/win10/merfish_scverse/data.zarr")
print(sd_merfish)

## ---- Xenium ---------------------------------------------------------------
xen_pts <- as.data.frame(spatialPoints(sd_xenium)[["transcripts"]])
colnames(xen_pts)[colnames(xen_pts) == "feature_name"] <- "gene"
xen_pts$gene <- as.character(xen_pts$gene)
xen_pts <- xen_pts[xen_pts$qv >= 20 & !grepl("^NegControl|^BLANK", xen_pts$gene), ]

xen_cells <- as.data.frame(shapes(sd_xenium)[["cell_circles"]])
xen_coords <- t(sapply(xen_cells$geometry, function(g) {
    rb <- as.raw(g)
    c(readBin(rb[6:13], "double", size = 8, endian = "little"),
      readBin(rb[14:21], "double", size = 8, endian = "little"))
}))
xen_cells$x <- xen_coords[, 1]; xen_cells$y <- xen_coords[, 2]

cat("Xenium:", format(nrow(xen_pts), big.mark = ","), "transcripts,",
    length(unique(xen_pts$gene)), "genes,",
    format(nrow(xen_cells), big.mark = ","), "cells\n")

## ---- MERFISH: extract anatomical polygons ---------------------------------
mer_pts <- as.data.frame(spatialPoints(sd_merfish)[["single_molecule"]])
mer_cells <- as.data.frame(shapes(sd_merfish)[["cells"]])
anat <- as.data.frame(shapes(sd_merfish)[["anatomical"]])

cat("MERFISH:", format(nrow(mer_pts), big.mark = ","), "transcripts,",
    format(nrow(mer_cells), big.mark = ","), "cells\n")
cat("Anatomical regions:", nrow(anat), "\n")

## Parse all 6 anatomical polygons
## The regions are ordered in the Parquet: rows 1-6
## We need to figure out which is which from spatial position
poly_list <- lapply(seq_len(nrow(anat)), function(i) {
    coords <- .parseWKBPolygon(as.raw(anat$geometry[[i]]))
    data.frame(x = coords[, 1], y = coords[, 2], region_id = i)
})
poly_df <- do.call(rbind, poly_list)

## Assign layer names by matching with transcript cell_type annotations
## Use centroid y-coordinate to order layers (higher y = more pial/superficial)
centroids <- sapply(poly_list, function(p) mean(p$y))
## In mouse cortex, Layer I is most superficial (highest y in this dataset)
## Order: most superficial to deepest
layer_order <- order(centroids, decreasing = TRUE)
layer_names <- c("VISp_I", "VISp_II/III", "VISp_IV",
                  "VISp_V", "VISp_VI", "VISp_wm")
poly_df$layer <- layer_names[match(poly_df$region_id, layer_order)]

cat("Layer assignment by centroid y (highest=superficial):\n")
for (i in seq_along(layer_order)) {
    cat("  ", layer_names[i], ": region", layer_order[i],
        "centroid_y =", round(centroids[layer_order[i]]), "\n")
}

## ---- Verify layer assignment against transcript cell_type ----
## Sample points from each polygon and check their cell_type
cat("\nVerifying layer assignment:\n")
for (lid in seq_len(6)) {
    p <- poly_list[[layer_order[lid]]]
    cx <- mean(p$x); cy <- mean(p$y)
    ## Find transcripts near centroid
    nearby <- mer_pts[abs(mer_pts$x - cx) < 30 & abs(mer_pts$y - cy) < 30, ]
    if (nrow(nearby) > 0) {
        dom <- names(sort(table(nearby$cell_type), decreasing = TRUE))[1]
        cat("  ", layer_names[lid], "-> dominant cell_type:", dom,
            "(n=", nrow(nearby), ")\n")
    }
}

## ======================================================================
## FIG 5: Build figure
## ======================================================================
cat("\n--- Fig 5 ---\n")

## Panel a: Xenium density
xen_bin <- 30
xen_pts$xb <- round(xen_pts$x / xen_bin) * xen_bin
xen_pts$yb <- round(xen_pts$y / xen_bin) * xen_bin
xen_density <- aggregate(gene ~ xb + yb, data = xen_pts, FUN = length)
colnames(xen_density)[3] <- "count"

p5a <- ggplot(xen_density, aes(x = xb, y = yb, fill = count)) +
    geom_raster() +
    scale_fill_viridis_c(option = "inferno", trans = "sqrt",
                          name = "Transcripts") +
    coord_equal(expand = FALSE) +
    annotate("segment",
             x = max(xen_pts$x) - 1100, xend = max(xen_pts$x) - 100,
             y = min(xen_pts$y) + 100, yend = min(xen_pts$y) + 100,
             linewidth = 1.5, colour = "white") +
    annotate("text",
             x = max(xen_pts$x) - 600, y = min(xen_pts$y) + 100,
             label = "1 mm", vjust = -0.7, size = 2.5,
             colour = "white", fontface = "bold") +
    labs(title = "10x Xenium",
         subtitle = paste0("Human breast cancer\n",
                           format(nrow(xen_pts), big.mark = ","),
                           " transcripts, ", length(unique(xen_pts$gene)),
                           " genes, ", format(nrow(xen_cells), big.mark = ","),
                           " cells"),
         x = expression(italic(x)~"("*mu*"m)"),
         y = expression(italic(y)~"("*mu*"m)")) +
    th(8.5)

## Panel b: MERFISH with anatomical polygons + transcript density overlay
layer_cols <- c(
    "VISp_I"       = "#D62728",
    "VISp_II/III"  = "#1F77B4",
    "VISp_IV"      = "#2CA02C",
    "VISp_V"       = "#9467BD",
    "VISp_VI"      = "#FF7F0E",
    "VISp_wm"      = "#8C564B"
)
layer_labels <- c(
    "VISp_I" = "I", "VISp_II/III" = "II/III", "VISp_IV" = "IV",
    "VISp_V" = "V", "VISp_VI" = "VI", "VISp_wm" = "WM"
)

poly_df$layer_f <- factor(poly_df$layer,
                           levels = names(layer_cols))

## Cell centroids for overlay
mer_cell_coords <- t(sapply(mer_cells$geometry, function(g) {
    rb <- as.raw(g)
    c(readBin(rb[6:13], "double", size = 8, endian = "little"),
      readBin(rb[14:21], "double", size = 8, endian = "little"))
}))
mer_cells$x <- mer_cell_coords[, 1]; mer_cells$y <- mer_cell_coords[, 2]

p5b <- ggplot() +
    ## Anatomical region polygons (filled)
    geom_polygon(data = poly_df,
                 aes(x = x, y = y, fill = layer_f, group = region_id),
                 colour = "white", linewidth = 0.4) +
    ## Cell centroids
    geom_point(data = mer_cells, aes(x = x, y = y),
               size = 0.3, alpha = 0.4, colour = "grey20") +
    scale_fill_manual(values = layer_cols, name = "Layer",
                      labels = layer_labels) +
    coord_equal(expand = FALSE) +
    annotate("segment",
             x = max(mer_pts$x) - 550, xend = max(mer_pts$x) - 50,
             y = min(mer_pts$y) + 50, yend = min(mer_pts$y) + 50,
             linewidth = 1.5, colour = "black") +
    annotate("text",
             x = max(mer_pts$x) - 300, y = min(mer_pts$y) + 50,
             label = "500 \u00b5m", vjust = -0.7, size = 2.5,
             fontface = "bold") +
    labs(title = "MERFISH",
         subtitle = paste0("Mouse primary visual cortex\n",
                           format(nrow(mer_pts), big.mark = ","),
                           " transcripts, ",
                           format(nrow(mer_cells), big.mark = ","),
                           " cells, 268 genes"),
         x = expression(italic(x)~"("*mu*"m)"),
         y = expression(italic(y)~"("*mu*"m)")) +
    th(8.5)

fig5 <- (p5a + labs(tag = "a")) + (p5b + labs(tag = "b")) +
    plot_layout(ncol = 2, widths = c(1, 0.85)) +
    plot_annotation(
        title = "Cross-platform compatibility: one API, multiple technologies",
        subtitle = paste0("Both datasets read with readSpatialData() \u2014 ",
                          "no platform-specific code, no Python dependency"),
        theme = theme(
            plot.title = element_text(size = 11, face = "bold"),
            plot.subtitle = element_text(size = 8.5, colour = "grey30"))) &
    theme(plot.tag = element_text(size = 10, face = "bold"))

ggsave(file.path(od, "fig5_crossplatform.png"), fig5,
       width = 220, height = 110, units = "mm", dpi = 300, bg = "white")
cat("  saved (", round(file.size(file.path(od, "fig5_crossplatform.png"))/1024),
    "KB)\n")
