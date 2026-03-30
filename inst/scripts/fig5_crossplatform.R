#!/usr/bin/env Rscript
# ============================================================================
# Fig 5: Cross-platform compatibility — same API, different technologies
# Xenium (10x) + MERFISH (Vizgen) read by the same readSpatialData() call
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
                                        hjust = 0, colour = "black",
                                        margin = margin(b = 1)),
        plot.subtitle    = element_text(size = rel(0.75), colour = "grey30",
                                        hjust = 0, margin = margin(b = 3)),
        plot.margin      = margin(4, 6, 4, 4),
        legend.title     = element_text(size = rel(0.85), face = "bold"),
        legend.text      = element_text(size = rel(0.75)),
        legend.key.size  = unit(0.3, "cm"),
        legend.background = element_blank()
    )
}

## ======================================================================
## READ BOTH DATASETS WITH THE SAME API
## ======================================================================

cat("=== Reading Xenium breast cancer ===\n")
sd_xenium <- readSpatialData("C:/Users/win10/xenium_breast/data.zarr")
print(sd_xenium)

cat("\n=== Reading MERFISH mouse brain ===\n")
sd_merfish <- readSpatialData("C:/Users/win10/merfish_scverse/data.zarr")
print(sd_merfish)

## ---- Xenium data ----------------------------------------------------------
xen_pts <- as.data.frame(spatialPoints(sd_xenium)[["transcripts"]])
colnames(xen_pts)[colnames(xen_pts) == "feature_name"] <- "gene"
xen_pts$gene <- as.character(xen_pts$gene)
xen_pts <- xen_pts[xen_pts$qv >= 20 & !grepl("^NegControl|^BLANK", xen_pts$gene), ]

xen_cells <- as.data.frame(shapes(sd_xenium)[["cell_circles"]])
xen_coords <- t(sapply(xen_cells$geometry, function(g) {
    raw_bytes <- as.raw(g)
    x <- readBin(raw_bytes[6:13], "double", size = 8, endian = "little")
    y <- readBin(raw_bytes[14:21], "double", size = 8, endian = "little")
    c(x = x, y = y)
}))
xen_cells$x <- xen_coords[, 1]
xen_cells$y <- xen_coords[, 2]

cat("Xenium: ", format(nrow(xen_pts), big.mark = ","), " transcripts, ",
    length(unique(xen_pts$gene)), " genes, ",
    format(nrow(xen_cells), big.mark = ","), " cells\n")

## ---- MERFISH data ---------------------------------------------------------
mer_pts <- as.data.frame(spatialPoints(sd_merfish)[["single_molecule"]])
mer_pts$gene <- mer_pts$cell_type  # In this dataset, cell_type column is the gene/region annotation

## Actually check what columns we have
cat("MERFISH columns:", paste(names(mer_pts), collapse = ", "), "\n")

## The MERFISH scverse dataset has cell_type as annotation
## Need to check if there's a gene column
## Let's see what unique cell_types there are
cat("Unique cell_type values:", length(unique(mer_pts$cell_type)), "\n")
cat("First 10:", paste(head(unique(mer_pts$cell_type), 10), collapse = ", "), "\n")

## Get cells
mer_cells <- as.data.frame(shapes(sd_merfish)[["cells"]])
mer_cell_coords <- t(sapply(mer_cells$geometry, function(g) {
    raw_bytes <- as.raw(g)
    x <- readBin(raw_bytes[6:13], "double", size = 8, endian = "little")
    y <- readBin(raw_bytes[14:21], "double", size = 8, endian = "little")
    c(x = x, y = y)
}))
mer_cells$x <- mer_cell_coords[, 1]
mer_cells$y <- mer_cell_coords[, 2]

cat("MERFISH: ", format(nrow(mer_pts), big.mark = ","), " transcripts, ",
    format(nrow(mer_cells), big.mark = ","), " cells\n")

## ======================================================================
## FIG 5: Side-by-side comparison
## ======================================================================
cat("\n--- Fig 5: Cross-platform ---\n")

## Panel a: Xenium tissue (transcript density)
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

## Panel b: MERFISH tissue (transcript density colored by cell type / region)
mer_bin <- 10
mer_pts$xb <- round(mer_pts$x / mer_bin) * mer_bin
mer_pts$yb <- round(mer_pts$y / mer_bin) * mer_bin

## Assign dominant region per bin
region_layers <- c("VISp_I", "VISp_II/III", "VISp_IV",
                   "VISp_V", "VISp_VI", "VISp_wm")
region_cols <- c(
    "VISp_I"       = "#D62728",
    "VISp_II/III"  = "#1F77B4",
    "VISp_IV"      = "#2CA02C",
    "VISp_V"       = "#9467BD",
    "VISp_VI"      = "#FF7F0E",
    "VISp_wm"      = "#8C564B",
    "Tissue"       = "#D9D9D9"
)
region_labels <- c(
    "VISp_I" = "I", "VISp_II/III" = "II/III", "VISp_IV" = "IV",
    "VISp_V" = "V", "VISp_VI" = "VI", "VISp_wm" = "WM",
    "Tissue" = "Other"
)

mer_bin_layer <- aggregate(cell_type ~ xb + yb, data = mer_pts,
                            FUN = function(x) {
                                tt <- table(x)
                                cx <- intersect(names(tt), region_layers)
                                if (length(cx) > 0) {
                                    cx_tt <- tt[cx]
                                    names(cx_tt)[which.max(cx_tt)]
                                } else {
                                    "Tissue"
                                }
                            })
mer_bin_layer$layer_f <- factor(mer_bin_layer$cell_type,
                                 levels = c(region_layers, "Tissue"))

p5b <- ggplot(mer_bin_layer, aes(x = xb, y = yb, fill = layer_f)) +
    geom_raster() +
    scale_fill_manual(values = region_cols, name = "Layer",
                      labels = region_labels,
                      guide = guide_legend(ncol = 1)) +
    coord_equal(expand = FALSE) +
    annotate("segment",
             x = max(mer_pts$x) - 550, xend = max(mer_pts$x) - 50,
             y = min(mer_pts$y) + 50, yend = min(mer_pts$y) + 50,
             linewidth = 1.5, colour = "black") +
    annotate("text",
             x = max(mer_pts$x) - 300, y = min(mer_pts$y) + 50,
             label = "500 \u00b5m", vjust = -0.7, size = 2.5,
             colour = "black", fontface = "bold") +
    labs(title = "Vizgen MERFISH",
         subtitle = paste0("Mouse primary visual cortex\n",
                           format(nrow(mer_pts), big.mark = ","),
                           " transcripts, ",
                           format(nrow(mer_cells), big.mark = ","),
                           " cells"),
         x = expression(italic(x)~"("*mu*"m)"),
         y = expression(italic(y)~"("*mu*"m)")) +
    th(8.5)

## Combine
fig5 <- (p5a + labs(tag = "a")) + (p5b + labs(tag = "b")) +
    plot_layout(ncol = 2, widths = c(1, 0.85)) +
    plot_annotation(
        title = "Cross-platform compatibility: one API, multiple technologies",
        subtitle = paste0("Both datasets read with readSpatialData() — ",
                          "no platform-specific code, no Python dependency"),
        theme = theme(
            plot.title = element_text(size = 11, face = "bold"),
            plot.subtitle = element_text(size = 8.5, colour = "grey30"))) &
    theme(plot.tag = element_text(size = 10, face = "bold"))

ggsave(file.path(od, "fig5_crossplatform.png"), fig5,
       width = 220, height = 110, units = "mm", dpi = 300, bg = "white")
cat("  saved\n")

sz <- round(file.size(file.path(od, "fig5_crossplatform.png")) / 1024)
cat("  fig5_crossplatform.png (", sz, "KB)\n")
