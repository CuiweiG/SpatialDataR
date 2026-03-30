#!/usr/bin/env Rscript
# ============================================================================
# SpatialDataR README figures — v4 (Nature Methods standard)
#
# 4 figures, all from real MERFISH VISp data (Moffitt 2018 Science)
# Descriptive titles, colorblind-safe palettes, clustered heatmaps
# ============================================================================

.libPaths(c("C:/Users/win10/R/win-library/4.4",
            "C:/Users/win10/AppData/Local/R/win-library/4.4",
            .libPaths()))

suppressPackageStartupMessages({
    library(SpatialDataR)
    library(S4Vectors)
    library(ggplot2)
    library(patchwork)
    library(viridis)
    library(ggrepel)
    library(reshape2)
    library(scales)
    library(RColorBrewer)
})

setwd("C:/Users/win10/SpatialDataR")
od <- "man/figures"
dir.create(od, showWarnings = FALSE, recursive = TRUE)

## ---- Global theme (Nature Methods style) ----------------------------------
## Nature Methods: 7pt min font, Helvetica/Arial, minimal chrome
th <- function(base = 9) {
    theme_minimal(base_size = base, base_family = "sans") +
    theme(
        panel.grid       = element_blank(),
        panel.border     = element_rect(colour = "grey30", fill = NA,
                                        linewidth = 0.4),
        axis.ticks       = element_line(linewidth = 0.25, colour = "grey30"),
        axis.title       = element_text(size = rel(1), colour = "grey10"),
        axis.text        = element_text(size = rel(0.9), colour = "grey30"),
        plot.title       = element_text(size = rel(1.15), face = "bold",
                                        hjust = 0, margin = margin(b = 2)),
        plot.subtitle    = element_text(size = rel(0.85), colour = "grey40",
                                        hjust = 0, margin = margin(b = 4)),
        plot.margin      = margin(4, 6, 4, 4),
        legend.title     = element_text(size = rel(0.9), face = "bold"),
        legend.text      = element_text(size = rel(0.8)),
        legend.key.size  = unit(0.35, "cm"),
        legend.background = element_blank(),
        strip.text       = element_text(face = "bold", size = rel(0.9))
    )
}

## Colorblind-safe cortical layer palette (Okabe-Ito derived)
layer_cols <- c(
    "VISp_I"       = "#E69F00",  # orange
    "VISp_II/III"  = "#56B4E9",  # sky blue
    "VISp_IV"      = "#009E73",  # teal
    "VISp_V"       = "#CC79A7",  # pink
    "VISp_VI"      = "#0072B2",  # dark blue
    "VISp_wm"      = "#D55E00"   # vermillion
)
layer_labels <- c(
    "VISp_I"       = "I",
    "VISp_II/III"  = "II/III",
    "VISp_IV"      = "IV",
    "VISp_V"       = "V",
    "VISp_VI"      = "VI",
    "VISp_wm"      = "WM"
)
cortex_layers <- c("VISp_I", "VISp_II/III", "VISp_IV",
                   "VISp_V", "VISp_VI", "VISp_wm")

## ---- Load data ------------------------------------------------------------
cat("=== Loading MERFISH data ===\n")
store <- "C:/Users/win10/merfish_spatialdata.zarr"
sd <- readSpatialData(store)
pts <- as.data.frame(spatialPoints(sd)[["transcripts"]])
shp <- as.data.frame(shapes(sd)[["cell_boundaries"]])
obs <- as.data.frame(tables(sd)[["table"]]$obs)

## Layer assignment via coordinate-key lookup
raw <- read.csv("C:/Users/win10/merfish_real.csv", stringsAsFactors = FALSE)
raw$key <- paste0(round(raw$x_um, 2), "_", round(raw$y_um, 2))
pts$key <- paste0(round(pts$x, 2), "_", round(pts$y, 2))
layer_lookup <- setNames(raw$layer, raw$key)
layer_lookup <- layer_lookup[!duplicated(names(layer_lookup))]
pts$layer <- layer_lookup[pts$key]
cat("  Layer match rate:", round(100 * mean(!is.na(pts$layer)), 1), "%\n")
cat("  Transcripts:", format(nrow(pts), big.mark = ","), "\n")
cat("  Genes:", length(unique(pts$gene)), "\n")
cat("  Cells:", nrow(shp), "\n\n")

x_rng <- range(pts$x); y_rng <- range(pts$y)

## ======================================================================
## FIG 1: Laminar architecture of mouse primary visual cortex
## ======================================================================
cat("--- Fig 1: Laminar architecture ---\n")

bin_size <- 10
pts$xbin <- round(pts$x / bin_size) * bin_size
pts$ybin <- round(pts$y / bin_size) * bin_size

## Binning: prefer cortex layer; non-cortex bins → "Tissue" background
bin_layer <- aggregate(layer ~ xbin + ybin, data = pts,
                       FUN = function(x) {
                           tt <- table(x)
                           cx_names <- intersect(names(tt), cortex_layers)
                           if (length(cx_names) > 0) {
                               cx_tt <- tt[cx_names]
                               names(cx_tt)[which.max(cx_tt)]
                           } else {
                               "Tissue"
                           }
                       })

fig1_cols <- c(layer_cols, "Tissue" = "#F0F0F0")
all_levels <- c(cortex_layers, "Tissue")
fig1_labels <- c(layer_labels, "Tissue" = "Non-cortex")
bin_layer$layer_f <- factor(bin_layer$layer, levels = all_levels)

fig1 <- ggplot(bin_layer, aes(x = xbin, y = ybin, fill = layer_f)) +
    geom_raster() +
    scale_fill_manual(values = fig1_cols, name = "Layer",
                      labels = fig1_labels,
                      guide = guide_legend(ncol = 1,
                          override.aes = list(colour = NA))) +
    coord_equal(expand = FALSE) +
    ## Scale bar (black bar, bottom-right)
    annotate("segment", x = x_rng[2] - 550, xend = x_rng[2] - 50,
             y = y_rng[1] + 50, yend = y_rng[1] + 50,
             linewidth = 1.8, colour = "black") +
    annotate("text", x = x_rng[2] - 300, y = y_rng[1] + 50,
             label = "500 \u00b5m", vjust = -0.8, size = 2.8,
             colour = "black") +
    labs(title = "Laminar architecture of mouse primary visual cortex",
         subtitle = paste0(format(nrow(pts), big.mark = ","),
                           " MERFISH transcripts across ",
                           length(unique(pts$gene)),
                           " genes (Moffitt et al. 2018)"),
         x = expression(italic(x)~"("*mu*"m)"),
         y = expression(italic(y)~"("*mu*"m)")) +
    th(9) +
    theme(legend.position = c(0.98, 0.5),
          legend.justification = c(1, 0.5),
          legend.background = element_rect(fill = alpha("white", 0.85),
                                           colour = NA),
          legend.key.size = unit(0.4, "cm"))

ggsave(file.path(od, "fig1_store_reading.png"), fig1,
       width = 180, height = 100, units = "mm", dpi = 300, bg = "white")
cat("  saved (", round(file.size(file.path(od, "fig1_store_reading.png"))/1024),
    "KB)\n")

## ======================================================================
## FIG 2: Spatial bounding-box query isolates region of interest
## ======================================================================
cat("--- Fig 2: Spatial query ---\n")

cx <- median(pts$x); cy <- median(pts$y); hw <- 200
qx <- c(cx - hw, cx + hw); qy <- c(cy - hw, cy + hw)

pts_roi <- as.data.frame(bboxQuery(
    spatialPoints(sd)[["transcripts"]],
    qx[1], qx[2], qy[1], qy[2]))
n_in <- nrow(pts_roi)

## Panel a: tissue overview with ROI box
## Use the same binned layer map as Fig 1 for consistency
p2a <- ggplot(bin_layer, aes(x = xbin, y = ybin, fill = layer_f)) +
    geom_raster() +
    scale_fill_manual(values = fig1_cols, guide = "none") +
    annotate("rect", xmin = qx[1], xmax = qx[2],
             ymin = qy[1], ymax = qy[2],
             fill = NA, colour = "black", linewidth = 0.8,
             linetype = "solid") +
    coord_equal(expand = FALSE) +
    labs(x = expression(italic(x)~"("*mu*"m)"),
         y = expression(italic(y)~"("*mu*"m)")) +
    th(8)

## Panel b: zoomed ROI — transcript density by gene
## Assign layer to ROI points
pts_roi$layer <- pts$layer[match(
    paste0(round(pts_roi$x, 2), "_", round(pts_roi$y, 2)),
    pts$key)]

bin_roi <- 5
pts_roi_cx <- pts_roi[!is.na(pts_roi$layer) &
                       pts_roi$layer %in% cortex_layers, ]
pts_roi_cx$xbr <- round(pts_roi_cx$x / bin_roi) * bin_roi
pts_roi_cx$ybr <- round(pts_roi_cx$y / bin_roi) * bin_roi
roi_layer <- aggregate(layer ~ xbr + ybr, data = pts_roi_cx,
                        FUN = function(x) {
                            tt <- table(x)
                            names(tt)[which.max(tt)]
                        })
roi_layer$layer_f <- factor(roi_layer$layer, levels = cortex_layers)

p2b <- ggplot(roi_layer, aes(x = xbr, y = ybr, fill = layer_f)) +
    geom_raster() +
    scale_fill_manual(values = layer_cols, name = "Layer",
                      labels = layer_labels,
                      na.value = "#F0F0F0") +
    coord_equal(xlim = qx, ylim = qy, expand = FALSE) +
    ## Scale bar
    annotate("segment", x = qx[2] - 110, xend = qx[2] - 10,
             y = qy[1] + 12, yend = qy[1] + 12,
             linewidth = 1.2, colour = "white") +
    annotate("text", x = qx[2] - 60, y = qy[1] + 12,
             label = "100 \u00b5m", vjust = -0.7, size = 2.2,
             colour = "white") +
    labs(x = expression(italic(x)~"("*mu*"m)"), y = NULL) +
    th(8) +
    theme(panel.background = element_rect(fill = "#F0F0F0", colour = NA),
          panel.border = element_rect(colour = "black", linewidth = 0.8,
                                       fill = NA))

fig2 <- (p2a + labs(tag = "a")) + (p2b + labs(tag = "b")) +
    plot_layout(ncol = 2, widths = c(1, 0.75)) +
    plot_annotation(
        title = "Spatial bounding-box query isolates region of interest",
        subtitle = paste0(format(n_in, big.mark = ","), " of ",
                          format(nrow(pts), big.mark = ","),
                          " transcripts retained in 400 \u00d7 400 \u00b5m ROI"),
        theme = theme(
            plot.title = element_text(size = 11, face = "bold"),
            plot.subtitle = element_text(size = 8.5, colour = "grey40",
                                          face = "plain"))) &
    theme(plot.tag = element_text(size = 10, face = "bold"))

ggsave(file.path(od, "fig2_spatial_query.png"), fig2,
       width = 180, height = 95, units = "mm", dpi = 300, bg = "white")
cat("  saved\n")

## ======================================================================
## FIG 3: Layer-specific gene expression from transcript aggregation
## Hierarchically clustered heatmap with cell and gene dendrograms
## ======================================================================
cat("--- Fig 3: Aggregation heatmap ---\n")

counts <- aggregatePoints(
    spatialPoints(sd)[["transcripts"]],
    shapes(sd)[["cell_boundaries"]],
    feature_col = "gene", region_col = "cell_id")
counts_df <- as.data.frame(counts)
cids <- counts_df$cell_id
counts_df$cell_id <- NULL
counts_mat <- as.matrix(counts_df)
rownames(counts_mat) <- cids

## Layer assignment per cell
pts_with_cell <- pts[pts$cell_id %in% as.integer(cids), ]
layer_votes <- tapply(pts_with_cell$layer, pts_with_cell$cell_id, function(x) {
    tt <- table(x)
    names(tt)[which.max(tt)]
})
cell_layers <- layer_votes[as.character(cids)]

## Filter: cortex layer + >5 total counts
keep <- !is.na(cell_layers) &
        cell_layers %in% cortex_layers &
        rowSums(counts_mat) > 5
counts_filt <- counts_mat[keep, ]
layers_filt <- cell_layers[keep]
cat("  Cells after filtering:", nrow(counts_filt), "\n")

## Top 25 genes by variance
gene_vars <- apply(counts_filt, 2, var)
top25 <- names(sort(gene_vars, decreasing = TRUE))[1:25]
mat25 <- counts_filt[, top25]

## Column-scale (z-score per gene)
mat_scaled <- scale(mat25)
mat_scaled[mat_scaled > 2.5] <- 2.5
mat_scaled[mat_scaled < -2.5] <- -2.5

## Hierarchical clustering — BOTH axes
cell_hc <- hclust(dist(mat_scaled), method = "ward.D2")
gene_hc <- hclust(dist(t(mat_scaled)), method = "ward.D2")

cell_order <- cell_hc$labels[cell_hc$order]
gene_order <- gene_hc$labels[gene_hc$order]

## Reorder layers to match cell dendrogram
layers_ordered <- layers_filt[match(cell_order, rownames(mat_scaled))]

## Melt for ggplot
hm_df <- reshape2::melt(mat_scaled)
colnames(hm_df) <- c("Cell", "Gene", "Zscore")
hm_df$Cell <- factor(hm_df$Cell, levels = cell_order)
hm_df$Gene <- factor(hm_df$Gene, levels = gene_order)

## Layer annotation strip
strip_df <- data.frame(
    Cell = factor(cell_order, levels = cell_order),
    Layer = layers_ordered
)

p3_hm <- ggplot(hm_df, aes(x = Gene, y = Cell, fill = Zscore)) +
    geom_raster() +
    scale_fill_distiller(palette = "RdBu", direction = -1,
                         limits = c(-2.5, 2.5),
                         name = "Z-score") +
    labs(x = NULL, y = NULL) +
    th(8) +
    theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 7.5,
                                      face = "italic"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.border = element_rect(colour = "grey30", fill = NA,
                                       linewidth = 0.3))

p3_strip <- ggplot(strip_df, aes(x = 1, y = Cell, fill = Layer)) +
    geom_raster() +
    scale_fill_manual(values = layer_cols, name = "Layer",
                      labels = layer_labels) +
    theme_void() +
    theme(legend.position = "left",
          legend.title = element_text(size = 7.5, face = "bold"),
          legend.text = element_text(size = 7),
          legend.key.size = unit(0.3, "cm"))

fig3 <- p3_strip + p3_hm +
    plot_layout(widths = c(0.03, 1), guides = "collect") +
    plot_annotation(
        title = "Layer-specific gene expression from transcript aggregation",
        subtitle = paste0(nrow(mat25), " cells \u00d7 ", ncol(mat25),
                          " genes (top variable), column-scaled, ",
                          "hierarchically clustered (Ward\u2019s D2)"),
        theme = theme(
            plot.title = element_text(size = 11, face = "bold"),
            plot.subtitle = element_text(size = 8.5, colour = "grey40")))

ggsave(file.path(od, "fig3_aggregation.png"), fig3,
       width = 180, height = 130, units = "mm", dpi = 300, bg = "white")
cat("  saved\n")

## ======================================================================
## FIG 4: Cortical layer identity recovered from aggregated expression
## PCA + layer-marker dot plot (replaces old transform + roundtrip figs)
## ======================================================================
cat("--- Fig 4: Downstream analysis ---\n")

## PCA on log-normalized counts
mat_norm <- log1p(mat25[, colSums(mat25) > 0])
pca <- prcomp(mat_norm, center = TRUE, scale. = TRUE)
pca_df <- data.frame(
    PC1 = pca$x[, 1], PC2 = pca$x[, 2],
    Layer = factor(layers_filt[match(rownames(mat_norm),
                                     names(layers_filt))],
                   levels = cortex_layers)
)
var1 <- round(100 * summary(pca)$importance[2, 1], 1)
var2 <- round(100 * summary(pca)$importance[2, 2], 1)

p4a <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = Layer)) +
    geom_point(size = 2, alpha = 0.85) +
    scale_colour_manual(values = layer_cols, name = "Layer",
                        labels = layer_labels) +
    stat_ellipse(level = 0.68, linewidth = 0.35, linetype = "dashed",
                 show.legend = FALSE) +
    labs(x = paste0("PC1 (", var1, "%)"),
         y = paste0("PC2 (", var2, "%)")) +
    th(8.5) +
    theme(legend.position = c(0.02, 0.02),
          legend.justification = c(0, 0),
          legend.background = element_rect(fill = alpha("white", 0.8),
                                           colour = NA))

## Dot plot: mean expression + fraction detected per layer × gene
## Top 12 genes by F-statistic (most layer-discriminating)
fstats <- sapply(colnames(mat25), function(g) {
    mod <- try(summary(aov(mat25[, g] ~ layers_filt))[[1]][1, 4],
               silent = TRUE)
    if (inherits(mod, "try-error")) return(0)
    mod
})
top12 <- names(sort(fstats, decreasing = TRUE))[1:12]

dot_data <- do.call(rbind, lapply(cortex_layers, function(lay) {
    idx <- which(layers_filt == lay)
    if (length(idx) == 0) return(NULL)
    sub <- counts_filt[idx, top12, drop = FALSE]
    do.call(rbind, lapply(top12, function(g) {
        vals <- sub[, g]
        data.frame(
            Layer = lay,
            Gene = g,
            MeanExpr = mean(vals),
            FracDetected = mean(vals > 0),
            stringsAsFactors = FALSE
        )
    }))
}))
dot_data$Layer <- factor(dot_data$Layer, levels = cortex_layers)
dot_data$Gene <- factor(dot_data$Gene, levels = rev(top12))
## Scale mean expression per gene to [0, 1] for colour
dot_data <- do.call(rbind, lapply(split(dot_data, dot_data$Gene), function(d) {
    rng <- range(d$MeanExpr)
    d$ScaledExpr <- if (diff(rng) == 0) 0.5 else
        (d$MeanExpr - rng[1]) / diff(rng)
    d
}))

p4b <- ggplot(dot_data, aes(x = Layer, y = Gene)) +
    geom_point(aes(size = FracDetected, colour = ScaledExpr)) +
    scale_size_continuous(range = c(0.5, 5), name = "Fraction\ndetected",
                          breaks = c(0.25, 0.5, 0.75, 1.0)) +
    scale_colour_viridis_c(option = "magma", direction = -1,
                           name = "Scaled\nmean expr.",
                           limits = c(0, 1)) +
    scale_x_discrete(labels = layer_labels) +
    labs(x = NULL, y = NULL) +
    th(8.5) +
    theme(axis.text.y = element_text(face = "italic", size = 7.5),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.2))

fig4 <- (p4a + labs(tag = "a")) + (p4b + labs(tag = "b")) +
    plot_layout(ncol = 2, widths = c(1, 0.85)) +
    plot_annotation(
        title = "Cortical layer identity recovered from aggregated expression",
        subtitle = paste0(nrow(pca_df), " cells, ",
                          length(unique(pts$gene)), " genes | ",
                          "PCA of log-counts; dot plot of top 12 ",
                          "layer-discriminating genes (one-way ANOVA F)"),
        theme = theme(
            plot.title = element_text(size = 11, face = "bold"),
            plot.subtitle = element_text(size = 8, colour = "grey40"))) &
    theme(plot.tag = element_text(size = 10, face = "bold"))

ggsave(file.path(od, "fig4_downstream.png"), fig4,
       width = 195, height = 110, units = "mm", dpi = 300, bg = "white")
cat("  saved\n")

## ---- Done ----------------------------------------------------------------
cat("\n=== All figures generated ===\n")
for (f in sort(list.files(od, "^fig[1-4].*\\.png$"))) {
    sz <- round(file.size(file.path(od, f)) / 1024)
    cat("  ", f, " (", sz, "KB)\n")
}
