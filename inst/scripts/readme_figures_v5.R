#!/usr/bin/env Rscript
# ============================================================================
# SpatialDataR README figures — v5
# Comprehensive rewrite: saturated colors, gap-filled tissue, clean labels
# All from real MERFISH VISp (Moffitt 2018), executed locally
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

## ---- Publication theme (high contrast, minimal chrome) --------------------
th <- function(base = 9) {
    theme_classic(base_size = base, base_family = "sans") +
    theme(
        axis.line        = element_line(linewidth = 0.4, colour = "black"),
        axis.ticks       = element_line(linewidth = 0.3, colour = "black"),
        axis.title       = element_text(size = rel(1), colour = "black"),
        axis.text        = element_text(size = rel(0.9), colour = "grey20"),
        plot.title       = element_text(size = rel(1.15), face = "bold",
                                        hjust = 0, colour = "black",
                                        margin = margin(b = 1)),
        plot.subtitle    = element_text(size = rel(0.8), colour = "grey30",
                                        hjust = 0, margin = margin(b = 4)),
        plot.margin      = margin(6, 8, 6, 6),
        legend.title     = element_text(size = rel(0.9), face = "bold"),
        legend.text      = element_text(size = rel(0.8)),
        legend.key.size  = unit(0.35, "cm"),
        legend.background = element_blank(),
        strip.text       = element_text(face = "bold", size = rel(0.9))
    )
}

## High-saturation cortical layer palette (distinguishable in print + screen)
layer_cols <- c(
    "VISp_I"       = "#D62728",  # strong red
    "VISp_II/III"  = "#1F77B4",  # strong blue
    "VISp_IV"      = "#2CA02C",  # strong green
    "VISp_V"       = "#9467BD",  # strong purple
    "VISp_VI"      = "#FF7F0E",  # strong orange
    "VISp_wm"      = "#8C564B"   # strong brown
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
## FIG 1: Laminar architecture — gap-filled tissue map
## ======================================================================
cat("--- Fig 1 ---\n")

bin_size <- 8  # smaller bin for smoother edges
pts$xbin <- round(pts$x / bin_size) * bin_size
pts$ybin <- round(pts$y / bin_size) * bin_size

## Step 1: assign dominant cortex layer per bin (cortex-preferring)
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

## Step 2: morphological gap-filling
## For each empty position surrounded by >= 3 filled neighbors,
## assign the majority neighbor layer. Repeat 3 times.
all_x <- sort(unique(bin_layer$xbin))
all_y <- sort(unique(bin_layer$ybin))
filled <- bin_layer
rownames(filled) <- paste(filled$xbin, filled$ybin)

for (iteration in 1:5) {
    ## Build grid of all possible bins within data range
    grid_x <- seq(min(all_x), max(all_x), by = bin_size)
    grid_y <- seq(min(all_y), max(all_y), by = bin_size)
    all_keys <- paste(rep(grid_x, each = length(grid_y)),
                      rep(grid_y, length(grid_x)))
    existing_keys <- paste(filled$xbin, filled$ybin)
    missing_keys <- setdiff(all_keys, existing_keys)

    new_rows <- list()
    for (mk in missing_keys) {
        coords <- as.numeric(strsplit(mk, " ")[[1]])
        mx <- coords[1]; my <- coords[2]
        ## Check 8 neighbors
        nbr_keys <- paste(
            rep(mx + c(-1, 0, 1) * bin_size, each = 3),
            rep(my + c(-1, 0, 1) * bin_size, 3))
        nbr_keys <- setdiff(nbr_keys, mk)
        nbr_layers <- filled$layer[match(nbr_keys, existing_keys)]
        nbr_layers <- nbr_layers[!is.na(nbr_layers)]
        ## Fill if >= 3 neighbors exist (aggressive gap-filling)
        if (length(nbr_layers) >= 3) {
            tt <- table(nbr_layers)
            new_rows[[mk]] <- data.frame(
                xbin = mx, ybin = my,
                layer = names(tt)[which.max(tt)],
                stringsAsFactors = FALSE)
        }
    }
    if (length(new_rows) > 0) {
        added <- do.call(rbind, new_rows)
        filled <- rbind(filled, added)
        rownames(filled) <- paste(filled$xbin, filled$ybin)
        existing_keys <- paste(filled$xbin, filled$ybin)
        cat("  Gap-fill iteration", iteration, ": added", nrow(added), "bins\n")
    } else {
        cat("  Gap-fill iteration", iteration, ": no gaps to fill\n")
        break
    }
}

bin_layer <- filled

fig1_cols <- c(layer_cols, "Tissue" = "#D9D9D9")  # visible grey, not nearly-white
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
    annotate("segment", x = x_rng[2] - 550, xend = x_rng[2] - 50,
             y = y_rng[1] + 50, yend = y_rng[1] + 50,
             linewidth = 2, colour = "black") +
    annotate("text", x = x_rng[2] - 300, y = y_rng[1] + 50,
             label = "500 \u00b5m", vjust = -0.8, size = 3,
             colour = "black", fontface = "bold") +
    labs(title = "Laminar architecture of mouse primary visual cortex",
         subtitle = paste0(format(nrow(pts), big.mark = ","),
                           " MERFISH transcripts, ",
                           length(unique(pts$gene)),
                           " genes (Moffitt et al. 2018)"),
         x = expression(italic(x)~"("*mu*"m)"),
         y = expression(italic(y)~"("*mu*"m)")) +
    th(10) +
    theme(legend.position = c(0.02, 0.5),
          legend.justification = c(0, 0.5),
          legend.background = element_rect(fill = alpha("white", 0.9),
                                           colour = "grey70",
                                           linewidth = 0.3),
          legend.key.size = unit(0.45, "cm"),
          legend.text = element_text(size = 9))

ggsave(file.path(od, "fig1_store_reading.png"), fig1,
       width = 200, height = 115, units = "mm", dpi = 300, bg = "white")
cat("  saved\n")

## ======================================================================
## FIG 2: Spatial query — saturated colormap, gap-filled ROI
## ======================================================================
cat("--- Fig 2 ---\n")

cx <- median(pts$x); cy <- median(pts$y); hw <- 200
qx <- c(cx - hw, cx + hw); qy <- c(cy - hw, cy + hw)

pts_roi <- as.data.frame(bboxQuery(
    spatialPoints(sd)[["transcripts"]],
    qx[1], qx[2], qy[1], qy[2]))
n_in <- nrow(pts_roi)

## Panel a: tissue overview with layer colors + ROI box
## Re-use fig1's gap-filled bin_layer for consistent appearance
p2a <- ggplot(bin_layer, aes(x = xbin, y = ybin, fill = layer_f)) +
    geom_raster() +
    scale_fill_manual(values = fig1_cols, guide = "none") +
    annotate("rect", xmin = qx[1], xmax = qx[2],
             ymin = qy[1], ymax = qy[2],
             fill = NA, colour = "black", linewidth = 1,
             linetype = "solid") +
    coord_equal(expand = FALSE) +
    labs(x = expression(italic(x)~"("*mu*"m)"),
         y = expression(italic(y)~"("*mu*"m)")) +
    th(8.5)

## Panel b: zoomed ROI with gap-filled layer raster
pts_roi$layer <- pts$layer[match(
    paste0(round(pts_roi$x, 2), "_", round(pts_roi$y, 2)),
    pts$key)]

## Panel b: just crop from the already gap-filled Fig 1 raster
## This guarantees no gaps (same data, same gap-fill)
roi_data <- bin_layer[bin_layer$xbin >= qx[1] & bin_layer$xbin <= qx[2] &
                       bin_layer$ybin >= qy[1] & bin_layer$ybin <= qy[2], ]

p2b <- ggplot(roi_data, aes(x = xbin, y = ybin, fill = layer_f)) +
    geom_raster() +
    scale_fill_manual(values = fig1_cols, name = "Layer",
                      labels = fig1_labels,
                      na.value = "#D9D9D9") +
    coord_equal(xlim = qx, ylim = qy, expand = FALSE) +
    annotate("segment", x = qx[2] - 110, xend = qx[2] - 10,
             y = qy[1] + 10, yend = qy[1] + 10,
             linewidth = 1.5, colour = "black") +
    annotate("text", x = qx[2] - 60, y = qy[1] + 10,
             label = "100 \u00b5m", vjust = -0.7, size = 2.5,
             colour = "black", fontface = "bold") +
    labs(x = expression(italic(x)~"("*mu*"m)"), y = NULL) +
    th(8.5) +
    theme(panel.border = element_rect(colour = "black", linewidth = 0.8,
                                       fill = NA))

fig2 <- (p2a + labs(tag = "a")) + (p2b + labs(tag = "b")) +
    plot_layout(ncol = 2, widths = c(1, 0.8)) +
    plot_annotation(
        title = "Spatial bounding-box query isolates region of interest",
        subtitle = paste0(format(n_in, big.mark = ","), " of ",
                          format(nrow(pts), big.mark = ","),
                          " transcripts retained in 400 \u00d7 400 \u00b5m ROI"),
        theme = theme(
            plot.title = element_text(size = 11, face = "bold"),
            plot.subtitle = element_text(size = 8.5, colour = "grey30"))) &
    theme(plot.tag = element_text(size = 10, face = "bold"))

ggsave(file.path(od, "fig2_spatial_query.png"), fig2,
       width = 195, height = 100, units = "mm", dpi = 300, bg = "white")
cat("  saved\n")

## ======================================================================
## FIG 3: Heatmap — high-contrast, clustered, clean labels
## ======================================================================
cat("--- Fig 3 ---\n")

counts <- aggregatePoints(
    spatialPoints(sd)[["transcripts"]],
    shapes(sd)[["cell_boundaries"]],
    feature_col = "gene", region_col = "cell_id")
counts_df <- as.data.frame(counts)
cids <- counts_df$cell_id
counts_df$cell_id <- NULL
counts_mat <- as.matrix(counts_df)
rownames(counts_mat) <- cids

pts_with_cell <- pts[pts$cell_id %in% as.integer(cids), ]
layer_votes <- tapply(pts_with_cell$layer, pts_with_cell$cell_id, function(x) {
    tt <- table(x)
    names(tt)[which.max(tt)]
})
cell_layers <- layer_votes[as.character(cids)]

keep <- !is.na(cell_layers) &
        cell_layers %in% cortex_layers &
        rowSums(counts_mat) > 5
counts_filt <- counts_mat[keep, ]
layers_filt <- cell_layers[keep]
cat("  Cells:", nrow(counts_filt), "\n")

## Top 25 by variance
gene_vars <- apply(counts_filt, 2, var)
top25 <- names(sort(gene_vars, decreasing = TRUE))[1:25]
mat25 <- counts_filt[, top25]

## Column z-score
mat_scaled <- scale(mat25)
mat_scaled[mat_scaled > 2] <- 2
mat_scaled[mat_scaled < -2] <- -2

## Hierarchical clustering both axes
cell_hc <- hclust(dist(mat_scaled), method = "ward.D2")
gene_hc <- hclust(dist(t(mat_scaled)), method = "ward.D2")

cell_order <- cell_hc$labels[cell_hc$order]
gene_order <- gene_hc$labels[gene_hc$order]
layers_ordered <- layers_filt[match(cell_order, rownames(mat_scaled))]

## Melt
hm_df <- reshape2::melt(mat_scaled)
colnames(hm_df) <- c("Cell", "Gene", "Zscore")
hm_df$Cell <- factor(hm_df$Cell, levels = cell_order)
hm_df$Gene <- factor(hm_df$Gene, levels = gene_order)

## Layer annotation strip
strip_df <- data.frame(
    Cell = factor(cell_order, levels = cell_order),
    Layer = factor(layers_ordered, levels = cortex_layers)
)

## HIGH-CONTRAST heatmap color scale
## Deep blue → white → deep red (not the washed-out RdBu)
hm_colors <- colorRampPalette(c("#08306B", "#2171B5", "#6BAED6",
                                 "#FFFFFF",
                                 "#FB6A4A", "#CB181D", "#67000D"))(100)

p3_hm <- ggplot(hm_df, aes(x = Gene, y = Cell, fill = Zscore)) +
    geom_raster() +
    scale_fill_gradientn(colours = hm_colors,
                         limits = c(-2, 2),
                         name = "Z-score") +
    labs(x = NULL, y = NULL) +
    th(8) +
    theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 8,
                                      face = "italic", colour = "black"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA,
                                       linewidth = 0.4))

p3_strip <- ggplot(strip_df, aes(x = 1, y = Cell, fill = Layer)) +
    geom_raster() +
    scale_fill_manual(values = layer_cols, name = "Layer",
                      labels = layer_labels) +
    theme_void() +
    theme(legend.position = "left",
          legend.title = element_text(size = 8, face = "bold"),
          legend.text = element_text(size = 7.5),
          legend.key.size = unit(0.35, "cm"))

fig3 <- p3_strip + p3_hm +
    plot_layout(widths = c(0.035, 1), guides = "collect") +
    plot_annotation(
        title = "Layer-resolved gene expression from transcript aggregation",
        subtitle = paste0(length(cell_order), " cells \u00d7 ",
                          length(gene_order),
                          " genes, column Z-scored (clipped \u00b12), ",
                          "hierarchically clustered (Ward\u2019s D2)"),
        theme = theme(
            plot.title = element_text(size = 11, face = "bold"),
            plot.subtitle = element_text(size = 8.5, colour = "grey30")))

ggsave(file.path(od, "fig3_aggregation.png"), fig3,
       width = 195, height = 140, units = "mm", dpi = 300, bg = "white")
cat("  saved\n")

## ======================================================================
## FIG 4: PCA + dot plot — high contrast, proper scaling
## ======================================================================
cat("--- Fig 4 ---\n")

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
    geom_point(size = 2.5, alpha = 0.9) +
    scale_colour_manual(values = layer_cols, name = "Layer",
                        labels = layer_labels) +
    stat_ellipse(level = 0.68, linewidth = 0.5, linetype = "dashed",
                 show.legend = FALSE) +
    labs(x = paste0("PC1 (", var1, "%)"),
         y = paste0("PC2 (", var2, "%)")) +
    th(9) +
    theme(legend.position = c(0.02, 0.02),
          legend.justification = c(0, 0),
          legend.background = element_rect(fill = alpha("white", 0.85),
                                           colour = "grey70",
                                           linewidth = 0.3))

## Dot plot: ANOVA F-statistic top 12 genes
fstats <- sapply(colnames(counts_filt), function(g) {
    mod <- try(summary(aov(counts_filt[, g] ~ layers_filt))[[1]][1, 4],
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
        data.frame(Layer = lay, Gene = g,
                   MeanExpr = mean(vals),
                   FracDetected = mean(vals > 0),
                   stringsAsFactors = FALSE)
    }))
}))
dot_data$Layer <- factor(dot_data$Layer, levels = cortex_layers)
dot_data$Gene <- factor(dot_data$Gene, levels = rev(top12))
## Scale mean expression per gene to [0, 1]
dot_data <- do.call(rbind, lapply(split(dot_data, dot_data$Gene), function(d) {
    rng <- range(d$MeanExpr)
    d$ScaledExpr <- if (diff(rng) == 0) 0.5 else
        (d$MeanExpr - rng[1]) / diff(rng)
    d
}))

## Dot plot with light grey background so all dots visible
p4b <- ggplot(dot_data, aes(x = Layer, y = Gene)) +
    geom_point(aes(size = FracDetected, fill = ScaledExpr),
               shape = 21, colour = "grey20", stroke = 0.5) +
    scale_size_continuous(range = c(1.5, 7), name = "Fraction\ndetected",
                          breaks = c(0.25, 0.5, 0.75, 1.0)) +
    scale_fill_gradientn(
        colours = c("#FEE0D2", "#FC9272", "#DE2D26", "#67000D"),
        name = "Scaled\nmean expr.",
        limits = c(0, 1)) +
    scale_x_discrete(labels = layer_labels) +
    labs(x = NULL, y = NULL) +
    th(9) +
    theme(axis.text.y = element_text(face = "italic", size = 8,
                                      colour = "black"),
          panel.background = element_rect(fill = "#F5F5F5"),
          panel.grid.major = element_line(colour = "white",
                                          linewidth = 0.3))

fig4 <- (p4a + labs(tag = "a")) + (p4b + labs(tag = "b")) +
    plot_layout(ncol = 2, widths = c(1, 0.9)) +
    plot_annotation(
        title = "Cortical layer identity recovered from aggregated counts",
        subtitle = paste0(nrow(pca_df), " cells, ",
                          length(unique(pts$gene)), " genes | ",
                          "dot plot: top 12 layer-discriminating ",
                          "genes (one-way ANOVA)"),
        theme = theme(
            plot.title = element_text(size = 11, face = "bold"),
            plot.subtitle = element_text(size = 8, colour = "grey30"))) &
    theme(plot.tag = element_text(size = 10, face = "bold"))

ggsave(file.path(od, "fig4_downstream.png"), fig4,
       width = 200, height = 115, units = "mm", dpi = 300, bg = "white")
cat("  saved\n")

## ---- Done ----------------------------------------------------------------
cat("\n=== All figures generated ===\n")
for (f in sort(list.files(od, "^fig[1-4].*\\.png$"))) {
    sz <- round(file.size(file.path(od, f)) / 1024)
    cat("  ", f, " (", sz, "KB)\n")
}
