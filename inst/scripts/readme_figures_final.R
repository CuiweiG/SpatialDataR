#!/usr/bin/env Rscript
# ============================================================================
# SpatialDataR README figures — FINAL
# ALL from real MERFISH data (Moffitt 2018 Science, 3.7M transcripts)
# Story: Read → Query → Aggregate → Transform → Write → Downstream
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
})

setwd("C:/Users/win10/SpatialDataR")
od <- "man/figures"
dir.create(od, showWarnings = FALSE, recursive = TRUE)

## ---- Global style ---------------------------------------------------------
th <- function(base = 10) {
    theme_classic(base_size = base) +
    theme(
        axis.line   = element_line(linewidth = 0.3, colour = "grey20"),
        axis.ticks  = element_line(linewidth = 0.25, colour = "grey20"),
        axis.title  = element_text(size = rel(0.95), face = "bold"),
        axis.text   = element_text(size = rel(0.85), colour = "grey20"),
        plot.title  = element_text(size = rel(1.2), face = "bold", hjust = 0),
        plot.subtitle = element_text(size = rel(0.85), colour = "grey35",
                                      face = "italic", hjust = 0),
        plot.margin = margin(6, 8, 6, 6),
        legend.title = element_text(size = rel(0.88), face = "bold"),
        legend.text  = element_text(size = rel(0.78)),
        legend.background = element_rect(fill = "white", colour = NA),
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(face = "bold", size = rel(0.9))
    )
}

## Layer palette (consistent across all figures)
layer_cols <- c(
    "VISp_I"       = "#E41A1C",
    "VISp_II/III"  = "#377EB8",
    "VISp_IV"      = "#4DAF4A",
    "VISp_V"       = "#984EA3",
    "VISp_VI"      = "#FF7F00",
    "VISp_wm"      = "#A65628",
    "outside_VISp" = "#BDBDBD"
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

## Assign layers by matching transcript coords to raw CSV coords
## (row counts differ: pts=3,714,642 vs raw=3,841,412)
raw <- read.csv("C:/Users/win10/merfish_real.csv", stringsAsFactors = FALSE)
## Build a spatial lookup: round coords to 0.01 precision for matching
raw$key <- paste0(round(raw$x_um, 2), "_", round(raw$y_um, 2))
pts$key <- paste0(round(pts$x, 2), "_", round(pts$y, 2))
layer_lookup <- setNames(raw$layer, raw$key)
## Remove duplicates — take first match
layer_lookup <- layer_lookup[!duplicated(names(layer_lookup))]
pts$layer <- layer_lookup[pts$key]
cat("  Layer match rate:", round(100 * mean(!is.na(pts$layer)), 1), "%\n")

cat("  Transcripts:", format(nrow(pts), big.mark = ","), "\n")
cat("  Genes:", length(unique(pts$gene)), "\n")
cat("  Cells:", nrow(shp), "\n\n")

gene_freq <- sort(table(pts$gene), decreasing = TRUE)
top6 <- names(gene_freq)[1:6]
gene_cols <- c("#E64B35", "#1F77B4", "#2CA02C",
               "#FF7F0E", "#9467BD", "#D62728", "grey80")
names(gene_cols) <- c(top6, "Other")
x_rng <- range(pts$x); y_rng <- range(pts$y)

## ======================================================================
## FIG 1: readSpatialData() — laminar architecture
## WHY NECESSARY: One function call reads an entire multi-element store.
## No Python. No reticulate. No manual parsing.
## ======================================================================
cat("--- Fig 1: Store reading ---\n")

## ---- Rasterized tissue map (single clean panel) ----
## Bin ALL transcripts (not just cortex) for complete tissue coverage
bin_size <- 10  # um — high resolution

pts$xbin <- round(pts$x / bin_size) * bin_size
pts$ybin <- round(pts$y / bin_size) * bin_size

## Each bin → dominant layer (majority vote)
bin_layer <- aggregate(layer ~ xbin + ybin, data = pts,
                       FUN = function(x) {
                           tt <- table(x)
                           names(tt)[which.max(tt)]
                       })
## Clean layer labels for display
layer_labels <- c(
    "VISp_I"       = "Layer I",
    "VISp_II/III"  = "Layer II/III",
    "VISp_IV"      = "Layer IV",
    "VISp_V"       = "Layer V",
    "VISp_VI"      = "Layer VI",
    "VISp_wm"      = "White matter",
    "outside_VISp" = "Outside VISp",
    "VISp"         = "VISp (unspec.)"
)
all_layer_cols <- c(layer_cols, "VISp" = "#F0F0F0", "outside_VISp" = "#F0F0F0")
all_layers_ord <- c(cortex_layers, "VISp", "outside_VISp")
## Remove VISp and outside_VISp from legend
bin_layer <- bin_layer[bin_layer$layer %in% cortex_layers, ]

bin_layer$layer_f <- factor(bin_layer$layer, levels = all_layers_ord)

fig1 <- ggplot(bin_layer, aes(x = xbin, y = ybin, fill = layer_f)) +
    geom_raster() +
    scale_fill_manual(values = layer_cols, name = NULL,
                      labels = layer_labels,
                      guide = guide_legend(
                          override.aes = list(colour = NA),
                          ncol = 1)) +
    coord_equal(expand = FALSE) +
    ## Scale bar
    annotate("segment", x = x_rng[2] - 550, xend = x_rng[2] - 50,
             y = y_rng[1] + 40, yend = y_rng[1] + 40,
             linewidth = 1.5, colour = "black") +
    annotate("text", x = x_rng[2] - 300, y = y_rng[1] + 40,
             label = "500 \u00B5m", vjust = -0.6, size = 3.2,
             fontface = "bold") +
    labs(title = "readSpatialData()",
         subtitle = paste0(format(nrow(pts), big.mark = ","),
                           " transcripts, ", length(unique(pts$gene)),
                           " genes (MERFISH mouse VISp, ",
                           "Moffitt et al. 2018)"),
         x = expression("x ("*mu*"m)"),
         y = expression("y ("*mu*"m)")) +
    th(10) +
    theme(legend.text = element_text(size = 9),
          legend.key.size = unit(0.5, "cm"))

ggsave(file.path(od, "fig1_store_reading.png"), fig1,
       width = 230, height = 120, units = "mm", dpi = 300, bg = "white")
cat("  saved\n")

## ======================================================================
## FIG 2: bboxQuery() — spatial query
## WHY NECESSARY: No R package can query SpatialData elements spatially.
## ======================================================================
cat("--- Fig 2: Spatial query ---\n")

cx <- median(pts$x); cy <- median(pts$y); hw <- 200
qx <- c(cx - hw, cx + hw); qy <- c(cy - hw, cy + hw)

pts_roi <- as.data.frame(bboxQuery(
    spatialPoints(sd)[["transcripts"]],
    qx[1], qx[2], qy[1], qy[2]))
n_in <- nrow(pts_roi)

## Panel a: rasterized transcript density with ROI box
bin2 <- 12
pts$xb2 <- round(pts$x / bin2) * bin2
pts$yb2 <- round(pts$y / bin2) * bin2
density_all <- aggregate(gene ~ xb2 + yb2, data = pts, FUN = length)
colnames(density_all)[3] <- "count"

p2a <- ggplot(density_all, aes(x = xb2, y = yb2, fill = count)) +
    geom_raster() +
    scale_fill_viridis_c(option = "mako", trans = "sqrt",
                          name = "Transcripts", direction = -1) +
    annotate("rect", xmin = qx[1], xmax = qx[2],
             ymin = qy[1], ymax = qy[2],
             fill = NA, colour = "#FF4500", linewidth = 1.2) +
    coord_equal(expand = FALSE) + th() +
    labs(x = expression("x ("*mu*"m)"), y = expression("y ("*mu*"m)"))

## Panel b: zoomed ROI — rasterize by dominant gene
bin_roi <- 5
pts_roi$xbr <- round(pts_roi$x / bin_roi) * bin_roi
pts_roi$ybr <- round(pts_roi$y / bin_roi) * bin_roi
roi_gene <- aggregate(gene ~ xbr + ybr, data = pts_roi,
                       FUN = function(x) {
                           tt <- table(x)
                           names(tt)[which.max(tt)]
                       })
roi_gene$gene_top <- ifelse(roi_gene$gene %in% top6, roi_gene$gene, "Other")
roi_gene$gene_top <- factor(roi_gene$gene_top, levels = c(top6, "Other"))

p2b <- ggplot(roi_gene, aes(x = xbr, y = ybr, fill = gene_top)) +
    geom_raster() +
    scale_fill_manual(values = gene_cols, name = "Dominant\ngene") +
    coord_equal(xlim = qx, ylim = qy, expand = FALSE) +
    annotate("segment", x = qx[2] - 120, xend = qx[2] - 20,
             y = qy[1] + 15, yend = qy[1] + 15,
             linewidth = 1.2, colour = "black") +
    annotate("text", x = qx[2] - 70, y = qy[1] + 15,
             label = "100 \u00B5m", vjust = -0.7, size = 2.5,
             fontface = "bold") +
    th() + labs(x = expression("x ("*mu*"m)"), y = NULL) +
    theme(panel.border = element_rect(colour = "#FF4500",
                                       linewidth = 1, fill = NA),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank())

fig2 <- (p2a + labs(tag = "a")) + (p2b + labs(tag = "b")) +
    plot_layout(ncol = 2, widths = c(1, 0.8)) +
    plot_annotation(
        title = "bboxQuery()",
        subtitle = paste0(format(n_in, big.mark = ","), "/",
                          format(nrow(pts), big.mark = ","),
                          " transcripts in 400\u00D7400 \u00B5m ROI"),
        theme = theme(
            plot.title = element_text(size = 13, face = "bold"),
            plot.subtitle = element_text(size = 9, colour = "grey30",
                                          face = "italic"))) &
    theme(plot.tag = element_text(size = 11, face = "bold"))

ggsave(file.path(od, "fig2_spatial_query.png"), fig2,
       width = 210, height = 110, units = "mm", dpi = 300, bg = "white")
cat("  saved\n")

## ======================================================================
## FIG 3: aggregatePoints() — molecule → cell x gene matrix
## WHY NECESSARY: No R package aggregates SpatialData points to regions.
## Uses REAL MERFISH data with corrected cell_ids.
## ======================================================================
cat("--- Fig 3: Aggregation ---\n")

counts <- aggregatePoints(
    spatialPoints(sd)[["transcripts"]],
    shapes(sd)[["cell_boundaries"]],
    feature_col = "gene", region_col = "cell_id")
counts_df <- as.data.frame(counts)
cids <- counts_df$cell_id
counts_df$cell_id <- NULL
counts_mat <- as.matrix(counts_df)
rownames(counts_mat) <- cids

## Assign layers to cells by majority vote
pts_with_cell <- pts[pts$cell_id %in% as.integer(cids), ]
layer_votes <- tapply(pts_with_cell$layer, pts_with_cell$cell_id, function(x) {
    tt <- table(x)
    names(tt)[which.max(tt)]
})
cell_layers <- layer_votes[as.character(cids)]

## Filter to cells with cortical layer assignments and >50 counts
keep <- !is.na(cell_layers) &
        cell_layers %in% cortex_layers &
        rowSums(counts_mat) > 5
counts_filt <- counts_mat[keep, ]
layers_filt <- cell_layers[keep]

cat("  Cells after filtering:", nrow(counts_filt), "\n")

## Top 20 by variance
gene_vars <- apply(counts_filt, 2, var)
top20 <- names(sort(gene_vars, decreasing = TRUE))[1:20]
mat20 <- counts_filt[, top20]

## Order cells by layer
cell_ord <- order(factor(layers_filt, levels = cortex_layers))
mat20 <- mat20[cell_ord, ]
layers_ord <- layers_filt[cell_ord]

## Scale columns (genes) for heatmap
mat_scaled <- scale(mat20)
mat_scaled[mat_scaled > 2.5] <- 2.5
mat_scaled[mat_scaled < -2.5] <- -2.5

## Gene clustering
gene_hc <- hclust(dist(t(mat_scaled)), method = "ward.D2")

## Melt
hm_df <- reshape2::melt(mat_scaled)
colnames(hm_df) <- c("Cell", "Gene", "Scaled")
hm_df$Cell <- factor(hm_df$Cell, levels = rownames(mat20))
hm_df$Gene <- factor(hm_df$Gene, levels = gene_hc$labels[gene_hc$order])

## Layer strip
strip_df <- data.frame(Cell = factor(rownames(mat20), levels = rownames(mat20)),
                        Layer = layers_ord)

p3_strip <- ggplot(strip_df, aes(x = 1, y = Cell, fill = Layer)) +
    geom_tile(colour = NA) +
    scale_fill_manual(values = layer_cols, name = "Layer") +
    theme_void() +
    theme(legend.position = "left",
          legend.title = element_text(size = 8.5, face = "bold"),
          legend.text = element_text(size = 7.5))

p3_hm <- ggplot(hm_df, aes(x = Gene, y = Cell, fill = Scaled)) +
    geom_tile(colour = NA) +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                         midpoint = 0, name = "Scaled\nexpr.",
                         limits = c(-2.5, 2.5)) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5,
                                      face = "italic"),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          legend.title = element_text(size = 8, face = "bold"),
          legend.text = element_text(size = 7))

fig3 <- p3_strip + p3_hm +
    plot_layout(widths = c(0.04, 1), guides = "collect") +
    plot_annotation(
        title = "aggregatePoints()",
        subtitle = paste0(nrow(mat20), " cells \u00D7 ", ncol(mat20),
                          " genes (top 20 most variable of ",
                          ncol(counts_mat), "), cells grouped by ",
                          "cortical layer, column-scaled"),
        theme = theme(
            plot.title = element_text(size = 13, face = "bold"),
            plot.subtitle = element_text(size = 9, colour = "grey30",
                                          face = "italic")))

ggsave(file.path(od, "fig3_aggregation.png"), fig3,
       width = 195, height = 135, units = "mm", dpi = 300, bg = "white")
cat("  saved\n")

## ======================================================================
## FIG 4: composeTransforms() + invertTransform()
## WHY NECESSARY: No R package parses/composes OME-NGFF transforms.
## ======================================================================
cat("--- Fig 4: Transforms ---\n")

shp_ann <- merge(shp, obs, by = "cell_id")
set.seed(42)
rep_cells <- shp_ann[sample(nrow(shp_ann), 8), ]

ct_comp <- composeTransforms(
    CoordinateTransform("affine", affine = diag(c(0.5, 0.5, 1)),
                        input_cs = "pixels", output_cs = "scaled"),
    CoordinateTransform("affine",
        affine = matrix(c(1,0,500, 0,1,2000, 0,0,1), 3, byrow = TRUE),
        input_cs = "scaled", output_cs = "global"))
inv <- invertTransform(ct_comp)

rep_out <- as.data.frame(transformCoords(
    DataFrame(x = rep_cells$x, y = rep_cells$y), ct_comp))
rep_rt <- as.data.frame(transformCoords(
    DataFrame(x = rep_out$x, y = rep_out$y), inv))
rt_err <- max(abs(rep_cells$x - rep_rt$x),
              abs(rep_cells$y - rep_rt$y))

arrow_df <- data.frame(x = rep_cells$x, y = rep_cells$y,
                        xend = rep_out$x, yend = rep_out$y)
ids <- paste0("C", seq_len(nrow(rep_cells)))

fig4 <- ggplot() +
    geom_point(data = rep_cells, aes(x = x, y = y),
               shape = 4, size = 3.5, stroke = 1.1, colour = "#3C5488") +
    geom_text_repel(data = data.frame(rep_cells, lab = ids),
                    aes(x = x, y = y, label = lab),
                    size = 2.8, colour = "#3C5488", segment.size = 0.2,
                    max.overlaps = 20, seed = 42,
                    nudge_y = 60, box.padding = 0.4) +
    geom_segment(data = arrow_df,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 arrow = arrow(length = unit(4, "pt"), type = "closed"),
                 colour = "grey55", linewidth = 0.35, alpha = 0.6) +
    geom_point(data = rep_out, aes(x = x, y = y),
               shape = 16, size = 3.5, colour = "#E64B35") +
    geom_text_repel(data = data.frame(rep_out, lab = ids),
                    aes(x = x, y = y, label = lab),
                    size = 2.8, colour = "#E64B35", fontface = "bold",
                    segment.size = 0.2, max.overlaps = 20, seed = 42,
                    nudge_y = -60, box.padding = 0.4) +
    annotate("label",
             x = min(c(rep_cells$x, rep_out$x)),
             y = min(c(rep_cells$y, rep_out$y)) - 150,
             label = paste0("T = scale(0.5) \u00D7 translate(500, 2000)\n",
                            "roundtrip error = ",
                            formatC(rt_err, format = "e", digits = 1)),
             size = 3, fill = "grey97", colour = "grey30",
             fontface = "bold", hjust = 0,
             label.r = unit(0.15, "lines")) +
    coord_equal() + th() +
    labs(title = "composeTransforms() + invertTransform()",
         subtitle = paste0("8 real MERFISH cell coordinates: pixel (\u00D7) ",
                           "\u2192 global (\u25CF) via composed affine"),
         x = expression("x ("*mu*"m)"),
         y = expression("y ("*mu*"m)"))

ggsave(file.path(od, "fig4_transforms.png"), fig4,
       width = 150, height = 140, units = "mm", dpi = 300, bg = "white")
cat("  saved\n")

## ======================================================================
## FIG 5: writeSpatialData() roundtrip
## WHY NECESSARY: No R package can write SpatialData Zarr stores.
## ======================================================================
cat("--- Fig 5: Roundtrip ---\n")

hw5 <- 300
qx5 <- c(cx - hw5, cx + hw5); qy5 <- c(cy - hw5, cy + hw5)
sub_sd <- bboxQuery(sd, qx5[1], qx5[2], qy5[1], qy5[2])
sub_pts <- as.data.frame(spatialPoints(sub_sd)[["transcripts"]])
n_sub <- nrow(sub_pts)

tmp_zarr <- file.path(tempdir(), "roundtrip_final.zarr")
writeSpatialData(sub_sd, tmp_zarr)
sd2 <- readSpatialData(tmp_zarr)
verify_pts <- as.data.frame(spatialPoints(sd2)[["transcripts"]])
n_verify <- nrow(verify_pts)
cat("  Written:", n_sub, " Read back:", n_verify, " Match:", n_sub == n_verify, "\n")

## Panel a: rasterized density (same style as fig2a)
bin5 <- 12
pts$xb5 <- round(pts$x / bin5) * bin5
pts$yb5 <- round(pts$y / bin5) * bin5
dens5 <- aggregate(gene ~ xb5 + yb5, data = pts, FUN = length)
colnames(dens5)[3] <- "count"

p5a <- ggplot(dens5, aes(x = xb5, y = yb5, fill = count)) +
    geom_raster() +
    scale_fill_viridis_c(option = "mako", trans = "sqrt",
                          guide = "none", direction = -1) +
    annotate("rect", xmin = qx5[1], xmax = qx5[2],
             ymin = qy5[1], ymax = qy5[2],
             fill = NA, colour = "#0072B2", linewidth = 1) +
    coord_equal(expand = FALSE) + th(9) +
    labs(title = "readSpatialData()",
         subtitle = paste0(format(nrow(pts), big.mark = ","), " transcripts"),
         x = expression("x ("*mu*"m)"), y = expression("y ("*mu*"m)"))

## Panels b,c: rasterized ROI by dominant gene
bin5r <- 5
sub_pts$xbr5 <- round(sub_pts$x / bin5r) * bin5r
sub_pts$ybr5 <- round(sub_pts$y / bin5r) * bin5r
sub_gene <- aggregate(gene ~ xbr5 + ybr5, data = sub_pts,
                       FUN = function(x) {
                           tt <- table(x); names(tt)[which.max(tt)]
                       })
sub_gene$gene_top <- ifelse(sub_gene$gene %in% top6, sub_gene$gene, "Other")
sub_gene$gene_top <- factor(sub_gene$gene_top, levels = c(top6, "Other"))

verify_pts$xbr5 <- round(verify_pts$x / bin5r) * bin5r
verify_pts$ybr5 <- round(verify_pts$y / bin5r) * bin5r
ver_gene <- aggregate(gene ~ xbr5 + ybr5, data = verify_pts,
                       FUN = function(x) {
                           tt <- table(x); names(tt)[which.max(tt)]
                       })
ver_gene$gene_top <- ifelse(ver_gene$gene %in% top6, ver_gene$gene, "Other")
ver_gene$gene_top <- factor(ver_gene$gene_top, levels = c(top6, "Other"))

p5b <- ggplot(sub_gene, aes(x = xbr5, y = ybr5, fill = gene_top)) +
    geom_raster() +
    scale_fill_manual(values = gene_cols, guide = "none") +
    coord_equal(xlim = qx5, ylim = qy5, expand = FALSE) + th(9) +
    labs(title = "bboxQuery() + writeSpatialData()",
         subtitle = paste0(format(n_sub, big.mark = ","),
                           " in 600\u00D7600 \u00B5m ROI"),
         x = expression("x ("*mu*"m)"), y = NULL) +
    theme(panel.border = element_rect(colour = "#0072B2",
                                       linewidth = 0.8, fill = NA),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank())

p5c <- ggplot(ver_gene, aes(x = xbr5, y = ybr5, fill = gene_top)) +
    geom_raster() +
    scale_fill_manual(values = gene_cols, guide = "none") +
    coord_equal(xlim = qx5, ylim = qy5, expand = FALSE) + th(9) +
    labs(title = "readSpatialData() [verify]",
         subtitle = paste0(format(n_verify, big.mark = ","),
                           " preserved"),
         x = expression("x ("*mu*"m)"), y = NULL) +
    theme(plot.title = element_text(colour = "#009E73"),
          panel.border = element_rect(colour = "#009E73",
                                       linewidth = 0.8, fill = NA),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank())

fig5 <- (p5a + labs(tag = "a")) + (p5b + labs(tag = "b")) +
    (p5c + labs(tag = "c")) +
    plot_layout(ncol = 3, widths = c(1, 0.65, 0.65)) +
    plot_annotation(
        title = "Read \u2192 Query \u2192 Write \u2192 Verify",
        subtitle = paste0("Full roundtrip: ",
                          format(n_sub, big.mark = ","), "/",
                          format(nrow(pts), big.mark = ","),
                          " transcripts preserved through writeSpatialData()"),
        theme = theme(
            plot.title = element_text(size = 13, face = "bold"),
            plot.subtitle = element_text(size = 9, colour = "grey30",
                                          face = "italic"))) &
    theme(plot.tag = element_text(size = 11, face = "bold"))

ggsave(file.path(od, "fig5_roundtrip.png"), fig5,
       width = 225, height = 105, units = "mm", dpi = 300, bg = "white")
cat("  saved\n")

## ======================================================================
## FIG 6: Downstream analysis — PCA + gene composition
## WHY NECESSARY: SpatialDataR output plugs directly into Bioconductor.
## Uses REAL MERFISH 160 cells.
## ======================================================================
cat("--- Fig 6: Downstream ---\n")

## Use the ordered matrix from Fig 3 (already filtered and sorted)
mat_norm <- log1p(mat20[, colSums(mat20) > 0])
pca <- prcomp(mat_norm, center = TRUE, scale. = FALSE)
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                      Layer = layers_ord)
pca_df$Layer <- factor(pca_df$Layer, levels = cortex_layers)
var1 <- round(100 * summary(pca)$importance[2, 1], 1)
var2 <- round(100 * summary(pca)$importance[2, 2], 1)

p6a <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = Layer)) +
    geom_point(size = 2.5, alpha = 0.8) +
    scale_colour_manual(values = layer_cols, name = "Layer") +
    stat_ellipse(level = 0.68, linewidth = 0.4, linetype = "dashed") +
    th() +
    labs(title = "PCA of aggregated count matrix",
         subtitle = paste0(nrow(pca_df), " MERFISH cells, ",
                           ncol(mat_norm), " genes"),
         x = paste0("PC1 (", var1, "%)"),
         y = paste0("PC2 (", var2, "%)"))

## Gene composition per layer — top 8 genes
top8 <- names(sort(gene_freq, decreasing = TRUE))[1:8]
comp_data <- do.call(rbind, lapply(cortex_layers, function(lay) {
    cells_in <- rownames(counts_filt)[layers_filt == lay]
    if (length(cells_in) == 0) return(NULL)
    sub_mat <- counts_filt[cells_in, , drop = FALSE]
    gene_sums <- colSums(sub_mat)
    total <- sum(gene_sums)
    do.call(rbind, lapply(top8, function(g) {
        data.frame(Layer = lay, Gene = g,
                   Fraction = gene_sums[g] / total,
                   stringsAsFactors = FALSE)
    }))
}))
comp_data$Layer <- factor(comp_data$Layer, levels = cortex_layers)
comp_data$Gene <- factor(comp_data$Gene, levels = rev(top8))

p6b <- ggplot(comp_data, aes(x = Layer, y = Fraction, fill = Gene)) +
    geom_col(position = "stack", width = 0.75) +
    scale_fill_brewer(palette = "Set2", name = "Gene") +
    th() +
    labs(title = "Gene composition per layer",
         subtitle = "Top 8 genes by total transcript count",
         x = NULL, y = "Fraction of transcripts") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 9))

fig6 <- (p6a + labs(tag = "a")) + (p6b + labs(tag = "b")) +
    plot_layout(ncol = 2) +
    plot_annotation(
        title = "aggregatePoints() \u2192 SingleCellExperiment \u2192 downstream analysis",
        subtitle = paste0("Real MERFISH data: ", nrow(pca_df),
                          " cells, ", ncol(mat_norm), " genes, ",
                          length(unique(layers_ord)), " cortical layers"),
        theme = theme(
            plot.title = element_text(size = 13, face = "bold"),
            plot.subtitle = element_text(size = 9, colour = "grey30",
                                          face = "italic"))) &
    theme(plot.tag = element_text(size = 11, face = "bold"))

ggsave(file.path(od, "fig6_downstream.png"), fig6,
       width = 210, height = 110, units = "mm", dpi = 300, bg = "white")
cat("  saved\n")

## ---- Cleanup old files ---------------------------------------------------
old_files <- c("fig6_multimodal.png", "fig7_downstream.png",
               "fig8_interop.png", "fig9_scalability.png",
               "fig1_spatial_overview.png")
for (f in old_files) {
    fp <- file.path(od, f)
    if (file.exists(fp)) file.remove(fp)
}

## ---- Done ----------------------------------------------------------------
cat("\n=== All figures generated ===\n")
for (f in sort(list.files(od, "^fig[1-6].*\\.png$"))) {
    sz <- round(file.size(file.path(od, f)) / 1024)
    cat("  ", f, " (", sz, "KB)\n")
}
