#!/usr/bin/env Rscript
# ============================================================================
# SpatialDataR README figures v3 — Technically correct, visually compelling
# Real MERFISH data (Moffitt et al. 2018 Science)
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
})

setwd("C:/Users/win10/SpatialDataR")
od <- "man/figures"
if (!dir.exists(od)) dir.create(od, recursive = TRUE)

## ---- Theme ----------------------------------------------------------------
th <- theme_classic(base_size = 9.5) +
    theme(
        axis.line        = element_line(linewidth = 0.3),
        axis.ticks       = element_line(linewidth = 0.25),
        axis.title       = element_text(size = 9.5, face = "bold"),
        axis.text        = element_text(size = 8.5, colour = "black"),
        plot.title       = element_text(size = 12.5, face = "bold", hjust = 0),
        plot.subtitle    = element_text(size = 8.5, colour = "grey30",
                                        face = "italic", hjust = 0),
        plot.margin      = margin(8, 10, 8, 8),
        legend.title     = element_text(size = 8.5, face = "bold"),
        legend.text      = element_text(size = 7.5),
        legend.background = element_rect(fill = "white", colour = NA)
    )

## ---- Load data ------------------------------------------------------------
cat("=== Loading data ===\n")
store <- "C:/Users/win10/merfish_spatialdata.zarr"
sd <- readSpatialData(store)
show(sd)

pts <- as.data.frame(spatialPoints(sd)[["transcripts"]])
shp <- as.data.frame(shapes(sd)[["cell_boundaries"]])
obs <- as.data.frame(tables(sd)[["table"]]$obs)

raw <- read.csv("C:/Users/win10/merfish_real.csv", stringsAsFactors = FALSE)
pts$layer <- raw$layer[seq_len(nrow(pts))]

cat("Transcripts:", format(nrow(pts), big.mark = ","), "\n")
cat("Genes:", length(unique(pts$gene)), "\n\n")

gene_freq <- sort(table(pts$gene), decreasing = TRUE)
x_rng <- range(pts$x); y_rng <- range(pts$y)

## ====================================================================
## FIG 1: readSpatialData() — cortical layer structure
## KEY FIX: bright saturated colours, outside_VISp nearly invisible,
## higher point density for in-cortex layers
## ====================================================================
cat("--- Fig 1 ---\n")

## Highly saturated, maximally distinct colours for cortical layers
## outside_VISp rendered nearly invisible to let structure dominate
layer_cols <- c(
    "VISp_I"       = "#E41A1C",   # bright red
    "VISp_II/III"  = "#377EB8",   # strong blue
    "VISp_IV"      = "#4DAF4A",   # green
    "VISp_V"       = "#984EA3",   # purple
    "VISp_VI"      = "#FF7F00",   # orange
    "VISp_wm"      = "#A65628",   # brown
    "VISp"         = "#999999",   # medium grey
    "outside_VISp" = "#E8E8E8"    # near-white
)
layer_order <- c("VISp_I", "VISp_II/III", "VISp_IV", "VISp_V",
                 "VISp_VI", "VISp_wm", "VISp", "outside_VISp")

set.seed(42)

## Plot cortical layers ONLY (exclude outside_VISp and generic VISp)
## to make laminar structure visible
cortex_layers <- c("VISp_I", "VISp_II/III", "VISp_IV",
                   "VISp_V", "VISp_VI", "VISp_wm")
pts_cortex <- pts[pts$layer %in% cortex_layers, ]
cat("  Cortical transcripts:", format(nrow(pts_cortex), big.mark = ","), "\n")

## Background: ALL points in very faint grey for tissue outline
set.seed(42)
pts_bg <- pts[sample(nrow(pts), 80000), ]

## Cortex: dense sample, plot ON TOP of background
set.seed(42)
pts_cx <- pts_cortex[sample(nrow(pts_cortex),
                             min(150000, nrow(pts_cortex))), ]
pts_cx$layer_f <- factor(pts_cx$layer, levels = cortex_layers)

fig1 <- ggplot() +
    ## Faint tissue outline
    geom_point(data = pts_bg, aes(x = x, y = y),
               size = 0.005, alpha = 0.03, colour = "grey70", stroke = 0) +
    ## CORTICAL LAYERS: bold, saturated, clearly visible
    geom_point(data = pts_cx, aes(x = x, y = y, colour = layer_f),
               size = 0.15, alpha = 1, stroke = 0) +
    scale_colour_manual(values = layer_cols, name = "Cortical\nlayer",
                        guide = guide_legend(
                            override.aes = list(size = 3.5, alpha = 1))) +
    coord_equal() +
    annotate("segment", x = x_rng[2] - 600, xend = x_rng[2] - 100,
             y = y_rng[1] + 60, yend = y_rng[1] + 60,
             linewidth = 1.3, colour = "black") +
    annotate("text", x = x_rng[2] - 350, y = y_rng[1] + 60,
             label = "500 \u00B5m", vjust = -0.7, size = 2.8,
             fontface = "bold") +
    labs(title = "readSpatialData()",
         subtitle = paste0(format(nrow(pts), big.mark = ","),
                           " transcripts, ", length(unique(pts$gene)),
                           " genes, 8 cortical layers ",
                           "(MERFISH mouse VISp, Moffitt et al. 2018)"),
         x = expression("x ("*mu*"m)"),
         y = expression("y ("*mu*"m)")) +
    th

ggsave(file.path(od, "fig1_store_reading.png"), fig1,
       width = 185, height = 125, units = "mm", dpi = 300, bg = "white")
cat("  saved\n")

## ====================================================================
## FIG 2: bboxQuery() — ROI query
## ====================================================================
cat("--- Fig 2 ---\n")

top6 <- names(gene_freq)[1:6]
gene_cols <- c("#E64B35", "#1F77B4", "#2CA02C",
               "#FF7F0E", "#9467BD", "#D62728", "grey82")
names(gene_cols) <- c(top6, "Other")

cx <- median(pts$x); cy <- median(pts$y); hw <- 200
qx <- c(cx - hw, cx + hw); qy <- c(cy - hw, cy + hw)

pts_roi <- as.data.frame(bboxQuery(
    spatialPoints(sd)[["transcripts"]],
    qx[1], qx[2], qy[1], qy[2]))
n_in <- nrow(pts_roi)

set.seed(42)
pts_s2 <- pts[sample(nrow(pts), 100000), ]
pts_s2$gene_top <- ifelse(pts_s2$gene %in% top6, pts_s2$gene, "Other")
pts_s2$gene_top <- factor(pts_s2$gene_top, levels = c(top6, "Other"))

p2a <- ggplot(pts_s2, aes(x = x, y = y, colour = gene_top)) +
    geom_point(size = 0.02, alpha = 0.3, stroke = 0) +
    scale_colour_manual(values = gene_cols, guide = "none") +
    annotate("rect", xmin = qx[1], xmax = qx[2],
             ymin = qy[1], ymax = qy[2],
             fill = NA, colour = "#D55E00",
             linewidth = 0.9, linetype = "dashed") +
    coord_equal() +
    labs(x = expression("x ("*mu*"m)"),
         y = expression("y ("*mu*"m)")) + th

set.seed(42)
pts_roi_sub <- pts_roi[sample(nrow(pts_roi), min(10000, nrow(pts_roi))), ]
pts_roi_sub$gene_top <- ifelse(pts_roi_sub$gene %in% top6,
                                pts_roi_sub$gene, "Other")
pts_roi_sub$gene_top <- factor(pts_roi_sub$gene_top,
                                levels = c(top6, "Other"))

p2b <- ggplot(pts_roi_sub, aes(x = x, y = y, colour = gene_top)) +
    geom_point(size = 0.5, alpha = 0.7, stroke = 0) +
    scale_colour_manual(values = gene_cols, name = "Gene") +
    coord_equal(xlim = qx, ylim = qy, expand = FALSE) +
    annotate("segment", x = qx[2] - 120, xend = qx[2] - 20,
             y = qy[1] + 15, yend = qy[1] + 15,
             linewidth = 1.2, colour = "black") +
    annotate("text", x = qx[2] - 70, y = qy[1] + 15,
             label = "100 \u00B5m", vjust = -0.7, size = 2.5,
             fontface = "bold") +
    labs(x = expression("x ("*mu*"m)"), y = NULL) +
    guides(colour = guide_legend(
        override.aes = list(size = 2.5, alpha = 1))) +
    th + theme(axis.text.y = element_blank(),
               axis.ticks.y = element_blank(),
               axis.line.y = element_blank(),
               panel.border = element_rect(
                   colour = "#D55E00", linewidth = 0.8, fill = NA))

fig2 <- (p2a + labs(tag = "a")) + (p2b + labs(tag = "b")) +
    plot_layout(ncol = 2, widths = c(1, 0.8)) +
    plot_annotation(
        title = "bboxQuery()",
        subtitle = paste0(format(n_in, big.mark = ","), "/",
                          format(nrow(pts), big.mark = ","),
                          " transcripts in 400\u00D7400 \u00B5m ROI"),
        theme = theme(
            plot.title = element_text(size = 12.5, face = "bold"),
            plot.subtitle = element_text(size = 8.5, colour = "grey30",
                                          face = "italic"))) &
    theme(plot.tag = element_text(size = 11, face = "bold"))

ggsave(file.path(od, "fig2_spatial_query.png"), fig2,
       width = 210, height = 110, units = "mm", dpi = 300, bg = "white")
cat("  saved\n")

## ====================================================================
## FIG 3: aggregatePoints() — show ACTUAL function output
## Use xenium_mini test data where cell_id assignments are correct.
## The MERFISH store has sequential unique IDs per transcript (not
## per cell), making aggregation produce near-empty matrices.
## ====================================================================
cat("--- Fig 3 ---\n")

## Load xenium_mini — this has correct cell_id to transcript mapping
xen_path <- system.file("extdata", "xenium_mini.zarr",
                         package = "SpatialDataR")
xen <- readSpatialData(xen_path)

xen_counts <- aggregatePoints(
    spatialPoints(xen)[["transcripts"]],
    shapes(xen)[["cell_boundaries"]],
    feature_col = "gene",
    region_col = "cell_id")
xen_df <- as.data.frame(xen_counts)
xen_cids <- xen_df$cell_id
xen_df$cell_id <- NULL
xen_mat <- as.matrix(xen_df)
rownames(xen_mat) <- xen_cids

## Cell type annotations
xen_obs <- as.data.frame(tables(xen)[["table"]]$obs)
xen_ct <- setNames(xen_obs$cell_type, xen_obs$cell_id)

## Order cells by cell type
cell_types <- xen_ct[rownames(xen_mat)]
ct_ord <- c("Epithelial", "Endothelial", "Immune", "Stromal")
cell_order <- order(factor(cell_types, levels = ct_ord))
xen_ord <- xen_mat[cell_order, ]
ct_ordered <- cell_types[cell_order]

## Use log1p(counts) — raw counts are more honest and informative
## than Z-score on this small sparse dataset
xen_log <- log1p(xen_ord)

## Gene clustering on log counts
gene_hc <- hclust(dist(t(xen_log)), method = "ward.D2")

## Melt
hm_df <- reshape2::melt(xen_log)
colnames(hm_df) <- c("Cell", "Gene", "Value")
hm_df$Cell <- factor(hm_df$Cell, levels = rownames(xen_ord))
hm_df$Gene <- factor(hm_df$Gene, levels = gene_hc$labels[gene_hc$order])

ct_cols <- c("Epithelial" = "#E64B35", "Endothelial" = "#3C5488",
             "Immune" = "#00A087", "Stromal" = "#4DBBD5")

strip_df <- data.frame(
    Cell = factor(rownames(xen_ord), levels = rownames(xen_ord)),
    CellType = ct_ordered)

p3_strip <- ggplot(strip_df, aes(x = 1, y = Cell, fill = CellType)) +
    geom_tile(colour = NA) +
    scale_fill_manual(values = ct_cols, name = "Cell type") +
    theme_void() +
    theme(legend.position = "left",
          legend.title = element_text(size = 8, face = "bold"),
          legend.text = element_text(size = 7))

p3_hm <- ggplot(hm_df, aes(x = Gene, y = Cell, fill = Value)) +
    geom_tile(colour = NA) +
    scale_fill_viridis_c(option = "inferno", name = "log(count+1)") +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9,
                                      face = "italic"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          legend.position = "right",
          legend.title = element_text(size = 8, face = "bold"),
          legend.text = element_text(size = 7))

fig3 <- p3_strip + p3_hm +
    plot_layout(widths = c(0.04, 1), guides = "collect") +
    plot_annotation(
        title = "aggregatePoints()",
        subtitle = paste0(nrow(xen_ord), " cells \u00D7 ",
                          ncol(xen_ord), " genes \u2192 count matrix, ",
                          "log(count+1), cells grouped by phenotype ",
                          "(Xenium mini, 500 transcripts)"),
        theme = theme(
            plot.title = element_text(size = 12.5, face = "bold"),
            plot.subtitle = element_text(size = 8.5, colour = "grey30",
                                          face = "italic")))

ggsave(file.path(od, "fig3_aggregation.png"), fig3,
       width = 190, height = 130, units = "mm", dpi = 300, bg = "white")
cat("  saved\n")

## ====================================================================
## FIG 4: composeTransforms() + invertTransform()
## ====================================================================
cat("--- Fig 4 ---\n")

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

fig4 <- ggplot() +
    geom_point(data = rep_cells, aes(x = x, y = y),
               shape = 4, size = 3.5, stroke = 1.1, colour = "#3C5488") +
    geom_text_repel(data = data.frame(rep_cells,
                        lab = paste0("C", seq_len(nrow(rep_cells)))),
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
    geom_text_repel(data = data.frame(rep_out,
                        lab = paste0("C", seq_len(nrow(rep_out)))),
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
    coord_equal() +
    labs(title = "composeTransforms() + invertTransform()",
         subtitle = paste0("8 real MERFISH cell coordinates: ",
                           "pixel (\u00D7) \u2192 global (\u25CF) ",
                           "via composed affine"),
         x = expression("x ("*mu*"m)"),
         y = expression("y ("*mu*"m)")) +
    th

ggsave(file.path(od, "fig4_transforms.png"), fig4,
       width = 150, height = 140, units = "mm", dpi = 300, bg = "white")
cat("  saved\n")

## ====================================================================
## FIG 5: Roundtrip
## ====================================================================
cat("--- Fig 5 ---\n")

hw5 <- 300
qx5 <- c(cx - hw5, cx + hw5); qy5 <- c(cy - hw5, cy + hw5)

sub_sd <- bboxQuery(sd, qx5[1], qx5[2], qy5[1], qy5[2])
sub_pts <- as.data.frame(spatialPoints(sub_sd)[["transcripts"]])
n_sub <- nrow(sub_pts)

tmp_zarr <- file.path(tempdir(), "roundtrip_test.zarr")
writeSpatialData(sub_sd, tmp_zarr)
sd2 <- readSpatialData(tmp_zarr)
verify_pts <- as.data.frame(spatialPoints(sd2)[["transcripts"]])
n_verify <- nrow(verify_pts)
cat("  Written:", n_sub, " Read back:", n_verify, " Match:", n_sub == n_verify, "\n")

set.seed(42)
pts_s5 <- pts[sample(nrow(pts), 100000), ]

p5a <- ggplot(pts_s5, aes(x = x, y = y)) +
    geom_point(size = 0.02, alpha = 0.25, colour = "grey40") +
    annotate("rect", xmin = qx5[1], xmax = qx5[2],
             ymin = qy5[1], ymax = qy5[2],
             fill = NA, colour = "#0072B2",
             linewidth = 0.8, linetype = "dashed") +
    coord_equal() +
    labs(title = "readSpatialData()",
         subtitle = paste0(format(nrow(pts), big.mark = ","), " transcripts"),
         x = expression("x ("*mu*"m)"),
         y = expression("y ("*mu*"m)")) +
    th + theme(plot.title = element_text(size = 10, face = "bold"),
               plot.subtitle = element_text(size = 7.5))

set.seed(42)
sub_plot <- sub_pts[sample(nrow(sub_pts), min(12000, nrow(sub_pts))), ]
sub_plot$gene_top <- ifelse(sub_plot$gene %in% top6, sub_plot$gene, "Other")

p5b <- ggplot(sub_plot, aes(x = x, y = y, colour = gene_top)) +
    geom_point(size = 0.2, alpha = 0.55, stroke = 0) +
    scale_colour_manual(values = gene_cols, guide = "none") +
    coord_equal(xlim = qx5, ylim = qy5) +
    labs(title = "bboxQuery() + writeSpatialData()",
         subtitle = paste0(format(n_sub, big.mark = ","),
                           " in 600\u00D7600 \u00B5m ROI"),
         x = expression("x ("*mu*"m)"), y = NULL) +
    th + theme(plot.title = element_text(size = 10, face = "bold"),
               plot.subtitle = element_text(size = 7.5),
               panel.border = element_rect(colour = "#0072B2",
                                            linewidth = 0.8, fill = NA),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank(),
               axis.line.y = element_blank())

set.seed(42)
ver_plot <- verify_pts[sample(nrow(verify_pts), min(12000, nrow(verify_pts))), ]
ver_plot$gene_top <- ifelse(ver_plot$gene %in% top6, ver_plot$gene, "Other")

p5c <- ggplot(ver_plot, aes(x = x, y = y, colour = gene_top)) +
    geom_point(size = 0.2, alpha = 0.55, stroke = 0) +
    scale_colour_manual(values = gene_cols, guide = "none") +
    coord_equal(xlim = qx5, ylim = qy5) +
    labs(title = "readSpatialData() [verify]",
         subtitle = paste0(format(n_verify, big.mark = ","),
                           " transcripts preserved"),
         x = expression("x ("*mu*"m)"), y = NULL) +
    th + theme(plot.title = element_text(size = 10, face = "bold",
                                          colour = "#009E73"),
               plot.subtitle = element_text(size = 7.5),
               panel.border = element_rect(colour = "#009E73",
                                            linewidth = 0.8, fill = NA),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank(),
               axis.line.y = element_blank())

fig5 <- (p5a + labs(tag = "a")) +
    (p5b + labs(tag = "b")) +
    (p5c + labs(tag = "c")) +
    plot_layout(ncol = 3, widths = c(1, 0.65, 0.65)) +
    plot_annotation(
        title = "Read \u2192 Query \u2192 Write \u2192 Verify",
        subtitle = paste0("Full roundtrip: ",
                          format(n_sub, big.mark = ","), "/",
                          format(nrow(pts), big.mark = ","),
                          " transcripts preserved through ",
                          "writeSpatialData() + re-read"),
        theme = theme(
            plot.title = element_text(size = 12.5, face = "bold"),
            plot.subtitle = element_text(size = 8.5, colour = "grey30",
                                          face = "italic"))) &
    theme(plot.tag = element_text(size = 11, face = "bold"))

ggsave(file.path(od, "fig5_roundtrip.png"), fig5,
       width = 225, height = 105, units = "mm", dpi = 300, bg = "white")
cat("  saved\n")

## ---- Done -----------------------------------------------------------------
cat("\n=== All figures regenerated ===\n")
for (f in list.files(od, "^fig[1-5]_.*\\.png$")) {
    sz <- round(file.size(file.path(od, f)) / 1024)
    cat("  ", f, " (", sz, "KB)\n")
}
