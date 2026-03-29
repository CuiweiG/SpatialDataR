#!/usr/bin/env Rscript
# ============================================================================
# Generate all README figures using the latest SpatialDataR
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
})

## ---- Global theme ---------------------------------------------------------
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

## ---- Output directory -----------------------------------------------------
od <- "man/figures"
if (!dir.exists(od)) dir.create(od, recursive = TRUE)
setwd("C:/Users/win10/SpatialDataR")

## ---- Load data ------------------------------------------------------------
cat("=== Loading real MERFISH SpatialData ===\n")
store <- "C:/Users/win10/merfish_spatialdata.zarr"
sd <- readSpatialData(store)
show(sd)

pts <- as.data.frame(spatialPoints(sd)[["transcripts"]])
shp <- as.data.frame(shapes(sd)[["cell_boundaries"]])
obs <- as.data.frame(tables(sd)[["table"]]$obs)

cat("Transcripts:", format(nrow(pts), big.mark = ","), "\n")
cat("Genes:", length(unique(pts$gene)), "\n")
cat("Cells:", nrow(shp), "\n")

## Load layer annotations from raw CSV
raw <- read.csv("C:/Users/win10/merfish_real.csv", stringsAsFactors = FALSE)
pts$layer <- raw$layer[seq_len(nrow(pts))]

## Gene rankings
gene_freq <- sort(table(pts$gene), decreasing = TRUE)
top6 <- names(gene_freq)[1:6]
## Colorblind-safe Okabe-Ito derived palette
gene_cols <- c("#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "grey80")
names(gene_cols) <- c(top6, "Other")

## ---- Fig 1: readSpatialData() — full-field spatial map --------------------
cat("\n--- Fig 1: Store reading ---\n")
set.seed(42)
pts_sample <- pts[sample(nrow(pts), 60000), ]
pts_sample$gene_top <- ifelse(pts_sample$gene %in% top6,
                               pts_sample$gene, "Other")
pts_sample$gene_top <- factor(pts_sample$gene_top,
                               levels = c(top6, "Other"))
x_rng <- range(pts$x)
y_rng <- range(pts$y)

fig1 <- ggplot(pts_sample, aes(x = x, y = y, colour = gene_top)) +
    geom_point(size = 0.05, alpha = 0.4, stroke = 0) +
    scale_colour_manual(values = gene_cols, name = "Gene") +
    coord_equal() +
    ## Scale bar 500 um
    annotate("segment", x = x_rng[2] - 600, xend = x_rng[2] - 100,
             y = y_rng[1] + 60, yend = y_rng[1] + 60,
             linewidth = 1.3, colour = "black") +
    annotate("text", x = x_rng[2] - 350, y = y_rng[1] + 60,
             label = "500 \u00B5m", vjust = -0.7, size = 2.8,
             fontface = "bold") +
    labs(title = "readSpatialData()",
         subtitle = paste0(format(nrow(pts), big.mark = ","),
                           " transcripts, ",
                           length(unique(pts$gene)),
                           " genes (MERFISH mouse VISp, Moffitt et al. 2018)"),
         x = expression("x ("*mu*"m)"),
         y = expression("y ("*mu*"m)")) +
    guides(colour = guide_legend(
        override.aes = list(size = 2.5, alpha = 1))) +
    th

ggsave(file.path(od, "fig1_store_reading.png"), fig1,
       width = 180, height = 115, units = "mm", dpi = 300, bg = "white")
cat("  fig1_store_reading.png saved\n")

## ---- Fig 2: bboxQuery() — ROI zoom into cortex ---------------------------
cat("--- Fig 2: Spatial query ---\n")
cx <- median(pts$x); cy <- median(pts$y); hw <- 200
qx <- c(cx - hw, cx + hw); qy <- c(cy - hw, cy + hw)

pts_roi <- as.data.frame(bboxQuery(
    spatialPoints(sd)[["transcripts"]],
    qx[1], qx[2], qy[1], qy[2]))
n_in <- nrow(pts_roi)

set.seed(42)
pts_roi_sub <- pts_roi[sample(nrow(pts_roi), min(8000, nrow(pts_roi))), ]
pts_roi_sub$gene_top <- ifelse(pts_roi_sub$gene %in% top6,
                                pts_roi_sub$gene, "Other")
pts_roi_sub$gene_top <- factor(pts_roi_sub$gene_top,
                                levels = c(top6, "Other"))

p2a <- ggplot(pts_sample, aes(x = x, y = y)) +
    geom_point(size = 0.03, alpha = 0.2, colour = "grey50") +
    annotate("rect", xmin = qx[1], xmax = qx[2],
             ymin = qy[1], ymax = qy[2],
             fill = NA, colour = "#D55E00",
             linewidth = 0.8, linetype = "dashed") +
    coord_equal() +
    labs(x = expression("x ("*mu*"m)"),
         y = expression("y ("*mu*"m)")) + th

p2b <- ggplot(pts_roi_sub, aes(x = x, y = y, colour = gene_top)) +
    geom_point(size = 0.5, alpha = 0.65, stroke = 0) +
    scale_colour_manual(values = gene_cols, name = NULL) +
    coord_equal(xlim = qx, ylim = qy, expand = FALSE) +
    ## Scale bar 100 um
    annotate("segment", x = qx[2] - 120, xend = qx[2] - 20,
             y = qy[1] + 15, yend = qy[1] + 15,
             linewidth = 1.2, colour = "black") +
    annotate("text", x = qx[2] - 70, y = qy[1] + 15,
             label = "100 \u00B5m", vjust = -0.7, size = 2.5,
             fontface = "bold") +
    labs(x = expression("x ("*mu*"m)"), y = NULL) +
    guides(colour = guide_legend(
        override.aes = list(size = 2, alpha = 1))) +
    th + theme(axis.text.y = element_blank(),
               axis.ticks.y = element_blank(),
               axis.line.y = element_blank(),
               panel.border = element_rect(
                   colour = "#D55E00", linewidth = 0.8, fill = NA))

fig2 <- (p2a + labs(tag = "a")) + (p2b + labs(tag = "b")) +
    plot_layout(ncol = 2) +
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
       width = 200, height = 105, units = "mm", dpi = 300, bg = "white")
cat("  fig2_spatial_query.png saved\n")

## ---- Fig 3: aggregatePoints() — dot plot by layer -------------------------
cat("--- Fig 3: Aggregation ---\n")
top12 <- names(gene_freq)[1:12]
layers_use <- c("VISp_I", "VISp_II/III", "VISp_IV",
                "VISp_V", "VISp_VI", "VISp_wm")

dot_data <- do.call(rbind, lapply(layers_use, function(lay) {
    sub_l <- pts[pts$layer == lay, ]
    total <- nrow(sub_l)
    do.call(rbind, lapply(top12, function(g) {
        n_g <- sum(sub_l$gene == g)
        data.frame(layer = lay, gene = g,
                   pct = 100 * n_g / total,
                   frac = n_g / total,
                   stringsAsFactors = FALSE)
    }))
}))
dot_data$layer <- factor(dot_data$layer, levels = rev(layers_use))
dot_data$gene <- factor(dot_data$gene, levels = top12)

fig3 <- ggplot(dot_data, aes(x = gene, y = layer)) +
    geom_point(aes(size = pct, colour = frac)) +
    scale_size_continuous(name = "% of\ntranscripts",
                          range = c(1, 10),
                          breaks = c(1, 3, 5, 8)) +
    scale_colour_viridis_c(name = "Fraction", option = "viridis") +
    labs(title = "aggregatePoints()",
         subtitle = paste0("Gene enrichment across ",
                           length(layers_use), " cortical layers (",
                           format(nrow(pts), big.mark = ","),
                           " transcripts)"),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 9.5) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9,
                                    face = "italic"),
        axis.text.y = element_text(size = 9),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 7),
        plot.title = element_text(size = 12.5, face = "bold"),
        plot.subtitle = element_text(size = 8.5, colour = "grey30",
                                      face = "italic"),
        plot.margin = margin(8, 10, 8, 8))

ggsave(file.path(od, "fig3_aggregation.png"), fig3,
       width = 165, height = 105, units = "mm", dpi = 300, bg = "white")
cat("  fig3_aggregation.png saved\n")

## ---- Fig 4: composeTransforms() — fixed with ggrepel ---------------------
cat("--- Fig 4: Transforms ---\n")
set.seed(42)
shp_ann <- merge(shp, obs, by = "cell_id")
rep_cells <- shp_ann[sample(nrow(shp_ann), 5), ]

ct_comp <- composeTransforms(
    CoordinateTransform("affine", affine = diag(c(0.5, 0.5, 1))),
    CoordinateTransform("affine",
        affine = matrix(c(1,0,500, 0,1,2000, 0,0,1),
                        nrow = 3, byrow = TRUE)))
inv <- invertTransform(ct_comp)

rep_out <- as.data.frame(transformCoords(
    DataFrame(x = rep_cells$x, y = rep_cells$y), ct_comp))
rep_rt <- as.data.frame(transformCoords(
    DataFrame(x = rep_out$x, y = rep_out$y), inv))

rt_err <- max(abs(rep_cells$x - rep_rt$x),
              abs(rep_cells$y - rep_rt$y))

arrow_df <- data.frame(x = rep_cells$x, y = rep_cells$y,
                        xend = rep_out$x, yend = rep_out$y)
ids <- LETTERS[seq_len(nrow(rep_cells))]

fig4 <- ggplot() +
    ## Original (pixels)
    geom_point(data = rep_cells, aes(x = x, y = y),
               shape = 4, size = 4, stroke = 1.1, colour = "#555555") +
    geom_text_repel(data = data.frame(rep_cells, lab = paste0(ids, " [px]")),
                    aes(x = x, y = y, label = lab),
                    size = 3, colour = "#555555", segment.size = 0.25,
                    max.overlaps = 20, seed = 42,
                    nudge_y = 80, box.padding = 0.5) +
    ## Arrows
    geom_segment(data = arrow_df,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 arrow = arrow(length = unit(5, "pt"), type = "closed"),
                 colour = "grey50", linewidth = 0.4, alpha = 0.7) +
    ## Transformed (global)
    geom_point(data = rep_out, aes(x = x, y = y),
               shape = 16, size = 4, colour = "#D55E00") +
    geom_text_repel(data = data.frame(rep_out, lab = paste0(ids, " [gl]")),
                    aes(x = x, y = y, label = lab),
                    size = 3, colour = "#D55E00", fontface = "bold",
                    segment.size = 0.25, max.overlaps = 20, seed = 42,
                    nudge_y = -80, box.padding = 0.5) +
    ## Transform label
    annotate("label",
             x = min(c(rep_cells$x, rep_out$x)),
             y = min(c(rep_cells$y, rep_out$y)) - 200,
             label = paste0("T = scale(0.5) * translate(500, 2000)\n",
                            "roundtrip error = ",
                            formatC(rt_err, format = "e", digits = 1)),
             size = 3, fill = "grey97", colour = "grey30",
             fontface = "bold", hjust = 0,
             label.r = unit(0.15, "lines")) +
    coord_equal() +
    labs(title = "composeTransforms()",
         subtitle = paste0("Real MERFISH coordinates: pixel [px] -> ",
                           "global [gl] via composed affine"),
         x = expression("x ("*mu*"m)"),
         y = expression("y ("*mu*"m)")) +
    th

ggsave(file.path(od, "fig4_transforms.png"), fig4,
       width = 145, height = 135, units = "mm", dpi = 300, bg = "white")
cat("  fig4_transforms.png saved\n")

## ---- Fig 5: Roundtrip (read -> query -> write -> verify) ------------------
cat("--- Fig 5: Roundtrip ---\n")

## Use a 600x600 um ROI
hw5 <- 300
qx5 <- c(cx - hw5, cx + hw5)
qy5 <- c(cy - hw5, cy + hw5)

sub_sd <- bboxQuery(sd, qx5[1], qx5[2], qy5[1], qy5[2])
sub_pts <- as.data.frame(spatialPoints(sub_sd)[["transcripts"]])
n_sub <- nrow(sub_pts)

## Write and read back
tmp_zarr <- file.path(tempdir(), "roundtrip_test.zarr")
writeSpatialData(sub_sd, tmp_zarr)
sd2 <- readSpatialData(tmp_zarr)
verify_pts <- as.data.frame(spatialPoints(sd2)[["transcripts"]])
n_verify <- nrow(verify_pts)

cat("  Written:", n_sub, "transcripts\n")
cat("  Read back:", n_verify, "transcripts\n")
cat("  Match:", n_sub == n_verify, "\n")

## Subsample for plotting
set.seed(42)
sub_pts_plot <- sub_pts[sample(nrow(sub_pts), min(10000, nrow(sub_pts))), ]
sub_pts_plot$gene_top <- ifelse(sub_pts_plot$gene %in% top6,
                                 sub_pts_plot$gene, "Other")
verify_plot <- verify_pts[sample(nrow(verify_pts),
                                  min(10000, nrow(verify_pts))), ]
verify_plot$gene_top <- ifelse(verify_plot$gene %in% top6,
                                verify_plot$gene, "Other")

## Panel a: Read original
p5a <- ggplot(pts_sample, aes(x = x, y = y)) +
    geom_point(size = 0.03, alpha = 0.2, colour = "grey50") +
    annotate("rect", xmin = qx5[1], xmax = qx5[2],
             ymin = qy5[1], ymax = qy5[2],
             fill = NA, colour = "#0072B2",
             linewidth = 0.8, linetype = "dashed") +
    coord_equal() +
    labs(title = "readSpatialData()",
         subtitle = paste0(format(nrow(pts), big.mark = ","),
                           " transcripts"),
         x = expression("x ("*mu*"m)"),
         y = expression("y ("*mu*"m)")) +
    th + theme(plot.title = element_text(size = 10, face = "bold"),
               plot.subtitle = element_text(size = 7.5))

## Panel b: Query + write
p5b <- ggplot(sub_pts_plot, aes(x = x, y = y, colour = gene_top)) +
    geom_point(size = 0.2, alpha = 0.5, stroke = 0) +
    scale_colour_manual(values = gene_cols, guide = "none") +
    coord_equal(xlim = qx5, ylim = qy5) +
    labs(title = "bboxQuery() + writeSpatialData()",
         subtitle = paste0(format(n_sub, big.mark = ","),
                           " in 600x600 \u00B5m ROI"),
         x = expression("x ("*mu*"m)"), y = NULL) +
    th + theme(plot.title = element_text(size = 10, face = "bold"),
               plot.subtitle = element_text(size = 7.5),
               panel.border = element_rect(colour = "#0072B2",
                                            linewidth = 0.8, fill = NA),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank(),
               axis.line.y = element_blank())

## Panel c: Verify read-back
p5c <- ggplot(verify_plot, aes(x = x, y = y, colour = gene_top)) +
    geom_point(size = 0.2, alpha = 0.5, stroke = 0) +
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
        title = "Read -> Query -> Write -> Verify",
        subtitle = paste0("Full roundtrip on real MERFISH data: ",
                          format(n_sub, big.mark = ","), "/",
                          format(nrow(pts), big.mark = ","),
                          " transcripts preserved through write + re-read"),
        theme = theme(
            plot.title = element_text(size = 12.5, face = "bold"),
            plot.subtitle = element_text(size = 8.5, colour = "grey30",
                                          face = "italic"))) &
    theme(plot.tag = element_text(size = 11, face = "bold"))

ggsave(file.path(od, "fig5_roundtrip.png"), fig5,
       width = 220, height = 100, units = "mm", dpi = 300, bg = "white")
cat("  fig5_roundtrip.png saved\n")

## ---- Summary --------------------------------------------------------------
cat("\n=== All README figures regenerated ===\n")
for (f in list.files(od, "^fig[1-5].*\\.png$")) {
    sz <- round(file.size(file.path(od, f)) / 1024)
    cat("  ", f, " (", sz, "KB)\n")
}
cat("\nPackage version:", as.character(packageVersion("SpatialDataR")), "\n")
