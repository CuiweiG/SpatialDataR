#!/usr/bin/env Rscript
# ============================================================================
# SpatialDataR: Publication Figure from Real MERFISH Data
# Moffitt et al. 2018, Science — Mouse brain VISp
# 3.7M transcripts, 268 genes, 160 cells, 8 cortical layers
# ============================================================================

## ---- Setup ----------------------------------------------------------------
.libPaths(c("C:/Users/win10/R/win-library/4.4",
            "C:/Users/win10/AppData/Local/R/win-library/4.4",
            .libPaths()))
suppressPackageStartupMessages({
    library(SpatialDataR)
    library(S4Vectors)
    library(ggplot2)
    library(patchwork)
    library(viridis)
    library(scales)
    library(reshape2)
    library(grid)
    library(gridExtra)
    library(ggrepel)
    library(RColorBrewer)
})

## ---- Theme ----------------------------------------------------------------
theme_pub <- function(base_size = 10) {
    theme_classic(base_size = base_size) %+replace%
    theme(
        text             = element_text(family = "sans", colour = "grey10"),
        plot.title       = element_text(face = "bold", size = rel(1.15),
                                        hjust = 0, margin = margin(b = 3)),
        plot.subtitle    = element_text(size = rel(0.82), colour = "grey40",
                                        hjust = 0, margin = margin(b = 5),
                                        face = "italic"),
        axis.title       = element_text(size = rel(0.9), face = "bold"),
        axis.text        = element_text(size = rel(0.82), colour = "grey30"),
        axis.line        = element_line(linewidth = 0.35, colour = "grey30"),
        axis.ticks       = element_line(linewidth = 0.25, colour = "grey30"),
        legend.title     = element_text(size = rel(0.88), face = "bold"),
        legend.text      = element_text(size = rel(0.78)),
        legend.key.size  = unit(0.4, "cm"),
        legend.background = element_rect(fill = "white", colour = NA),
        plot.margin      = margin(6, 10, 6, 8),
        panel.grid       = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text       = element_text(face = "bold", size = rel(0.85))
    )
}

## ---- Layer colours (anatomically ordered) ---------------------------------
layer_cols <- c(
    "outside_VISp" = "#BDBDBD",
    "VISp"         = "#969696",
    "VISp_I"       = "#FBB4AE",
    "VISp_II/III"  = "#B3CDE3",
    "VISp_IV"      = "#CCEBC5",
    "VISp_V"       = "#DECBE4",
    "VISp_VI"      = "#FED9A6",
    "VISp_wm"      = "#E5D8BD"
)
layer_order <- c("VISp_I", "VISp_II/III", "VISp_IV", "VISp_V",
                 "VISp_VI", "VISp_wm", "VISp", "outside_VISp")

## ---- Load data ------------------------------------------------------------
cat("=== Reading real MERFISH SpatialData store ===\n")
store <- "C:/Users/win10/merfish_spatialdata.zarr"
sd <- readSpatialData(store)
show(sd)

pts_df <- as.data.frame(spatialPoints(sd)[["transcripts"]])
shp_df <- as.data.frame(shapes(sd)[["cell_boundaries"]])
obs_df <- as.data.frame(tables(sd)[["table"]]$obs)

cat("\nTranscripts:", format(nrow(pts_df), big.mark = ","), "\n")
cat("Genes:", length(unique(pts_df$gene)), "\n")
cat("Cells:", nrow(shp_df), "\n")

raw <- read.csv("C:/Users/win10/merfish_real.csv", stringsAsFactors = FALSE)
pts_df$layer <- raw$layer[seq_len(nrow(pts_df))]
pts_df$depth_um <- raw$depth_um[seq_len(nrow(pts_df))]
shp_ann <- merge(shp_df, obs_df, by = "cell_id")

val <- validateSpatialData(store)
cat("Validation: valid =", val$valid, "\n\n")

## ---- Gene rankings --------------------------------------------------------
gene_freq <- sort(table(pts_df$gene), decreasing = TRUE)
top6 <- names(gene_freq)[1:6]
top6_cols <- c("#E64B35", "#4DBBD5", "#00A087",
               "#3C5488", "#F39B7F", "#8491B4")
names(top6_cols) <- top6
all_cols <- c(top6_cols, "Other" = "grey80")

x_rng <- range(pts_df$x)
y_rng <- range(pts_df$y)

## ---- Panel A: Full-field transcript spatial map ---------------------------
cat("Panel A...\n")
set.seed(42)
pts_sub <- pts_df[sample(nrow(pts_df), 80000), ]
pts_sub$gene_top <- ifelse(pts_sub$gene %in% top6, pts_sub$gene, "Other")
pts_sub$gene_top <- factor(pts_sub$gene_top, levels = c(top6, "Other"))

pa <- ggplot(pts_sub, aes(x = x, y = y, colour = gene_top)) +
    geom_point(size = 0.04, alpha = 0.35, stroke = 0) +
    scale_colour_manual(values = all_cols, name = "Gene",
                        guide = guide_legend(
                            override.aes = list(size = 2.5, alpha = 1))) +
    coord_equal() +
    annotate("segment",
             x = x_rng[2] - 600, xend = x_rng[2] - 100,
             y = y_rng[1] + 80, yend = y_rng[1] + 80,
             linewidth = 1.5, colour = "black") +
    annotate("text", x = x_rng[2] - 350, y = y_rng[1] + 80,
             label = "500 \u00B5m", vjust = -0.8, size = 2.8,
             fontface = "bold") +
    labs(title = "readSpatialData()",
         subtitle = paste0(format(nrow(pts_df), big.mark = ","),
                           " transcripts, ", length(unique(pts_df$gene)),
                           " genes (MERFISH mouse VISp, Moffitt 2018)"),
         x = expression("x ("*mu*"m)"),
         y = expression("y ("*mu*"m)")) +
    theme_pub()

## ---- Panel B: Cortical layer map ------------------------------------------
cat("Panel B...\n")
pts_sub$layer_f <- factor(pts_sub$layer, levels = layer_order)

pb <- ggplot(pts_sub[!is.na(pts_sub$layer_f), ],
             aes(x = x, y = y, colour = layer_f)) +
    geom_point(size = 0.04, alpha = 0.3, stroke = 0) +
    scale_colour_manual(values = layer_cols, name = "Layer",
                        guide = guide_legend(
                            override.aes = list(size = 3, alpha = 1))) +
    coord_equal() +
    annotate("segment",
             x = x_rng[2] - 600, xend = x_rng[2] - 100,
             y = y_rng[1] + 80, yend = y_rng[1] + 80,
             linewidth = 1.5, colour = "black") +
    annotate("text", x = x_rng[2] - 350, y = y_rng[1] + 80,
             label = "500 \u00B5m", vjust = -0.8, size = 2.8,
             fontface = "bold") +
    labs(title = "Cortical layer annotation",
         subtitle = "8 layers: VISp I-VI + white matter + outside",
         x = expression("x ("*mu*"m)"),
         y = expression("y ("*mu*"m)")) +
    theme_pub()

## ---- Panel C: bboxQuery — no inset (avoids extra tag) ---------------------
cat("Panel C...\n")
cx <- median(pts_df$x); cy <- median(pts_df$y); hw <- 200
bb <- list(xmin = cx - hw, xmax = cx + hw, ymin = cy - hw, ymax = cy + hw)

pts_roi <- as.data.frame(bboxQuery(
    spatialPoints(sd)[["transcripts"]],
    bb$xmin, bb$xmax, bb$ymin, bb$ymax))
n_roi <- nrow(pts_roi)

roi_idx <- which(pts_df$x >= bb$xmin & pts_df$x <= bb$xmax &
                 pts_df$y >= bb$ymin & pts_df$y <= bb$ymax)

## Mark inside/outside on subsample
pts_sub$in_bbox <- pts_sub$x >= bb$xmin & pts_sub$x <= bb$xmax &
                   pts_sub$y >= bb$ymin & pts_sub$y <= bb$ymax

## Subsample ROI for colour
set.seed(42)
pts_roi_sub <- pts_roi[sample(nrow(pts_roi), min(12000, nrow(pts_roi))), ]
pts_roi_sub$gene_top <- ifelse(pts_roi_sub$gene %in% top6,
                               pts_roi_sub$gene, "Other")

pc <- ggplot() +
    ## Grey outside
    geom_point(data = pts_sub[!pts_sub$in_bbox, ],
               aes(x = x, y = y),
               size = 0.03, alpha = 0.12, colour = "grey65") +
    ## Coloured inside
    geom_point(data = pts_roi_sub,
               aes(x = x, y = y, colour = gene_top),
               size = 0.08, alpha = 0.5) +
    ## Bounding box
    annotate("rect", xmin = bb$xmin, xmax = bb$xmax,
             ymin = bb$ymin, ymax = bb$ymax,
             fill = NA, colour = "#E64B35", linewidth = 0.9,
             linetype = "dashed") +
    annotate("label", x = bb$xmax + 30, y = bb$ymax,
             label = paste0(format(n_roi, big.mark = ","), "/",
                            format(nrow(pts_df), big.mark = ","),
                            "\nretained"),
             size = 2.5, colour = "#E64B35", fontface = "bold",
             fill = "white", hjust = 0, label.r = unit(0.1, "lines")) +
    scale_colour_manual(values = all_cols, guide = "none") +
    coord_equal() +
    annotate("segment",
             x = x_rng[2] - 600, xend = x_rng[2] - 100,
             y = y_rng[1] + 80, yend = y_rng[1] + 80,
             linewidth = 1.5, colour = "black") +
    annotate("text", x = x_rng[2] - 350, y = y_rng[1] + 80,
             label = "500 \u00B5m", vjust = -0.8, size = 2.8,
             fontface = "bold") +
    labs(title = "bboxQuery()",
         subtitle = paste0("400x400 \u00B5m ROI, ",
                           format(n_roi, big.mark = ","), " transcripts"),
         x = expression("x ("*mu*"m)"),
         y = expression("y ("*mu*"m)")) +
    theme_pub()

## ---- Panel D: Dot plot — gene enrichment by layer -------------------------
cat("Panel D...\n")
top12 <- names(gene_freq)[1:12]
layers_use <- c("VISp_I", "VISp_II/III", "VISp_IV", "VISp_V",
                "VISp_VI", "VISp_wm")

dot_data <- do.call(rbind, lapply(layers_use, function(lay) {
    sub_l <- pts_df[pts_df$layer == lay, ]
    total <- nrow(sub_l)
    do.call(rbind, lapply(top12, function(g) {
        n_g <- sum(sub_l$gene == g)
        data.frame(layer = lay, gene = g,
                   count = n_g,
                   pct = 100 * n_g / total,
                   frac = n_g / total,
                   stringsAsFactors = FALSE)
    }))
}))
dot_data$layer <- factor(dot_data$layer, levels = rev(layers_use))
dot_data$gene <- factor(dot_data$gene, levels = top12)

pd <- ggplot(dot_data, aes(x = gene, y = layer)) +
    geom_point(aes(size = pct, colour = frac)) +
    scale_size_continuous(name = "% of layer\ntranscripts",
                          range = c(1, 10),
                          breaks = c(1, 3, 5, 8)) +
    scale_colour_viridis_c(name = "Fraction", option = "viridis",
                            limits = c(0, NA)) +
    labs(title = "Gene enrichment across cortical layers",
         subtitle = paste0("aggregatePoints(): top 12/",
                           length(unique(pts_df$gene)),
                           " genes, 6 VISp layers"),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9,
                                    face = "italic"),
        axis.text.y = element_text(size = 9),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 8.5, face = "bold"),
        legend.text = element_text(size = 7.5),
        plot.title = element_text(face = "bold", size = 11.5, hjust = 0),
        plot.subtitle = element_text(size = 9, colour = "grey40",
                                      hjust = 0, face = "italic"),
        plot.margin = margin(6, 10, 6, 8)
    )

## ---- Panel E: Cell x gene heatmap ----------------------------------------
cat("Panel E...\n")
counts <- aggregatePoints(
    spatialPoints(sd)[["transcripts"]], shapes(sd)[["cell_boundaries"]],
    feature_col = "gene", region_col = "cell_id")
counts_df <- as.data.frame(counts)
rownames(counts_df) <- counts_df$cell_id
counts_df$cell_id <- NULL

## Top 25 genes by total count across cells
gene_totals <- colSums(counts_df)
top25 <- names(sort(gene_totals, decreasing = TRUE))[1:25]
counts_top <- as.matrix(counts_df[, top25])

## Cell type ordering
ct_map <- setNames(obs_df$cell_type, obs_df$cell_id)
cell_order_df <- data.frame(
    cell_id = rownames(counts_top),
    layer = ct_map[rownames(counts_top)],
    stringsAsFactors = FALSE)
cell_order_df <- cell_order_df[order(cell_order_df$layer), ]
counts_mat <- counts_top[cell_order_df$cell_id, ]

## Row-normalize: scale each cell to Z-score across genes for visibility
## This highlights relative enrichment patterns per cell
counts_scaled <- t(scale(t(counts_mat)))
counts_scaled[is.nan(counts_scaled)] <- 0
## Clamp extremes for better colour range
counts_scaled[counts_scaled > 3]  <- 3
counts_scaled[counts_scaled < -3] <- -3

## ggplot heatmap with diverging colour scale
hm_df <- reshape2::melt(counts_scaled)
colnames(hm_df) <- c("Cell", "Gene", "Zscore")
hm_df$Cell <- factor(hm_df$Cell, levels = cell_order_df$cell_id)
hm_df$Layer <- ct_map[as.character(hm_df$Cell)]

## Gene clustering
gene_hc <- hclust(dist(t(counts_scaled)), method = "ward.D2")
hm_df$Gene <- factor(hm_df$Gene, levels = gene_hc$labels[gene_hc$order])

pe_main <- ggplot(hm_df, aes(x = Gene, y = Cell, fill = Zscore)) +
    geom_tile(colour = NA) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "Z-score",
                         limits = c(-3, 3)) +
    labs(x = NULL, y = "Cells (grouped by layer)") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 10,
                                      face = "italic"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "right")

## Layer colour strip
pe_strip <- ggplot(data.frame(Cell = factor(cell_order_df$cell_id,
                                             levels = cell_order_df$cell_id),
                               Layer = cell_order_df$layer),
                    aes(x = 1, y = Cell, fill = Layer)) +
    geom_tile(colour = NA) +
    scale_fill_manual(values = layer_cols, name = "Layer") +
    theme_void() +
    theme(legend.position = "right",
          legend.title = element_text(size = 8.5, face = "bold"),
          legend.text = element_text(size = 7.5))

## Combine strip + heatmap into single grob to avoid extra tag
pe_combined <- pe_strip + pe_main +
    plot_layout(widths = c(0.03, 1), guides = "collect") +
    plot_annotation(
        title = "Cell x gene expression matrix",
        subtitle = paste0(nrow(counts_mat), " cells x ", ncol(counts_mat),
                          " genes | row Z-score normalized"),
        theme = theme(
            plot.title = element_text(face = "bold", size = 11.5, hjust = 0),
            plot.subtitle = element_text(size = 9, colour = "grey40",
                                          hjust = 0, face = "italic")
        )
    )
pe_wrapped <- wrap_elements(full = patchworkGrob(pe_combined))

## ---- Panel F: Transform on real coords ------------------------------------
cat("Panel F...\n")
set.seed(42)
rep_cells <- shp_ann[sample(nrow(shp_ann), 8), ]

scale_tf <- CoordinateTransform("affine",
    affine = diag(c(0.5, 0.5, 1)),
    input_cs = "pixels", output_cs = "scaled")
translate_tf <- CoordinateTransform("affine",
    affine = matrix(c(1,0,500, 0,1,2000, 0,0,1), nrow = 3, byrow = TRUE),
    input_cs = "scaled", output_cs = "global")
composed <- composeTransforms(scale_tf, translate_tf)
inv <- invertTransform(composed)

cell_pts <- DataFrame(x = rep_cells$x, y = rep_cells$y)
cell_global <- as.data.frame(transformCoords(cell_pts, composed))
cell_rt <- as.data.frame(transformCoords(
    DataFrame(x = cell_global$x, y = cell_global$y), inv))
rt_err <- max(abs(rep_cells$x - cell_rt$x),
              abs(rep_cells$y - cell_rt$y))

arrow_df <- data.frame(
    x = rep_cells$x, y = rep_cells$y,
    xend = cell_global$x, yend = cell_global$y)
labels_orig <- data.frame(x = rep_cells$x, y = rep_cells$y,
                           label = paste0("C", seq_len(nrow(rep_cells))))
labels_gl   <- data.frame(x = cell_global$x, y = cell_global$y,
                           label = paste0("C", seq_len(nrow(rep_cells))))

pf <- ggplot() +
    geom_point(data = rep_cells, aes(x = x, y = y),
               shape = 4, size = 3.5, stroke = 1.2, colour = "#3C5488") +
    geom_text_repel(data = labels_orig,
                    aes(x = x, y = y, label = label),
                    size = 2.5, colour = "#3C5488", segment.size = 0.2,
                    max.overlaps = 20, nudge_y = 50) +
    geom_segment(data = arrow_df,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 arrow = arrow(length = unit(4, "pt"), type = "closed"),
                 colour = "grey50", linewidth = 0.3, alpha = 0.6) +
    geom_point(data = cell_global, aes(x = x, y = y),
               shape = 16, size = 3.5, colour = "#E64B35") +
    geom_text_repel(data = labels_gl,
                    aes(x = x, y = y, label = label),
                    size = 2.5, colour = "#E64B35", fontface = "bold",
                    segment.size = 0.2, max.overlaps = 20, nudge_y = -50) +
    annotate("label",
             x = min(c(rep_cells$x, cell_global$x)) - 20,
             y = max(c(rep_cells$y, cell_global$y)) + 150,
             label = "T = scale(0.5) * translate(500,2000)",
             size = 3, fill = "grey97", colour = "grey30",
             fontface = "bold", hjust = 0, label.r = unit(0.15, "lines")) +
    coord_equal() +
    labs(title = "composeTransforms() + invertTransform()",
         subtitle = paste0("Real MERFISH coordinates | roundtrip error = ",
                           formatC(rt_err, format = "e", digits = 1)),
         x = expression("x ("*mu*"m)"),
         y = expression("y ("*mu*"m)")) +
    theme_pub()

## ---- Assemble 6-panel figure ----------------------------------------------
cat("\nAssembling figure...\n")

composite <- (pa | pb) /
             (pc | pd) /
             (pe_wrapped | pf) +
    plot_layout(heights = c(1, 1, 1)) +
    plot_annotation(
        title = "SpatialDataR: Real MERFISH data analysis pipeline",
        subtitle = paste0("Moffitt et al. 2018 Science  |  Mouse VISp  |  ",
                          format(nrow(pts_df), big.mark = ","),
                          " transcripts  |  ",
                          length(unique(pts_df$gene)),
                          " genes  |  ",
                          nrow(shp_df), " cells  |  ",
                          length(unique(raw$layer)), " cortical layers"),
        tag_levels = "a",
        theme = theme(
            plot.title = element_text(face = "bold", size = 16,
                                      hjust = 0, family = "sans"),
            plot.subtitle = element_text(size = 11, colour = "grey40",
                                         hjust = 0, margin = margin(b = 8),
                                         face = "italic")
        )
    ) &
    theme(plot.tag = element_text(face = "bold", size = 14))

## ---- Save -----------------------------------------------------------------
out_dir <- "C:/Users/win10/SpatialDataR/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_png <- file.path(out_dir, "MERFISH_real_data_figure.png")
out_pdf <- file.path(out_dir, "MERFISH_real_data_figure.pdf")

cat("Saving PNG (300 DPI)...\n")
png(out_png, width = 16, height = 17, units = "in", res = 300,
    type = "windows")
print(composite)
dev.off()

cat("Saving PDF...\n")
pdf(out_pdf, width = 16, height = 17)
print(composite)
dev.off()

cat("\n=== COMPLETE ===\n")
cat("PNG:", normalizePath(out_png), "\n")
cat("PDF:", normalizePath(out_pdf), "\n")
cat("6 panels (a-f):\n")
cat("  a: Full-field transcript map (80k/3.7M subsampled, top 6 genes)\n")
cat("  b: Cortical layer spatial annotation (8 layers)\n")
cat("  c: bboxQuery 400x400 um ROI\n")
cat("  d: Gene enrichment dot plot (top 12 genes x 6 layers)\n")
cat("  e: Cell x gene heatmap (160 cells x 25 genes, Z-score)\n")
cat("  f: composeTransforms on real MERFISH coordinates\n")
