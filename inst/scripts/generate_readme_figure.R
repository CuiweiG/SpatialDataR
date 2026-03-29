#!/usr/bin/env Rscript
## Publication figure for SpatialDataR README
## Standard: Nature Methods / Genome Biology
.libPaths("C:/Users/win10/R/win-library/4.4")
library(SpatialDataR)
library(S4Vectors)
library(ggplot2)

store <- system.file("extdata", "xenium_mini.zarr",
    package = "SpatialDataR")
sd <- readSpatialData(store)
pts <- as.data.frame(spatialPoints(sd)[["transcripts"]])
shp <- as.data.frame(shapes(sd)[["cell_boundaries"]])
obs <- as.data.frame(tables(sd)[["table"]]$obs)
shp_ann <- merge(shp, obs, by = "cell_id")

## Okabe-Ito colorblind-safe (4 types)
pal <- c(Epithelial = "#0072B2", Stromal = "#D55E00",
    Immune = "#009E73", Endothelial = "#CC79A7")

## Shared theme
th <- theme_classic(base_size = 7, base_family = "sans") +
    theme(
        plot.margin = margin(4, 6, 4, 4),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6, color = "black"))

## ==========================================================
## Panel a: Spatial map — cell outlines + transcripts
## ==========================================================
pa <- ggplot() +
    ## Cell outlines (NOT filled circles)
    geom_point(data = shp_ann,
        aes(x = x, y = y, color = cell_type),
        shape = 1, size = shp_ann$radius * 18,
        stroke = 0.5, alpha = 0.7) +
    ## Transcripts as tiny dots
    geom_point(data = pts, aes(x = x, y = y),
        size = 0.08, alpha = 0.35, color = "grey25") +
    scale_color_manual(values = pal, name = NULL) +
    coord_equal(xlim = c(-0.1, 4.5),
        ylim = c(-0.1, 4.5)) +
    labs(x = expression("x ("*mu*"m)"),
        y = expression("y ("*mu*"m)")) +
    th +
    theme(
        legend.position = c(0.82, 0.18),
        legend.background = element_blank(),
        legend.key.size = unit(7, "pt"),
        legend.text = element_text(size = 6),
        legend.spacing.y = unit(1, "pt"))

## ==========================================================
## Panel b: Bounding box query — before/after
## ==========================================================
qx <- c(1, 3); qy <- c(1, 3)
pts$inside <- pts$x >= qx[1] & pts$x <= qx[2] &
    pts$y >= qy[1] & pts$y <= qy[2]
n_in <- sum(pts$inside)

## Show only top 4 genes for clarity
top4 <- names(sort(table(pts$gene[pts$inside]),
    decreasing = TRUE))[1:4]
pts_sub <- pts[pts$inside & pts$gene %in% top4, ]

gene_pal <- c("#0072B2", "#D55E00", "#009E73", "#E69F00")
names(gene_pal) <- top4

pb <- ggplot() +
    geom_point(data = pts[!pts$inside, ],
        aes(x = x, y = y),
        size = 0.1, alpha = 0.15, color = "grey75") +
    annotate("rect", xmin = qx[1], xmax = qx[2],
        ymin = qy[1], ymax = qy[2],
        fill = NA, color = "#CC79A7",
        linewidth = 0.6, linetype = "solid") +
    geom_point(data = pts_sub,
        aes(x = x, y = y, color = gene),
        size = 0.7, alpha = 0.85, shape = 16) +
    scale_color_manual(values = gene_pal, name = NULL) +
    coord_equal(xlim = c(-0.1, 4.5),
        ylim = c(-0.1, 4.5)) +
    annotate("text", x = 3.05, y = 3.15,
        label = paste0(n_in, "/", nrow(pts)),
        size = 2.2, hjust = 0, color = "#CC79A7") +
    labs(x = expression("x ("*mu*"m)"),
        y = expression("y ("*mu*"m)")) +
    th +
    theme(
        legend.position = c(0.88, 0.18),
        legend.background = element_blank(),
        legend.key.size = unit(6, "pt"),
        legend.text = element_text(size = 6))

## ==========================================================
## Panel c: Aggregation heatmap — wider, shorter
## ==========================================================
counts_df <- aggregatePoints(
    spatialPoints(sd)[["transcripts"]],
    shapes(sd)[["cell_boundaries"]])
counts_r <- as.data.frame(counts_df)
gene_cols <- setdiff(colnames(counts_r), "cell_id")
counts_ann <- merge(counts_r,
    obs[, c("cell_id", "cell_type")], by = "cell_id")
counts_ann <- counts_ann[order(counts_ann$cell_type), ]

## Normalize per cell (fraction)
row_sums <- rowSums(counts_ann[, gene_cols])
row_sums[row_sums == 0] <- 1
counts_norm <- counts_ann
counts_norm[, gene_cols] <- counts_ann[, gene_cols] /
    row_sums

## Melt
hm <- do.call(rbind, lapply(seq_len(nrow(counts_norm)),
    function(i) {
    data.frame(
        cell = factor(i),
        cell_type = counts_norm$cell_type[i],
        gene = factor(gene_cols, levels = gene_cols),
        frac = as.numeric(counts_norm[i, gene_cols]))
}))

## Reorder cell_type factor
hm$cell_type <- factor(hm$cell_type,
    levels = c("Epithelial", "Stromal",
        "Immune", "Endothelial"))

pc <- ggplot(hm, aes(x = gene, y = cell, fill = frac)) +
    geom_tile(color = NA) +
    scale_fill_gradient(low = "grey97", high = "#0072B2",
        name = "Fraction", limits = c(0, 0.5),
        oob = scales::squish,
        breaks = c(0, 0.25, 0.5)) +
    facet_grid(cell_type ~ ., scales = "free_y",
        space = "free_y", switch = "y") +
    labs(x = NULL, y = NULL) +
    scale_x_discrete(position = "bottom") +
    theme_minimal(base_size = 7) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1,
            size = 6, color = "black"),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.text.y.left = element_text(size = 6,
            angle = 0, hjust = 1, face = "italic"),
        strip.placement = "outside",
        legend.key.height = unit(10, "pt"),
        legend.key.width = unit(5, "pt"),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 6),
        panel.grid = element_blank(),
        panel.spacing = unit(1, "pt"),
        plot.margin = margin(4, 6, 4, 4))

## ==========================================================
## Panel d: Transform before/after (side by side effect)
## ==========================================================
## 5 representative points
rep_pts <- data.frame(
    x = c(2, 8, 14, 5, 11),
    y = c(3, 15, 7, 10, 18),
    id = factor(1:5))

ct_composed <- composeTransforms(
    CoordinateTransform("affine",
        affine = diag(c(0.2125, 0.2125, 1)),
        input_cs = "pixels", output_cs = "um"),
    CoordinateTransform("affine",
        affine = matrix(c(1, 0, 1, 0, 1, 0.5, 0, 0, 1),
            nrow = 3, byrow = TRUE),
        input_cs = "um", output_cs = "global"))

rep_df <- DataFrame(x = rep_pts$x, y = rep_pts$y)
rep_out <- as.data.frame(
    transformCoords(rep_df, ct_composed))

arrow_df <- data.frame(
    x = rep_pts$x, y = rep_pts$y,
    xend = rep_out$x, yend = rep_out$y,
    id = factor(1:5))

pd <- ggplot() +
    ## Before (pixel grid background)
    geom_point(data = rep_pts,
        aes(x = x, y = y),
        shape = 4, size = 2.5, stroke = 0.6,
        color = "grey55") +
    ## Arrows
    geom_segment(data = arrow_df,
        aes(x = x, y = y, xend = xend, yend = yend),
        arrow = arrow(length = unit(3, "pt"),
            type = "closed"),
        color = "grey45", linewidth = 0.35) +
    ## After (global)
    geom_point(data = rep_out,
        aes(x = x, y = y),
        shape = 16, size = 2.5,
        color = "#D55E00") +
    ## Labels
    annotate("label", x = 14, y = 2,
        label = "pixels",
        size = 2.3, color = "grey50",
        fill = "white", label.size = 0,
        fontface = "italic") +
    annotate("label", x = 2.5, y = 1,
        label = "global",
        size = 2.3, color = "#D55E00",
        fill = "white", label.size = 0,
        fontface = "bold.italic") +
    annotate("text", x = 10, y = 19.5,
        label = expression(
            "scale" %*% "0.2125 + translate"),
        size = 2, color = "grey45") +
    coord_equal(xlim = c(0, 20), ylim = c(0, 20)) +
    labs(x = "x", y = "y") +
    th

## ==========================================================
## Compose with patchwork
## ==========================================================
library(patchwork)

fig <- (
    (pa + labs(tag = "a")) |
    (pb + labs(tag = "b"))
) / (
    (pc + labs(tag = "c")) |
    (pd + labs(tag = "d"))
) + plot_layout(heights = c(1, 0.9)) &
    theme(plot.tag = element_text(size = 9,
        face = "bold"))

od <- "man/figures"
if (!dir.exists(od)) dir.create(od, recursive = TRUE)

ggsave(file.path(od, "fig1_spatial_overview.png"),
    fig, width = 180, height = 155, units = "mm",
    dpi = 300, bg = "white")

cat("Saved:", file.path(od, "fig1_spatial_overview.png"),
    "\n")
cat("Size:", round(file.size(file.path(od,
    "fig1_spatial_overview.png")) / 1024), "KB\n")
