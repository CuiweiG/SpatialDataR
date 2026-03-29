#!/usr/bin/env Rscript
## Publication-grade README figure for SpatialDataR
## 4-panel: store overview | bbox query | aggregation | transforms
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

## Wong colorblind-safe palette (Nature standard)
pal4 <- c(Epithelial = "#0072B2", Stromal = "#D55E00",
    Immune = "#009E73", Endothelial = "#E69F00")

## ============================================================
## Panel a: Full spatial map (transcripts + cell boundaries)
## Shows: readSpatialData + spatialPoints + shapes
## ============================================================
pa <- ggplot() +
    geom_point(data = shp_ann,
        aes(x = x, y = y, color = cell_type,
            size = radius),
        alpha = 0.25, show.legend = FALSE) +
    geom_point(data = pts,
        aes(x = x, y = y),
        size = 0.15, alpha = 0.5, color = "grey30") +
    scale_color_manual(values = pal4) +
    scale_size_continuous(range = c(3, 12)) +
    coord_equal(xlim = c(0, 4.5), ylim = c(0, 4.5)) +
    annotate("text", x = 0.15, y = 4.35, label = "a",
        fontface = "bold", size = 4, hjust = 0) +
    annotate("text", x = 2.25, y = 4.35,
        label = "readSpatialData()",
        fontface = "italic", size = 2.5,
        color = "grey40") +
    labs(x = expression("x ("*mu*"m)"),
        y = expression("y ("*mu*"m)")) +
    theme_classic(base_size = 8) +
    theme(plot.margin = margin(5, 5, 5, 5))

## ============================================================
## Panel b: Bounding box query
## Shows: bboxQuery() — highlight queried region
## ============================================================
qx <- c(1, 3)
qy <- c(1, 3)
sub_pts <- pts[pts$x >= qx[1] & pts$x <= qx[2] &
    pts$y >= qy[1] & pts$y <= qy[2], ]

pb <- ggplot() +
    geom_point(data = pts,
        aes(x = x, y = y), size = 0.15,
        alpha = 0.2, color = "grey70") +
    annotate("rect", xmin = qx[1], xmax = qx[2],
        ymin = qy[1], ymax = qy[2],
        fill = NA, color = "#CC79A7",
        linewidth = 0.7, linetype = "dashed") +
    geom_point(data = sub_pts,
        aes(x = x, y = y, color = gene),
        size = 0.6, alpha = 0.8) +
    coord_equal(xlim = c(0, 4.5), ylim = c(0, 4.5)) +
    annotate("text", x = 0.15, y = 4.35, label = "b",
        fontface = "bold", size = 4, hjust = 0) +
    annotate("text", x = 2.25, y = 4.35,
        label = paste0("bboxQuery()
", nrow(sub_pts), "/", nrow(pts), " points"),
        fontface = "italic", size = 2.5,
        color = "grey40") +
    labs(x = expression("x ("*mu*"m)"),
        y = expression("y ("*mu*"m)"),
        color = NULL) +
    theme_classic(base_size = 8) +
    theme(legend.key.size = unit(5, "pt"),
        legend.text = element_text(size = 6),
        legend.position = c(0.88, 0.3),
        legend.background = element_blank(),
        plot.margin = margin(5, 5, 5, 5))

## ============================================================
## Panel c: Region aggregation heatmap
## Shows: aggregatePoints() → cell-gene count matrix
## ============================================================
counts_df <- aggregatePoints(
    spatialPoints(sd)[["transcripts"]],
    shapes(sd)[["cell_boundaries"]])
counts_r <- as.data.frame(counts_df)
gene_cols <- setdiff(colnames(counts_r), "cell_id")

## Merge cell type for row annotation
counts_ann <- merge(counts_r, obs[, c("cell_id", "cell_type")],
    by = "cell_id")
counts_ann <- counts_ann[order(counts_ann$cell_type), ]

## Melt for heatmap
hm_data <- do.call(rbind, lapply(seq_len(nrow(counts_ann)),
    function(i) {
    data.frame(
        cell = factor(i, levels = seq_len(nrow(counts_ann))),
        cell_type = counts_ann$cell_type[i],
        gene = gene_cols,
        count = as.numeric(counts_ann[i, gene_cols]),
        stringsAsFactors = FALSE)
}))

pc <- ggplot(hm_data, aes(x = gene, y = cell,
    fill = count)) +
    geom_tile(color = "white", linewidth = 0.1) +
    scale_fill_gradient(low = "white", high = "#0072B2",
        name = "Count") +
    facet_grid(cell_type ~ ., scales = "free_y",
        space = "free_y") +
    labs(x = NULL, y = "Cells",
        title = expression(bold("c") ~
            ~ italic("aggregatePoints()"))) +
    theme_minimal(base_size = 8) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1,
            size = 6),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(size = 6, angle = 0),
        legend.key.size = unit(6, "pt"),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 6),
        plot.title = element_text(size = 7,
            face = "italic", color = "grey40"),
        panel.grid = element_blank(),
        plot.margin = margin(5, 5, 5, 5))

## ============================================================
## Panel d: Coordinate transform visualization
## Shows: composeTransforms + transformCoords
## ============================================================
## Create a grid of points in pixel space
grid_pts <- expand.grid(
    x = seq(0, 20, by = 2),
    y = seq(0, 20, by = 2))

## Scale transform (pixels -> microns)
ct_scale <- CoordinateTransform("affine",
    affine = diag(c(0.2125, 0.2125, 1)),
    input_cs = "pixels", output_cs = "microns")

## Translate transform
ct_shift <- CoordinateTransform("affine",
    affine = matrix(c(1, 0, 1, 0, 1, 1, 0, 0, 1),
        nrow = 3, byrow = TRUE),
    input_cs = "microns", output_cs = "global")

## Composed
ct_composed <- composeTransforms(ct_scale, ct_shift)

grid_df_px <- DataFrame(
    x = grid_pts$x, y = grid_pts$y)
grid_df_um <- transformCoords(grid_df_px, ct_scale)
grid_df_gl <- transformCoords(grid_df_px, ct_composed)

all_grid <- rbind(
    data.frame(x = grid_df_px$x, y = grid_df_px$y,
        cs = "pixels (raw)"),
    data.frame(x = grid_df_um$x, y = grid_df_um$y,
        cs = "microns (scaled)"),
    data.frame(x = grid_df_gl$x, y = grid_df_gl$y,
        cs = "global (composed)")
)
all_grid$cs <- factor(all_grid$cs,
    levels = c("pixels (raw)", "microns (scaled)",
        "global (composed)"))

## Draw arrows between corresponding points
arrow_data <- data.frame(
    x = grid_df_px$x, y = grid_df_px$y,
    xend = grid_df_gl$x, yend = grid_df_gl$y)

pd <- ggplot() +
    geom_point(data = all_grid[all_grid$cs ==
        "pixels (raw)", ],
        aes(x = x, y = y),
        color = "grey70", size = 0.5, shape = 3) +
    geom_point(data = all_grid[all_grid$cs ==
        "global (composed)", ],
        aes(x = x, y = y),
        color = "#D55E00", size = 0.8) +
    geom_segment(data = arrow_data[seq(1, nrow(arrow_data),
        by = 3), ],
        aes(x = x, y = y, xend = xend, yend = yend),
        arrow = arrow(length = unit(1.5, "pt"),
            type = "closed"),
        color = "#999999", linewidth = 0.2, alpha = 0.5) +
    annotate("text", x = 0.5, y = 21, label = "d",
        fontface = "bold", size = 4, hjust = 0) +
    annotate("text", x = 10, y = 21,
        label = "composeTransforms()",
        fontface = "italic", size = 2.5,
        color = "grey40") +
    annotate("text", x = 16, y = 17,
        label = "pixels", color = "grey60",
        size = 2.5) +
    annotate("text", x = 3, y = 2.5,
        label = "global",
        color = "#D55E00", size = 2.5,
        fontface = "bold") +
    annotate("text", x = 10, y = 18.5,
        label = "scale \u00D7 0.2125\n+ translate (1, 1)",
        color = "grey50", size = 2, lineheight = 0.9) +
    coord_equal() +
    labs(x = "x", y = "y") +
    theme_classic(base_size = 8) +
    theme(plot.margin = margin(5, 5, 5, 5))

## ============================================================
## Compose final figure
## ============================================================
library(patchwork)
fig <- (pa | pb) / (pc | pd) +
    plot_layout(heights = c(1, 1))

od <- "man/figures"
if (!dir.exists(od)) dir.create(od, recursive = TRUE)

ggsave(file.path(od, "fig1_spatial_overview.png"),
    fig, width = 183, height = 160, units = "mm",
    dpi = 300, bg = "white")

cat("Saved:", file.path(od, "fig1_spatial_overview.png"), "\n")
cat("Size:", round(file.size(
    file.path(od, "fig1_spatial_overview.png")) / 1024),
    "KB\n")
