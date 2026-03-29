#!/usr/bin/env Rscript
.libPaths("C:/Users/win10/R/win-library/4.4")
library(SpatialDataR)
library(S4Vectors)
library(ggplot2)
library(patchwork)

store <- system.file("extdata", "xenium_mini.zarr",
    package = "SpatialDataR")
sd <- readSpatialData(store)
pts <- as.data.frame(spatialPoints(sd)[["transcripts"]])
shp <- as.data.frame(shapes(sd)[["cell_boundaries"]])
obs <- as.data.frame(tables(sd)[["table"]]$obs)
shp_ann <- merge(shp, obs, by = "cell_id")

od <- "man/figures"
if (!dir.exists(od)) dir.create(od, recursive = TRUE)

cell_pal <- c(Epithelial = "#0072B2",
    Stromal = "#D55E00", Immune = "#F0E442",
    Endothelial = "#000000")
top5 <- names(sort(table(pts$gene), TRUE))[1:5]
gene_pal <- c("#0072B2", "#D55E00", "#009E73",
    "#CC79A7", "#E69F00")
names(gene_pal) <- top5
gene_shp <- c(16, 17, 15, 18, 3)
names(gene_shp) <- top5

th <- theme_classic(base_size = 8) +
    theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7,
            color = "black"),
        plot.margin = margin(8, 10, 8, 8))

pts$gene5 <- ifelse(pts$gene %in% top5, pts$gene, NA)
pts_top <- pts[!is.na(pts$gene5), ]

## ==========================================================
## Fig 1: 3-panel — cells | transcripts | combined
## ==========================================================
cat("fig1 ...\n")

p1a <- ggplot(shp_ann, aes(x = x, y = y,
    color = cell_type)) +
    geom_point(shape = 1, size = 6, stroke = 0.8) +
    scale_color_manual(values = cell_pal, name = NULL) +
    coord_equal(xlim = c(-0.1, 4.5),
        ylim = c(-0.1, 4.5)) +
    labs(subtitle = "Shapes: 50 cells",
        x = expression(mu*"m"),
        y = expression(mu*"m")) +
    th + theme(legend.position = "bottom",
        legend.text = element_text(size = 7),
        plot.subtitle = element_text(size = 8,
            face = "bold"))

p1b <- ggplot(pts_top, aes(x = x, y = y,
    color = gene5, shape = gene5)) +
    geom_point(size = 1.8, alpha = 0.8) +
    scale_color_manual(values = gene_pal, name = NULL) +
    scale_shape_manual(values = gene_shp, name = NULL) +
    coord_equal(xlim = c(-0.1, 4.5),
        ylim = c(-0.1, 4.5)) +
    labs(subtitle = "Points: 500 transcripts",
        x = expression(mu*"m"), y = NULL) +
    th + theme(legend.position = "bottom",
        legend.text = element_text(size = 7),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.subtitle = element_text(size = 8,
            face = "bold"))

p1c <- ggplot() +
    geom_point(data = shp_ann,
        aes(x = x, y = y, fill = cell_type),
        shape = 21, size = 5, stroke = 0.3,
        alpha = 0.2, color = "grey50") +
    geom_point(data = pts_top,
        aes(x = x, y = y, color = gene5,
            shape = gene5),
        size = 1.5, alpha = 0.75) +
    scale_fill_manual(values = cell_pal,
        guide = "none") +
    scale_color_manual(values = gene_pal,
        guide = "none") +
    scale_shape_manual(values = gene_shp,
        guide = "none") +
    coord_equal(xlim = c(-0.1, 4.5),
        ylim = c(-0.1, 4.5)) +
    labs(subtitle = "Combined overlay",
        x = expression(mu*"m"), y = NULL) +
    th + theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.subtitle = element_text(size = 8,
            face = "bold"))

fig1 <- p1a + p1b + p1c +
    plot_layout(ncol = 3) +
    plot_annotation(
        title = "Native SpatialData Zarr Reading",
        subtitle = paste0(
            "readSpatialData(): one call reads ",
            "shapes, points, tables, images, labels"),
        theme = theme(
            plot.title = element_text(size = 11,
                face = "bold"),
            plot.subtitle = element_text(size = 8,
                color = "grey35",
                face = "italic")))

ggsave(file.path(od, "fig1_store_reading.png"), fig1,
    width = 220, height = 95, units = "mm",
    dpi = 300, bg = "white")

## ==========================================================
## Fig 2: Bbox — left=full, right=queried
## ==========================================================
cat("fig2 ...\n")

qx <- c(1, 3); qy <- c(1, 3)
pts$inside <- pts$x >= qx[1] & pts$x <= qx[2] &
    pts$y >= qy[1] & pts$y <= qy[2]
n_in <- sum(pts$inside)

## Left: all points with box outline
p2a <- ggplot(pts, aes(x = x, y = y)) +
    geom_point(size = 0.6, alpha = 0.4,
        color = "grey40") +
    annotate("rect", xmin = qx[1], xmax = qx[2],
        ymin = qy[1], ymax = qy[2],
        fill = NA, color = "#D55E00",
        linewidth = 0.8, linetype = "dashed") +
    coord_equal(xlim = c(-0.1, 4.5),
        ylim = c(-0.1, 4.5)) +
    labs(subtitle = paste0("All ", nrow(pts),
        " transcripts"),
        x = expression(mu*"m"),
        y = expression(mu*"m")) +
    th + theme(plot.subtitle = element_text(size = 8,
        face = "bold"))

## Right: only inside points, colored by gene
pts_sel <- pts[pts$inside & !is.na(pts$gene5), ]

p2b <- ggplot(pts_sel, aes(x = x, y = y,
    color = gene5, shape = gene5)) +
    geom_point(size = 2.5, alpha = 0.9) +
    scale_color_manual(values = gene_pal,
        name = NULL) +
    scale_shape_manual(values = gene_shp,
        name = NULL) +
    coord_equal(xlim = c(qx[1] - 0.1, qx[2] + 0.1),
        ylim = c(qy[1] - 0.1, qy[2] + 0.1)) +
    labs(subtitle = paste0(n_in,
        " selected in ROI"),
        x = expression(mu*"m"), y = NULL) +
    th + theme(legend.position = "right",
        legend.text = element_text(size = 7),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.subtitle = element_text(size = 8,
            face = "bold", color = "#D55E00"))

fig2 <- p2a + p2b +
    plot_layout(ncol = 2, widths = c(1, 1)) +
    plot_annotation(
        title = "Bounding Box Spatial Query",
        subtitle = paste0(
            "bboxQuery(): ", n_in, "/", nrow(pts),
            " transcripts in [1,3] x [1,3]"),
        theme = theme(
            plot.title = element_text(size = 11,
                face = "bold"),
            plot.subtitle = element_text(size = 8,
                color = "grey35",
                face = "italic")))

ggsave(file.path(od, "fig2_spatial_query.png"), fig2,
    width = 200, height = 95, units = "mm",
    dpi = 300, bg = "white")

## ==========================================================
## Fig 3: Aggregation heatmap
## ==========================================================
cat("fig3 ...\n")

counts_df <- aggregatePoints(
    spatialPoints(sd)[["transcripts"]],
    shapes(sd)[["cell_boundaries"]])
counts_r <- as.data.frame(counts_df)
gene_cols <- setdiff(colnames(counts_r), "cell_id")
counts_ann <- merge(counts_r,
    obs[, c("cell_id", "cell_type")], by = "cell_id")
counts_ann <- counts_ann[order(counts_ann$cell_type), ]
rs <- rowSums(counts_ann[, gene_cols])
rs[rs == 0] <- 1
cn <- counts_ann
cn[, gene_cols] <- counts_ann[, gene_cols] / rs

hm <- do.call(rbind, lapply(seq_len(nrow(cn)),
    function(i) {
    data.frame(cell = factor(i),
        cell_type = cn$cell_type[i],
        gene = factor(gene_cols, levels = gene_cols),
        frac = as.numeric(cn[i, gene_cols]))
}))
hm$cell_type <- factor(hm$cell_type,
    levels = c("Epithelial", "Stromal",
        "Immune", "Endothelial"))

fig3 <- ggplot(hm, aes(x = gene, y = cell,
    fill = frac)) +
    geom_tile(color = NA) +
    scale_fill_gradient(low = "grey98",
        high = "#0072B2", name = "Fraction",
        limits = c(0, 0.5), oob = scales::squish,
        breaks = c(0, 0.25, 0.5)) +
    facet_grid(cell_type ~ ., scales = "free_y",
        space = "free_y", switch = "y") +
    labs(title = "Region Aggregation",
        subtitle = paste0("aggregatePoints(): ",
            nrow(cn), " cells x ",
            length(gene_cols), " genes"),
        x = NULL, y = NULL) +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45,
            hjust = 1, size = 7, color = "black"),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.text.y.left = element_text(size = 7,
            angle = 0, hjust = 1, face = "italic"),
        strip.placement = "outside",
        legend.key.height = unit(14, "pt"),
        legend.key.width = unit(6, "pt"),
        legend.text = element_text(size = 6.5),
        legend.title = element_text(size = 7,
            face = "bold"),
        panel.grid = element_blank(),
        panel.spacing = unit(1, "pt"),
        plot.title = element_text(size = 11,
            face = "bold"),
        plot.subtitle = element_text(size = 8,
            color = "grey35", face = "italic"),
        plot.margin = margin(8, 10, 8, 8))

ggsave(file.path(od, "fig3_aggregation.png"), fig3,
    width = 160, height = 120, units = "mm",
    dpi = 300, bg = "white")

## ==========================================================
## Fig 4: Transform composition
## ==========================================================
cat("fig4 ...\n")

rep_pts <- data.frame(
    x = c(2, 8, 14, 5, 17),
    y = c(3, 15, 7, 10, 18),
    id = LETTERS[1:5])

ct_comp <- composeTransforms(
    CoordinateTransform("affine",
        affine = diag(c(0.2125, 0.2125, 1))),
    CoordinateTransform("affine",
        affine = matrix(c(1, 0, 1, 0, 1, 0.5, 0, 0, 1),
            nrow = 3, byrow = TRUE)))
rep_out <- as.data.frame(transformCoords(
    DataFrame(x = rep_pts$x, y = rep_pts$y), ct_comp))

arr <- data.frame(x = rep_pts$x, y = rep_pts$y,
    xend = rep_out$x, yend = rep_out$y,
    id = rep_pts$id)

fig4 <- ggplot() +
    geom_point(data = rep_pts, aes(x = x, y = y),
        shape = 4, size = 4, stroke = 0.9,
        color = "grey50") +
    geom_text(data = rep_pts,
        aes(x = x - 0.9, y = y, label = id),
        size = 3, color = "grey50", fontface = "bold") +
    geom_segment(data = arr,
        aes(x = x, y = y, xend = xend, yend = yend),
        arrow = arrow(length = unit(5, "pt"),
            type = "closed"),
        color = "grey40", linewidth = 0.4) +
    geom_point(data = rep_out, aes(x = x, y = y),
        shape = 16, size = 4, color = "#D55E00") +
    geom_text(data = data.frame(
        x = rep_out$x + 0.9, y = rep_out$y,
        id = rep_pts$id),
        aes(x = x, y = y, label = id),
        size = 3, color = "#D55E00",
        fontface = "bold") +
    annotate("text", x = 16, y = 19, label = "pixels",
        size = 3.5, color = "grey50",
        fontface = "italic") +
    annotate("text", x = 3, y = 0.8, label = "global",
        size = 3.5, color = "#D55E00",
        fontface = "bold.italic") +
    coord_equal(xlim = c(-1, 20), ylim = c(-1, 20)) +
    labs(title = "Transform Composition",
        subtitle = "composeTransforms(): scale(0.2125) + translate(1, 0.5)",
        x = "x", y = "y") +
    th + theme(
        plot.title = element_text(size = 11,
            face = "bold"),
        plot.subtitle = element_text(size = 8,
            color = "grey35", face = "italic"))

ggsave(file.path(od, "fig4_transforms.png"), fig4,
    width = 140, height = 120, units = "mm",
    dpi = 300, bg = "white")

## Combined
cat("combined ...\n")
combined <- (
    (p1a + labs(subtitle = NULL, tag = "a") +
        theme(legend.position = "none")) |
    (p1b + labs(subtitle = NULL, tag = "b") +
        theme(legend.position = "none")) |
    (p1c + labs(subtitle = NULL, tag = "c"))
) / (
    (p2a + labs(subtitle = NULL, tag = "d")) |
    (p2b + labs(subtitle = NULL, tag = "e") +
        theme(legend.position = "none"))
) + plot_layout(heights = c(1, 1)) &
    theme(plot.tag = element_text(size = 10,
        face = "bold"))

ggsave(file.path(od, "fig1_spatial_overview.png"),
    combined, width = 220, height = 160, units = "mm",
    dpi = 300, bg = "white")

cat("\nDone:\n")
for (f in list.files(od, "^fig")) {
    cat("  ", f, ":",
        round(file.size(file.path(od, f)) / 1024),
        "KB\n")
}
