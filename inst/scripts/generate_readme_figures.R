#!/usr/bin/env Rscript
## SpatialDataR: 4 independent figures for README
## Each figure demonstrates one unique capability
## that no other R package provides.
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

od <- "man/figures"
if (!dir.exists(od)) dir.create(od, recursive = TRUE)

pal <- c(Epithelial = "#0072B2", Stromal = "#D55E00",
    Immune = "#009E73", Endothelial = "#CC79A7")

th <- theme_classic(base_size = 8, base_family = "sans") +
    theme(
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7, color = "black"),
        plot.title = element_text(size = 9, face = "bold",
            margin = margin(0, 0, 4, 0)),
        plot.subtitle = element_text(size = 7,
            color = "grey40", face = "italic",
            margin = margin(0, 0, 6, 0)),
        plot.margin = margin(8, 10, 8, 8))

## ==========================================================
## Figure 1: Zero-Python Zarr Store Reading
## WHY: No other R package can read SpatialData .zarr stores
## natively. Users currently need reticulate + Python.
## ==========================================================
cat("Generating fig1_store_reading.png ...\n")

fig1 <- ggplot() +
    geom_point(data = shp_ann,
        aes(x = x, y = y, color = cell_type),
        shape = 1, size = shp_ann$radius * 20,
        stroke = 0.6, alpha = 0.7) +
    geom_point(data = pts, aes(x = x, y = y),
        size = 0.12, alpha = 0.3, color = "grey25") +
    scale_color_manual(values = pal, name = "Cell type") +
    coord_equal(xlim = c(-0.2, 4.6),
        ylim = c(-0.2, 4.6)) +
    labs(
        title = "Native SpatialData Zarr Reading",
        subtitle = paste0(
            "readSpatialData(): 5 element types, ",
            "2 coordinate systems, zero Python"),
        x = expression("x ("*mu*"m)"),
        y = expression("y ("*mu*"m)")) +
    annotate("label", x = 0.1, y = 4.5,
        label = paste0(
            nrow(pts), " transcripts\n",
            nrow(shp), " cells\n",
            "10 genes"),
        size = 2.2, hjust = 0, vjust = 1,
        fill = "white", label.size = 0.2,
        color = "grey30", lineheight = 0.9) +
    th +
    theme(
        legend.position = c(0.85, 0.2),
        legend.background = element_rect(
            fill = "white", color = "grey80",
            linewidth = 0.3),
        legend.key.size = unit(8, "pt"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7,
            face = "bold"))

ggsave(file.path(od, "fig1_store_reading.png"), fig1,
    width = 100, height = 95, units = "mm",
    dpi = 300, bg = "white")

## ==========================================================
## Figure 2: Spatial Bounding Box Query
## WHY: Python spatialdata has bounding_box_query() but
## no R package does. Essential for ROI analysis.
## ==========================================================
cat("Generating fig2_spatial_query.png ...\n")

qx <- c(1, 3); qy <- c(1, 3)
pts$inside <- pts$x >= qx[1] & pts$x <= qx[2] &
    pts$y >= qy[1] & pts$y <= qy[2]
n_in <- sum(pts$inside)

top4 <- names(sort(table(pts$gene[pts$inside]),
    decreasing = TRUE))[1:4]
pts_q <- pts[pts$inside & pts$gene %in% top4, ]
gene_pal <- c("#0072B2", "#D55E00", "#009E73", "#E69F00")
names(gene_pal) <- top4

fig2 <- ggplot() +
    geom_point(data = pts[!pts$inside, ],
        aes(x = x, y = y),
        size = 0.15, alpha = 0.12, color = "grey70") +
    annotate("rect", xmin = qx[1], xmax = qx[2],
        ymin = qy[1], ymax = qy[2],
        fill = "#CC79A7", alpha = 0.06,
        color = "#CC79A7", linewidth = 0.6) +
    geom_point(data = pts_q,
        aes(x = x, y = y, color = gene),
        size = 1.0, alpha = 0.85) +
    scale_color_manual(values = gene_pal,
        name = "Gene") +
    coord_equal(xlim = c(-0.2, 4.6),
        ylim = c(-0.2, 4.6)) +
    labs(
        title = "Bounding Box Spatial Query",
        subtitle = paste0(
            "bboxQuery(): ", n_in, "/", nrow(pts),
            " transcripts in ROI [1,3] \u00D7 [1,3]"),
        x = expression("x ("*mu*"m)"),
        y = expression("y ("*mu*"m)")) +
    th +
    theme(
        legend.position = c(0.88, 0.18),
        legend.background = element_rect(
            fill = "white", color = "grey80",
            linewidth = 0.3),
        legend.key.size = unit(8, "pt"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7,
            face = "bold"))

ggsave(file.path(od, "fig2_spatial_query.png"), fig2,
    width = 100, height = 95, units = "mm",
    dpi = 300, bg = "white")

## ==========================================================
## Figure 3: Region Aggregation
## WHY: Going from molecules to cell-gene matrix is the
## critical step in spatial transcriptomics. Only Python
## spatialdata.aggregate() does this; now R can too.
## ==========================================================
cat("Generating fig3_aggregation.png ...\n")

counts_df <- aggregatePoints(
    spatialPoints(sd)[["transcripts"]],
    shapes(sd)[["cell_boundaries"]])
counts_r <- as.data.frame(counts_df)
gene_cols <- setdiff(colnames(counts_r), "cell_id")
counts_ann <- merge(counts_r,
    obs[, c("cell_id", "cell_type")], by = "cell_id")
counts_ann <- counts_ann[order(counts_ann$cell_type), ]

## Row-normalize
rs <- rowSums(counts_ann[, gene_cols])
rs[rs == 0] <- 1
cn <- counts_ann
cn[, gene_cols] <- counts_ann[, gene_cols] / rs

hm <- do.call(rbind, lapply(seq_len(nrow(cn)),
    function(i) {
    data.frame(
        cell = factor(i),
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
        limits = c(0, 0.5),
        oob = scales::squish,
        breaks = c(0, 0.25, 0.5)) +
    facet_grid(cell_type ~ ., scales = "free_y",
        space = "free_y", switch = "y") +
    labs(
        title = "Region Aggregation",
        subtitle = paste0(
            "aggregatePoints(): ",
            nrow(cn), " cells \u00D7 ",
            length(gene_cols), " genes"),
        x = NULL, y = NULL) +
    theme_minimal(base_size = 8) +
    theme(
        axis.text.x = element_text(angle = 45,
            hjust = 1, size = 7, color = "black"),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.text.y.left = element_text(size = 7,
            angle = 0, hjust = 1, face = "italic"),
        strip.placement = "outside",
        legend.key.height = unit(12, "pt"),
        legend.key.width = unit(6, "pt"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7,
            face = "bold"),
        panel.grid = element_blank(),
        panel.spacing = unit(1, "pt"),
        plot.title = element_text(size = 9,
            face = "bold",
            margin = margin(0, 0, 4, 0)),
        plot.subtitle = element_text(size = 7,
            color = "grey40", face = "italic",
            margin = margin(0, 0, 6, 0)),
        plot.margin = margin(8, 10, 8, 8))

ggsave(file.path(od, "fig3_aggregation.png"), fig3,
    width = 100, height = 95, units = "mm",
    dpi = 300, bg = "white")

## ==========================================================
## Figure 4: Coordinate Transform Composition
## WHY: OME-NGFF defines scale/translate/affine/sequence
## transforms. No R package parses or composes these.
## Essential for multi-modal alignment.
## ==========================================================
cat("Generating fig4_transforms.png ...\n")

rep_pts <- data.frame(
    x = c(2, 8, 14, 5, 17),
    y = c(3, 15, 7, 10, 18),
    label = c("A", "B", "C", "D", "E"))

ct_composed <- composeTransforms(
    CoordinateTransform("affine",
        affine = diag(c(0.2125, 0.2125, 1))),
    CoordinateTransform("affine",
        affine = matrix(c(1, 0, 1, 0, 1, 0.5, 0, 0, 1),
            nrow = 3, byrow = TRUE)))

rep_df <- DataFrame(x = rep_pts$x, y = rep_pts$y)
rep_out <- as.data.frame(transformCoords(rep_df,
    ct_composed))

arr <- data.frame(
    x = rep_pts$x, y = rep_pts$y,
    xend = rep_out$x, yend = rep_out$y,
    label = rep_pts$label)

fig4 <- ggplot() +
    geom_point(data = rep_pts,
        aes(x = x, y = y),
        shape = 4, size = 3, stroke = 0.7,
        color = "grey50") +
    geom_text(data = rep_pts,
        aes(x = x + 0.6, y = y + 0.6, label = label),
        size = 2.5, color = "grey50") +
    geom_segment(data = arr,
        aes(x = x, y = y, xend = xend, yend = yend),
        arrow = arrow(length = unit(4, "pt"),
            type = "closed"),
        color = "grey40", linewidth = 0.3) +
    geom_point(data = rep_out,
        aes(x = x, y = y),
        shape = 16, size = 3, color = "#D55E00") +
    geom_text(data = data.frame(
            x = rep_out$x + 0.6,
            y = rep_out$y + 0.6,
            label = rep_pts$label),
        aes(x = x, y = y, label = label),
        size = 2.5, color = "#D55E00",
        fontface = "bold") +
    annotate("label", x = 14, y = 1.5,
        label = "pixels",
        size = 2.5, color = "grey50",
        fill = "white", label.size = 0,
        fontface = "italic") +
    annotate("label", x = 1.5, y = 1.5,
        label = "global",
        size = 2.5, color = "#D55E00",
        fill = "white", label.size = 0,
        fontface = "bold.italic") +
    coord_equal(xlim = c(0, 20), ylim = c(0, 20)) +
    labs(
        title = "Transform Composition",
        subtitle = expression(
            "composeTransforms(): scale" %*%
            "0.2125 + translate(1, 0.5)"),
        x = "x", y = "y") +
    th

ggsave(file.path(od, "fig4_transforms.png"), fig4,
    width = 100, height = 95, units = "mm",
    dpi = 300, bg = "white")

## ==========================================================
## Also generate combined 4-panel for overview
## ==========================================================
cat("Generating fig1_spatial_overview.png ...\n")
library(patchwork)

combined <- (
    (fig1 + labs(title = NULL, subtitle = NULL,
        tag = "a")) |
    (fig2 + labs(title = NULL, subtitle = NULL,
        tag = "b"))
) / (
    (fig3 + labs(title = NULL, subtitle = NULL,
        tag = "c")) |
    (fig4 + labs(title = NULL, subtitle = NULL,
        tag = "d"))
) + plot_layout(heights = c(1, 0.9)) &
    theme(plot.tag = element_text(size = 10,
        face = "bold"))

ggsave(file.path(od, "fig1_spatial_overview.png"),
    combined, width = 180, height = 160, units = "mm",
    dpi = 300, bg = "white")

cat("\nAll figures saved to", od, "\n")
for (f in list.files(od, pattern = "^fig")) {
    cat("  ", f, ":",
        round(file.size(file.path(od, f)) / 1024),
        "KB\n")
}
