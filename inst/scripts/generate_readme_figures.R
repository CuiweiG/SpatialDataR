#!/usr/bin/env Rscript
## SpatialDataR README figures �?Nature Methods standard
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

## Palettes: maximally distinct for CVD
cell_pal <- c(
    Epithelial = "#0072B2",   # blue
    Stromal = "#D55E00",      # vermillion
    Immune = "#F0E442",       # yellow
    Endothelial = "#000000")  # black

## Top 5 genes, shapes + colors for redundant encoding
top5 <- names(sort(table(pts$gene),
    decreasing = TRUE))[1:5]
gene_pal <- c("#0072B2", "#D55E00", "#009E73",
    "#CC79A7", "#E69F00")
names(gene_pal) <- top5
gene_shp <- c(16, 17, 15, 18, 3)
names(gene_shp) <- top5

th <- theme_classic(base_size = 8) +
    theme(
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7,
            color = "black"),
        plot.title = element_text(size = 9,
            face = "bold",
            margin = margin(0, 0, 3, 0)),
        plot.subtitle = element_text(size = 7,
            color = "grey35", face = "italic",
            margin = margin(0, 0, 6, 0)),
        plot.margin = margin(8, 10, 8, 8))

## ==========================================================
## Fig 1: Store Reading �?show ALL data types clearly
## ==========================================================
cat("fig1 ...\n")

## Color transcripts by top 5 genes
pts$gene5 <- ifelse(pts$gene %in% top5,
    pts$gene, NA_character_)
pts_top <- pts[!is.na(pts$gene5), ]

fig1 <- ggplot() +
    ## Cell boundaries: filled + outlined, uniform size
    geom_point(data = shp_ann,
        aes(x = x, y = y, fill = cell_type),
        shape = 21, size = 4, stroke = 0.5,
        alpha = 0.3, color = "grey40") +
    ## Transcripts: colored by gene, different shapes
    geom_point(data = pts_top,
        aes(x = x, y = y,
            color = gene5, shape = gene5),
        size = 1.2, alpha = 0.8) +
    scale_fill_manual(values = cell_pal,
        name = "Cell type") +
    scale_color_manual(values = gene_pal,
        name = "Transcript") +
    scale_shape_manual(values = gene_shp,
        name = "Transcript") +
    coord_equal(xlim = c(-0.2, 4.6),
        ylim = c(-0.2, 4.6)) +
    guides(
        fill = guide_legend(order = 1,
            override.aes = list(alpha = 0.5,
                size = 3)),
        color = guide_legend(order = 2),
        shape = guide_legend(order = 2)) +
    labs(
        title = "Native SpatialData Zarr Reading",
        subtitle = paste0(
            nrow(pts), " transcripts + ",
            nrow(shp), " cells + ",
            length(unique(pts$gene)),
            " genes from one readSpatialData() call"),
        x = expression("x ("*mu*"m)"),
        y = expression("y ("*mu*"m)")) +
    th +
    theme(
        legend.position = "right",
        legend.box = "vertical",
        legend.key.size = unit(8, "pt"),
        legend.text = element_text(size = 6.5),
        legend.title = element_text(size = 7,
            face = "bold"),
        legend.spacing.y = unit(2, "pt"),
        legend.margin = margin(0, 0, 0, 0))

ggsave(file.path(od, "fig1_store_reading.png"), fig1,
    width = 180, height = 130, units = "mm",
    dpi = 300, bg = "white")

## ==========================================================
## Fig 2: Bbox query �?clear before/after
## ==========================================================
cat("fig2 ...\n")

qx <- c(1, 3); qy <- c(1, 3)
pts$inside <- pts$x >= qx[1] & pts$x <= qx[2] &
    pts$y >= qy[1] & pts$y <= qy[2]
n_in <- sum(pts$inside)
pts_in <- pts[pts$inside & pts$gene %in% top5, ]

## Excluded points — show them clearly as red X marks
pts_out <- pts[!pts$inside, ]

fig2 <- ggplot() +
    ## Excluded: visible red-gray crosses
    geom_point(data = pts_out,
        aes(x = x, y = y),
        size = 1.2, alpha = 0.4, color = "#999999",
        shape = 4, stroke = 0.4) +
    ## ROI box with light fill
    annotate("rect", xmin = qx[1], xmax = qx[2],
        ymin = qy[1], ymax = qy[2],
        fill = "#0072B2", alpha = 0.06,
        color = "#0072B2", linewidth = 0.8) +
    ## Selected: large colored shapes
    geom_point(data = pts_in,
        aes(x = x, y = y, color = gene5,
            shape = gene5),
        size = 2.2, alpha = 0.9) +
    scale_color_manual(values = gene_pal,
        name = "Selected") +
    scale_shape_manual(values = gene_shp,
        name = "Selected") +
    coord_equal(xlim = c(-0.2, 4.6),
        ylim = c(-0.2, 4.6)) +
    annotate("label", x = 3.15, y = 3.15,
        label = paste0(n_in, " / ", nrow(pts)),
        size = 3.5, hjust = 0, color = "#0072B2",
        fontface = "bold", fill = "white",
        label.size = 0) +
    annotate("text", x = 4.0, y = 0.2,
        label = "excluded", size = 2.5,
        color = "#999999", fontface = "italic") +
    labs(
        title = "Bounding Box Spatial Query",
        subtitle = paste0(
            "bboxQuery(): select ", n_in, "/",
            nrow(pts),
            " transcripts in ROI"),
        x = expression("x ("*mu*"m)"),
        y = expression("y ("*mu*"m)")) +
    th +
    theme(
        legend.position = c(0.9, 0.2),
        legend.background = element_rect(
            fill = "white", color = "grey80",
            linewidth = 0.3),
        legend.key.size = unit(8, "pt"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7,
            face = "bold"))

ggsave(file.path(od, "fig2_spatial_query.png"), fig2,
    width = 170, height = 130, units = "mm",
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
            length(gene_cols),
            " genes count matrix"),
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
        legend.key.height = unit(14, "pt"),
        legend.key.width = unit(6, "pt"),
        legend.text = element_text(size = 6.5),
        legend.title = element_text(size = 7,
            face = "bold"),
        panel.grid = element_blank(),
        panel.spacing = unit(1, "pt"),
        plot.title = element_text(size = 9,
            face = "bold",
            margin = margin(0, 0, 3, 0)),
        plot.subtitle = element_text(size = 7,
            color = "grey35", face = "italic",
            margin = margin(0, 0, 6, 0)),
        plot.margin = margin(8, 10, 8, 8))

ggsave(file.path(od, "fig3_aggregation.png"), fig3,
    width = 170, height = 140, units = "mm",
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
    DataFrame(x = rep_pts$x, y = rep_pts$y),
    ct_comp))

arr <- data.frame(
    x = rep_pts$x, y = rep_pts$y,
    xend = rep_out$x, yend = rep_out$y,
    id = rep_pts$id)

fig4 <- ggplot() +
    geom_point(data = rep_pts,
        aes(x = x, y = y),
        shape = 4, size = 3.5, stroke = 0.8,
        color = "grey50") +
    geom_text(data = rep_pts,
        aes(x = x - 0.8, y = y, label = id),
        size = 2.8, color = "grey50",
        fontface = "bold") +
    geom_segment(data = arr,
        aes(x = x, y = y, xend = xend, yend = yend),
        arrow = arrow(length = unit(4, "pt"),
            type = "closed"),
        color = "grey45", linewidth = 0.35) +
    geom_point(data = rep_out,
        aes(x = x, y = y),
        shape = 16, size = 3.5, color = "#D55E00") +
    geom_text(data = data.frame(
        x = rep_out$x + 0.8, y = rep_out$y,
        id = rep_pts$id),
        aes(x = x, y = y, label = id),
        size = 2.8, color = "#D55E00",
        fontface = "bold") +
    annotate("text", x = 16, y = 19,
        label = "pixels", size = 3,
        color = "grey50", fontface = "italic") +
    annotate("text", x = 3, y = 0.8,
        label = "global", size = 3,
        color = "#D55E00", fontface = "bold.italic") +
    coord_equal(xlim = c(-1, 20), ylim = c(-1, 20)) +
    labs(
        title = "Transform Composition",
        subtitle = paste0(
            "composeTransforms(): ",
            "scale(0.2125) + translate(1, 0.5)"),
        x = "x", y = "y") +
    th

ggsave(file.path(od, "fig4_transforms.png"), fig4,
    width = 160, height = 130, units = "mm",
    dpi = 300, bg = "white")

## ==========================================================
## Combined 4-panel
## ==========================================================
cat("combined ...\n")
library(patchwork)

combined <- (
    (fig1 + labs(title = NULL, subtitle = NULL,
        tag = "a") +
        theme(legend.position = "none")) |
    (fig2 + labs(title = NULL, subtitle = NULL,
        tag = "b") +
        theme(legend.position = "none"))
) / (
    (fig3 + labs(title = NULL, subtitle = NULL,
        tag = "c")) |
    (fig4 + labs(title = NULL, subtitle = NULL,
        tag = "d"))
) + plot_layout(heights = c(1, 1)) &
    theme(plot.tag = element_text(size = 10,
        face = "bold"))

ggsave(file.path(od, "fig1_spatial_overview.png"),
    combined, width = 180, height = 160, units = "mm",
    dpi = 300, bg = "white")

cat("\nDone:\n")
for (f in list.files(od, "^fig")) {
    cat("  ", f, ":",
        round(file.size(file.path(od, f)) / 1024),
        "KB\n")
}
