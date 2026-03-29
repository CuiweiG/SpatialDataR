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

th <- theme_classic(base_size = 9) +
    theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8,
            color = "black"),
        plot.margin = margin(8, 10, 8, 8))

pts$gene5 <- ifelse(pts$gene %in% top5, pts$gene, NA)
pts_top <- pts[!is.na(pts$gene5), ]

## ==========================================================
## Fig 1: cells | transcripts (2-panel, no overlay)
## ==========================================================
cat("fig1 ...\n")

p1a <- ggplot(shp_ann, aes(x = x, y = y,
    color = cell_type)) +
    geom_point(shape = 21, size = 7, stroke = 1.0,
        aes(fill = cell_type), alpha = 0.15) +
    geom_point(shape = 1, size = 7, stroke = 0.8) +
    scale_color_manual(values = cell_pal, name = NULL) +
    scale_fill_manual(values = cell_pal, guide = "none") +
    coord_equal(xlim = c(-0.1, 4.5),
        ylim = c(-0.1, 4.5)) +
    labs(x = expression(mu*"m"),
        y = expression(mu*"m")) +
    th + theme(legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.key.size = unit(10, "pt"))

p1b <- ggplot(pts_top, aes(x = x, y = y,
    color = gene5)) +
    geom_point(size = 1.8, alpha = 0.75) +
    scale_color_manual(values = gene_pal, name = NULL) +
    coord_equal(xlim = c(-0.1, 4.5),
        ylim = c(-0.1, 4.5)) +
    labs(x = expression(mu*"m"), y = NULL) +
    th + theme(legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.key.size = unit(10, "pt"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

fig1 <- (p1a + labs(tag = "a")) +
    (p1b + labs(tag = "b")) +
    plot_layout(ncol = 2) +
    plot_annotation(
        title = "readSpatialData()",
        subtitle = "50 cell boundaries + 500 transcripts from one Zarr store",
        theme = theme(
            plot.title = element_text(size = 12,
                face = "bold"),
            plot.subtitle = element_text(size = 9,
                color = "grey30",
                face = "italic"))) &
    theme(plot.tag = element_text(size = 11,
        face = "bold"))

ggsave(file.path(od, "fig1_store_reading.png"), fig1,
    width = 180, height = 100, units = "mm",
    dpi = 300, bg = "white")

## ==========================================================
## Fig 2: Hard-clipped bbox query
## ==========================================================
cat("fig2 ...\n")

qx <- c(1, 3); qy <- c(1, 3)
pts$inside <- pts$x >= qx[1] & pts$x <= qx[2] &
    pts$y >= qy[1] & pts$y <= qy[2]
n_in <- sum(pts$inside)
pts_sel <- pts[pts$inside & !is.na(pts$gene5), ]

p2a <- ggplot(pts, aes(x = x, y = y)) +
    geom_point(size = 0.8, alpha = 0.4,
        color = "grey40") +
    annotate("rect", xmin = qx[1], xmax = qx[2],
        ymin = qy[1], ymax = qy[2],
        fill = NA, color = "#D55E00",
        linewidth = 1.0, linetype = "dashed") +
    coord_equal(xlim = c(-0.1, 4.5),
        ylim = c(-0.1, 4.5)) +
    labs(x = expression(mu*"m"),
        y = expression(mu*"m")) + th

p2b <- ggplot(pts_sel, aes(x = x, y = y,
    color = gene5)) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_manual(values = gene_pal,
        name = NULL) +
    coord_equal(xlim = c(1, 3), ylim = c(1, 3),
        expand = FALSE) +
    labs(x = expression(mu*"m"), y = NULL) +
    th + theme(legend.position = "right",
        legend.text = element_text(size = 8),
        legend.key.size = unit(10, "pt"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        panel.border = element_rect(
            color = "#D55E00", linewidth = 0.8,
            fill = NA))

fig2 <- (p2a + labs(tag = "a")) +
    (p2b + labs(tag = "b")) +
    plot_layout(ncol = 2, widths = c(1, 1)) +
    plot_annotation(
        title = "bboxQuery()",
        subtitle = paste0(n_in, "/", nrow(pts),
            " transcripts in [1, 3] \u00D7 [1, 3] \u00B5m"),
        theme = theme(
            plot.title = element_text(size = 12,
                face = "bold"),
            plot.subtitle = element_text(size = 9,
                color = "grey30",
                face = "italic"))) &
    theme(plot.tag = element_text(size = 11,
        face = "bold"))

ggsave(file.path(od, "fig2_spatial_query.png"), fig2,
    width = 190, height = 95, units = "mm",
    dpi = 300, bg = "white")

## ==========================================================
## Fig 3: DOT PLOT — mean expression + % expressing
## Scanpy/Seurat-style, aggregated by cell type
## ==========================================================
cat("fig3 ...\n")

counts_df <- aggregatePoints(
    spatialPoints(sd)[["transcripts"]],
    shapes(sd)[["cell_boundaries"]])
counts_r <- as.data.frame(counts_df)
gene_cols <- setdiff(colnames(counts_r), "cell_id")
counts_ann <- merge(counts_r,
    obs[, c("cell_id", "cell_type")], by = "cell_id")

## Compute per-cell-type: mean fraction + pct expressing
ct_levels <- c("Epithelial", "Stromal",
    "Immune", "Endothelial")

dot_data <- do.call(rbind, lapply(ct_levels,
    function(ct) {
    sub <- counts_ann[counts_ann$cell_type == ct,
        gene_cols]
    rs <- rowSums(sub); rs[rs == 0] <- 1
    frac <- sweep(sub, 1, rs, "/")
    do.call(rbind, lapply(gene_cols, function(g) {
        vals <- frac[[g]]
        data.frame(
            cell_type = ct, gene = g,
            mean_frac = mean(vals),
            pct_expr = 100 * mean(vals > 0))
    }))
}))

dot_data$cell_type <- factor(dot_data$cell_type,
    levels = rev(ct_levels))
## Order genes: epithelial > immune > stromal > other
gene_order <- c("EPCAM", "KRT18", "ESR1", "PGR",
    "HER2", "CD45", "MKI67",
    "VIM", "ACTB", "ERBB2")
gene_order <- intersect(gene_order, gene_cols)
dot_data$gene <- factor(dot_data$gene,
    levels = gene_order)

## Show 3 clean cell types
dot_data <- dot_data[dot_data$cell_type %in%
    c("Epithelial", "Immune", "Stromal"), ]
dot_data$cell_type <- factor(dot_data$cell_type,
    levels = rev(c("Epithelial", "Immune", "Stromal")))

fig3 <- ggplot(dot_data, aes(x = gene,
    y = cell_type)) +
    geom_point(aes(size = pct_expr,
        color = mean_frac)) +
    scale_size_continuous(name = "% expressing",
        range = c(0.5, 8),
        breaks = c(20, 50, 80)) +
    scale_color_viridis_c(name = "Mean fraction",
        option = "viridis",
        limits = c(0, 0.45),
        oob = scales::squish) +
    labs(title = "aggregatePoints()",
        subtitle = paste0(
            "Molecule counts \u2192 cell-type ",
            "expression profiles (", nrow(counts_ann),
            " cells, ", length(gene_cols), " genes)"),
        x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(
        axis.text.x = element_text(angle = 45,
            hjust = 1, size = 8, color = "black",
            face = "italic"),
        axis.text.y = element_text(size = 8,
            color = "black"),
        panel.grid.major = element_line(
            color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7,
            face = "bold"),
        plot.title = element_text(size = 12,
            face = "bold"),
        plot.subtitle = element_text(size = 9,
            color = "grey30", face = "italic"),
        plot.margin = margin(8, 10, 8, 8))

ggsave(file.path(od, "fig3_aggregation.png"), fig3,
    width = 160, height = 90, units = "mm",
    dpi = 300, bg = "white")

## ==========================================================
## Fig 4: Transform composition
## ==========================================================
cat("fig4 ...\n")

rep_pts <- data.frame(
    x = c(2, 10, 16, 5, 14),
    y = c(4, 16, 8, 12, 18),
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
    xend = rep_out$x, yend = rep_out$y)

fig4 <- ggplot() +
    geom_point(data = rep_pts, aes(x = x, y = y),
        shape = 4, size = 5, stroke = 1.2,
        color = "#555555") +
    geom_text(data = rep_pts,
        aes(x = x, y = y + 1.2,
            label = paste0(id, "[px]")),
        size = 3, color = "#555555") +
    geom_segment(data = arr,
        aes(x = x, y = y, xend = xend, yend = yend),
        arrow = arrow(length = unit(6, "pt"),
            type = "closed"),
        color = "grey50", linewidth = 0.5) +
    geom_point(data = rep_out, aes(x = x, y = y),
        shape = 16, size = 5, color = "#D55E00") +
    geom_text(data = data.frame(
        x = rep_out$x, y = rep_out$y - 0.8,
        id = rep_pts$id),
        aes(x = x, y = y,
            label = paste0(id, "[gl]")),
        size = 3, color = "#D55E00",
        fontface = "bold") +
    annotate("label", x = 12, y = 1,
        label = expression(
            bold("T") == "scale(0.2125)"
            %*% "translate(1, 0.5)"),
        size = 3, fill = "grey95",
        label.size = 0.3, color = "grey30") +
    coord_equal(xlim = c(-1, 20), ylim = c(-1, 20)) +
    labs(title = "composeTransforms()",
        subtitle = paste0("5 landmarks: ",
            "pixel [px] \u2192 global [gl] ",
            "via composed affine"),
        x = "x", y = "y") +
    th + theme(
        plot.title = element_text(size = 12,
            face = "bold"),
        plot.subtitle = element_text(size = 9,
            color = "grey30", face = "italic"))

ggsave(file.path(od, "fig4_transforms.png"), fig4,
    width = 140, height = 130, units = "mm",
    dpi = 300, bg = "white")

cat("\nDone:\n")
for (f in list.files(od, "^fig")) {
    cat("  ", f, ":",
        round(file.size(file.path(od, f)) / 1024),
        "KB\n")
}
