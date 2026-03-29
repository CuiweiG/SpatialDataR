#!/usr/bin/env Rscript
## ============================================================
## SpatialDataR — Publication-Quality Figures
## Real data from xenium_mini.zarr test store
## ============================================================
## Style: Nature Methods / Genome Biology figure conventions
##   - Clean, minimal axes; no grey background
##   - Viridis-family colormaps (colourblind-safe)
##   - Panel labels (a–f) in bold 14pt
##   - Common theme_classic() base with fine adjustments

suppressPackageStartupMessages({
    library(SpatialDataR)
    library(S4Vectors)
    library(ggplot2)
    library(viridis)
    library(patchwork)
    library(scales)
    library(grid)
})

## ---- Output directory ------------------------------------------
out_dir <- file.path(
    system.file(package = "SpatialDataR"), "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
## Also save in man/figures for README
man_fig <- file.path(
    system.file(package = "SpatialDataR"), "..", "man", "figures")
if (dir.exists(file.path(
    system.file(package = "SpatialDataR"), ".."))) {
    man_fig <- normalizePath(man_fig, mustWork = FALSE)
    if (!dir.exists(man_fig)) dir.create(man_fig, recursive = TRUE)
}

## ---- Global theme (Nature style) -------------------------------
theme_pub <- function(base_size = 11) {
    theme_classic(base_size = base_size) %+replace%
        theme(
            text = element_text(family = ""),
            plot.title = element_text(
                size = rel(1.15), face = "bold",
                hjust = 0, margin = margin(b = 8)),
            axis.title = element_text(size = rel(0.95)),
            axis.text = element_text(
                size = rel(0.85), colour = "grey30"),
            legend.title = element_text(
                size = rel(0.9), face = "bold"),
            legend.text = element_text(size = rel(0.8)),
            legend.key.size = unit(0.4, "cm"),
            strip.text = element_text(
                size = rel(0.9), face = "bold"),
            strip.background = element_blank(),
            panel.grid = element_blank(),
            plot.margin = margin(8, 12, 8, 8)
        )
}
theme_set(theme_pub())

## ---- 1. Load real data -----------------------------------------
store <- system.file("extdata", "xenium_mini.zarr",
    package = "SpatialDataR")
sd <- readSpatialData(store)
cat("Loaded SpatialData:\n")
show(sd)

## Extract points DataFrame
pts_df <- as.data.frame(spatialPoints(sd)[["transcripts"]])
cat("\nTranscripts:", nrow(pts_df), "molecules,",
    length(unique(pts_df$gene)), "genes,",
    length(unique(pts_df$cell_id)), "cells\n")

## Shapes (cell boundaries as circles)
shapes_df <- as.data.frame(shapes(sd)[["cell_boundaries"]])

## Validation
val <- validateSpatialData(store)
cat("\nValidation: valid =", val$valid, "\n")

## ============================================================
## Panel (a): Spatial transcript map — all genes
## ============================================================
## Colourblind-safe palette (Okabe-Ito extended, 10 distinct)
gene_pal <- c(
    "ACTB"  = "#E69F00", "CD45"  = "#56B4E9",
    "EPCAM" = "#009E73", "ERBB2" = "#0072B2",
    "ESR1"  = "#D55E00", "HER2"  = "#CC79A7",
    "KRT18" = "#F0E442", "MKI67" = "#000000",
    "PGR"   = "#882255", "VIM"   = "#999999"
)

p_a <- ggplot(pts_df, aes(x = x, y = y, colour = gene)) +
    geom_point(size = 1.0, alpha = 0.75, stroke = 0) +
    scale_colour_manual(
        values = gene_pal,
        name = "Gene") +
    coord_fixed(ratio = 1) +
    labs(
        title = "Spatial transcript map",
        subtitle = paste0(
            nrow(pts_df), " molecules \u00b7 ",
            length(unique(pts_df$gene)), " genes"),
        x = expression(italic(x) ~ "(spatial units)"),
        y = expression(italic(y) ~ "(spatial units)")) +
    guides(colour = guide_legend(
        override.aes = list(size = 2.5, alpha = 1),
        ncol = 2)) +
    theme(legend.position = "right")

## ============================================================
## Panel (b): Bounding-box query demonstration
## ============================================================
bbox_xmin <- 1.5; bbox_xmax <- 3.5
bbox_ymin <- 1.0; bbox_ymax <- 3.0

pts_bbox <- bboxQuery(
    spatialPoints(sd)[["transcripts"]],
    xmin = bbox_xmin, xmax = bbox_xmax,
    ymin = bbox_ymin, ymax = bbox_ymax)
pts_bbox_df <- as.data.frame(pts_bbox)

cat("\nbboxQuery: ", nrow(pts_df), " -> ",
    nrow(pts_bbox_df), " molecules\n")

## Build rect annotation
bbox_rect <- data.frame(
    xmin = bbox_xmin, xmax = bbox_xmax,
    ymin = bbox_ymin, ymax = bbox_ymax)

p_b <- ggplot(pts_df, aes(x = x, y = y)) +
    geom_point(colour = "grey75", size = 0.6,
        alpha = 0.4, stroke = 0) +
    geom_rect(
        data = bbox_rect,
        aes(xmin = xmin, xmax = xmax,
            ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        fill = NA, colour = "#E64B35",
        linewidth = 0.7, linetype = "dashed") +
    geom_point(
        data = pts_bbox_df,
        aes(colour = gene),
        size = 1.2, alpha = 0.85, stroke = 0) +
    scale_colour_manual(values = gene_pal, name = "Gene") +
    coord_fixed(ratio = 1) +
    labs(
        title = "Bounding-box query",
        subtitle = paste0(
            nrow(pts_bbox_df), "/", nrow(pts_df),
            " molecules retained"),
        x = expression(italic(x) ~ "(spatial units)"),
        y = expression(italic(y) ~ "(spatial units)")) +
    guides(colour = guide_legend(
        override.aes = list(size = 2.5, alpha = 1),
        ncol = 2)) +
    theme(legend.position = "right")

## ============================================================
## Panel (c): Cell-by-gene aggregation heatmap
## ============================================================
regions_df <- DataFrame(
    cell_id = sort(unique(pts_df$cell_id)),
    x = as.numeric(tapply(pts_df$x, pts_df$cell_id, mean)),
    y = as.numeric(tapply(pts_df$y, pts_df$cell_id, mean)))

## Use aggregatePoints (the package's own function)
counts <- aggregatePoints(
    spatialPoints(sd)[["transcripts"]],
    regions_df,
    feature_col = "gene",
    region_col = "cell_id")
count_cols <- setdiff(colnames(counts), "cell_id")
count_mat <- as.matrix(as.data.frame(counts[, count_cols]))
rownames(count_mat) <- counts$cell_id
cat("\nAggregation: ", nrow(count_mat), " cells x ",
    ncol(count_mat), " genes\n")

## Top 30 cells by total counts (readable heatmap)
total_per_cell <- rowSums(count_mat)
top_cells <- names(sort(total_per_cell, decreasing = TRUE))[
    seq_len(min(30, length(total_per_cell)))]
count_sub <- count_mat[top_cells, , drop = FALSE]

## Long form
hm_long <- expand.grid(
    cell = rownames(count_sub),
    gene = colnames(count_sub),
    stringsAsFactors = FALSE)
hm_long$count <- as.vector(count_sub)
hm_long$cell <- factor(hm_long$cell,
    levels = rev(rownames(count_sub)))
hm_long$gene <- factor(hm_long$gene,
    levels = colnames(count_sub))

p_c <- ggplot(hm_long, aes(x = gene, y = cell,
    fill = count)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_viridis(
        option = "inferno",
        name = "Counts",
        breaks = scales::pretty_breaks(4)) +
    labs(
        title = "Cell \u00d7 gene counts",
        subtitle = paste0(
            nrow(count_sub), " cells \u00b7 ",
            ncol(count_sub), " genes"),
        x = NULL, y = "Cell ID") +
    theme(
        axis.text.x = element_text(
            angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 7),
        legend.position = "right")

## ============================================================
## Panel (d): Gene density per cell (violin/dot)
## ============================================================
gene_per_cell <- data.frame(
    cell_id = rownames(count_mat),
    total = rowSums(count_mat),
    n_genes = rowSums(count_mat > 0))

## Gene frequency table (used here and in panel f)
gene_counts <- sort(table(pts_df$gene), decreasing = TRUE)

gene_totals <- data.frame(
    gene = names(gene_counts),
    count = as.integer(gene_counts),
    stringsAsFactors = FALSE)
gene_totals$gene <- factor(gene_totals$gene,
    levels = gene_totals$gene[order(gene_totals$count,
        decreasing = TRUE)])

p_d <- ggplot(gene_totals, aes(
    x = gene, y = count, fill = gene)) +
    geom_col(alpha = 0.85, width = 0.7) +
    geom_text(aes(label = count),
        vjust = -0.4, size = 3, colour = "grey30") +
    scale_fill_manual(values = gene_pal) +
    labs(
        title = "Transcript abundance",
        subtitle = paste0(
            sum(gene_totals$count), " total molecules"),
        x = NULL,
        y = "Molecule count") +
    theme(
        axis.text.x = element_text(
            angle = 45, hjust = 1, size = 9),
        legend.position = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12)))

## ============================================================
## Panel (e): Transform composition demonstration
## ============================================================
## Create a realistic transform chain: pixels → scaled → microns
scale_tf <- CoordinateTransform("affine",
    affine = diag(c(0.2125, 0.2125, 1)),
    input_cs = "pixels", output_cs = "scaled")
translate_tf <- CoordinateTransform("affine",
    affine = matrix(c(1, 0, 100, 0, 1, 200, 0, 0, 1),
        nrow = 3, byrow = TRUE),
    input_cs = "scaled", output_cs = "microns")
composed <- composeTransforms(scale_tf, translate_tf)
inv <- invertTransform(composed)

## Demo: transform a grid of points
grid_pts <- expand.grid(
    x = seq(0, 20, by = 4),
    y = seq(0, 20, by = 4))
grid_df <- DataFrame(x = grid_pts$x, y = grid_pts$y)

fwd <- as.data.frame(transformCoords(grid_df, composed))
back <- as.data.frame(
    transformCoords(DataFrame(x = fwd$x, y = fwd$y), inv))

## Combine for plotting
orig <- as.data.frame(grid_df)
orig$space <- "Pixels (input)"
fwd$space <- "Microns (output)"
back$space <- "Roundtrip (inv)"

tf_all <- rbind(orig, fwd, back)
tf_all$space <- factor(tf_all$space,
    levels = c("Pixels (input)", "Microns (output)",
        "Roundtrip (inv)"))

## Roundtrip error
rt_err <- max(abs(orig$x - back$x), abs(orig$y - back$y))

p_e <- ggplot(tf_all, aes(x = x, y = y, colour = space)) +
    geom_point(size = 2.5, alpha = 0.8) +
    facet_wrap(~ space, scales = "free", nrow = 1) +
    scale_colour_manual(values = c(
        "Pixels (input)" = "#3C5488",
        "Microns (output)" = "#E64B35",
        "Roundtrip (inv)" = "#00A087")) +
    labs(
        title = "Affine transform composition",
        subtitle = paste0(
            "pixels \u2192 microns | roundtrip error = ",
            format(rt_err, scientific = TRUE, digits = 2)),
        x = "x", y = "y") +
    theme(legend.position = "none",
        strip.text = element_text(size = 9))

## ============================================================
## Panel (f): Spatial gene expression map (per-gene facets)
## ============================================================
## Show top 6 genes by abundance
gene_counts <- sort(table(pts_df$gene), decreasing = TRUE)
top6 <- names(gene_counts)[seq_len(min(6, length(gene_counts)))]
pts_top6 <- pts_df[pts_df$gene %in% top6, ]
pts_top6$gene <- factor(pts_top6$gene, levels = top6)

p_f <- ggplot(pts_top6, aes(x = x, y = y)) +
    geom_point(
        data = pts_df[, c("x", "y")],
        colour = "grey85", size = 0.4, alpha = 0.25) +
    geom_point(aes(colour = gene),
        size = 1.3, alpha = 0.8, stroke = 0) +
    facet_wrap(~ gene, nrow = 2) +
    scale_colour_manual(values = gene_pal) +
    coord_fixed(ratio = 1) +
    labs(
        title = "Per-gene spatial distribution",
        subtitle = "Top 6 genes by abundance",
        x = expression(italic(x) ~ "(spatial units)"),
        y = expression(italic(y) ~ "(spatial units)")) +
    theme(legend.position = "none")

## ============================================================
## Assemble composite figure
## ============================================================
## Layout: 3 rows × 2 columns
## Row 1: (a) transcript map  |  (b) bbox query
## Row 2: (c) heatmap         |  (d) violin
## Row 3: (e) transforms      |  (f) per-gene facets

composite <- (
    (p_a + p_b) /
    (p_c + p_d) /
    (p_e + p_f)
) +
    plot_annotation(
        tag_levels = "a",
        title = "SpatialDataR: R-native analysis of SpatialData stores",
        subtitle = paste0(
            "Xenium mini dataset \u2014 ",
            nrow(pts_df), " transcripts, ",
            length(unique(pts_df$gene)), " genes, ",
            length(unique(pts_df$cell_id)), " cells"),
        theme = theme(
            plot.title = element_text(
                size = 16, face = "bold", hjust = 0),
            plot.subtitle = element_text(
                size = 11, colour = "grey40", hjust = 0))
    ) &
    theme(plot.tag = element_text(
        size = 14, face = "bold"))

## ---- Save figures ----------------------------------------------
fig_path <- file.path(out_dir, "SpatialDataR_publication_figure.png")
ggsave(fig_path, composite,
    width = 16, height = 18, dpi = 300, bg = "white")
cat("\nSaved:", fig_path, "\n")

## Also save as PDF (vector, preferred for publications)
pdf_path <- file.path(out_dir, "SpatialDataR_publication_figure.pdf")
ggsave(pdf_path, composite,
    width = 16, height = 18, device = cairo_pdf)
cat("Saved:", pdf_path, "\n")

## Copy PNG to man/figures if possible
if (dir.exists(man_fig)) {
    man_png <- file.path(man_fig,
        "SpatialDataR_publication_figure.png")
    file.copy(fig_path, man_png, overwrite = TRUE)
    cat("Copied to:", man_png, "\n")
}

## ---- Individual panels (high-res) ------------------------------
for (panel_name in c("a", "b", "c", "d", "e", "f")) {
    p <- switch(panel_name,
        "a" = p_a, "b" = p_b, "c" = p_c,
        "d" = p_d, "e" = p_e, "f" = p_f)
    fname <- paste0("panel_", panel_name, ".png")
    ggsave(file.path(out_dir, fname), p,
        width = 8, height = 6, dpi = 300, bg = "white")
}
cat("\nIndividual panels saved to:", out_dir, "\n")

cat("\n=== Done! ===\n")
