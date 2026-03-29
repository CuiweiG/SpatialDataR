#!/usr/bin/env Rscript
# ============================================================================
# SpatialDataR: Publication-Quality 6-Panel Figure
# Full analytical pipeline on real Xenium mini data
# Output: 300 DPI PNG + PDF, Nature Methods style
# ============================================================================

## ---- Setup ----------------------------------------------------------------
suppressPackageStartupMessages({
    library(SpatialDataR)
    library(S4Vectors)
    library(ggplot2)
    library(patchwork)
    library(viridis)
    library(scales)
    library(pheatmap)
    library(reshape2)
    library(grid)
    library(gridExtra)
    library(ggrepel)
    library(RColorBrewer)
})

## ---- Publication theme ----------------------------------------------------
theme_pub <- function(base_size = 10.5) {
    theme_classic(base_size = base_size) %+replace%
    theme(
        text             = element_text(family = "sans", colour = "grey10"),
        plot.title       = element_text(face = "bold", size = rel(1.1),
                                        hjust = 0, margin = margin(b = 3)),
        plot.subtitle    = element_text(size = rel(0.82), colour = "grey40",
                                        hjust = 0, margin = margin(b = 5)),
        axis.title       = element_text(size = rel(0.9), face = "bold"),
        axis.text        = element_text(size = rel(0.82), colour = "grey30"),
        axis.line        = element_line(linewidth = 0.4, colour = "grey30"),
        axis.ticks       = element_line(linewidth = 0.3, colour = "grey30"),
        legend.title     = element_text(size = rel(0.88), face = "bold"),
        legend.text      = element_text(size = rel(0.78)),
        legend.key.size  = unit(0.4, "cm"),
        legend.background = element_rect(fill = NA),
        plot.margin      = margin(6, 10, 6, 8),
        panel.grid       = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text       = element_text(face = "bold", size = rel(0.88))
    )
}

## ---- Colour palettes ------------------------------------------------------
cell_cols <- c(
    "Epithelial"   = "#E64B35",
    "Stromal"      = "#4DBBD5",
    "Immune"       = "#00A087",
    "Endothelial"  = "#3C5488"
)

## 10 maximally perceptually distinct colours (colorblind-aware)
gene_cols <- c(
    "EPCAM" = "#E64B35",   # vermillion
    "KRT18" = "#FF7F0E",   # orange
    "VIM"   = "#1F77B4",   # blue
    "CD45"  = "#2CA02C",   # green
    "HER2"  = "#9467BD",   # purple
    "ESR1"  = "#D62728",   # crimson
    "PGR"   = "#BCBD22",   # olive/yellow-green
    "MKI67" = "#17BECF",   # cyan
    "ERBB2" = "#8C564B",   # brown
    "ACTB"  = "#E377C2"    # pink
)

## ---- Scale bar helper -----------------------------------------------------
add_scalebar <- function(xrng, yrng, bar_um = 1, label = NULL) {
    x_end   <- xrng[2] - diff(xrng) * 0.04
    x_start <- x_end - bar_um
    y_pos   <- yrng[1] + diff(yrng) * 0.06
    if (is.null(label)) label <- paste0(bar_um, " \u00B5m")
    list(
        annotate("segment", x = x_start, xend = x_end,
                 y = y_pos, yend = y_pos,
                 linewidth = 1.5, colour = "black", lineend = "butt"),
        annotate("text", x = (x_start + x_end) / 2, y = y_pos,
                 label = label, vjust = -0.9, size = 3, fontface = "bold")
    )
}

## ---- Load data ------------------------------------------------------------
store_path <- system.file("extdata", "xenium_mini.zarr",
    package = "SpatialDataR")
cat("Reading SpatialData store:", store_path, "\n")
sd <- readSpatialData(store_path)
show(sd)
cat("\n")

pts_df   <- as.data.frame(spatialPoints(sd)[[1]])
shp_df   <- as.data.frame(shapes(sd)[[1]])
tbl      <- tables(sd)[[1]]
obs_df   <- as.data.frame(tbl$obs)
shp_ann  <- merge(shp_df, obs_df, by = "cell_id")
x_rng    <- range(pts_df$x)
y_rng    <- range(pts_df$y)

val <- validateSpatialData(store_path)
cat("Validation: valid=", val$valid, "\n\n")

## ---- Panel A: Transcript spatial map by gene ------------------------------
pa <- ggplot(pts_df, aes(x = x, y = y, colour = gene)) +
    geom_point(size = 0.8, alpha = 0.78, stroke = 0) +
    scale_colour_manual(values = gene_cols, name = "Gene") +
    add_scalebar(x_rng, y_rng, bar_um = 1) +
    coord_fixed(ratio = 1) +
    labs(title = "Transcript spatial map",
         subtitle = paste0(nrow(pts_df), " molecules, ",
                           length(unique(pts_df$gene)), " genes"),
         x = expression("x ("*mu*"m)"),
         y = expression("y ("*mu*"m)")) +
    theme_pub() +
    guides(colour = guide_legend(
        override.aes = list(size = 3, alpha = 1), ncol = 2))

## ---- Panel B: Cell segmentation & phenotype -------------------------------
pb <- ggplot() +
    geom_point(data = pts_df, aes(x = x, y = y),
               size = 0.12, alpha = 0.10, colour = "grey50") +
    geom_point(data = shp_ann,
               aes(x = x, y = y, fill = cell_type, size = radius),
               shape = 21, colour = "grey25", stroke = 0.35, alpha = 0.6) +
    scale_fill_manual(values = cell_cols, name = "Cell type") +
    scale_size_continuous(range = c(2.5, 9), guide = "none") +
    add_scalebar(x_rng, y_rng, bar_um = 1) +
    coord_fixed(ratio = 1) +
    labs(title = "Cell segmentation & phenotype",
         subtitle = paste0(nrow(shp_ann), " cells, ",
                           length(unique(shp_ann$cell_type)), " types"),
         x = expression("x ("*mu*"m)"),
         y = expression("y ("*mu*"m)")) +
    theme_pub() +
    guides(fill = guide_legend(override.aes = list(size = 5)))

## ---- Panel C: bboxQuery ---------------------------------------------------
bb <- list(xmin = 1.0, xmax = 3.5, ymin = 1.0, ymax = 3.5)
pts_sub <- as.data.frame(bboxQuery(
    spatialPoints(sd)[[1]], bb$xmin, bb$xmax, bb$ymin, bb$ymax))
pts_df$in_bbox <- pts_df$x >= bb$xmin & pts_df$x <= bb$xmax &
                  pts_df$y >= bb$ymin & pts_df$y <= bb$ymax

pc <- ggplot(pts_df, aes(x = x, y = y)) +
    geom_point(data = pts_df[!pts_df$in_bbox, ],
               size = 0.4, alpha = 0.15, colour = "grey70") +
    geom_point(data = pts_df[pts_df$in_bbox, ],
               aes(colour = gene), size = 1.0, alpha = 0.85) +
    annotate("rect", xmin = bb$xmin, xmax = bb$xmax,
             ymin = bb$ymin, ymax = bb$ymax,
             fill = NA, colour = "#E64B35", linewidth = 1,
             linetype = "dashed") +
    annotate("label", x = (bb$xmin + bb$xmax) / 2, y = bb$ymax + 0.2,
             label = paste0(nrow(pts_sub), "/", nrow(pts_df), " retained"),
             size = 2.8, colour = "#E64B35", fontface = "bold",
             fill = "white", label.r = unit(0.15, "lines")) +
    scale_colour_manual(values = gene_cols, guide = "none") +
    add_scalebar(x_rng, y_rng, bar_um = 1) +
    coord_fixed(ratio = 1) +
    labs(title = "Bounding box spatial query",
         subtitle = "bboxQuery() on SpatialData object",
         x = expression("x ("*mu*"m)"),
         y = expression("y ("*mu*"m)")) +
    theme_pub()

## ---- Panel D: ggplot2 heatmap (avoid pheatmap grob issues) ----------------
counts <- aggregatePoints(
    spatialPoints(sd)[[1]], shapes(sd)[[1]],
    feature_col = "gene", region_col = "cell_id")
counts_df <- as.data.frame(counts)
rownames(counts_df) <- counts_df$cell_id
counts_df$cell_id <- NULL

ct_map <- setNames(obs_df$cell_type, obs_df$cell_id)
row_ann <- data.frame(cell_id = rownames(counts_df),
                       CellType = ct_map[rownames(counts_df)],
                       stringsAsFactors = FALSE)
row_ann <- row_ann[order(row_ann$CellType), ]
counts_mat <- as.matrix(counts_df[row_ann$cell_id, ])

## Melt for ggplot heatmap
hm_df <- reshape2::melt(counts_mat)
colnames(hm_df) <- c("Cell", "Gene", "Count")
hm_df$Cell <- factor(hm_df$Cell, levels = row_ann$cell_id)
hm_df$CellType <- ct_map[as.character(hm_df$Cell)]

## Gene clustering order via hclust
gene_dist <- dist(t(counts_mat))
gene_hc <- hclust(gene_dist, method = "complete")
hm_df$Gene <- factor(hm_df$Gene, levels = gene_hc$labels[gene_hc$order])

pd <- ggplot(hm_df, aes(x = Gene, y = Cell, fill = Count)) +
    geom_tile(colour = NA) +
    scale_fill_viridis_c(option = "inferno", name = "Counts") +
    ## Cell type side bar using facet trick: annotate via colour strip
    labs(title = "Cell x gene count matrix",
         subtitle = paste0("aggregatePoints(): ", nrow(counts_mat),
                           " cells x ", ncol(counts_mat), " genes"),
         x = NULL, y = "Cell ID (grouped by type)") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.9),
                                      face = "bold"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "right")

## Add cell type colour strip as a separate narrow panel
ct_strip <- ggplot(data.frame(Cell = factor(row_ann$cell_id,
                                            levels = row_ann$cell_id),
                              CellType = row_ann$CellType),
                   aes(x = 1, y = Cell, fill = CellType)) +
    geom_tile(colour = NA) +
    scale_fill_manual(values = cell_cols, name = "Cell type") +
    theme_void() +
    theme(legend.position = "right",
          legend.title = element_text(size = 9, face = "bold"),
          legend.text = element_text(size = 8)) +
    guides(fill = guide_legend(override.aes = list(colour = NA)))

pd_combined <- ct_strip + pd + plot_layout(widths = c(0.04, 1),
                                            guides = "collect")

## ---- Panel E: Transform pipeline -----------------------------------------
scale_tf <- CoordinateTransform("affine",
    affine = diag(c(0.2125, 0.2125, 1)),
    input_cs = "pixels", output_cs = "microns")
translate_tf <- CoordinateTransform("affine",
    affine = matrix(c(1,0,100, 0,1,200, 0,0,1), nrow = 3, byrow = TRUE),
    input_cs = "microns", output_cs = "global")
composed <- composeTransforms(scale_tf, translate_tf)
inv <- invertTransform(composed)

pts_s4 <- spatialPoints(sd)[[1]]
pts_micro <- transformCoords(pts_s4, scale_tf)
pts_global <- transformCoords(pts_s4, composed)
pts_roundtrip <- transformCoords(pts_global, inv)

tf_data <- rbind(
    data.frame(x = pts_s4$x, y = pts_s4$y,
               stage = "Original (px)"),
    data.frame(x = as.numeric(pts_micro$x),
               y = as.numeric(pts_micro$y),
               stage = "Scaled (um)"),
    data.frame(x = as.numeric(pts_global$x),
               y = as.numeric(pts_global$y),
               stage = "Composed (global)")
)
tf_data$stage <- factor(tf_data$stage,
    levels = c("Original (px)", "Scaled (um)", "Composed (global)"))
rterr <- max(abs(pts_s4$x - pts_roundtrip$x),
             abs(pts_s4$y - pts_roundtrip$y))

pe <- ggplot(tf_data, aes(x = x, y = y)) +
    geom_point(size = 0.5, alpha = 0.55, colour = "#3C5488") +
    facet_wrap(~ stage, scales = "free", nrow = 1) +
    labs(title = "Coordinate transform pipeline",
         subtitle = paste0("composeTransforms() + invertTransform()  |  ",
                           "roundtrip error = ",
                           formatC(rterr, format = "e", digits = 1)),
         x = "x", y = "y") +
    theme_pub() +
    theme(strip.text = element_text(size = rel(0.88)))

## ---- Panel F: Violin + boxplot by cell type --------------------------------
top3 <- c("ERBB2", "VIM", "CD45")
total_per_cell <- data.frame(
    cell_id = as.integer(rownames(counts_mat)),
    total = rowSums(counts_mat))
total_per_cell <- merge(total_per_cell, obs_df, by = "cell_id")
for (g in top3) {
    if (g %in% colnames(counts_mat))
        total_per_cell[[g]] <- counts_mat[
            as.character(total_per_cell$cell_id), g]
}
gene_long <- reshape2::melt(total_per_cell,
    id.vars = c("cell_id", "cell_type", "region", "total"),
    measure.vars = top3,
    variable.name = "Gene", value.name = "Count")

pf <- ggplot(gene_long, aes(x = cell_type, y = Count, fill = cell_type)) +
    geom_violin(alpha = 0.5, linewidth = 0.3, scale = "width") +
    geom_boxplot(width = 0.15, outlier.size = 0.6, alpha = 0.85,
                 linewidth = 0.3) +
    facet_wrap(~ Gene, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = cell_cols, guide = "none") +
    labs(title = "Gene expression by cell type",
         subtitle = "aggregatePoints() -> cell-by-gene count matrix",
         x = NULL, y = "Transcript count") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = rel(0.85)))

## ---- Assemble composite figure --------------------------------------------
## Wrap the combined heatmap as a single element to avoid extra tag
pd_wrapped <- wrap_elements(full = patchworkGrob(pd_combined))

composite <- (pa | pb | pc) /
             (pd_wrapped | pe) /
             pf +
    plot_layout(heights = c(1, 0.9, 0.55)) +
    plot_annotation(
        title = "SpatialDataR: R-native SpatialData analysis pipeline",
        subtitle = paste0("Xenium mini dataset  |  ",
                          nrow(pts_df), " transcripts  |  ",
                          nrow(shp_ann), " cells  |  ",
                          length(unique(pts_df$gene)), " genes"),
        tag_levels = "a",
        theme = theme(
            plot.title = element_text(face = "bold", size = 16,
                                      hjust = 0, family = "sans"),
            plot.subtitle = element_text(size = 11, colour = "grey40",
                                         hjust = 0, margin = margin(b = 8))
        )
    ) &
    theme(plot.tag = element_text(face = "bold", size = 14))

## ---- Save -----------------------------------------------------------------
out_dir <- "C:/Users/win10/SpatialDataR/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_png <- file.path(out_dir, "SpatialDataR_publication_figure.png")
out_pdf <- file.path(out_dir, "SpatialDataR_publication_figure.pdf")

cat("Saving PNG...\n")
png(out_png, width = 17, height = 13.5, units = "in", res = 300,
    type = "windows")
print(composite)
dev.off()

cat("Saving PDF...\n")
pdf(out_pdf, width = 17, height = 13.5)
print(composite)
dev.off()

cat("\n=== COMPLETE ===\n")
cat("PNG:", normalizePath(out_png), "\n")
cat("PDF:", normalizePath(out_pdf), "\n")
cat("Dimensions: 17 x 13.5 in @ 300 DPI\n")
cat("6 panels (a-f): transcript map, cell phenotype, bboxQuery,\n")
cat("  heatmap, coordinate transforms, violin expression\n")
