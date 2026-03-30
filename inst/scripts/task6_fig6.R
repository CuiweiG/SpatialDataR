## Task 6: Multimodal visualization — publication quality rewrite
## Output: man/figures/fig6_multimodal.png

.libPaths(c("C:/Users/win10/R/win-library/4.4",
            "C:/Users/win10/AppData/Local/R/win-library/4.4",
            .libPaths()))

library(SpatialDataR)
library(S4Vectors)
library(ggplot2)
library(viridis)
library(patchwork)
library(Rarr)
library(scales)

cat("=== Task 6: Fig 6 — Multimodal visualization ===\n")

## ---- Xenium mini data ----
store <- system.file("extdata", "xenium_mini.zarr",
    package = "SpatialDataR")

img_path <- file.path(store, "images", "morphology", "scale0")
img <- readZarrArray(img_path)
cat("Image dims:", dim(img), "\n")

ny <- dim(img)[2]; nx <- dim(img)[3]
img_df <- expand.grid(y = seq_len(ny), x = seq_len(nx))
img_df$gray <- (as.numeric(img[1,,]) + as.numeric(img[2,,]) +
                as.numeric(img[3,,])) / (3 * 255)

x_scale <- 4.0 / nx; y_scale <- 4.0 / ny
img_df$x_s <- (img_df$x - 0.5) * x_scale
img_df$y_s <- (img_df$y - 0.5) * y_scale

pts <- read.csv(file.path(store, "points/transcripts/transcripts.csv"),
    stringsAsFactors = FALSE)
cells <- read.csv(file.path(store, "shapes/cell_boundaries/circles.csv"),
    stringsAsFactors = FALSE)

## Top 6 genes
gene_counts <- sort(table(pts$gene), decreasing = TRUE)
top6 <- names(gene_counts)[1:min(6, length(gene_counts))]
pts$gene_group <- ifelse(pts$gene %in% top6, pts$gene, "Other")
pts$gene_group <- factor(pts$gene_group, levels = c(top6, "Other"))

gene_pal <- c("#E64B35", "#4DBBD5", "#00A087",
              "#3C5488", "#F39B7F", "#8491B4", "#B2B2B2")
names(gene_pal) <- c(top6, "Other")

## ---- Panel A: Greyscale raster + transcript overlay ----
p_a <- ggplot() +
    geom_tile(data = img_df,
        aes(x = x_s, y = y_s),
        fill = grey(img_df$gray)) +
    geom_point(data = pts,
        aes(x = x, y = y, colour = gene_group),
        size = 1.4, alpha = 0.85) +
    scale_colour_manual(values = gene_pal, name = "Gene") +
    coord_fixed(expand = FALSE) +
    theme_classic(base_size = 10) +
    labs(x = expression("x ("*mu*"m)"),
         y = expression("y ("*mu*"m)")) +
    theme(legend.position = "right",
          legend.key.size = unit(3.5, "mm"),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9, face = "bold"),
          plot.margin = margin(5, 5, 5, 5))

## ---- Panel B: Cell segmentation mask (no ID legend — use continuous fill) ----
label_mat <- matrix(0L, nrow = ny, ncol = nx)
for (i in seq_len(nrow(cells))) {
    cx_px <- cells$x[i] / x_scale + 0.5
    cy_px <- cells$y[i] / y_scale + 0.5
    r_px  <- cells$radius[i] / x_scale
    for (py in seq_len(ny)) {
        for (px in seq_len(nx)) {
            if (sqrt((px - cx_px)^2 + (py - cy_px)^2) <= r_px)
                label_mat[py, px] <- cells$cell_id[i]
        }
    }
}
lbl_df <- expand.grid(y = seq_len(ny), x = seq_len(nx))
lbl_df$cell_id <- as.vector(label_mat)
lbl_df$x_s <- (lbl_df$x - 0.5) * x_scale
lbl_df$y_s <- (lbl_df$y - 0.5) * y_scale

active <- lbl_df[lbl_df$cell_id > 0, ]

p_b <- ggplot() +
    geom_tile(data = img_df,
        aes(x = x_s, y = y_s),
        fill = grey(img_df$gray)) +
    geom_tile(data = active,
        aes(x = x_s, y = y_s, fill = cell_id),
        alpha = 0.6) +
    scale_fill_viridis(option = "turbo", name = "Cell ID",
        breaks = c(min(active$cell_id),
                   round(median(active$cell_id)),
                   max(active$cell_id))) +
    geom_point(data = cells,
        aes(x = x, y = y), colour = "white",
        size = 1, shape = 4, stroke = 0.7) +
    coord_fixed(expand = FALSE) +
    theme_classic(base_size = 10) +
    labs(x = expression("x ("*mu*"m)"),
         y = expression("y ("*mu*"m)")) +
    theme(legend.position = "right",
          legend.key.size = unit(3.5, "mm"),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9, face = "bold"),
          plot.margin = margin(5, 5, 5, 5))

## ---- Panel C: cropImage() — full + cropped inset ----
source("C:/Users/win10/SpatialDataR/R/image-ops.R")

crop_xmin <- 6; crop_xmax <- 15
crop_ymin <- 6; crop_ymax <- 15

crop_arr <- cropImage(img_path,
    xmin = crop_xmin, xmax = crop_xmax,
    ymin = crop_ymin, ymax = crop_ymax)

box_x1 <- (crop_xmin - 0.5) * x_scale
box_x2 <- (crop_xmax - 0.5) * x_scale
box_y1 <- (crop_ymin - 0.5) * y_scale
box_y2 <- (crop_ymax - 0.5) * y_scale

cn <- dim(crop_arr)
crop_df <- expand.grid(y = seq_len(cn[2]), x = seq_len(cn[3]))
crop_df$gray <- (as.numeric(crop_arr[1,,]) + as.numeric(crop_arr[2,,]) +
                 as.numeric(crop_arr[3,,])) / (3 * 255)
## Place cropped region in an inset at top-right
## Map crop to a small region offset to the right
crop_df$x_s <- 4.3 + (crop_df$x - 1) * 0.15
crop_df$y_s <- 2.5 + (crop_df$y - 1) * 0.15

pts_crop <- pts[pts$x >= box_x1 & pts$x <= box_x2 &
                pts$y >= box_y1 & pts$y <= box_y2, ]

p_c <- ggplot() +
    ## Full image
    geom_tile(data = img_df,
        aes(x = x_s, y = y_s),
        fill = grey(img_df$gray)) +
    ## Red crop box on full image
    annotate("rect",
        xmin = box_x1, xmax = box_x2,
        ymin = box_y1, ymax = box_y2,
        colour = "#E41A1C", fill = NA,
        linewidth = 0.9, linetype = "solid") +
    ## Cropped inset
    geom_tile(data = crop_df,
        aes(x = x_s, y = y_s),
        fill = grey(crop_df$gray)) +
    annotate("rect",
        xmin = min(crop_df$x_s) - 0.08,
        xmax = max(crop_df$x_s) + 0.08,
        ymin = min(crop_df$y_s) - 0.08,
        ymax = max(crop_df$y_s) + 0.08,
        colour = "#E41A1C", fill = NA,
        linewidth = 0.6) +
    ## Arrow from box to inset
    annotate("segment",
        x = box_x2, y = (box_y1 + box_y2) / 2,
        xend = min(crop_df$x_s) - 0.08,
        yend = mean(crop_df$y_s),
        colour = "#E41A1C", linewidth = 0.4,
        arrow = arrow(length = unit(2, "mm"), type = "closed")) +
    ## Transcripts in crop region
    geom_point(data = pts_crop,
        aes(x = x, y = y, colour = gene_group),
        size = 1.2, alpha = 0.8) +
    scale_colour_manual(values = gene_pal, guide = "none") +
    annotate("text", x = mean(crop_df$x_s), y = max(crop_df$y_s) + 0.25,
        label = "cropImage()", size = 3, fontface = "italic",
        colour = "#E41A1C") +
    coord_fixed(expand = FALSE, xlim = c(-0.2, 6), ylim = c(-0.2, 4.3)) +
    theme_classic(base_size = 10) +
    labs(x = expression("x ("*mu*"m)"),
         y = expression("y ("*mu*"m)")) +
    theme(plot.margin = margin(5, 5, 5, 5))

## ---- Panel D: MERFISH transcript density + cell centroids ----
cat("Loading MERFISH data for panel D...\n")
merfish_path <- "C:/Users/win10/merfish_spatialdata.zarr"
sd_m <- readSpatialData(merfish_path)
pts_m <- as.data.frame(spatialPoints(sd_m)[["transcripts"]])
cells_m <- as.data.frame(shapes(sd_m)[["cell_boundaries"]])

set.seed(42)
pts_sub <- pts_m[sample(nrow(pts_m), 200000), ]

p_d <- ggplot() +
    stat_density_2d(data = pts_sub,
        aes(x = x, y = y, fill = after_stat(density)),
        geom = "raster", contour = FALSE, n = 150) +
    scale_fill_viridis(option = "magma", name = "Transcript\ndensity",
        labels = function(x) {
            ifelse(x == 0, "0", sprintf("%.1f\u00d710\u207b\u2076",
                x * 1e6))
        }) +
    geom_point(data = cells_m,
        aes(x = x, y = y),
        colour = "#00FF99", size = 1.8, shape = 1, stroke = 0.8) +
    coord_fixed(expand = FALSE) +
    theme_classic(base_size = 10) +
    labs(x = expression("x ("*mu*"m)"),
         y = expression("y ("*mu*"m)")) +
    theme(legend.position = "right",
          legend.key.size = unit(3.5, "mm"),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 9, face = "bold"),
          plot.margin = margin(5, 5, 5, 5))

## ---- Assemble 2x2 ----
p_final <- (p_a + p_b) / (p_c + p_d) +
    plot_annotation(
        tag_levels = "a",
        title = "SpatialDataR: multimodal spatial data integration",
        subtitle = "Xenium mini (a\u2013c) and MERFISH (d) datasets",
        theme = theme(
            plot.title = element_text(face = "bold", size = 13),
            plot.subtitle = element_text(size = 10, colour = "grey40")
        )
    ) & theme(
        plot.tag = element_text(face = "bold", size = 13)
    )

outfile <- "C:/Users/win10/SpatialDataR/man/figures/fig6_multimodal.png"
png(outfile, width = 3200, height = 2600, res = 300, type = "windows")
print(p_final)
dev.off()
cat("Saved:", outfile, "\n")
cat("=== Task 6 complete ===\n")
