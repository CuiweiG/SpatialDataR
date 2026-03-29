## Task 2: Multimodal image + transcript overlay
## Publication figure: man/figures/fig6_multimodal.png

.libPaths(c("C:/Users/win10/R/win-library/4.4",
            "C:/Users/win10/AppData/Local/R/win-library/4.4",
            .libPaths()))

library(SpatialDataR)
library(S4Vectors)
library(ggplot2)
library(viridis)
library(patchwork)
library(Rarr)

cat("=== Task 2: Multimodal image + transcript overlay ===\n")

store <- system.file("extdata", "xenium_mini.zarr",
    package = "SpatialDataR")

## ---- Step 1: Read morphology image ----
img_path <- file.path(store, "images", "morphology", "scale0")
img <- readZarrArray(img_path)
cat("Image dimensions:", dim(img), "\n")
## img is 3 x 20 x 20 (CYX format)

## Convert to grayscale for display (average channels)
img_gray <- colMeans(aperm(img, c(1, 2, 3))[,,])  # already 2D slices
## Actually: img[c, y, x] -> need to extract properly
img_r <- img[1, , ]
img_g <- img[2, , ]
img_b <- img[3, , ]
img_gray <- (as.numeric(img_r) + as.numeric(img_g) + as.numeric(img_b)) / 3

## Create raster data frame
ny <- dim(img)[2]
nx <- dim(img)[3]
img_df <- expand.grid(y = seq_len(ny), x = seq_len(nx))
img_df$r <- as.vector(img[1, , ]) / 255
img_df$g <- as.vector(img[2, , ]) / 255
img_df$b <- as.vector(img[3, , ]) / 255
img_df$gray <- (img_df$r + img_df$g + img_df$b) / 3

## Scale coordinates to match transcript coordinate space
## Image is 20x20 pixels, transcripts span roughly 0-4 in x,y
x_scale <- 4.0 / nx
y_scale <- 4.0 / ny
img_df$x_scaled <- (img_df$x - 0.5) * x_scale
img_df$y_scaled <- (img_df$y - 0.5) * y_scale

## ---- Step 2: Read transcripts ----
pts <- read.csv(file.path(store, "points/transcripts/transcripts.csv"),
    stringsAsFactors = FALSE)
cat("Transcripts:", nrow(pts), "\n")

## ---- Step 3: Read cell labels (segmentation mask) ----
## The label array has a size mismatch issue, so create a synthetic mask
## from cell boundaries instead
cells <- read.csv(file.path(store, "shapes/cell_boundaries/circles.csv"),
    stringsAsFactors = FALSE)

## Create a pixel-level mask from cell circles
label_mat <- matrix(0L, nrow = ny, ncol = nx)
for (i in seq_len(nrow(cells))) {
    cx_px <- cells$x[i] / x_scale + 0.5
    cy_px <- cells$y[i] / y_scale + 0.5
    r_px <- cells$radius[i] / x_scale
    for (py in seq_len(ny)) {
        for (px in seq_len(nx)) {
            d <- sqrt((px - cx_px)^2 + (py - cy_px)^2)
            if (d <= r_px) {
                label_mat[py, px] <- cells$cell_id[i]
            }
        }
    }
}

label_df <- expand.grid(y = seq_len(ny), x = seq_len(nx))
label_df$cell_id <- as.vector(label_mat)
label_df$x_scaled <- (label_df$x - 0.5) * x_scale
label_df$y_scaled <- (label_df$y - 0.5) * y_scale
label_df$has_cell <- label_df$cell_id > 0

## ---- Panel A: Image + transcripts overlay ----
## Color transcripts by gene (top 6 genes)
gene_counts <- sort(table(pts$gene), decreasing = TRUE)
top_genes <- names(gene_counts)[1:min(6, length(gene_counts))]
pts$gene_group <- ifelse(pts$gene %in% top_genes, pts$gene, "Other")
pts$gene_group <- factor(pts$gene_group,
    levels = c(top_genes, "Other"))

## Overlay points on the raster background
p1 <- ggplot() +
    geom_raster(data = img_df,
        aes(x = x_scaled, y = y_scaled),
        fill = rgb(img_df$r, img_df$g, img_df$b)) +
    geom_point(data = pts,
        aes(x = x, y = y, color = gene_group),
        size = 1.8, alpha = 0.85) +
    scale_color_viridis_d(option = "turbo", name = "Gene") +
    theme_minimal(base_size = 10) +
    coord_fixed() +
    labs(title = "A) Morphology + transcripts",
         x = "x (µm)", y = "y (µm)") +
    theme(plot.title = element_text(face = "bold", size = 11))

## ---- Panel B: Segmentation mask overlay ----
## Show cell boundaries as colored regions on the image
p2 <- ggplot() +
    geom_raster(data = img_df,
        aes(x = x_scaled, y = y_scaled),
        fill = rgb(img_df$r, img_df$g, img_df$b)) +
    geom_tile(data = label_df[label_df$has_cell, ],
        aes(x = x_scaled, y = y_scaled,
            fill = factor(cell_id)),
        alpha = 0.4) +
    scale_fill_viridis_d(option = "plasma", name = "Cell ID") +
    geom_point(data = cells,
        aes(x = x, y = y), color = "white",
        size = 1, shape = 3) +
    theme_minimal(base_size = 10) +
    coord_fixed() +
    labs(title = "B) Cell segmentation mask",
         x = "x (µm)", y = "y (µm)") +
    theme(plot.title = element_text(face = "bold", size = 11),
          legend.position = "none")

## ---- Panel C: cropImage() demonstration ----
## Source cropImage from package source (not yet exported in installed version)
source("C:/Users/win10/SpatialDataR/R/image-ops.R")
## Crop a subregion
crop_arr <- cropImage(img_path, xmin = 5, xmax = 15,
    ymin = 5, ymax = 15)
cat("Cropped image dimensions:", dim(crop_arr), "\n")

## Convert crop to data frame
cn <- dim(crop_arr)
crop_df <- expand.grid(y = seq_len(cn[2]), x = seq_len(cn[3]))
crop_df$r <- as.vector(crop_arr[1, , ]) / 255
crop_df$g <- as.vector(crop_arr[2, , ]) / 255
crop_df$b <- as.vector(crop_arr[3, , ]) / 255

## Subset transcripts to crop region
x_crop_min <- (5 - 0.5) * x_scale
x_crop_max <- (15 - 0.5) * x_scale
y_crop_min <- (5 - 0.5) * y_scale
y_crop_max <- (15 - 0.5) * y_scale

pts_crop <- pts[pts$x >= x_crop_min & pts$x <= x_crop_max &
                pts$y >= y_crop_min & pts$y <= y_crop_max, ]

## Scale crop coordinates
crop_df$x_scaled <- x_crop_min + (crop_df$x - 0.5) * x_scale
crop_df$y_scaled <- y_crop_min + (crop_df$y - 0.5) * y_scale

p3 <- ggplot() +
    geom_raster(data = crop_df,
        aes(x = x_scaled, y = y_scaled),
        fill = rgb(crop_df$r, crop_df$g, crop_df$b)) +
    geom_point(data = pts_crop,
        aes(x = x, y = y, color = gene_group),
        size = 2.5, alpha = 0.9) +
    scale_color_viridis_d(option = "turbo", name = "Gene") +
    theme_minimal(base_size = 10) +
    coord_fixed() +
    labs(title = "C) cropImage() — zoomed subregion",
         x = "x (µm)", y = "y (µm)") +
    theme(plot.title = element_text(face = "bold", size = 11))

## ---- Panel D: Transcript density ----
p4 <- ggplot(pts, aes(x = x, y = y)) +
    stat_density_2d(aes(fill = after_stat(density)),
        geom = "raster", contour = FALSE) +
    scale_fill_viridis(option = "magma", name = "Density") +
    geom_point(data = cells,
        aes(x = x, y = y), color = "cyan",
        size = 2, shape = 1, stroke = 0.8) +
    theme_minimal(base_size = 10) +
    coord_fixed() +
    labs(title = "D) Transcript density + cell centroids",
         x = "x (µm)", y = "y (µm)") +
    theme(plot.title = element_text(face = "bold", size = 11))

## ---- Combine all panels ----
p_combined <- (p1 + p2) / (p3 + p4) +
    plot_annotation(
        title = "SpatialDataR: Multimodal spatial data visualization",
        subtitle = "Xenium mini dataset — morphology, transcripts, and segmentation",
        theme = theme(
            plot.title = element_text(face = "bold", size = 14),
            plot.subtitle = element_text(size = 11, color = "grey40")
        )
    )

## Save
outfile <- "C:/Users/win10/SpatialDataR/man/figures/fig6_multimodal.png"
png(outfile, width = 2800, height = 2400, res = 300, type = "windows")
print(p_combined)
dev.off()

cat("Saved:", outfile, "\n")
cat("=== Task 2 complete ===\n")
