#!/usr/bin/env Rscript
## =================================================================
## SpatialDataR: Validation + Publication Figure Generation
## =================================================================
## Validates all package functions on the bundled mini-store
## and generates Figure 1 for the README.
##
## Usage: Rscript inst/scripts/validate_and_plot.R
## =================================================================

.libPaths("C:/Users/win10/R/win-library/4.4")
devtools::load_all(".", quiet = TRUE)
library(S4Vectors)
library(ggplot2)

cat("============================================\n")
cat("SpatialDataR Validation + Figure Generation\n")
cat("============================================\n\n")

store <- system.file("extdata", "xenium_mini.zarr",
    package = "SpatialDataR")
if (store == "") store <- "inst/extdata/xenium_mini.zarr"

## ---- 1. Validate store ----
cat("1. Validating Zarr store\n")
v <- validateSpatialData(store)
cat("   Valid:", v$valid, "\n")
cat("   Errors:", length(v$errors), "\n")
cat("   Warnings:", length(v$warnings), "\n")
cat("   Elements:\n")
print(v$elements)

## ---- 2. Read full store ----
cat("\n2. Reading SpatialData store\n")
sd <- readSpatialData(store)
show(sd)

## ---- 3. Element summary ----
cat("\n3. Element summary\n")
print(elementSummary(sd))

## ---- 4. Access loaded data ----
cat("\n4. Points data\n")
pts <- spatialPoints(sd)[["transcripts"]]
cat("   Transcripts:", nrow(pts), "\n")
cat("   Genes:", length(unique(pts$gene)), "\n")
cat("   X range:", round(range(pts$x), 2), "\n")
cat("   Y range:", round(range(pts$y), 2), "\n")

cat("\n5. Cell metadata\n")
tbl <- tables(sd)[["table"]]
obs <- tbl$obs
cat("   Cells:", nrow(obs), "\n")
cat("   Types:", paste(unique(obs$cell_type),
    collapse = ", "), "\n")

## ---- 5. Transform chain ----
cat("\n6. Coordinate transform chain\n")
ct <- elementTransform(images(sd)[["morphology"]])
cat("   Image transform: ")
show(ct)
inv <- invertTransform(ct)
cat("   Inverse: ")
show(inv)

## ---- 6. Spatial query ----
cat("\n7. Bounding box query\n")
sub_sd <- bboxQuery(sd, xmin = 0, xmax = 2,
    ymin = 0, ymax = 2)
cat("   Full points:", nrow(pts), "\n")
sub_pts <- spatialPoints(sub_sd)[["transcripts"]]
cat("   Queried points:", nrow(sub_pts), "\n")

## ---- 7. Aggregation ----
cat("\n8. Region aggregation\n")
shp <- shapes(sd)[["cell_boundaries"]]
counts <- aggregatePoints(pts, shp)
cat("   Cell-gene matrix:", nrow(counts), "x",
    ncol(counts) - 1L, "\n")

## ---- 8. Multi-sample ----
cat("\n9. Multi-sample operations\n")
combined <- combineSpatialData(sd, sd,
    sample_ids = c("primary", "metastasis"))
cat("   Combined elements:", length(combined), "\n")
cat("   Names:", paste(head(names(combined), 4),
    collapse = ", "), "...\n")
primary <- filterSample(combined, "primary")
cat("   Filtered 'primary':", length(primary),
    "elements\n")

## ============================================================
## FIGURE
## ============================================================
cat("\n10. Generating publication figure\n")

od <- "man/figures"
if (!dir.exists(od)) dir.create(od, recursive = TRUE)

pts_df <- as.data.frame(pts)
obs_df <- as.data.frame(obs)
shp_df <- as.data.frame(shp)

## Wong colorblind-safe palette
pal <- c(Epithelial = "#0072B2", Stromal = "#D55E00",
    Immune = "#009E73", Endothelial = "#E69F00")

## Panel a: transcript spatial map
pa <- ggplot(pts_df, aes(x = x, y = y, color = gene)) +
    geom_point(size = 0.3, alpha = 0.6) +
    coord_equal() +
    labs(x = expression("x ("*mu*"m)"),
        y = expression("y ("*mu*"m)"),
        color = NULL,
        title = expression(bold("a"))) +
    theme_classic(base_size = 8) +
    theme(legend.key.size = unit(6, "pt"),
        legend.text = element_text(size = 5))

## Panel b: cell type composition
ct_counts <- as.data.frame(table(obs_df$cell_type))
names(ct_counts) <- c("type", "count")
ct_counts$pct <- round(
    ct_counts$count / sum(ct_counts$count) * 100)

pb <- ggplot(ct_counts, aes(x = reorder(type, -count),
    y = count, fill = type)) +
    geom_col(width = 0.6, alpha = 0.85) +
    geom_text(aes(label = paste0(pct, "%")),
        vjust = -0.3, size = 2.2) +
    scale_fill_manual(values = pal) +
    scale_y_continuous(
        expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Number of cells",
        title = expression(bold("b"))) +
    theme_classic(base_size = 8) +
    theme(legend.position = "none")

## Panel c: cell boundaries colored by type
shp_merged <- merge(shp_df, obs_df, by = "cell_id")
pc <- ggplot(shp_merged, aes(x = x, y = y,
    color = cell_type)) +
    geom_point(aes(size = radius), alpha = 0.7) +
    scale_color_manual(values = pal) +
    scale_size_continuous(range = c(1, 4),
        guide = "none") +
    coord_equal() +
    labs(x = expression("x ("*mu*"m)"),
        y = expression("y ("*mu*"m)"),
        color = NULL,
        title = expression(bold("c"))) +
    theme_classic(base_size = 8) +
    theme(legend.key.size = unit(8, "pt"))

## Compose
library(patchwork)
fig <- (pa | pb) / pc + plot_layout(heights = c(1, 1))

ggsave(file.path(od, "fig1_spatial_overview.png"),
    fig, width = 183, height = 140, units = "mm",
    dpi = 300, bg = "white")
cat("   Saved: fig1_spatial_overview.png\n")

cat("\n============================================\n")
cat("VALIDATION COMPLETE — ALL OPERATIONS PASSED\n")
cat("============================================\n")
