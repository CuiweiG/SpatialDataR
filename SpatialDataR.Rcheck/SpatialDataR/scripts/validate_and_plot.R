.libPaths("C:/Users/win10/R/win-library/4.4")
devtools::load_all(".", quiet = TRUE)
library(S4Vectors)
library(ggplot2)

cat("============================================\n")
cat("SpatialDataR Validation on Xenium-mini Store\n")
cat("============================================\n\n")

store <- system.file("extdata", "xenium_mini.zarr",
                      package = "SpatialDataR")
if (store == "") store <- "inst/extdata/xenium_mini.zarr"

## ---- 1. Read full store ----
cat("1. Reading SpatialData store\n")
sd <- readSpatialData(store)
cat("   "); show(sd)

## ---- 2. Verify accessors ----
cat("\n2. Accessors\n")
cat("   images:", length(images(sd)), names(images(sd)), "\n")
cat("   labels:", length(spatialLabels(sd)), names(spatialLabels(sd)), "\n")
cat("   points:", length(spatialPoints(sd)), names(spatialPoints(sd)), "\n")
cat("   shapes:", length(shapes(sd)), names(shapes(sd)), "\n")
cat("   tables:", length(tables(sd)), names(tables(sd)), "\n")

## ---- 3. Coordinate systems ----
cat("\n3. Coordinate systems\n")
cs <- coordinateSystems(sd)
cat("   Found:", length(cs), "systems:",
    paste(names(cs), collapse = ", "), "\n")

## ---- 4. Read points (CSV fallback) ----
cat("\n4. Reading transcript points\n")
pts_path <- file.path(store, "points", "transcripts",
                       "transcripts.csv")
pts <- read.csv(pts_path, stringsAsFactors = FALSE)
pts_df <- DataFrame(pts)
cat("   Transcripts:", nrow(pts_df), "\n")
cat("   Genes:", length(unique(pts$gene)), "\n")
cat("   X range:", round(range(pts$x), 2), "\n")
cat("   Y range:", round(range(pts$y), 2), "\n")

## ---- 5. Read cell metadata ----
cat("\n5. Reading cell metadata\n")
obs <- read.csv(file.path(store, "tables", "table", "obs",
                           "obs.csv"), stringsAsFactors = FALSE)
cat("   Cells:", nrow(obs), "\n")
cat("   Cell types:", paste(unique(obs$cell_type), collapse = ", "), "\n")
cat("   Type distribution:\n")
print(table(obs$cell_type))

## ---- 6. Coordinate transform ----
cat("\n6. Coordinate transformation\n")
pixel_to_um <- CoordinateTransform("affine",
    affine = matrix(c(0.2125, 0, 0, 0, 0.2125, 0, 0, 0, 1),
                    nrow = 3, byrow = TRUE),
    input_cs = "pixels", output_cs = "global")

test_pts <- DataFrame(x = c(0, 100, 200), y = c(0, 100, 200))
transformed <- transformCoords(test_pts, pixel_to_um)
cat("   (0,0) px -> (", transformed$x[1], ",",
    transformed$y[1], ") um\n")
cat("   (100,100) px -> (", transformed$x[2], ",",
    transformed$y[2], ") um\n")
cat("   (200,200) px -> (", transformed$x[3], ",",
    transformed$y[3], ") um\n")

## ---- 7. Read shapes ----
cat("\n7. Reading cell shapes\n")
shp <- read.csv(file.path(store, "shapes", "cell_boundaries",
                           "circles.csv"), stringsAsFactors = FALSE)
cat("   Cell boundaries:", nrow(shp), "\n")
cat("   Mean radius:", round(mean(shp$radius), 3), "um\n")

## ============================================================
## FIGURES
## ============================================================
cat("\n8. Generating figures\n")

od <- "man/figures"
if (!dir.exists(od)) dir.create(od, recursive = TRUE)

## Wong palette
pal <- c(Epithelial = "#0072B2", Stromal = "#D55E00",
         Immune = "#009E73", Endothelial = "#E69F00")

## Panel A: Transcript spatial map colored by gene
pa <- ggplot(pts, aes(x = x, y = y, color = gene)) +
    geom_point(size = 0.3, alpha = 0.6) +
    scale_color_manual(values = c(
        EPCAM = "#0072B2", KRT18 = "#56B4E9", VIM = "#D55E00",
        CD45 = "#009E73", HER2 = "#E69F00", ESR1 = "#CC79A7",
        PGR = "#F0E442", MKI67 = "#000000", ERBB2 = "#999999",
        ACTB = "#0072B2")) +
    coord_equal() +
    labs(x = expression("x ("*mu*"m)"),
         y = expression("y ("*mu*"m)"),
         color = NULL,
         title = expression(bold("a"))) +
    theme_classic(base_size = 8, base_family = "sans") +
    theme(legend.position = "right",
          legend.key.size = unit(6, "pt"),
          legend.text = element_text(size = 5),
          panel.grid = element_blank())

## Panel B: Cell type composition
ct_counts <- as.data.frame(table(obs$cell_type))
names(ct_counts) <- c("type", "count")
ct_counts$pct <- round(ct_counts$count / sum(ct_counts$count) * 100)

pb <- ggplot(ct_counts, aes(x = reorder(type, -count),
                             y = count, fill = type)) +
    geom_col(width = 0.6, alpha = 0.85) +
    geom_text(aes(label = paste0(pct, "%")),
              vjust = -0.3, size = 2.2) +
    scale_fill_manual(values = pal) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Number of cells",
         title = expression(bold("b"))) +
    theme_classic(base_size = 8, base_family = "sans") +
    theme(legend.position = "none",
          panel.grid = element_blank())

## Panel C: Cell boundaries with type coloring
shp_merged <- merge(shp, obs, by = "cell_id")
pc <- ggplot(shp_merged, aes(x = x, y = y,
                              color = cell_type)) +
    geom_point(aes(size = radius), alpha = 0.7) +
    scale_color_manual(values = pal) +
    scale_size_continuous(range = c(1, 4), guide = "none") +
    coord_equal() +
    labs(x = expression("x ("*mu*"m)"),
         y = expression("y ("*mu*"m)"),
         color = NULL,
         title = expression(bold("c"))) +
    theme_classic(base_size = 8, base_family = "sans") +
    theme(legend.position = "right",
          legend.key.size = unit(8, "pt"),
          panel.grid = element_blank())

## Compose
library(patchwork)
fig <- (pa | pb) / pc + plot_layout(heights = c(1, 1))

ggsave(file.path(od, "fig1_spatial_overview.png"),
       fig, width = 183, height = 140, units = "mm",
       dpi = 300, bg = "white")
cat("   Saved: fig1_spatial_overview.png\n")

cat("\n============================================\n")
cat("VALIDATION COMPLETE\n")
cat("============================================\n")
