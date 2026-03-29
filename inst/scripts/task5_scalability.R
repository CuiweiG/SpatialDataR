## Task 5: Scalability architecture demonstration
## Publication figure: man/figures/fig9_scalability.png

.libPaths(c("C:/Users/win10/R/win-library/4.4",
            "C:/Users/win10/AppData/Local/R/win-library/4.4",
            .libPaths()))

library(SpatialDataR)
library(S4Vectors)
library(DelayedArray)
library(Rarr)
library(ggplot2)
library(patchwork)
library(viridis)

cat("=== Task 5: Scalability architecture demonstration ===\n")

## ---- Step 1: Show lazy path references ----
store <- system.file("extdata", "xenium_mini.zarr",
    package = "SpatialDataR")
sd <- readSpatialData(store)

## Images and labels are stored as path references, not loaded
img_ref <- images(sd)[["morphology"]]
lbl_ref <- spatialLabels(sd)[["cell_labels"]]

cat("Image element class:", class(img_ref), "\n")
cat("Image is path reference:", is.list(img_ref) && "path" %in% names(img_ref), "\n")
cat("Image path:", img_ref$path, "\n")
cat("Label element class:", class(lbl_ref), "\n")
cat("Label is path reference:", is.list(lbl_ref) && "path" %in% names(lbl_ref), "\n")

## Memory footprint of the SpatialData object
sd_size <- as.numeric(object.size(sd))
cat(sprintf("\nSpatialData object size: %s bytes (%.1f KB)\n",
    format(sd_size, big.mark = ","), sd_size / 1024))

## ---- Step 2: DelayedArray demonstration ----
img_path <- file.path(store, "images", "morphology", "scale0")
da <- readZarrDelayed(img_path)
cat("\nDelayedArray class:", class(da), "\n")
cat("DelayedArray dims:", dim(da), "\n")
cat("Is DelayedArray:", is(da, "DelayedArray"), "\n")

## Memory comparison: DelayedArray seed vs realized array
da_size <- as.numeric(object.size(da))
realized <- as.array(da)
realized_size <- as.numeric(object.size(realized))
cat(sprintf("DelayedArray object: %s bytes\n",
    format(da_size, big.mark = ",")))
cat(sprintf("Realized array: %s bytes\n",
    format(realized_size, big.mark = ",")))

## ---- Step 3: MERFISH scalability ----
merfish_path <- "C:/Users/win10/merfish_spatialdata.zarr"

## Time reading the full store
cat("\nReading MERFISH store...\n")
t0 <- proc.time()
sd_merfish <- readSpatialData(merfish_path)
t_read <- (proc.time() - t0)[3]
cat(sprintf("Full store read time: %.2f seconds\n", t_read))

merfish_size <- as.numeric(object.size(sd_merfish))
cat(sprintf("MERFISH SpatialData object size: %s bytes (%.1f MB)\n",
    format(merfish_size, big.mark = ","),
    merfish_size / 1024^2))

## Count total transcripts
pts_merfish <- spatialPoints(sd_merfish)[["transcripts"]]
n_total <- nrow(pts_merfish)
cat(sprintf("Total transcripts: %s\n", format(n_total, big.mark = ",")))

## ---- Step 4: Chunked bboxQuery benchmark ----
## Define several bbox sizes to demonstrate scaling
cat("\nBounding box query benchmarks:\n")

## Get coordinate ranges
x_range <- range(pts_merfish$x)
y_range <- range(pts_merfish$y)
x_mid <- mean(x_range)
y_mid <- mean(y_range)
x_span <- diff(x_range)
y_span <- diff(y_range)

## Different bbox sizes: 10%, 25%, 50%, 75%, 100%
fractions <- c(0.1, 0.25, 0.5, 0.75, 1.0)
bbox_results <- data.frame(
    fraction = numeric(),
    n_points = integer(),
    time_sec = numeric()
)

for (frac in fractions) {
    half_x <- x_span * frac / 2
    half_y <- y_span * frac / 2
    xmin <- x_mid - half_x
    xmax <- x_mid + half_x
    ymin <- y_mid - half_y
    ymax <- y_mid + half_y

    t0 <- proc.time()
    sub <- bboxQuery(sd_merfish,
        xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax)
    t_query <- (proc.time() - t0)[3]

    sub_pts <- spatialPoints(sub)[["transcripts"]]
    n_sub <- nrow(sub_pts)

    cat(sprintf("  %3.0f%% area: %s points in %.3f sec\n",
        frac * 100, format(n_sub, big.mark = ","), t_query))

    bbox_results <- rbind(bbox_results, data.frame(
        fraction = frac,
        n_points = n_sub,
        time_sec = t_query
    ))
}

## ---- Step 5: Generate figure ----

## Panel A: Memory comparison bar chart
mem_data <- data.frame(
    Object = c(
        "SpatialData\n(xenium mini)",
        "DelayedArray\n(image ref)",
        "Realized array\n(in memory)",
        "SpatialData\n(MERFISH 3.7M)"
    ),
    Size_KB = c(
        sd_size / 1024,
        da_size / 1024,
        realized_size / 1024,
        merfish_size / 1024
    ),
    Type = c("Lazy", "Lazy", "In-memory", "Lazy")
)
mem_data$Object <- factor(mem_data$Object,
    levels = mem_data$Object)

p1 <- ggplot(mem_data, aes(x = Object, y = Size_KB, fill = Type)) +
    geom_col(alpha = 0.85, width = 0.7) +
    geom_text(aes(label = sprintf("%.0f KB", Size_KB)),
        vjust = -0.3, size = 3.2) +
    scale_fill_manual(values = c("Lazy" = "#4DAF4A", "In-memory" = "#E41A1C")) +
    scale_y_log10(labels = function(x) paste0(round(x), " KB")) +
    theme_minimal(base_size = 11) +
    labs(
        title = "A) Memory footprint: lazy vs in-memory",
        x = "", y = "Size (KB, log scale)"
    ) +
    theme(
        plot.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(size = 9),
        legend.position = "top"
    )

## Panel B: bboxQuery timing benchmark
bbox_results$pct_label <- paste0(bbox_results$fraction * 100, "%")

p2 <- ggplot(bbox_results,
    aes(x = fraction * 100, y = time_sec)) +
    geom_line(color = "#377EB8", linewidth = 1.2) +
    geom_point(aes(size = n_points / 1000),
        color = "#377EB8", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.0fms\n%sk pts",
        time_sec * 1000,
        format(round(n_points / 1000), big.mark = ","))),
        vjust = -0.8, size = 2.8, color = "grey30") +
    scale_size_continuous(name = "Points (×1000)",
        range = c(3, 8)) +
    theme_minimal(base_size = 11) +
    labs(
        title = "B) bboxQuery() scaling on MERFISH (3.7M transcripts)",
        x = "Bounding box area (%)",
        y = "Query time (seconds)"
    ) +
    theme(
        plot.title = element_text(face = "bold", size = 12),
        legend.position = "right"
    )

## Panel C: Architecture summary
arch_lines <- c(
    "SpatialDataR Scalability Architecture",
    "",
    "1. Lazy loading: images/labels stored as path refs",
    paste0("   → SpatialData object: ",
        format(round(sd_size / 1024), big.mark = ","), " KB"),
    paste0("   → Full MERFISH: ",
        format(round(merfish_size / 1024), big.mark = ","), " KB"),
    "",
    "2. DelayedArray: out-of-memory array access via Rarr",
    paste0("   → Seed object: ", round(da_size / 1024), " KB vs ",
        "realized: ", round(realized_size / 1024), " KB"),
    "",
    "3. Vectorised spatial queries: bboxQuery() on DataFrames",
    paste0("   → 3.7M transcripts queried in <",
        sprintf("%.0f", max(bbox_results$time_sec) * 1000),
        " ms (full dataset)"),
    paste0("   → 10% bbox: ",
        sprintf("%.0f", bbox_results$time_sec[1] * 1000),
        " ms for ",
        format(bbox_results$n_points[1], big.mark = ","),
        " points")
)

arch_df <- data.frame(
    x = 0,
    y = rev(seq_along(arch_lines)),
    label = arch_lines
)

p3 <- ggplot(arch_df, aes(x = x, y = y, label = label)) +
    geom_text(hjust = 0, family = "mono", size = 3, color = "grey20") +
    xlim(-0.5, 12) +
    theme_void() +
    labs(title = "C) Architecture summary") +
    theme(plot.title = element_text(face = "bold", size = 12,
        margin = margin(b = 10)))

## Combine
p_combined <- (p1 | p2) / p3 +
    plot_layout(heights = c(3, 2)) +
    plot_annotation(
        title = "SpatialDataR: Scalability and out-of-memory architecture",
        subtitle = paste0("MERFISH: ", format(n_total, big.mark = ","),
            " transcripts | Lazy image refs | DelayedArray integration"),
        theme = theme(
            plot.title = element_text(face = "bold", size = 14),
            plot.subtitle = element_text(size = 11, color = "grey40")
        )
    )

## Save
outfile <- "C:/Users/win10/SpatialDataR/man/figures/fig9_scalability.png"
dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
png(outfile, width = 3000, height = 2200, res = 300, type = "windows")
print(p_combined)
dev.off()

cat("\nSaved:", outfile, "\n")
cat("=== Task 5 complete ===\n")
