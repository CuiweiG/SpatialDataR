#!/usr/bin/env Rscript
## Fig 5: writeSpatialData() roundtrip demonstration
## Uses REAL MERFISH data
.libPaths("C:/Users/win10/R/win-library/4.4")
library(SpatialDataR)
library(S4Vectors)
library(ggplot2)
library(patchwork)

cat("=== Roundtrip: read -> query -> write -> read ===\n")

## Step 1: Read real MERFISH store
store <- "C:/Users/win10/merfish_spatialdata.zarr"
sd <- readSpatialData(store)
pts <- spatialPoints(sd)[["transcripts"]]
cat("Original:", nrow(pts), "transcripts\n")

## Step 2: Spatial query — select cortex region
cx <- median(pts$x); cy <- median(pts$y)
sub <- bboxQuery(sd, xmin = cx - 300, xmax = cx + 300,
    ymin = cy - 300, ymax = cy + 300)
pts_sub <- spatialPoints(sub)[["transcripts"]]
cat("After bboxQuery:", nrow(pts_sub), "transcripts\n")

## Step 3: Write to new Zarr store
out_path <- "C:/Users/win10/merfish_subset.zarr"
if (dir.exists(out_path)) unlink(out_path, TRUE)
writeSpatialData(sub, out_path)
cat("Written to:", out_path, "\n")

## Step 4: Read back and validate
sd2 <- readSpatialData(out_path)
pts2 <- spatialPoints(sd2)[["transcripts"]]
cat("Read back:", nrow(pts2), "transcripts\n")
cat("Match:", nrow(pts_sub) == nrow(pts2), "\n")

v <- validateSpatialData(out_path)
cat("Validates:", v$valid, "\n\n")

## Step 5: Figure
od <- "man/figures"
pts_full <- as.data.frame(pts[sample(nrow(pts),
    min(50000, nrow(pts))), ])
pts_sub_df <- as.data.frame(pts_sub[sample(
    nrow(pts_sub), min(20000, nrow(pts_sub))), ])
pts_back <- as.data.frame(pts2[sample(nrow(pts2),
    min(20000, nrow(pts2))), ])

th <- theme_classic(base_size = 9) +
    theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8,
            color = "black"),
        plot.margin = margin(6, 8, 6, 6))

## Panel a: full dataset
pa <- ggplot(pts_full, aes(x = x, y = y)) +
    geom_point(size = 0.02, alpha = 0.15,
        color = "grey40") +
    annotate("rect",
        xmin = cx - 300, xmax = cx + 300,
        ymin = cy - 300, ymax = cy + 300,
        fill = NA, color = "#D55E00",
        linewidth = 0.6, linetype = "dashed") +
    coord_equal() +
    labs(subtitle = paste0("readSpatialData()\n",
        format(nrow(pts), big.mark = ","),
        " transcripts"),
        x = expression(mu*"m"),
        y = expression(mu*"m")) +
    th + theme(
        plot.subtitle = element_text(size = 7,
            face = "bold", lineheight = 1.1))

## Panel b: after query
pb <- ggplot(pts_sub_df, aes(x = x, y = y)) +
    geom_point(size = 0.05, alpha = 0.3,
        color = "#0072B2") +
    coord_equal() +
    labs(subtitle = paste0("bboxQuery()\n",
        format(nrow(pts_sub), big.mark = ","),
        " selected"),
        x = expression(mu*"m"), y = NULL) +
    th + theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.subtitle = element_text(size = 7,
            face = "bold", color = "#0072B2",
            lineheight = 1.1))

## Panel c: arrow showing write
pc <- ggplot() +
    annotate("text", x = 0.5, y = 0.65,
        label = "writeSpatialData()",
        size = 3.5, fontface = "bold",
        color = "#009E73") +
    annotate("segment", x = 0.2, xend = 0.8,
        y = 0.4, yend = 0.4,
        arrow = arrow(length = unit(6, "pt"),
            type = "closed"),
        linewidth = 0.8, color = "#009E73") +
    annotate("text", x = 0.5, y = 0.2,
        label = ".zarr", size = 3,
        color = "grey40", fontface = "italic") +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(plot.margin = margin(6, 2, 6, 2))

## Panel d: read back
pd <- ggplot(pts_back, aes(x = x, y = y)) +
    geom_point(size = 0.05, alpha = 0.3,
        color = "#D55E00") +
    coord_equal() +
    labs(subtitle = paste0("readSpatialData()\n",
        format(nrow(pts2), big.mark = ","),
        " verified"),
        x = expression(mu*"m"), y = NULL) +
    th + theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.subtitle = element_text(size = 7,
            face = "bold", color = "#D55E00",
            lineheight = 1.1))

fig5 <- (pa + labs(tag = "a")) +
    (pb + labs(tag = "b")) +
    (pc + labs(tag = "")) +
    (pd + labs(tag = "c")) +
    plot_layout(ncol = 4,
        widths = c(1, 0.8, 0.3, 0.8)) +
    plot_annotation(
        title = "Read \u2192 Query \u2192 Write \u2192 Verify",
        subtitle = paste0(
            "Full roundtrip on MERFISH data ",
            "(Moffitt 2018). ",
            "writeSpatialData() produces ",
            "spec-compliant Zarr."),
        theme = theme(
            plot.title = element_text(size = 12,
                face = "bold"),
            plot.subtitle = element_text(size = 8,
                color = "grey30",
                face = "italic"))) &
    theme(plot.tag = element_text(size = 10,
        face = "bold"))

ggsave(file.path(od, "fig5_roundtrip.png"), fig5,
    width = 220, height = 80, units = "mm",
    dpi = 300, bg = "white")

cat("Saved:", file.path(od, "fig5_roundtrip.png"),
    "\n")
cat("Size:", round(file.size(
    file.path(od, "fig5_roundtrip.png")) / 1024),
    "KB\n")
