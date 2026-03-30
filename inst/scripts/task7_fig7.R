## Task 7: Downstream analysis with MERFISH real data
## Output: man/figures/fig7_downstream.png

.libPaths(c("C:/Users/win10/R/win-library/4.4",
            "C:/Users/win10/AppData/Local/R/win-library/4.4",
            .libPaths()))

library(SpatialDataR)
library(S4Vectors)
library(SingleCellExperiment)
library(ggplot2)
library(viridis)
library(patchwork)
library(scales)
library(reshape2)

cat("=== Task 7: Fig 7 — Downstream analysis (MERFISH) ===\n")

## ---- Step 1: Read MERFISH store ----
merfish_path <- "C:/Users/win10/merfish_spatialdata.zarr"
sd <- readSpatialData(merfish_path)
pts <- spatialPoints(sd)[["transcripts"]]
regions <- shapes(sd)[["cell_boundaries"]]

cat("Transcripts:", nrow(pts), "\n")
cat("Cells:", nrow(regions), "\n")

## ---- Step 2: Aggregate points to count matrix ----
cat("Running aggregatePoints()...\n")
counts <- aggregatePoints(pts, regions,
    feature_col = "gene", region_col = "cell_id")

cid <- counts$cell_id
counts$cell_id <- NULL
count_mat <- as.matrix(as.data.frame(counts))
rownames(count_mat) <- as.character(cid)

cat(sprintf("Count matrix: %d cells x %d genes\n",
    nrow(count_mat), ncol(count_mat)))

## ---- Step 3: Assign cortical layers to cells by majority vote ----
cat("Assigning cortical layers via majority vote...\n")
layer_data <- read.csv("C:/Users/win10/merfish_real.csv",
    stringsAsFactors = FALSE)

## Read transcripts with cell_id for majority vote
trans <- read.csv(file.path(merfish_path,
    "points/transcripts/transcripts.csv"),
    stringsAsFactors = FALSE)
trans_assigned <- trans[trans$cell_id > 0, ]

## For each transcript, find its layer from layer_data by matching x,y
## layer_data has x_um, y_um; trans has x, y — same coordinate system
## Since there are millions of rows, use cell centroid-based assignment
## For each cell, use the nearest 100 layer_data points to the centroid
cells_df <- as.data.frame(regions)

## Sample layer_data for speed
set.seed(42)
n_sample <- min(500000, nrow(layer_data))
ld <- layer_data[sample(nrow(layer_data), n_sample), ]

cell_layers <- character(nrow(cells_df))
for (i in seq_len(nrow(cells_df))) {
    dx <- ld$x_um - cells_df$x[i]
    dy <- ld$y_um - cells_df$y[i]
    dists <- dx^2 + dy^2
    top_idx <- order(dists)[1:min(100, length(dists))]
    tbl <- table(ld$layer[top_idx])
    cell_layers[i] <- names(which.max(tbl))
}
names(cell_layers) <- as.character(cells_df$cell_id)

cat("Layer distribution:\n")
print(table(cell_layers))

## ---- Step 4: Build SingleCellExperiment ----
layer_vec <- cell_layers[rownames(count_mat)]

sce <- SingleCellExperiment(
    assays = list(counts = t(count_mat)),
    colData = DataFrame(
        cell_id = rownames(count_mat),
        layer = layer_vec
    )
)

## Library-size normalise + log
lib_size <- colSums(counts(sce))
norm_mat <- t(t(counts(sce)) / lib_size) * median(lib_size)
logcounts(sce) <- log1p(norm_mat)

## ---- Step 5: PCA ----
cat("Running PCA...\n")
gene_var <- apply(logcounts(sce), 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(50, length(gene_var))]
pca_input <- t(logcounts(sce)[top_genes, ])
pca_res <- prcomp(pca_input, center = TRUE, scale. = TRUE)
var_exp <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100

cat(sprintf("PC1: %.1f%%, PC2: %.1f%%\n", var_exp[1], var_exp[2]))

## ---- Layer colours (consistent with fig1) ----
layer_cols <- c(
    "VISp_I"      = "#E41A1C",
    "VISp_II/III" = "#377EB8",
    "VISp_IV"     = "#4DAF4A",
    "VISp_V"      = "#984EA3",
    "VISp_VI"     = "#FF7F00",
    "VISp_wm"     = "#A65628",
    "VISp"        = "#F781BF",
    "outside_VISp"= "#999999"
)

## ---- Panel A: PCA coloured by cortical layer ----
pca_df <- data.frame(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    Layer = layer_vec
)
pca_df$Layer <- factor(pca_df$Layer,
    levels = c("VISp_I", "VISp_II/III", "VISp_IV", "VISp_V",
               "VISp_VI", "VISp_wm", "VISp", "outside_VISp"))

p_a <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = Layer)) +
    geom_point(size = 3.5, alpha = 0.85) +
    scale_colour_manual(values = layer_cols, name = "Cortical layer",
        drop = FALSE) +
    theme_classic(base_size = 11) +
    labs(x = sprintf("PC1 (%.1f%%)", var_exp[1]),
         y = sprintf("PC2 (%.1f%%)", var_exp[2])) +
    theme(
        legend.position = "right",
        legend.key.size = unit(4, "mm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 11),
        plot.margin = margin(10, 10, 10, 10)
    )

## ---- Panel B: Stacked bar — gene composition per layer ----
## Get top 8 genes by total expression
gene_totals <- colSums(count_mat)
top8 <- names(sort(gene_totals, decreasing = TRUE))[1:8]

## Build per-layer gene counts
layers_all <- layer_vec
bar_data <- data.frame(
    cell = rep(rownames(count_mat), each = 8),
    gene = rep(top8, nrow(count_mat)),
    count = as.vector(t(count_mat[, top8])),
    layer = rep(layers_all, each = 8)
)

## Aggregate by layer x gene
agg <- aggregate(count ~ layer + gene, data = bar_data, sum)

## Calculate proportions within each layer
layer_totals <- aggregate(count ~ layer, data = agg, sum)
names(layer_totals)[2] <- "total"
agg <- merge(agg, layer_totals, by = "layer")
agg$proportion <- agg$count / agg$total

agg$layer <- factor(agg$layer,
    levels = c("VISp_I", "VISp_II/III", "VISp_IV", "VISp_V",
               "VISp_VI", "VISp_wm", "VISp", "outside_VISp"))
agg$gene <- factor(agg$gene, levels = rev(top8))

## Nature-style gene palette
gene_cols <- c("#264653", "#2A9D8F", "#E9C46A", "#F4A261",
               "#E76F51", "#606C38", "#283618", "#BC6C25")
names(gene_cols) <- top8

p_b <- ggplot(agg, aes(x = layer, y = proportion, fill = gene)) +
    geom_col(width = 0.75, colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = gene_cols, name = "Gene") +
    scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
    theme_classic(base_size = 11) +
    labs(x = "Cortical layer", y = "Proportion of counts (top 8 genes)") +
    theme(
        axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
        legend.position = "right",
        legend.key.size = unit(4, "mm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 11),
        plot.margin = margin(10, 10, 10, 10)
    )

## ---- Assemble ----
p_final <- p_a + p_b +
    plot_annotation(
        tag_levels = "a",
        title = "SpatialDataR \u2192 SingleCellExperiment: MERFISH downstream analysis",
        subtitle = sprintf("160 cells \u00d7 268 genes | aggregatePoints() \u2192 PCA + layer composition",
            nrow(count_mat)),
        theme = theme(
            plot.title = element_text(face = "bold", size = 13),
            plot.subtitle = element_text(size = 10, colour = "grey40")
        )
    ) & theme(
        plot.tag = element_text(face = "bold", size = 13)
    )

## Save
outfile <- "C:/Users/win10/SpatialDataR/man/figures/fig7_downstream.png"
png(outfile, width = 3200, height = 1500, res = 300, type = "windows")
print(p_final)
dev.off()
cat("Saved:", outfile, "\n")
cat("=== Task 7 complete ===\n")
