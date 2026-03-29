## Task 3: Downstream analysis integration
## Publication figure: man/figures/fig7_downstream.png

.libPaths(c("C:/Users/win10/R/win-library/4.4",
            "C:/Users/win10/AppData/Local/R/win-library/4.4",
            .libPaths()))

library(SpatialDataR)
library(S4Vectors)
library(SingleCellExperiment)
library(ggplot2)
library(viridis)
library(patchwork)

cat("=== Task 3: Downstream analysis integration ===\n")

## ---- Step 1: Load xenium_mini.zarr ----
store <- system.file("extdata", "xenium_mini.zarr",
    package = "SpatialDataR")
sd <- readSpatialData(store)
cat("Loaded SpatialData object\n")

## ---- Step 2: Aggregate points to count matrix ----
pts <- spatialPoints(sd)[["transcripts"]]
regions <- shapes(sd)[["cell_boundaries"]]
counts <- aggregatePoints(pts, regions,
    feature_col = "gene", region_col = "cell_id")

cid <- counts$cell_id
counts$cell_id <- NULL
count_mat <- as.matrix(as.data.frame(counts))
rownames(count_mat) <- as.character(cid)

cat(sprintf("Count matrix: %d cells x %d genes\n",
    nrow(count_mat), ncol(count_mat)))

## ---- Step 3: Convert to SingleCellExperiment ----
## Get cell_type from obs table
obs <- tables(sd)[["table"]]$obs
obs_df <- as.data.frame(obs)

## Match cell types to count matrix rows
ct_map <- setNames(obs_df$cell_type, as.character(obs_df$cell_id))
cell_types <- ct_map[rownames(count_mat)]
cell_types[is.na(cell_types)] <- "Unknown"

## Get spatial coordinates from shapes
coord_map_x <- setNames(as.data.frame(regions)$x,
    as.character(as.data.frame(regions)$cell_id))
coord_map_y <- setNames(as.data.frame(regions)$y,
    as.character(as.data.frame(regions)$cell_id))

spatial_x <- coord_map_x[rownames(count_mat)]
spatial_y <- coord_map_y[rownames(count_mat)]

## Build SCE
sce <- SingleCellExperiment(
    assays = list(counts = t(count_mat)),
    colData = DataFrame(
        cell_id = rownames(count_mat),
        cell_type = cell_types,
        x = spatial_x,
        y = spatial_y
    )
)

cat(sprintf("SingleCellExperiment: %d genes x %d cells\n",
    nrow(sce), ncol(sce)))

## ---- Step 4: Basic downstream analysis ----
## Library size normalization
lib_size <- colSums(counts(sce))
norm_mat <- t(t(counts(sce)) / lib_size) * median(lib_size)
logcounts(sce) <- log1p(norm_mat)

## PCA using base R prcomp
cat("Running PCA...\n")
## Use top variable genes
gene_var <- apply(logcounts(sce), 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[
    1:min(20, length(gene_var))]
pca_input <- t(logcounts(sce)[top_genes, ])
pca_res <- prcomp(pca_input, center = TRUE, scale. = TRUE)

## Store PCA results
reducedDim(sce, "PCA") <- pca_res$x[, 1:min(5, ncol(pca_res$x))]

## Variance explained
var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100

cat(sprintf("PC1: %.1f%% variance, PC2: %.1f%% variance\n",
    var_explained[1], var_explained[2]))

## ---- Step 5: Generate 2-panel figure ----
pca_df <- data.frame(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    cell_type = colData(sce)$cell_type
)

## Panel A: PCA coloured by cell type
ct_colors <- c(
    "Epithelial" = "#E41A1C",
    "Stromal" = "#377EB8",
    "Endothelial" = "#4DAF4A",
    "Immune" = "#FF7F00",
    "Unknown" = "#999999"
)
## Use available colors for present types
present_types <- unique(pca_df$cell_type)
if (!all(present_types %in% names(ct_colors))) {
    extra <- setdiff(present_types, names(ct_colors))
    extra_cols <- viridis(length(extra))
    names(extra_cols) <- extra
    ct_colors <- c(ct_colors, extra_cols)
}

p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cell_type)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = ct_colors, name = "Cell type") +
    theme_minimal(base_size = 12) +
    labs(
        title = "A) PCA of aggregated counts",
        x = sprintf("PC1 (%.1f%%)", var_explained[1]),
        y = sprintf("PC2 (%.1f%%)", var_explained[2])
    ) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "right"
    )

## Panel B: Cell type proportions bar chart
ct_counts <- as.data.frame(table(cell_type = pca_df$cell_type))
ct_counts <- ct_counts[order(-ct_counts$Freq), ]
ct_counts$cell_type <- factor(ct_counts$cell_type,
    levels = ct_counts$cell_type)

p2 <- ggplot(ct_counts, aes(x = cell_type, y = Freq, fill = cell_type)) +
    geom_col(alpha = 0.85) +
    scale_fill_manual(values = ct_colors, guide = "none") +
    theme_minimal(base_size = 12) +
    labs(
        title = "B) Cell type composition",
        x = "Cell type",
        y = "Number of cells"
    ) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        axis.text.x = element_text(angle = 30, hjust = 1)
    )

## Combine
p_combined <- p1 + p2 +
    plot_annotation(
        title = "SpatialDataR → SingleCellExperiment downstream analysis",
        subtitle = "Xenium mini: aggregatePoints() → normalization → PCA",
        theme = theme(
            plot.title = element_text(face = "bold", size = 14),
            plot.subtitle = element_text(size = 11, color = "grey40")
        )
    )

## Save
outfile <- "C:/Users/win10/SpatialDataR/man/figures/fig7_downstream.png"
dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
png(outfile, width = 2800, height = 1400, res = 300, type = "windows")
print(p_combined)
dev.off()

cat("Saved:", outfile, "\n")
cat("=== Task 3 complete ===\n")
