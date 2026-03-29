## Task 1: Fix MERFISH cell_id and generate real aggregation heatmap
## Publication figure: man/figures/fig3_aggregation.png

.libPaths(c("C:/Users/win10/R/win-library/4.4",
            "C:/Users/win10/AppData/Local/R/win-library/4.4",
            .libPaths()))

library(SpatialDataR)
library(S4Vectors)
library(ggplot2)
library(viridis)
library(reshape2)

cat("=== Task 1: MERFISH cell_id reassignment and aggregation ===\n")

## ---- Step 1: Read transcripts and cell centroids ----
merfish_path <- "C:/Users/win10/merfish_spatialdata.zarr"
transcripts <- read.csv(file.path(merfish_path,
    "points/transcripts/transcripts.csv"),
    stringsAsFactors = FALSE)
cells <- read.csv(file.path(merfish_path,
    "shapes/cell_boundaries/circles.csv"),
    stringsAsFactors = FALSE)

cat("Transcripts:", nrow(transcripts), "rows,",
    length(unique(transcripts$gene)), "genes\n")
cat("Cells:", nrow(cells), "centroids\n")

## ---- Step 2: Assign transcripts to nearest cell (within 2x radius) ----
## For efficiency, process in chunks
cat("Assigning transcripts to cells (nearest centroid within 2x radius)...\n")

cell_x <- cells$x
cell_y <- cells$y
cell_r <- cells$radius
cell_ids <- cells$cell_id
n_cells <- length(cell_x)

## Process in chunks of 100k to manage memory
chunk_size <- 100000L
n_trans <- nrow(transcripts)
new_cell_id <- integer(n_trans)

t0 <- proc.time()
for (start in seq(1L, n_trans, by = chunk_size)) {
    end <- min(start + chunk_size - 1L, n_trans)
    tx <- transcripts$x[start:end]
    ty <- transcripts$y[start:end]
    n_chunk <- length(tx)

    ## Compute distance to each cell centroid
    ## Use matrix operations: n_chunk x n_cells distance matrix
    ## But 100k x 160 is manageable (~128MB for doubles)
    dx <- outer(tx, cell_x, "-")
    dy <- outer(ty, cell_y, "-")
    dist_mat <- sqrt(dx^2 + dy^2)

    ## Find nearest cell for each transcript
    min_idx <- max.col(-dist_mat, ties.method = "first")
    min_dist <- dist_mat[cbind(seq_len(n_chunk), min_idx)]

    ## Check within 2x radius
    max_dist <- 2 * cell_r[min_idx]
    assigned <- ifelse(min_dist <= max_dist, cell_ids[min_idx], 0L)
    new_cell_id[start:end] <- assigned

    if (start %% 500000 < chunk_size) {
        cat(sprintf("  Processed %d / %d (%.0f%%)\n",
            end, n_trans, 100 * end / n_trans))
    }
}
elapsed <- (proc.time() - t0)[3]
cat(sprintf("Assignment complete in %.1f seconds\n", elapsed))

## Stats
assigned_count <- sum(new_cell_id > 0)
cat(sprintf("Assigned: %d / %d transcripts (%.1f%%)\n",
    assigned_count, n_trans, 100 * assigned_count / n_trans))

## ---- Step 3: Rewrite transcripts CSV with corrected cell_id ----
transcripts$cell_id <- new_cell_id
write.csv(transcripts,
    file.path(merfish_path, "points/transcripts/transcripts.csv"),
    row.names = FALSE)
cat("Rewrote transcripts.csv with corrected cell_id\n")

## ---- Step 4: Re-read store and aggregate ----
cat("Reading updated store...\n")
sd <- readSpatialData(merfish_path)
pts <- spatialPoints(sd)[["transcripts"]]
regions <- shapes(sd)[["cell_boundaries"]]

cat("Running aggregatePoints()...\n")
counts <- aggregatePoints(pts, regions,
    feature_col = "gene", region_col = "cell_id")

## Remove cell_id column for matrix
cid <- counts$cell_id
counts$cell_id <- NULL
count_mat <- as.matrix(as.data.frame(counts))
rownames(count_mat) <- as.character(cid)

cat(sprintf("Count matrix: %d cells x %d genes\n",
    nrow(count_mat), ncol(count_mat)))
cat(sprintf("Total counts: %d\n", sum(count_mat)))
cat(sprintf("Median counts per cell: %.0f\n", median(rowSums(count_mat))))

## ---- Step 5: Heatmap of top 20 variable genes, grouped by layer ----
## Get layer info from real CSV
layer_data <- read.csv("C:/Users/win10/merfish_real.csv",
    stringsAsFactors = FALSE)

## Assign layers to cells based on majority vote of assigned transcripts
## (use the corrected transcripts)
trans_with_layer <- transcripts[transcripts$cell_id > 0, ]
## Match each transcript to its layer from layer_data (by index/row order)
## Actually, layer_data has different nrow. Match by x,y coordinates.
## Simpler: for each cell, find the nearest transcript in layer_data
## and use its layer. Or: for each cell, majority-vote the layer
## of transcripts assigned to it.

## layer_data has x_um, y_um matching transcripts x, y
## The original transcripts are in the same order as layer_data rows
## (layer_data has 3841412 rows vs transcripts 3714642 - different size)
## Better approach: for each cell centroid, look up the layer from layer_data
## using nearest point

## Use cell centroids to assign layers
cell_layers <- character(nrow(cells))
for (i in seq_len(nrow(cells))) {
    dx <- layer_data$x_um[1:min(100000, nrow(layer_data))] - cells$x[i]
    dy <- layer_data$y_um[1:min(100000, nrow(layer_data))] - cells$y[i]
    dists <- dx^2 + dy^2
    nearest_idx <- which.min(dists)
    cell_layers[i] <- layer_data$layer[nearest_idx]
}

## For cells not in first 100k, use depth-based assignment
## Actually, let's use ALL transcripts assigned to each cell to majority-vote
cat("Assigning cortical layers to cells...\n")

## Match transcripts to layer_data by position
## Since layer_data has x_um,y_um and transcripts has x,y with same coords:
## Build a lookup: for each assigned transcript, find its layer
## layer_data may have more rows. Use the transcript coordinates to match.
## Fastest: just use the cell centroid depth to assign layer
## cells depth = y coordinate approximately

## Actually simplest: use the cell centroids and find nearest point in layer_data
## Do it properly with a sample of layer_data for speed
set.seed(42)
sample_idx <- sample(nrow(layer_data), min(500000, nrow(layer_data)))
ld_sample <- layer_data[sample_idx, ]

cell_layers <- character(nrow(cells))
for (i in seq_len(nrow(cells))) {
    dx <- ld_sample$x_um - cells$x[i]
    dy <- ld_sample$y_um - cells$y[i]
    dists <- dx^2 + dy^2
    ## Get top 50 nearest
    top_idx <- order(dists)[1:min(50, length(dists))]
    ## Majority vote
    tbl <- table(ld_sample$layer[top_idx])
    cell_layers[i] <- names(which.max(tbl))
}

cat("Layer distribution:\n")
print(table(cell_layers))

## Select top 20 most variable genes
gene_var <- apply(count_mat, 2, var)
top20 <- names(sort(gene_var, decreasing = TRUE))[1:20]
mat_top <- count_mat[, top20]

## Log-transform for visualization
mat_log <- log1p(mat_top)

## Scale per gene (z-score columns)
mat_scaled <- scale(mat_log)
mat_scaled[mat_scaled > 3] <- 3
mat_scaled[mat_scaled < -3] <- -3

## Order cells by layer
layer_order <- c("VISp_I", "VISp_II/III", "VISp_IV", "VISp_V",
                 "VISp_VI", "VISp_wm", "VISp", "outside_VISp")
layer_fac <- factor(cell_layers, levels = layer_order)
cell_order <- order(layer_fac)

mat_ordered <- mat_scaled[cell_order, ]
layers_ordered <- cell_layers[cell_order]

## Melt for ggplot
df_hm <- melt(mat_ordered)
colnames(df_hm) <- c("Cell", "Gene", "Expression")
df_hm$Layer <- rep(layers_ordered, ncol(mat_ordered))
df_hm$Cell <- factor(df_hm$Cell, levels = unique(df_hm$Cell))
df_hm$Gene <- factor(df_hm$Gene, levels = rev(top20))

## Create layer annotation bar data
layer_colors <- c(
    "VISp_I" = "#E41A1C", "VISp_II/III" = "#377EB8",
    "VISp_IV" = "#4DAF4A", "VISp_V" = "#984EA3",
    "VISp_VI" = "#FF7F00", "VISp_wm" = "#A65628",
    "VISp" = "#F781BF", "outside_VISp" = "#999999"
)

## Build heatmap
p_hm <- ggplot(df_hm, aes(x = Cell, y = Gene, fill = Expression)) +
    geom_tile() +
    scale_fill_viridis(option = "inferno", name = "Scaled\nexpression") +
    theme_minimal(base_size = 10) +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 12)
    ) +
    labs(y = "Gene", title = "MERFISH: Top 20 variable genes across 160 cells")

## Layer annotation bar
df_layer <- data.frame(
    Cell = factor(seq_along(layers_ordered), levels = seq_along(layers_ordered)),
    Layer = layers_ordered,
    y = 1
)

p_layer <- ggplot(df_layer, aes(x = Cell, y = y, fill = Layer)) +
    geom_tile() +
    scale_fill_manual(values = layer_colors, name = "Cortical\nlayer") +
    theme_void() +
    theme(
        legend.position = "right",
        plot.margin = margin(0, 5, 0, 5)
    )

## Summary stats panel
stats_text <- sprintf(
    paste0("MERFISH Dataset Summary\n\n",
    "Transcripts: %s\n",
    "Assigned to cells: %s (%.1f%%)\n",
    "Cells: %d\n",
    "Genes: %d\n",
    "Total counts: %s\n",
    "Median counts/cell: %.0f"),
    format(n_trans, big.mark = ","),
    format(assigned_count, big.mark = ","),
    100 * assigned_count / n_trans,
    nrow(count_mat), ncol(count_mat),
    format(sum(count_mat), big.mark = ","),
    median(rowSums(count_mat))
)

## Use patchwork to combine
library(patchwork)

p_combined <- p_layer / p_hm + plot_layout(heights = c(1, 15))

## Save
outfile <- "C:/Users/win10/SpatialDataR/man/figures/fig3_aggregation.png"
dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
png(outfile, width = 2400, height = 1600, res = 300, type = "windows")
print(p_combined)
dev.off()

cat("Saved:", outfile, "\n")
cat("=== Task 1 complete ===\n")
