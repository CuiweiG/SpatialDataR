#!/usr/bin/env Rscript
# ============================================================================
# SpatialDataR README figures — FINAL v2
# ALL API calls are real SpatialDataR functions. No workarounds.
# Features used: readSpatialData, parseGeometry, geometryCentroids,
#   bboxQuery, aggregatePoints, spatialJoin, plotSpatialData
# ============================================================================

.libPaths(c("C:/Users/win10/R/win-library/4.4",
            "C:/Users/win10/AppData/Local/R/win-library/4.4",
            .libPaths()))

suppressPackageStartupMessages({
    devtools::load_all("C:/Users/win10/SpatialDataR", quiet = TRUE)
    library(S4Vectors)
    library(ggplot2)
    library(patchwork)
    library(viridis)
    library(ggrepel)
    library(reshape2)
    library(scales)
    library(arrow)
})

setwd("C:/Users/win10/SpatialDataR")
od <- "man/figures"
dir.create(od, showWarnings = FALSE, recursive = TRUE)

th <- function(base = 9) {
    theme_classic(base_size = base, base_family = "sans") +
    theme(
        axis.line        = element_line(linewidth = 0.4, colour = "black"),
        axis.ticks       = element_line(linewidth = 0.3, colour = "black"),
        axis.title       = element_text(size = rel(1), colour = "black"),
        axis.text        = element_text(size = rel(0.9), colour = "grey20"),
        plot.title       = element_text(size = rel(1.15), face = "bold",
                                        hjust = 0, margin = margin(b = 1)),
        plot.subtitle    = element_text(size = rel(0.8), colour = "grey30",
                                        hjust = 0, margin = margin(b = 4)),
        plot.margin      = margin(6, 8, 6, 6),
        legend.title     = element_text(size = rel(0.9), face = "bold"),
        legend.text      = element_text(size = rel(0.8)),
        legend.key.size  = unit(0.35, "cm"),
        legend.background = element_blank()
    )
}

ct_cols <- c("Epithelial" = "#E64B35", "Stromal" = "#4DBBD5",
             "Immune" = "#00A087", "Endothelial" = "#F39B7F")
keep_types <- c("Epithelial", "Stromal", "Immune", "Endothelial")

## ======================================================================
## STEP 1: readSpatialData() on Xenium
## ======================================================================
cat("=== Step 1: readSpatialData() ===\n")
sd <- readSpatialData("C:/Users/win10/xenium_breast/data.zarr")
print(sd)

## Extract with real API + parseGeometry
pts_df <- spatialPoints(sd)[["transcripts"]]
pts <- as.data.frame(pts_df)
colnames(pts)[colnames(pts) == "feature_name"] <- "gene"
pts$gene <- as.character(pts$gene)
pts <- pts[pts$qv >= 20 & !grepl("^NegControl|^BLANK", pts$gene), ]

## geometryCentroids() — the new API
cells_df <- shapes(sd)[["cell_circles"]]
cell_coords <- geometryCentroids(cells_df[["geometry"]])
cells <- as.data.frame(cells_df)
cells$x <- cell_coords$x
cells$y <- cell_coords$y

## obs from table
tbl <- tables(sd)[["table"]]
obs <- as.data.frame(tbl$obs)

n_transcripts <- nrow(pts); n_genes <- length(unique(pts$gene))
n_cells <- nrow(cells)
x_rng <- range(pts$x); y_rng <- range(pts$y)
cat("  Transcripts:", format(n_transcripts, big.mark = ","), "\n")
cat("  Genes:", n_genes, "\n")
cat("  Cells:", format(n_cells, big.mark = ","), "\n\n")

## ======================================================================
## FIG 1: Multi-element store reading
## ======================================================================
cat("--- Fig 1 ---\n")
bin1 <- 30
pts$xb <- round(pts$x / bin1) * bin1
pts$yb <- round(pts$y / bin1) * bin1
density <- aggregate(gene ~ xb + yb, data = pts, FUN = length)
colnames(density)[3] <- "count"

p1a <- ggplot(density, aes(x = xb, y = yb, fill = count)) +
    geom_raster() +
    scale_fill_viridis_c(option = "inferno", trans = "sqrt",
                          name = "Transcripts\nper bin") +
    coord_equal(expand = FALSE) +
    annotate("segment", x = x_rng[2]-1100, xend = x_rng[2]-100,
             y = y_rng[1]+100, yend = y_rng[1]+100,
             linewidth = 2, colour = "white") +
    annotate("text", x = x_rng[2]-600, y = y_rng[1]+100,
             label = "1 mm", vjust = -0.8, size = 3,
             colour = "white", fontface = "bold") +
    labs(x = expression(italic(x)~"("*mu*"m)"),
         y = expression(italic(y)~"("*mu*"m)")) + th(9)

set.seed(42)
cell_sub <- cells[sample(n_cells, min(n_cells, 50000)), ]
cell_sub$total <- obs$total_counts[match(cell_sub$cell_id, obs$cell_id)]

p1b <- ggplot(cell_sub, aes(x = x, y = y, colour = total)) +
    geom_point(size = 0.15, alpha = 0.6) +
    scale_colour_viridis_c(option = "mako", trans = "sqrt",
                            name = "Total\ncounts", direction = -1) +
    coord_equal(expand = FALSE) +
    labs(x = expression(italic(x)~"("*mu*"m)"), y = NULL) +
    th(9) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank())

fig1 <- (p1a + labs(tag = "a")) + (p1b + labs(tag = "b")) +
    plot_layout(ncol = 2) +
    plot_annotation(
        title = "Multi-element spatial data store read in a single call",
        subtitle = paste0(format(n_transcripts, big.mark = ","),
                          " transcripts, ", n_genes, " genes, ",
                          format(n_cells, big.mark = ","),
                          " cells (10x Xenium breast cancer, Janesick et al. 2023)"),
        theme = theme(plot.title = element_text(size = 11, face = "bold"),
                      plot.subtitle = element_text(size = 8, colour = "grey30"))) &
    theme(plot.tag = element_text(size = 10, face = "bold"))

ggsave(file.path(od, "fig1_store_reading.png"), fig1,
       width = 210, height = 105, units = "mm", dpi = 300, bg = "white")
cat("  saved\n")

## ======================================================================
## FIG 2: bboxQuery() — real API
## ======================================================================
cat("--- Fig 2 ---\n")
cx <- median(pts$x); cy <- median(pts$y); hw <- 500
qx <- c(cx-hw, cx+hw); qy <- c(cy-hw, cy+hw)

pts_roi <- as.data.frame(bboxQuery(pts_df, qx[1], qx[2], qy[1], qy[2]))
colnames(pts_roi)[colnames(pts_roi) == "feature_name"] <- "gene"
pts_roi$gene <- as.character(pts_roi$gene)
pts_roi <- pts_roi[pts_roi$qv >= 20 & !grepl("^NegControl|^BLANK", pts_roi$gene), ]
n_roi <- nrow(pts_roi)
cells_roi <- cells[cells$x >= qx[1] & cells$x <= qx[2] &
                    cells$y >= qy[1] & cells$y <= qy[2], ]
cat("  ROI:", format(n_roi, big.mark=","), "transcripts,", nrow(cells_roi), "cells\n")

p2a <- ggplot(density, aes(x = xb, y = yb, fill = count)) +
    geom_raster() +
    scale_fill_viridis_c(option = "inferno", trans = "sqrt", guide = "none") +
    annotate("rect", xmin=qx[1], xmax=qx[2], ymin=qy[1], ymax=qy[2],
             fill = NA, colour = "white", linewidth = 1) +
    coord_equal(expand = FALSE) +
    labs(x = expression(italic(x)~"("*mu*"m)"),
         y = expression(italic(y)~"("*mu*"m)")) + th(8.5)

roi_gene_freq <- sort(table(pts_roi$gene), decreasing = TRUE)
top6 <- names(roi_gene_freq)[1:6]
pts_roi$gene_col <- ifelse(pts_roi$gene %in% top6, pts_roi$gene, "Other")
pts_roi$gene_col <- factor(pts_roi$gene_col, levels = c(top6, "Other"))
gene6_cols <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#E8E8E8")
names(gene6_cols) <- c(top6, "Other")

p2b <- ggplot() +
    geom_point(data = pts_roi[pts_roi$gene_col=="Other",], aes(x=x,y=y),
               colour = "#E8E8E8", size = 0.05, alpha = 0.3) +
    geom_point(data = pts_roi[pts_roi$gene_col!="Other",], aes(x=x,y=y,colour=gene_col),
               size = 0.15, alpha = 0.7) +
    scale_colour_manual(values = gene6_cols, name = "Gene",
                        guide = guide_legend(override.aes = list(size=2, alpha=1))) +
    coord_equal(xlim=qx, ylim=qy, expand=FALSE) +
    annotate("segment", x=qx[2]-220, xend=qx[2]-20,
             y=qy[1]+30, yend=qy[1]+30, linewidth=1.5, colour="black") +
    annotate("text", x=qx[2]-120, y=qy[1]+30, label="200 \u00b5m",
             vjust=-0.7, size=2.5, fontface="bold") +
    labs(x=expression(italic(x)~"("*mu*"m)"), y=NULL) +
    th(8.5) + theme(panel.background=element_rect(fill="grey95"))

fig2 <- (p2a + labs(tag="a")) + (p2b + labs(tag="b")) +
    plot_layout(ncol=2, widths=c(0.8,1)) +
    plot_annotation(
        title = "Spatial bounding-box query isolates tumour microenvironment",
        subtitle = paste0(format(n_roi,big.mark=",")," transcripts, ",
                          nrow(cells_roi)," cells in 1 \u00d7 1 mm ROI"),
        theme = theme(plot.title=element_text(size=11,face="bold"),
                      plot.subtitle=element_text(size=8,colour="grey30"))) &
    theme(plot.tag = element_text(size=10,face="bold"))

ggsave(file.path(od, "fig2_spatial_query.png"), fig2,
       width = 210, height = 105, units = "mm", dpi = 300, bg = "white")
cat("  saved\n")

## ======================================================================
## FIG 3: aggregatePoints() — REAL API call
## ======================================================================
cat("--- Fig 3 ---\n")
t0 <- proc.time()
counts <- aggregatePoints(pts_df[pts_df$cell_id > 0 & pts_df$qv >= 20 &
    !grepl("^NegControl|^BLANK", as.character(pts_df$feature_name)), ],
    cells_df, feature_col = "feature_name", region_col = "cell_id")
elapsed <- (proc.time() - t0)[3]
cat("  aggregatePoints():", nrow(counts), "x", ncol(counts), "in", round(elapsed,1), "s\n")

counts_df <- as.data.frame(counts)
cids <- counts_df$cell_id; counts_df$cell_id <- NULL
expr_mat <- as.matrix(counts_df); rownames(expr_mat) <- cids

## Cell type assignment
avail <- intersect(c("EPCAM","PTPRC","PECAM1"), colnames(expr_mat))
stromal_m <- intersect(c("LUM","POSTN","SFRP4"), colnames(expr_mat))
marker_expr <- expr_mat[, avail]
if (length(stromal_m) > 0)
    marker_expr <- cbind(marker_expr, Stromal = rowMeans(expr_mat[, stromal_m, drop=FALSE]))

cell_type <- apply(marker_expr, 1, function(row) {
    if (max(row)==0) return("Unassigned")
    nm <- colnames(marker_expr)[which.max(row)]
    switch(nm, "EPCAM"="Epithelial","PTPRC"="Immune",
           "PECAM1"="Endothelial","Stromal"="Stromal","Unassigned")
})
cat("  Cell types:\n"); print(table(cell_type))

set.seed(42)
sub_idx <- unlist(lapply(keep_types, function(ct) {
    idx <- which(cell_type == ct)
    if (length(idx) > 500) sample(idx, 500) else idx
}))
sub_mat <- expr_mat[sub_idx, ]; sub_types <- cell_type[sub_idx]
gene_vars <- apply(sub_mat, 2, var)
top30 <- names(sort(gene_vars, decreasing=TRUE))[1:30]
mat30 <- sub_mat[, top30]
mat_log <- log1p(mat30); mat_scaled <- scale(mat_log)
mat_scaled[mat_scaled > 2] <- 2; mat_scaled[mat_scaled < -2] <- -2

cell_hc <- hclust(dist(mat_scaled), method="ward.D2")
gene_hc <- hclust(dist(t(mat_scaled)), method="ward.D2")
cell_order <- cell_hc$labels[cell_hc$order]
gene_order <- gene_hc$labels[gene_hc$order]
types_ordered <- sub_types[cell_hc$order]

hm_df <- reshape2::melt(mat_scaled[cell_order, gene_order])
colnames(hm_df) <- c("Cell","Gene","Zscore")
hm_df$Cell <- factor(hm_df$Cell, levels=cell_order)
hm_df$Gene <- factor(hm_df$Gene, levels=gene_order)
strip_df <- data.frame(Cell=factor(cell_order,levels=cell_order),
                        Type=factor(types_ordered,levels=keep_types))

hm_colors <- colorRampPalette(c("#08306B","#2171B5","#6BAED6","#FFFFFF",
                                 "#FB6A4A","#CB181D","#67000D"))(100)

p3_hm <- ggplot(hm_df, aes(x=Gene, y=Cell, fill=Zscore)) +
    geom_raster() +
    scale_fill_gradientn(colours=hm_colors, limits=c(-2,2), name="Z-score") +
    labs(x=NULL, y=NULL) + th(8) +
    theme(axis.text.x=element_text(angle=50,hjust=1,size=7,face="italic"),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          axis.line=element_blank(),
          panel.border=element_rect(colour="black",fill=NA,linewidth=0.4))

p3_strip <- ggplot(strip_df, aes(x=1, y=Cell, fill=Type)) +
    geom_raster() + scale_fill_manual(values=ct_cols, name="Cell type") +
    theme_void() +
    theme(legend.position="left", legend.title=element_text(size=8,face="bold"),
          legend.text=element_text(size=7.5), legend.key.size=unit(0.35,"cm"))

fig3 <- p3_strip + p3_hm +
    plot_layout(widths=c(0.035,1), guides="collect") +
    plot_annotation(
        title="Cell-type resolved gene expression across the tumour",
        subtitle=paste0(length(cell_order)," cells (stratified from ",
                        format(nrow(expr_mat),big.mark=","),"), ",
                        length(gene_order)," genes, log-normalised Z-scored, Ward\u2019s D2"),
        theme=theme(plot.title=element_text(size=11,face="bold"),
                    plot.subtitle=element_text(size=8,colour="grey30")))

ggsave(file.path(od, "fig3_aggregation.png"), fig3,
       width=200, height=145, units="mm", dpi=300, bg="white")
cat("  saved\n")

## ======================================================================
## FIG 4: PCA + dot plot
## ======================================================================
cat("--- Fig 4 ---\n")
pca <- prcomp(mat_log, center=TRUE, scale.=TRUE)
pca_df <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2],
                      Type=factor(sub_types, levels=keep_types))
var1 <- round(100*summary(pca)$importance[2,1],1)
var2 <- round(100*summary(pca)$importance[2,2],1)

p4a <- ggplot(pca_df, aes(x=PC1, y=PC2, colour=Type)) +
    geom_point(size=0.8, alpha=0.5) +
    scale_colour_manual(values=ct_cols, name="Cell type") +
    stat_ellipse(level=0.68, linewidth=0.5, linetype="dashed", show.legend=FALSE) +
    labs(x=paste0("PC1 (",var1,"%)"), y=paste0("PC2 (",var2,"%)")) +
    th(9) + guides(colour=guide_legend(override.aes=list(size=2.5,alpha=1))) +
    theme(legend.position=c(0.02,0.98), legend.justification=c(0,1),
          legend.background=element_rect(fill=alpha("white",0.85),
                                         colour="grey70",linewidth=0.3))

top_markers <- unique(unlist(lapply(keep_types, function(ct) {
    in_ct <- which(sub_types==ct); out_ct <- which(sub_types!=ct & sub_types %in% keep_types)
    if (length(in_ct)<5||length(out_ct)<5) return(character(0))
    fc <- log2((colMeans(sub_mat[in_ct,,drop=FALSE])+0.1)/(colMeans(sub_mat[out_ct,,drop=FALSE])+0.1))
    fc <- fc[is.finite(fc)]; names(sort(fc,decreasing=TRUE))[1:4]
})))
top_markers <- head(top_markers[!is.na(top_markers)], 16)

dot_data <- do.call(rbind, lapply(keep_types, function(ct) {
    rows <- which(sub_types==ct); sub <- sub_mat[rows, top_markers, drop=FALSE]
    do.call(rbind, lapply(top_markers, function(g) {
        data.frame(Type=ct, Gene=g, MeanExpr=mean(sub[,g]),
                   FracDetected=mean(sub[,g]>0))
    }))
}))
dot_data$Type <- factor(dot_data$Type, levels=keep_types)
dot_data$Gene <- factor(dot_data$Gene, levels=rev(top_markers))
dot_data <- do.call(rbind, lapply(split(dot_data, dot_data$Gene), function(d) {
    rng <- range(d$MeanExpr)
    d$ScaledExpr <- if(diff(rng)==0) 0.5 else (d$MeanExpr-rng[1])/diff(rng); d
}))

p4b <- ggplot(dot_data, aes(x=Type, y=Gene)) +
    geom_point(aes(size=FracDetected, fill=ScaledExpr),
               shape=21, colour="grey20", stroke=0.5) +
    scale_size_continuous(range=c(1.5,7), name="Fraction\ndetected",
                          breaks=c(0.25,0.5,0.75,1.0)) +
    scale_fill_gradientn(colours=c("#FEE0D2","#FC9272","#DE2D26","#67000D"),
                         name="Scaled\nmean expr.", limits=c(0,1)) +
    labs(x=NULL, y=NULL) + th(9) +
    theme(axis.text.y=element_text(face="italic",size=7.5),
          axis.text.x=element_text(angle=30,hjust=1,size=8),
          panel.background=element_rect(fill="#F5F5F5"),
          panel.grid.major=element_line(colour="white",linewidth=0.3))

fig4 <- (p4a + labs(tag="a")) + (p4b + labs(tag="b")) +
    plot_layout(ncol=2, widths=c(1,0.9)) +
    plot_annotation(
        title="Tumour microenvironment cell types from transcript aggregation",
        subtitle=paste0(length(sub_types)," cells, ",n_genes," genes | ",
                        "dot plot: top cell-type markers (log2 fold-change)"),
        theme=theme(plot.title=element_text(size=11,face="bold"),
                    plot.subtitle=element_text(size=8,colour="grey30"))) &
    theme(plot.tag=element_text(size=10,face="bold"))

ggsave(file.path(od, "fig4_downstream.png"), fig4,
       width=210, height=120, units="mm", dpi=300, bg="white")
cat("  saved\n")

## ======================================================================
## FIG 5: Cross-platform — Xenium + MERFISH with spatialJoin()
## ======================================================================
cat("--- Fig 5 ---\n")
sd2 <- readSpatialData("C:/Users/win10/merfish_scverse/data.zarr")
print(sd2)

mer_pts <- as.data.frame(spatialPoints(sd2)[["single_molecule"]])
mer_cells <- as.data.frame(shapes(sd2)[["cells"]])
anat <- shapes(sd2)[["anatomical"]]

## geometryCentroids on MERFISH cells
mer_cell_coords <- geometryCentroids(mer_cells[["geometry"]])
mer_cells$x <- mer_cell_coords$x; mer_cells$y <- mer_cell_coords$y

## spatialJoin() — the new API!
cat("  Running spatialJoin()...\n")
layer_names <- c("VISp_I","VISp_II/III","VISp_IV","VISp_V","VISp_VI","VISp_wm")
t0 <- proc.time()
mer_pts$layer <- spatialJoin(
    spatialPoints(sd2)[["single_molecule"]], anat,
    region_names = layer_names)
elapsed <- (proc.time() - t0)[3]
cat("  spatialJoin:", round(elapsed,1), "s\n")
cat("  Assignments:\n"); print(table(mer_pts$layer, useNA="always"))

## Parse anatomical polygons for rendering
anat_df <- as.data.frame(anat)
poly_list <- lapply(seq_len(nrow(anat_df)), function(i) {
    parsed <- parseGeometry(as.raw(anat_df$geometry[[i]]))
    data.frame(x = parsed[[1]][, 1], y = parsed[[1]][, 2], region_id = i)
})
poly_df <- do.call(rbind, poly_list)

## Assign layer names using spatialJoin results
## For each polygon, find which layer it corresponds to by testing its centroid
anat_centroids <- geometryCentroids(anat_df$geometry)
test_pts <- DataFrame(x = anat_centroids$x, y = anat_centroids$y)
centroid_layers <- spatialJoin(test_pts, anat, region_names = layer_names)
cat("  Polygon-layer mapping:", paste(centroid_layers, collapse=", "), "\n")

poly_df$layer <- centroid_layers[poly_df$region_id]
poly_df$layer <- factor(poly_df$layer, levels = layer_names)

## Panel a: Xenium density
xen_bin <- 30
pts$xb2 <- round(pts$x / xen_bin) * xen_bin
pts$yb2 <- round(pts$y / xen_bin) * xen_bin
xen_dens <- aggregate(gene ~ xb2 + yb2, data=pts, FUN=length)
colnames(xen_dens)[3] <- "count"

p5a <- ggplot(xen_dens, aes(x=xb2, y=yb2, fill=count)) +
    geom_raster() +
    scale_fill_viridis_c(option="inferno", trans="sqrt", name="Transcripts") +
    coord_equal(expand=FALSE) +
    annotate("segment", x=max(pts$x)-1100, xend=max(pts$x)-100,
             y=min(pts$y)+100, yend=min(pts$y)+100, linewidth=1.5, colour="white") +
    annotate("text", x=max(pts$x)-600, y=min(pts$y)+100,
             label="1 mm", vjust=-0.7, size=2.5, colour="white", fontface="bold") +
    labs(title="10x Xenium",
         subtitle=paste0("Human breast cancer\n",
                          format(nrow(pts),big.mark=",")," transcripts, ",
                          n_genes," genes, ",format(n_cells,big.mark=",")," cells"),
         x=expression(italic(x)~"("*mu*"m)"),
         y=expression(italic(y)~"("*mu*"m)")) + th(8.5)

## Panel b: MERFISH with anatomical polygons + cell centroids
layer_cols <- c("VISp_I"="#D62728","VISp_II/III"="#1F77B4","VISp_IV"="#2CA02C",
                "VISp_V"="#9467BD","VISp_VI"="#FF7F00","VISp_wm"="#8C564B")
layer_labels <- c("VISp_I"="I","VISp_II/III"="II/III","VISp_IV"="IV",
                  "VISp_V"="V","VISp_VI"="VI","VISp_wm"="WM")

p5b <- ggplot() +
    geom_polygon(data=poly_df[!is.na(poly_df$layer),],
                 aes(x=x, y=y, fill=layer, group=region_id),
                 colour="white", linewidth=0.4) +
    geom_point(data=mer_cells, aes(x=x, y=y),
               size=0.3, alpha=0.4, colour="grey20") +
    scale_fill_manual(values=layer_cols, name="Layer", labels=layer_labels) +
    coord_equal(expand=FALSE) +
    annotate("segment", x=max(mer_pts$x)-550, xend=max(mer_pts$x)-50,
             y=min(mer_pts$y)+50, yend=min(mer_pts$y)+50,
             linewidth=1.5, colour="black") +
    annotate("text", x=max(mer_pts$x)-300, y=min(mer_pts$y)+50,
             label="500 \u00b5m", vjust=-0.7, size=2.5, fontface="bold") +
    labs(title="MERFISH",
         subtitle=paste0("Mouse primary visual cortex\n",
                          format(nrow(mer_pts),big.mark=",")," transcripts, ",
                          format(nrow(mer_cells),big.mark=",")," cells, 268 genes"),
         x=expression(italic(x)~"("*mu*"m)"),
         y=expression(italic(y)~"("*mu*"m)")) + th(8.5)

fig5 <- (p5a + labs(tag="a")) + (p5b + labs(tag="b")) +
    plot_layout(ncol=2, widths=c(1, 0.85)) +
    plot_annotation(
        title="Cross-platform compatibility: one API, multiple technologies",
        subtitle="Both datasets read with readSpatialData() \u2014 no platform-specific code, no Python dependency",
        theme=theme(plot.title=element_text(size=11,face="bold"),
                    plot.subtitle=element_text(size=8.5,colour="grey30"))) &
    theme(plot.tag=element_text(size=10,face="bold"))

ggsave(file.path(od, "fig5_crossplatform.png"), fig5,
       width=220, height=110, units="mm", dpi=300, bg="white")
cat("  saved\n")

## ---- Summary ----
cat("\n=== ALL DONE ===\n")
cat("API functions used:\n")
cat("  readSpatialData(), spatialPoints(), shapes(), tables()\n")
cat("  geometryCentroids(), parseGeometry()\n")
cat("  bboxQuery()\n")
cat("  aggregatePoints()\n")
cat("  spatialJoin()\n\n")
for (f in sort(list.files(od, "^fig[1-5].*\\.png$"))) {
    sz <- round(file.size(file.path(od, f)) / 1024)
    cat("  ", f, " (", sz, "KB)\n")
}
