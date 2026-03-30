#!/usr/bin/env Rscript
# ============================================================================
# SpatialDataR README — Final story: R as the best tool for spatial omics
# Story: Read → Bridge → Spatial Stats → Diff Expression → Cross-platform
# ============================================================================

.libPaths(c("C:/Users/win10/R/win-library/4.4",
            "C:/Users/win10/AppData/Local/R/win-library/4.4",
            .libPaths()))

suppressPackageStartupMessages({
    devtools::load_all("C:/Users/win10/SpatialDataR", quiet = TRUE)
    library(S4Vectors)
    library(SingleCellExperiment)
    library(SummarizedExperiment)
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

ct_cols <- c("Epithelial"="#E64B35", "Stromal"="#4DBBD5",
             "Immune"="#00A087", "Endothelial"="#F39B7F")
keep_types <- c("Epithelial","Stromal","Immune","Endothelial")

## ======================================================================
## STEP 1: Read + Convert to Bioconductor
## ======================================================================
cat("=== Step 1: Read + Bioconductor Bridge ===\n")
sd <- readSpatialData("C:/Users/win10/xenium_breast/data.zarr")
print(sd)

sce <- toSingleCellExperiment(sd)
cat("SCE:", nrow(sce), "genes x", ncol(sce), "cells\n")

coords <- geometryCentroids(shapes(sd)[["cell_circles"]][["geometry"]])
cells <- as.data.frame(shapes(sd)[["cell_circles"]])
cells$x <- coords$x; cells$y <- coords$y
obs <- as.data.frame(colData(sce))

pts_df <- spatialPoints(sd)[["transcripts"]]
pts <- as.data.frame(pts_df)
colnames(pts)[colnames(pts) == "feature_name"] <- "gene"
pts$gene <- as.character(pts$gene)
pts <- pts[pts$qv >= 20 & !grepl("^NegControl|^BLANK", pts$gene), ]

n_transcripts <- nrow(pts); n_genes <- nrow(sce); n_cells <- ncol(sce)
x_rng <- range(pts$x); y_rng <- range(pts$y)
cat("  ", format(n_transcripts, big.mark=","), "transcripts,",
    n_genes, "genes,", format(n_cells, big.mark=","), "cells\n\n")

## ======================================================================
## FIG 1: Read → SingleCellExperiment in two lines
## ======================================================================
cat("--- Fig 1 ---\n")
bin1 <- 30
pts$xb <- round(pts$x/bin1)*bin1; pts$yb <- round(pts$y/bin1)*bin1
density <- aggregate(gene ~ xb+yb, data=pts, FUN=length)
colnames(density)[3] <- "count"

p1a <- ggplot(density, aes(x=xb, y=yb, fill=count)) +
    geom_raster() +
    scale_fill_viridis_c(option="inferno", trans="sqrt", name="Transcripts\nper bin") +
    coord_equal(expand=FALSE) +
    annotate("segment", x=x_rng[2]-1100, xend=x_rng[2]-100,
             y=y_rng[1]+100, yend=y_rng[1]+100, linewidth=2, colour="white") +
    annotate("text", x=x_rng[2]-600, y=y_rng[1]+100,
             label="1 mm", vjust=-0.8, size=3, colour="white", fontface="bold") +
    labs(x=expression(italic(x)~"("*mu*"m)"),
         y=expression(italic(y)~"("*mu*"m)")) + th(9)

set.seed(42)
cell_sub <- cells[sample(nrow(cells), 50000), ]
cell_sub$total <- obs$total_counts[match(cell_sub$cell_id, obs$cell_id)]

p1b <- ggplot(cell_sub, aes(x=x, y=y, colour=total)) +
    geom_point(size=0.15, alpha=0.6) +
    scale_colour_viridis_c(option="mako", trans="sqrt", name="Total\ncounts", direction=-1) +
    coord_equal(expand=FALSE) +
    labs(x=expression(italic(x)~"("*mu*"m)"), y=NULL) +
    th(9) + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
                  axis.line.y=element_blank())

fig1 <- (p1a + labs(tag="a")) + (p1b + labs(tag="b")) +
    plot_layout(ncol=2) +
    plot_annotation(
        title="SpatialData Zarr store to Bioconductor in two function calls",
        subtitle=paste0("readSpatialData() \u2192 toSingleCellExperiment(): ",
                        format(n_transcripts, big.mark=","), " transcripts, ",
                        n_genes, " genes, ", format(n_cells, big.mark=","),
                        " cells (Xenium breast cancer, Janesick et al. 2023)"),
        theme=theme(plot.title=element_text(size=11, face="bold"),
                    plot.subtitle=element_text(size=8, colour="grey30"))) &
    theme(plot.tag=element_text(size=10, face="bold"))

ggsave(file.path(od, "fig1_store_reading.png"), fig1,
       width=210, height=105, units="mm", dpi=300, bg="white")
cat("  saved\n")

## ======================================================================
## FIG 2: Spatial query + cell type assignment
## ======================================================================
cat("--- Fig 2 ---\n")
cx <- median(pts$x); cy <- median(pts$y); hw <- 500
qx <- c(cx-hw, cx+hw); qy <- c(cy-hw, cy+hw)
pts_roi <- as.data.frame(bboxQuery(pts_df, qx[1], qx[2], qy[1], qy[2]))
colnames(pts_roi)[colnames(pts_roi)=="feature_name"] <- "gene"
pts_roi$gene <- as.character(pts_roi$gene)
pts_roi <- pts_roi[pts_roi$qv>=20 & !grepl("^NegControl|^BLANK", pts_roi$gene), ]
n_roi <- nrow(pts_roi)

p2a <- ggplot(density, aes(x=xb,y=yb,fill=count)) + geom_raster() +
    scale_fill_viridis_c(option="inferno",trans="sqrt",guide="none") +
    annotate("rect",xmin=qx[1],xmax=qx[2],ymin=qy[1],ymax=qy[2],
             fill=NA,colour="white",linewidth=1) +
    coord_equal(expand=FALSE) +
    labs(x=expression(italic(x)~"("*mu*"m)"),y=expression(italic(y)~"("*mu*"m)")) +
    th(8.5)

top6 <- names(sort(table(pts_roi$gene),decreasing=TRUE))[1:6]
pts_roi$gc <- ifelse(pts_roi$gene %in% top6, pts_roi$gene, "Other")
pts_roi$gc <- factor(pts_roi$gc, levels=c(top6,"Other"))
g6c <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#E8E8E8")
names(g6c) <- c(top6,"Other")

p2b <- ggplot() +
    geom_point(data=pts_roi[pts_roi$gc=="Other",],aes(x=x,y=y),
               colour="#E8E8E8",size=0.05,alpha=0.3) +
    geom_point(data=pts_roi[pts_roi$gc!="Other",],aes(x=x,y=y,colour=gc),
               size=0.15,alpha=0.7) +
    scale_colour_manual(values=g6c,name="Gene",
                        guide=guide_legend(override.aes=list(size=2,alpha=1))) +
    coord_equal(xlim=qx,ylim=qy,expand=FALSE) +
    annotate("segment",x=qx[2]-220,xend=qx[2]-20,y=qy[1]+30,yend=qy[1]+30,
             linewidth=1.5,colour="black") +
    annotate("text",x=qx[2]-120,y=qy[1]+30,label="200 \u00b5m",vjust=-0.7,
             size=2.5,fontface="bold") +
    labs(x=expression(italic(x)~"("*mu*"m)"),y=NULL) +
    th(8.5) + theme(panel.background=element_rect(fill="grey95"))

fig2 <- (p2a+labs(tag="a")) + (p2b+labs(tag="b")) +
    plot_layout(ncol=2,widths=c(0.8,1)) +
    plot_annotation(
        title="Spatial bounding-box query isolates tumour microenvironment",
        subtitle=paste0(format(n_roi,big.mark=",")," transcripts in 1 \u00d7 1 mm ROI"),
        theme=theme(plot.title=element_text(size=11,face="bold"),
                    plot.subtitle=element_text(size=8,colour="grey30"))) &
    theme(plot.tag=element_text(size=10,face="bold"))

ggsave(file.path(od, "fig2_spatial_query.png"), fig2,
       width=210, height=105, units="mm", dpi=300, bg="white")
cat("  saved\n")

## ======================================================================
## FIG 3: Moran's I — spatially variable genes (R's unique strength)
## ======================================================================
cat("--- Fig 3: Moran's I ---\n")
set.seed(42); idx <- sample(n_cells, 3000)
sub_mat <- t(as.matrix(assay(sce[, idx], "counts")))
sub_coords <- coords[idx, ]

t0 <- proc.time()
morans <- spatialAutocorrelation(sub_mat, sub_coords, k = 15L)
elapsed <- (proc.time()-t0)[3]
cat("  Moran's I on", nrow(morans), "genes:", round(elapsed,1), "s\n")

morans_sig <- morans[order(morans$observed, decreasing=TRUE), ]
top_spatial <- head(morans_sig, 20)
top_spatial$rank <- seq_len(20)
top_spatial$sig <- ifelse(top_spatial$adjusted.p < 0.001, "***",
                   ifelse(top_spatial$adjusted.p < 0.01, "**",
                   ifelse(top_spatial$adjusted.p < 0.05, "*", "ns")))

## Panel a: Moran's I ranked bar plot
p3a <- ggplot(top_spatial, aes(x=reorder(gene, observed), y=observed,
                                fill=observed)) +
    geom_col(width=0.7) +
    geom_text(aes(label=sig), hjust=-0.3, size=2.5, colour="grey30") +
    scale_fill_viridis_c(option="plasma", name="Moran's I", direction=-1) +
    coord_flip(ylim=c(0, max(top_spatial$observed)*1.15)) +
    labs(x=NULL, y="Moran's I (spatial autocorrelation)") +
    th(8.5) +
    theme(axis.text.y=element_text(face="italic", size=7.5))

## Panel b: spatial map of top gene using ALL 167K cells
top_gene <- top_spatial$gene[1]
all_expr <- as.numeric(assay(sce, "counts")[top_gene, ])
plot_df <- data.frame(x=coords$x, y=coords$y, expr=all_expr)
plot_df <- plot_df[order(plot_df$expr), ]  # low values first, high on top

p3b <- ggplot(plot_df, aes(x=x, y=y, colour=expr)) +
    geom_point(size=0.08, alpha=0.7, shape=16) +
    scale_colour_viridis_c(option="inferno", name="Counts",
                            trans="sqrt") +
    coord_equal(expand=FALSE) +
    labs(title=paste0(top_gene, " (I = ",
                       round(top_spatial$observed[1], 3), ")"),
         x=expression(italic(x)~"("*mu*"m)"),
         y=expression(italic(y)~"("*mu*"m)")) +
    th(8.5)

fig3 <- (p3a + labs(tag="a")) + (p3b + labs(tag="b")) +
    plot_layout(ncol=2, widths=c(0.9, 1)) +
    plot_annotation(
        title="Spatial autocorrelation identifies tissue-patterned genes",
        subtitle=paste0("Moran\u2019s I on ", nrow(morans), " genes (",
                        nrow(sub_mat), " cells, k=15 neighbours) \u2014 ",
                        "leverages R\u2019s statistical computing ecosystem"),
        theme=theme(plot.title=element_text(size=11,face="bold"),
                    plot.subtitle=element_text(size=8,colour="grey30"))) &
    theme(plot.tag=element_text(size=10,face="bold"))

ggsave(file.path(od, "fig3_spatial_autocorrelation.png"), fig3,
       width=210, height=110, units="mm", dpi=300, bg="white")
cat("  saved\n")

## ======================================================================
## FIG 4: Spatial differential expression — tumour vs stroma
## ======================================================================
cat("--- Fig 4: Spatial DE ---\n")
## Cell type assignment on subsampled cells
avail <- intersect(c("EPCAM","PTPRC","PECAM1"), colnames(sub_mat))
stromal_m <- intersect(c("LUM","POSTN","SFRP4"), colnames(sub_mat))
mk <- sub_mat[, avail]
if (length(stromal_m)>0) mk <- cbind(mk, Stromal=rowMeans(sub_mat[,stromal_m,drop=FALSE]))
sub_ct <- apply(mk, 1, function(row) {
    if(max(row)==0) return("Unassigned")
    nm <- colnames(mk)[which.max(row)]
    switch(nm,"EPCAM"="Epithelial","PTPRC"="Immune","PECAM1"="Endothelial",
           "Stromal"="Stromal","Unassigned")
})

epi_mask <- sub_ct == "Epithelial"
stro_mask <- sub_ct == "Stromal"
cat("  Epithelial:", sum(epi_mask), " Stromal:", sum(stro_mask), "\n")

de <- spatialDiffExpression(sub_mat, epi_mask, stro_mask)
de <- de[order(de$p.value), ]

## Panel a: Volcano plot
de$neglog10p <- -log10(pmax(de$adjusted.p, 1e-300))
de$sig_cat <- ifelse(abs(de$log2FC) > 1 & de$adjusted.p < 0.01, "Significant", "Not sig.")
de$label <- ifelse(de$sig_cat == "Significant" &
                    (rank(-abs(de$log2FC)) <= 8 | rank(de$adjusted.p) <= 5),
                   de$gene, "")

p4a <- ggplot(de, aes(x=log2FC, y=neglog10p, colour=sig_cat)) +
    geom_point(size=1, alpha=0.7) +
    geom_text_repel(aes(label=label), size=2.5, max.overlaps=20,
                    colour="black", fontface="italic", segment.size=0.2) +
    scale_colour_manual(values=c("Significant"="#E64B35","Not sig."="#D9D9D9"),
                        guide="none") +
    geom_vline(xintercept=c(-1,1), linetype="dashed", colour="grey60", linewidth=0.3) +
    geom_hline(yintercept=-log10(0.01), linetype="dashed", colour="grey60", linewidth=0.3) +
    labs(x="log2 fold-change (Epithelial / Stromal)",
         y=expression(-log[10]~"adjusted p-value")) +
    th(8.5)

## Panel b: spatial cell type map using ALL 167K cells
## Assign cell types on full data
full_mat <- t(as.matrix(assay(sce, "counts")))
full_avail <- intersect(c("EPCAM","PTPRC","PECAM1"), colnames(full_mat))
full_stro <- intersect(c("LUM","POSTN","SFRP4"), colnames(full_mat))
full_mk <- full_mat[, full_avail]
if (length(full_stro) > 0)
    full_mk <- cbind(full_mk, Stromal=rowMeans(full_mat[, full_stro, drop=FALSE]))
full_ct <- apply(full_mk, 1, function(row) {
    if(max(row)==0) return("Unassigned")
    nm <- colnames(full_mk)[which.max(row)]
    switch(nm,"EPCAM"="Epithelial","PTPRC"="Immune","PECAM1"="Endothelial",
           "Stromal"="Stromal","Unassigned")
})

ct_df <- data.frame(x=coords$x, y=coords$y, type=full_ct)
ct_df <- ct_df[ct_df$type %in% keep_types, ]
ct_df$type <- factor(ct_df$type, levels=keep_types)
## Shuffle to avoid overplotting bias
set.seed(42); ct_df <- ct_df[sample(nrow(ct_df)), ]

p4b <- ggplot(ct_df, aes(x=x, y=y, colour=type)) +
    geom_point(size=0.08, alpha=0.5, shape=16) +
    scale_colour_manual(values=ct_cols, name="Cell type") +
    coord_equal(expand=FALSE) +
    labs(x=expression(italic(x)~"("*mu*"m)"),
         y=expression(italic(y)~"("*mu*"m)")) +
    th(8.5) +
    guides(colour=guide_legend(override.aes=list(size=2.5,alpha=1)))

fig4 <- (p4a + labs(tag="a")) + (p4b + labs(tag="b")) +
    plot_layout(ncol=2, widths=c(1, 0.9)) +
    plot_annotation(
        title="Spatial differential expression between tumour and stroma",
        subtitle=paste0(sum(de$sig_cat=="Significant"), " genes with |log2FC| > 1 ",
                        "and adjusted p < 0.01 (Wilcoxon rank-sum test)"),
        theme=theme(plot.title=element_text(size=11,face="bold"),
                    plot.subtitle=element_text(size=8,colour="grey30"))) &
    theme(plot.tag=element_text(size=10,face="bold"))

ggsave(file.path(od, "fig4_spatial_de.png"), fig4,
       width=210, height=110, units="mm", dpi=300, bg="white")
cat("  saved\n")

## ======================================================================
## FIG 5: Cross-platform + spatstat integration
## ======================================================================
cat("--- Fig 5: Cross-platform ---\n")
sd2 <- readSpatialData("C:/Users/win10/merfish_scverse/data.zarr")

mer_pts <- as.data.frame(spatialPoints(sd2)[["single_molecule"]])
mer_cells <- as.data.frame(shapes(sd2)[["cells"]])
mer_coords <- geometryCentroids(mer_cells[["geometry"]])
mer_cells$x <- mer_coords$x; mer_cells$y <- mer_coords$y
anat <- shapes(sd2)[["anatomical"]]

## spatialJoin
layer_names <- c("VISp_I","VISp_II/III","VISp_IV","VISp_V","VISp_VI","VISp_wm")
mer_pts$layer <- spatialJoin(spatialPoints(sd2)[["single_molecule"]], anat,
                              region_names = layer_names)

## Parse anatomical polygons
anat_df <- as.data.frame(anat)
poly_list <- lapply(seq_len(nrow(anat_df)), function(i) {
    parsed <- parseGeometry(as.raw(anat_df$geometry[[i]]))
    data.frame(x=parsed[[1]][,1], y=parsed[[1]][,2], region_id=i)
})
poly_df <- do.call(rbind, poly_list)

## Map polygons to cortical layers using transcript cell_type annotations
## (centroid-based mapping fails for thin layers like VISp_I)
cat("  Mapping polygons to layers via transcript annotations...\n")
poly_idx <- spatialJoin(spatialPoints(sd2)[["single_molecule"]], anat,
                         region_names = paste0("poly", 1:6))
pts_with_poly <- data.frame(cell_type = mer_pts$cell_type, polygon = poly_idx)
pts_with_poly <- pts_with_poly[!is.na(pts_with_poly$polygon), ]

## For each polygon, find the layer it best represents
## Use exclusive assignment: which layer has the highest FRACTION inside this polygon?
cortex_types <- c("VISp_I","VISp_II/III","VISp_IV","VISp_V","VISp_VI","VISp_wm")
poly_layer_map <- character(6)
assigned_layers <- character(0)

for (pass in 1:6) {
    best_score <- -1; best_poly <- 0; best_layer <- ""
    for (pid in 1:6) {
        if (poly_layer_map[pid] != "") next
        in_poly <- pts_with_poly[pts_with_poly$polygon == paste0("poly", pid), ]
        for (layer in setdiff(cortex_types, assigned_layers)) {
            n_layer_in_poly <- sum(in_poly$cell_type == layer)
            n_layer_total <- sum(pts_with_poly$cell_type == layer)
            if (n_layer_total == 0) next
            frac <- n_layer_in_poly / n_layer_total
            if (frac > best_score) {
                best_score <- frac; best_poly <- pid; best_layer <- layer
            }
        }
    }
    if (best_poly > 0) {
        poly_layer_map[best_poly] <- best_layer
        assigned_layers <- c(assigned_layers, best_layer)
        cat("    poly", best_poly, "->", best_layer,
            "(", round(best_score*100,1), "% of layer's transcripts)\n")
    }
}

poly_df$layer <- poly_layer_map[poly_df$region_id]
poly_df$layer <- factor(poly_df$layer, levels=layer_names)

layer_cols <- c("VISp_I"="#D62728","VISp_II/III"="#1F77B4","VISp_IV"="#2CA02C",
                "VISp_V"="#9467BD","VISp_VI"="#FF7F00","VISp_wm"="#8C564B")
layer_labels <- c("VISp_I"="I","VISp_II/III"="II/III","VISp_IV"="IV",
                  "VISp_V"="V","VISp_VI"="VI","VISp_wm"="WM")

## Panel a: Xenium density
p5a <- ggplot(density, aes(x=xb,y=yb,fill=count)) + geom_raster() +
    scale_fill_viridis_c(option="inferno",trans="sqrt",name="Transcripts") +
    coord_equal(expand=FALSE) +
    annotate("segment",x=x_rng[2]-1100,xend=x_rng[2]-100,
             y=y_rng[1]+100,yend=y_rng[1]+100,linewidth=1.5,colour="white") +
    annotate("text",x=x_rng[2]-600,y=y_rng[1]+100,label="1 mm",
             vjust=-0.7,size=2.5,colour="white",fontface="bold") +
    labs(title="10x Xenium",
         subtitle=paste0("Human breast cancer\n",format(n_transcripts,big.mark=","),
                          " transcripts, ",n_genes," genes"),
         x=expression(italic(x)~"("*mu*"m)"),y=expression(italic(y)~"("*mu*"m)")) +
    th(8.5)

## Panel b: MERFISH polygons + cells
p5b <- ggplot() +
    geom_polygon(data=poly_df[!is.na(poly_df$layer),],
                 aes(x=x,y=y,fill=layer,group=region_id),
                 colour="white",linewidth=0.4) +
    geom_point(data=mer_cells, aes(x=x,y=y), size=0.3, alpha=0.4, colour="grey20") +
    scale_fill_manual(values=layer_cols, name="Layer", labels=layer_labels) +
    coord_equal(expand=FALSE) +
    annotate("segment",x=max(mer_pts$x)-550,xend=max(mer_pts$x)-50,
             y=min(mer_pts$y)+50,yend=min(mer_pts$y)+50,linewidth=1.5,colour="black") +
    annotate("text",x=max(mer_pts$x)-300,y=min(mer_pts$y)+50,
             label="500 \u00b5m",vjust=-0.7,size=2.5,fontface="bold") +
    labs(title="MERFISH",
         subtitle=paste0("Mouse primary visual cortex\n",
                          format(nrow(mer_pts),big.mark=",")," transcripts, ",
                          format(nrow(mer_cells),big.mark=",")," cells"),
         x=expression(italic(x)~"("*mu*"m)"),y=expression(italic(y)~"("*mu*"m)")) +
    th(8.5)

fig5 <- (p5a+labs(tag="a")) + (p5b+labs(tag="b")) +
    plot_layout(ncol=2,widths=c(1,0.85)) +
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
cat("API functions demonstrated:\n")
cat("  readSpatialData, toSingleCellExperiment, geometryCentroids\n")
cat("  bboxQuery, aggregatePoints, spatialJoin, parseGeometry\n")
cat("  spatialAutocorrelation, spatialDiffExpression, toPointPattern\n\n")
for (f in sort(list.files(od, "^fig[1-5].*\\.png$"))) {
    sz <- round(file.size(file.path(od,f))/1024)
    cat("  ",f," (",sz,"KB)\n")
}
