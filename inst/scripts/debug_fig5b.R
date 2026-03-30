.libPaths(c("C:/Users/win10/R/win-library/4.4", .libPaths()))
devtools::load_all("C:/Users/win10/SpatialDataR", quiet = TRUE)

sd <- readSpatialData("C:/Users/win10/merfish_scverse/data.zarr")
pts <- as.data.frame(spatialPoints(sd)[["single_molecule"]])

region_layers <- c("VISp_I", "VISp_II/III", "VISp_IV",
                   "VISp_V", "VISp_VI", "VISp_wm")

## Check: how many bins get each layer BEFORE gap-fill?
mer_bin <- 8
pts$xb <- round(pts$x / mer_bin) * mer_bin
pts$yb <- round(pts$y / mer_bin) * mer_bin

bin_layer <- aggregate(cell_type ~ xb + yb, data = pts,
                        FUN = function(x) {
                            tt <- table(x)
                            cx <- intersect(names(tt), region_layers)
                            if (length(cx) > 0) {
                                cx_tt <- tt[cx]
                                names(cx_tt)[which.max(cx_tt)]
                            } else {
                                "Tissue"
                            }
                        })

cat("BEFORE gap-fill:\n")
print(table(bin_layer$cell_type))
cat("VISp_I bins:", sum(bin_layer$cell_type == "VISp_I"), "\n")

## Now gap-fill and check
filled <- bin_layer
rownames(filled) <- paste(filled$xb, filled$yb)
for (iteration in 1:5) {
    existing_keys <- paste(filled$xb, filled$yb)
    grid_x <- seq(min(filled$xb), max(filled$xb), by = mer_bin)
    grid_y <- seq(min(filled$yb), max(filled$yb), by = mer_bin)
    all_keys <- paste(rep(grid_x, each = length(grid_y)),
                      rep(grid_y, length(grid_x)))
    missing_keys <- setdiff(all_keys, existing_keys)
    new_rows <- list()
    for (mk in missing_keys) {
        coords <- as.numeric(strsplit(mk, " ")[[1]])
        mx <- coords[1]; my <- coords[2]
        nbr_keys <- paste(
            rep(mx + c(-1, 0, 1) * mer_bin, each = 3),
            rep(my + c(-1, 0, 1) * mer_bin, 3))
        nbr_keys <- setdiff(nbr_keys, mk)
        nbr_layers <- filled$cell_type[match(nbr_keys, existing_keys)]
        nbr_layers <- nbr_layers[!is.na(nbr_layers)]
        if (length(nbr_layers) >= 3) {
            tt <- table(nbr_layers)
            new_rows[[mk]] <- data.frame(
                xb = mx, yb = my,
                cell_type = names(tt)[which.max(tt)])
        }
    }
    if (length(new_rows) > 0) {
        added <- do.call(rbind, new_rows)
        filled <- rbind(filled, added)
        rownames(filled) <- paste(filled$xb, filled$yb)
        cat("Iter", iteration, ": added", nrow(added), "bins\n")
    }
}

cat("\nAFTER gap-fill:\n")
print(table(filled$cell_type))
cat("VISp_I bins:", sum(filled$cell_type == "VISp_I"), "\n")
cat("Tissue bins:", sum(filled$cell_type == "Tissue"), "\n")
