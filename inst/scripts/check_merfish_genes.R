.libPaths(c("C:/Users/win10/R/win-library/4.4", .libPaths()))
devtools::load_all("C:/Users/win10/SpatialDataR", quiet = TRUE)

sd2 <- readSpatialData("C:/Users/win10/merfish_scverse/data.zarr")

## Try to get expression matrix via toSingleCellExperiment
tryCatch({
    sce2 <- toSingleCellExperiment(sd2)
    cat("SCE:", nrow(sce2), "genes x", ncol(sce2), "cells\n")
    cat("First 20 genes:", paste(head(rownames(sce2), 20), collapse=", "), "\n")
    
    ## Check for known cortical layer markers
    markers <- c("Rorb","Ctgf","Cux2","Foxp2","Bcl11b","Fezf2",
                 "Slc17a7","Gad1","Gad2","Pvalb","Sst","Vip","Aqp4")
    found <- intersect(markers, rownames(sce2))
    cat("\nKnown markers found:", paste(found, collapse=", "), "\n")
    
    ## Top expressed genes
    gene_sums <- rowSums(SummarizedExperiment::assay(sce2, "counts"))
    top20 <- head(sort(gene_sums, decreasing=TRUE), 20)
    cat("\nTop 20 by total counts:\n")
    print(top20)
    
    ## Check spatial coords
    cd <- SummarizedExperiment::colData(sce2)
    cat("\ncolData cols:", paste(names(cd), collapse=", "), "\n")
    if ("x_centroid" %in% names(cd)) {
        cat("x range:", range(cd$x_centroid, na.rm=TRUE), "\n")
        cat("y range:", range(cd$y_centroid, na.rm=TRUE), "\n")
    }
}, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
})
