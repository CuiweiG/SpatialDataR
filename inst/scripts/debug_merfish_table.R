.libPaths(c("C:/Users/win10/R/win-library/4.4", .libPaths()))
devtools::load_all("C:/Users/win10/SpatialDataR", quiet = TRUE)
sd <- readSpatialData("C:/Users/win10/merfish_scverse/data.zarr")
tbl <- tables(sd)[["table"]]

cat("Table class:", class(tbl), "\n")
cat("Table names:", paste(names(tbl), collapse = ", "), "\n")

obs <- tbl[["obs"]]
cat("\nobs cols:", paste(names(obs), collapse = ", "), "\n")
cat("obs rows:", nrow(obs), "\n")
if (nrow(obs) > 0) print(head(as.data.frame(obs), 3))

var <- tbl[["var"]]
cat("\nvar cols:", paste(names(var), collapse = ", "), "\n")
cat("var rows:", nrow(var), "\n")
if (nrow(var) > 0) {
    vdf <- as.data.frame(var)
    print(head(vdf, 10))
    cat("All gene names (first 20):\n")
    cat(paste(head(vdf[[1]], 20), collapse = ", "), "\n")
}

## Also check the scverse reference for this dataset
## Look for any attribution file
attr_files <- list.files("C:/Users/win10/merfish_scverse", 
                          pattern = "README|LICENSE|CITATION|reference",
                          recursive = TRUE, full.names = TRUE)
cat("\nAttribution files:", paste(attr_files, collapse = ", "), "\n")
