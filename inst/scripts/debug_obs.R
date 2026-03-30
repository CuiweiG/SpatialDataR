.libPaths(c("C:/Users/win10/R/win-library/4.4", .libPaths()))
devtools::load_all("C:/Users/win10/SpatialDataR", quiet = TRUE)
sd <- readSpatialData("C:/Users/win10/xenium_breast/data.zarr")
tbl <- tables(sd)[["table"]]
obs <- tbl[["obs"]]
cat("obs class:", class(obs), "\n")
cat("obs names:", paste(names(obs), collapse=", "), "\n")
cat("obs nrow:", nrow(obs), "\n")
for (nm in names(obs)) {
    v <- obs[[nm]]
    cat(nm, "=> class:", class(v)[1], " length:", length(v),
        " first:", head(v, 2), "\n")
}
