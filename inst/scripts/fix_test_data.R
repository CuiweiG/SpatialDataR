#!/usr/bin/env Rscript
## Regenerate transcripts.csv with biologically plausible
## cell-type-specific gene expression patterns
.libPaths("C:/Users/win10/R/win-library/4.4")

store <- "inst/extdata/xenium_mini.zarr"
obs <- read.csv(file.path(store, "tables", "table",
    "obs", "obs.csv"))
shp <- read.csv(file.path(store, "shapes",
    "cell_boundaries", "circles.csv"))

set.seed(2024)
n_tx <- 500

## Define cell-type-specific gene probabilities
## Based on known breast cancer markers
## Biologically accurate breast cancer Xenium panel
## Based on Janesick et al. 2022 bioRxiv (10x Xenium)
## Epithelial: EPCAM, KRT18, ESR1, PGR, HER2/ERBB2
## Immune (CD45+): CD45, MKI67 (proliferating)
## Stromal: VIM, ACTB
gene_probs <- list(
    Epithelial = c(EPCAM = 0.25, KRT18 = 0.25,
        ESR1 = 0.12, PGR = 0.10, HER2 = 0.10,
        ERBB2 = 0.08, MKI67 = 0.04, CD45 = 0.01,
        VIM = 0.02, ACTB = 0.03),
    Stromal = c(VIM = 0.40, ACTB = 0.30,
        MKI67 = 0.05, EPCAM = 0.01, KRT18 = 0.01,
        ESR1 = 0.01, PGR = 0.01, HER2 = 0.01,
        CD45 = 0.01, ERBB2 = 0.19),
    Immune = c(CD45 = 0.50, MKI67 = 0.20,
        VIM = 0.05, ACTB = 0.05, EPCAM = 0.01,
        KRT18 = 0.01, ESR1 = 0.01, PGR = 0.01,
        HER2 = 0.01, ERBB2 = 0.15),
    Endothelial = c(VIM = 0.25, ACTB = 0.25,
        CD45 = 0.05, MKI67 = 0.10, EPCAM = 0.01,
        KRT18 = 0.01, ESR1 = 0.01, PGR = 0.01,
        HER2 = 0.01, ERBB2 = 0.30)
)

## Assign transcripts to cells proportionally
cell_ids <- rep(obs$cell_id,
    length.out = n_tx)[sample(n_tx)]

## Get cell type for each transcript
cell_type_map <- setNames(obs$cell_type, obs$cell_id)
tx_types <- cell_type_map[as.character(cell_ids)]

## Sample genes based on cell type
genes <- character(n_tx)
for (i in seq_len(n_tx)) {
    ct <- tx_types[i]
    probs <- gene_probs[[ct]]
    genes[i] <- sample(names(probs), 1, prob = probs)
}

## Generate coordinates near cell centers
tx_x <- numeric(n_tx)
tx_y <- numeric(n_tx)
for (i in seq_len(n_tx)) {
    cid <- cell_ids[i]
    row <- shp[shp$cell_id == cid, ]
    if (nrow(row) > 0) {
        tx_x[i] <- row$x[1] + rnorm(1, 0, row$radius[1])
        tx_y[i] <- row$y[1] + rnorm(1, 0, row$radius[1])
    } else {
        tx_x[i] <- runif(1, 0, 4.25)
        tx_y[i] <- runif(1, 0, 4.25)
    }
}

tx_df <- data.frame(
    x = pmax(0, tx_x),
    y = pmax(0, tx_y),
    gene = genes,
    cell_id = cell_ids,
    stringsAsFactors = FALSE)

write.csv(tx_df,
    file.path(store, "points", "transcripts",
        "transcripts.csv"),
    row.names = FALSE)

cat("Regenerated", nrow(tx_df), "transcripts\n")
cat("Gene distribution by cell type:\n")
print(table(genes, tx_types))
