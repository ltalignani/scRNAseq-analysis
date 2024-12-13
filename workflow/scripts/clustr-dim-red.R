log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


# Fonction pour installer un package si non présent
install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}

# Liste des packages à vérifier et installer si nécessaire
packages <- c("scCustomize", "ggpubr", "hdf5r", "rliger")

# Boucle pour vérifier et installer chaque package
for (pkg in packages) {
  install_if_missing(pkg)
}


library(dplyr)
library(scCustomize)
library(ggplot2)

seurat_obj <- readRDS(snakemake@input[["seurat_object"]])


do_clustering <- function(obj) {
  obj <- RunUMAP(obj,
    dims = 1:snakemake@params[["dims"]],
    verbose = FALSE
  )
  obj <- FindNeighbors(obj,
    dims = 1:snakemake@params[["dims"]], verbose = FALSE
  )
  obj <- FindClusters(obj,
    verbose = FALSE,
    resolution = snakemake@params[["resolution"]]
  )
}

for (i in 1:length(seurat_obj)) {
  seurat_obj[[i]] <- do_clustering(seurat_obj[[i]])
}

saveRDS(seurat_obj, file = snakemake@output[["seurat_object"]])
