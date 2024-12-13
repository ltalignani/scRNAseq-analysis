log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

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

seurat_obj <- readRDS(snakemake@input[["sleuth_object"]])

heatmaps <- list()

for (i in 1: length(seurat_obj)){
  
  pos_markers[[i]] %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10_markers[[i]]

  write.csv(top10_markers[[i]],
            file = snakemake@output[["top_10_markers"]], quote = FALSE)
  pdf(file = snakemake@output[["heatmap"]])
  heatmaps[[i]] <-
    DoHeatmap(seurat_obj[[i]], features = top10_markers$gene) + NoLegend()
  dev.off()
}
