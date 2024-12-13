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

# Création d'une liste pour stocker les graphiques DimPlot avec des ajustements de tailles
dim_plot <- list()
for (i in 1:length(seurat_obj)) {
  dim_plot[[i]] <- DimPlot(
    seurat_obj[[i]],
    label = TRUE, # Affichage des étiquettes
    label.size = 3 # Taille réduite des étiquettes
  ) +
    labs(x = levels((seurat_obj[[i]]@meta.data$orig.ident))) +
    theme(
      axis.text = element_text(size = 8), # Réduction de la taille des textes des axes
      legend.text = element_text(size = 6), # Taille réduite des textes de la légende
      legend.title = element_text(size = 8) # Taille réduite du titre de la légende
    )
}

# Sauvegarde des graphiques en PDF
pdf(file = snakemake@output[["dim_plot"]])
for (plot in dim_plot) {
  print(plot) # Impression de chaque graphique dans le PDF
}
dev.off()
