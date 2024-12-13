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

groups <- snakemake@params[["column_name"]]

# Personnalisation de DimPlot avec des tailles ajustées pour les labels et la légende
merge_dim_plot <- DimPlot(
    seurat_obj,
    reduction = "umap.unintegrated",
    group.by = groups,
    label = TRUE,
    label.size = 3 # Réduction de la taille des étiquettes
) +
    theme(
        axis.text = element_text(size = 8), # Réduction de la taille du texte des axes
        legend.text = element_text(size = 6), # Réduction de la taille du texte de la légende
        legend.title = element_text(size = 8) # Réduction de la taille du titre de la légende
    )

# Sauvegarde du graphique en PDF
pdf(file = snakemake@output[["merge_dim_plot"]])
merge_dim_plot
dev.off()
