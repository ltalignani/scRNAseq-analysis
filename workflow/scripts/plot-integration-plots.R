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

seurat_obj <- readRDS(snakemake@input[["intergrated_seurat_object"]])

column_name <- snakemake@params[["column_name"]]

# Création du DimPlot avec des ajustements de tailles
merge_dim_plot <- DimPlot(
    seurat_obj,
    reduction = "umap.unintegrated",
    group.by = column_name,
    label = TRUE, # Ajout des étiquettes
    label.size = 3 # Taille des étiquettes réduite
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
