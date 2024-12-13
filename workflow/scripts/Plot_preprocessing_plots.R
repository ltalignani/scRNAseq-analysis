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
model <- snakemake@params[["model"]]

seurat_obj <- readRDS(snakemake@input[["seurat_object"]])

vln_plot <- list()
variable_features <- list()
elbow_plot <- list()

for (i in 1:length(seurat_obj)) {
  vln_plot[[i]] <-
    VlnPlot(seurat_obj[[i]],
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3
    )
  plot_nC_nF <-
    FeatureScatter(seurat_obj[[i]],
      feature1 = "nCount_RNA", feature2 = "nFeature_RNA"
    )

  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(seurat_obj[[i]]), 10)

  # plot variable features with and without labels
  plot1 <- VariableFeaturePlot(seurat_obj[[i]])
  plot2 <- LabelPoints(
    plot = plot1,
    points = top10, repel = TRUE
  ) +
    labs(x = levels(seurat_obj[[i]])) +
    theme(
      legend.text = element_text(size = 6), # Ajuste la taille de la légende
      legend.position = "right", # Déplace la légende pour éviter les chevauchements
      axis.text.x = element_text(size = 8), # Ajuste la taille du texte des axes
      axis.text.y = element_text(size = 8) # Ajuste la taille du texte des axes
    )
  variable_features[[i]] <- plot1 + plot2
  elbow_plot[[i]] <- ElbowPlot(seurat_obj[[i]]) +
    labs(x = levels(seurat_obj[[i]]))
}

# Sauvegarde des graphiques
pdf(file = snakemake@output[["QC_vln_plot"]])
vln_plot
dev.off()
pdf(file = snakemake@output[["QC_nCount_nFeatures"]])
plot_nC_nF
dev.off()
pdf(file = snakemake@output[["variable_features"]])
variable_features
dev.off()
pdf(file = snakemake@output[["elbow_plot"]])
elbow_plot
dev.off()