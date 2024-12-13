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
library(devtools)

#devtools::install_github('immunogenomics/presto')
# Installation conditionnelle de Presto
if (!requireNamespace("presto", quietly = TRUE)) {
  devtools::install_github("immunogenomics/presto")
}

seurat_obj <- readRDS(snakemake@input[["seurat_object"]])

genes_of_interest <- c(snakemake@params[["genes_of_interest"]])
sample <- snakemake@wildcards$sample

all_markers <- FindAllMarkers(seurat_obj[[sample]])
write.csv(all_markers,
          file = snakemake@output[["all_markers"]], quote = FALSE)
pos_markers <- FindAllMarkers(seurat_obj[[sample]], only.pos = TRUE)
pos_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_markers
write.csv(top10_markers,
          file = snakemake@output[["top_10_markers"]], quote = FALSE)

heatmap <-
  DoHeatmap(seurat_obj[[sample]], features = top10_markers$gene) + NoLegend()

pdf(file = snakemake@output[["heatmap"]], height = 20, width = 30)
heatmap
dev.off()

# Vérifier si une liste de gènes d'intérêt est fournie
if (length(genes_of_interest) > 0 && any(genes_of_interest %in% rownames(seurat_obj))) {
  # Utiliser les gènes d'intérêt
  feature_plot <- FeaturePlot(seurat_obj[[sample]],
    features = genes_of_interest, label = TRUE
  )
  pdf(file = snakemake@output[["feature_plot"]])
  feature_plot
  dev.off()

  vln_plot <- VlnPlot(seurat_obj[[sample]], features = genes_of_interest)
  pdf(file = snakemake@output[["Vln_plot"]])
  vln_plot
  dev.off()
} else if (nrow(top10_markers) > 0) {
  # Utiliser les top_10_markers
  feature_plot <- FeaturePlot(seurat_obj[[sample]],
    features = unique(top10_markers$gene), label = TRUE
  )
  pdf(file = snakemake@output[["feature_plot"]])
  feature_plot
  dev.off()

  vln_plot <- VlnPlot(seurat_obj[[sample]], features = unique(top10_markers$gene))
  pdf(file = snakemake@output[["Vln_plot"]])
  vln_plot
  dev.off()
} else {
  cat("Aucun gène d'intérêt ou top 10 markers trouvé pour générer les FeaturePlot et VlnPlot.\n")
}

# feature_plot <-
#   FeaturePlot(seurat_obj[[sample]],
#               features = genes_of_interest, label = TRUE)
# pdf(file = snakemake@output[["feature_plot"]])
# feature_plot
# dev.off()
# vln_plot <- VlnPlot(seurat_obj[[sample]], features = genes_of_interest)
# pdf(file = snakemake@output[["Vln_plot"]])
# vln_plot
# dev.off()