library(Seurat)
library(dplyr)
library(scCustomize)
library(ggplot2)
library(ComplexHeatmap)

# Charger les objets Seurat
seurat_obj <- readRDS(snakemake@input[["seurat_object"]])

# Fonction pour effectuer le clustering principal
do_clustering <- function(obj) {
    obj <- RunUMAP(obj, dims = 1:snakemake@params[["dims"]], verbose = FALSE)
    obj <- FindNeighbors(obj, dims = 1:snakemake@params[["dims"]], graph.name = "integrated_snn", verbose = FALSE)
    obj <- FindClusters(obj, resolution = snakemake@params[["resolution"]], verbose = FALSE)
    return(obj)
}

# Appliquer le clustering principal sur chaque objet Seurat
for (i in 1:length(seurat_obj)) {
    seurat_obj[[i]] <- do_clustering(seurat_obj[[i]])
}

# Fonction pour effectuer le subclustering sur tous les clusters
do_subclustering <- function(obj, resolution = snakemake@params[["resolution"]]) {
    clusters <- unique(Idents(obj)) # Obtenir tous les clusters
    for (cluster in clusters) {
        cat("Subclustering cluster:", cluster, "\n")
        obj <- FindSubCluster(
            obj,
            cluster = cluster,
            graph.name = "integrated_snn",
            subcluster.name = "sub.cluster", # Nom pour stocker les sous-clusters
            resolution = resolution,
            algorithm = 1 # Louvain par défaut
        )
    }
    return(obj)
}

# Appliquer le subclustering sur chaque objet Seurat
for (i in 1:length(seurat_obj)) {
    cat("Processing Seurat object:", names(seurat_obj)[i], "\n")
    seurat_obj[[i]] <- do_subclustering(seurat_obj[[i]])
}

# Créer des graphiques DimPlot pour les subclusters
dim_plot <- list()
for (i in 1:length(seurat_obj)) {
    dim_plot[[i]] <- DimPlot(
        seurat_obj[[i]],
        group.by = "sub.cluster", # Regrouper par sous-clusters
        label = TRUE, # Afficher les étiquettes
        label.size = 3 # Taille des étiquettes
    ) +
        labs(title = paste("Sample:", names(seurat_obj)[i])) +
        theme(
            axis.text = element_text(size = 10), # Texte des axes
            legend.text = element_text(size = 8), # Texte des légendes
            legend.title = element_text(size = 10) # Titre des légendes
        )
}

# Sauvegarder les graphiques en PDF
pdf(file = snakemake@output[["subdim_plot"]])
for (plot in dim_plot) {
    print(plot) # Impression de chaque graphique
}
dev.off()

# Sauvegarder les objets Seurat mis à jour
saveRDS(seurat_obj, file = snakemake@input[["seurat_object"]])
