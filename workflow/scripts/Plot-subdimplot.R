log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dplyr)
library(Seurat)
library(ggplot2)

# Charger la liste d'objets Seurat (P39, P40, etc.)
seurat_obj_list <- readRDS(snakemake@input[["seurat_object"]])

# Créer une liste pour stocker les DimPlots
dim_plots <- list()

for (obj_name in names(seurat_obj_list)) {
    cat("Generating DimPlots for:", obj_name, "\n")

    # Récupérer l'objet Seurat courant
    obj <- seurat_obj_list[[obj_name]]

    # Identifier les colonnes des subclusters dans les métadonnées
    subcluster_columns <- grep("Cluster_.*_Subcluster", colnames(obj@meta.data), value = TRUE)

    for (subcluster_column in subcluster_columns) {
        cat("Generating DimPlot for subcluster:", subcluster_column, "\n")

        # Définir les identités des cellules sur le sous-cluster courant
        Idents(obj) <- obj@meta.data[[subcluster_column]]

        # Générer le DimPlot
        dim_plots[[paste0(obj_name, "_", subcluster_column)]] <- DimPlot(
            obj,
            label = TRUE,
            label.size = 3,
            reduction = "umap"
        ) +
            ggtitle(paste("DimPlot for", obj_name, "-", subcluster_column)) +
            theme(
                axis.text = element_text(size = 10),
                legend.text = element_text(size = 8),
                legend.title = element_text(size = 10)
            )
    }
}

# Sauvegarder tous les DimPlots dans un fichier PDF
pdf(file = snakemake@output[["dim_plot"]])
for (plot in dim_plots) {
    print(plot)
}
dev.off()
