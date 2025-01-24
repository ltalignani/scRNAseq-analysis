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

# Liste des packages nécessaires
packages <- c("scCustomize", "ggpubr", "hdf5r", "rliger", "ggplot2", "dplyr", "tidyr", "tibble")

# Installation des packages
for (pkg in packages) {
    install_if_missing(pkg)
}

library(dplyr)
library(stringr)
library(scCustomize)
library(ggplot2)
library(tidyr)
library(tibble)
library(gridExtra)
library(grid)

# Chargement des données
seurat_obj <- readRDS(snakemake@input[["integrated_seurat_object"]])
file_celltype <- read.csv(snakemake@input[["file_celltype"]], sep = "\t", header = FALSE)
top_10_markers <- read.csv(snakemake@input[["top_10_celltype_markers"]], sep = ",", header = TRUE)

# Prétraitement des données
names(file_celltype) <- c("gene", "cell_type")

file_celltype$gene <- str_to_title(file_celltype$gene)

# Correction pour les noms de gènes contenant un tiret
file_celltype$gene <- sapply(file_celltype$gene, function(x) {
    parts <- unlist(strsplit(x, "-"))
    if (length(parts) == 2) {
        paste0(str_to_title(parts[1]), "-", tolower(parts[2]))
    } else {
        str_to_title(x) # Si pas de tiret, appliquer str_to_title
    }
})

top_10_markers <- left_join(top_10_markers, file_celltype, by = "gene")

top_10_markers$cell_type <- top_10_markers$cell_type %>% replace_na("unknown")

# Attribution des types cellulaires aux clusters
cell_type_markers <- top_10_markers %>%
    group_by(cluster, cell_type) %>%
    summarise(count = n()) %>%
    filter(count == max(count)) %>%
    distinct(cluster, .keep_all = TRUE)

new_cluster_ids <- c(cell_type_markers$cell_type)

# Renommer les clusters avec les types cellulaires dans l'objet Seurat
names(new_cluster_ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new_cluster_ids)
seurat_obj <- AddMetaData(seurat_obj, Idents(seurat_obj),
    col.name = "Cell_type"
)

# Mettre à jour seurat_clusters également dans meta.data
seurat_obj@meta.data$seurat_clusters <- Idents(seurat_obj)

# # Calcul du nombre de cellules par cluster
cluster_counts <- seurat_obj@meta.data %>%
    count(seurat_clusters) %>%
    mutate(label = paste0(seurat_clusters, " (", n, " cells)"))

# Récupérer les labels sous forme de texte
labels_text <- paste(cluster_counts$label, collapse = "\n")

# Génération du DimPlot
dim_plot <- DimPlot(seurat_obj, label = TRUE, label.size = 3) +
    theme(
        axis.text = element_text(size = 6), # Taille des textes des axes
        legend.text = element_text(size = 4), # Taille du texte de la légende
        legend.title = element_text(size = 4), # Taille du titre de la légende
        plot.margin = margin(10, 10, 30, 10) # Ajout d'espace pour le texte
    )

# Créer une annotation contenant le texte
annotation <- textGrob(
    labels_text,
    gp = gpar(fontsize = 8),
    x = unit(0.5, "npc"), hjust = 0.5
)

# Combiner le DimPlot et l'annotation
final_plot <- arrangeGrob(dim_plot, bottom = annotation)

pdf(file = snakemake@output[["dim_plot"]])
grid.draw(final_plot)
dev.off()

saveRDS(seurat_obj,
    file = snakemake@output[["celltype_seurat_object"]]
)

# Génération des top 10 marqueurs pour chaque cluster
top_10_markers_by_cluster <- top_10_markers %>%
    distinct(cluster, gene, .keep_all = TRUE) %>% # Supprime les doublons
    group_by(cluster) %>%
    slice_head(n = 10) %>%
    select(cluster, gene)

# Récupération des couleurs des clusters depuis DimPlot
cluster_colors <- scales::hue_pal()(length(levels(seurat_obj)))
names(cluster_colors) <- levels(seurat_obj)

# Construction des données pour le tableau
table_data <- top_10_markers_by_cluster %>%
    left_join(cell_type_markers, by = c("cluster")) %>%
    mutate(color = cluster_colors[match(cell_type, names(cluster_colors))]) %>%
    group_by(cell_type) %>%
    summarise(
        markers = paste(gene, collapse = "\n"), # Liste des marqueurs par type cellulaire
        color = unique(color),
        .groups = "drop"
    )

# Graphique tableau avec ggplot2
table_plot <- ggplot(table_data) +
    geom_text(aes(x = cell_type, y = 0, label = cell_type), color = table_data$color, size = 2, fontface = "bold") + # Taille réduite
    geom_text(aes(x = cell_type, y = -1, label = markers), color = "black", size = 3) + # Taille réduite
    theme_void() +
    theme(
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10) # Taille réduite du titre
    ) +
    labs(title = "Top 10 Markers per Cell Type")

# Sauvegarde du graphique
pdf(file = snakemake@output[["top_markers_plot"]], width = 10, height = 5)
print(table_plot)
dev.off()
