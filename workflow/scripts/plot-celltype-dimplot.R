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
library(stringr)
library(scCustomize)
library(ggplot2)
library(tidyr)

seurat_obj <- readRDS(snakemake@input[["integrated_seurat_object"]])

file_celltype <- read.csv(snakemake@input[["file_celltype"]], sep = "\t", header = FALSE)

top_10_markers <- read.csv(snakemake@input[["top_10_celltype_markers"]], sep = ",", header = TRUE)

names(file_celltype) <- c("gene", "cell_type")

file_celltype$gene <- str_to_title(file_celltype$gene)
top_10_markers <- left_join(top_10_markers, file_celltype, by = "gene")

top_10_markers$cell_type <- top_10_markers$cell_type %>% replace_na("unknown")
cell_type_markers <- top_10_markers %>%
  group_by(cluster, cell_type) %>%
  summarise(count = n()) %>%
  filter(count == max(count)) %>%
  distinct(cluster, .keep_all = TRUE)

new_cluster_ids <- c(cell_type_markers$cell_type)

names(new_cluster_ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new_cluster_ids)
seurat_obj <- AddMetaData(seurat_obj, Idents(seurat_obj), col.name = "Cell_type")

# Génération du DimPlot avec ajustements des tailles
dim_plot <- DimPlot(
  seurat_obj, 
  label = TRUE,         # Affichage des étiquettes
  label.size = 3        # Taille réduite des étiquettes
) +
  labs(x = levels((seurat_obj@meta.data$orig.ident))) +
  theme(
    axis.text = element_text(size = 10),     # Réduction de la taille des textes des axes
    legend.text = element_text(size = 8),   # Taille réduite des textes de la légende
    legend.title = element_text(size = 10)  # Taille réduite du titre de la légende
  )

# Sauvegarde du DimPlot en PDF
pdf(file = snakemake@output[["dim_plot"]])
print(dim_plot)
dev.off()

# Sauvegarde de l'objet Seurat avec les identifiants mis à jour
saveRDS(seurat_obj, file = snakemake@output[["celltype_seurat_object"]])
