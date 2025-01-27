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

# Chargement des données
seurat_obj <- readRDS(snakemake@input[["seurat_object"]])
sample <- snakemake@wildcards$sample

file_celltype <- read.csv(snakemake@input[["file_celltype"]], sep = "\t", header = FALSE)
top_10_markers <- read.csv(snakemake@input[["top_10_markers"]], sep = ",", header = TRUE)

# Préparation des données
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

cell_type_markers <- top_10_markers %>%
  group_by(cluster, cell_type) %>%
  summarise(count = n()) %>%
  filter(count == max(count)) %>%
  distinct(cluster, .keep_all = TRUE)

new_cluster_ids <- c(cell_type_markers$cell_type)
names(new_cluster_ids) <- levels(seurat_obj[[sample]])

# Mise à jour des identifiants des clusters
seurat_obj[[sample]] <- RenameIdents(seurat_obj[[sample]], new_cluster_ids)

# Génération du DimPlot avec ajustement des tailles
dim_plot <- DimPlot(
  seurat_obj[[sample]], 
  label = TRUE,         # Affichage des étiquettes
  label.size = 2        # Réduction de la taille des étiquettes
) +
  labs(x = levels((seurat_obj[[sample]]@meta.data$orig.ident))) +
  theme(
    axis.text = element_text(size = 10),     # Réduction de la taille des textes des axes
    legend.text = element_text(size = 8),   # Taille réduite des textes de la légende
    legend.title = element_text(size = 10)  # Taille réduite du titre de la légende
  )

# Sauvegarde du graphique dans un fichier PDF
pdf(file = snakemake@output[["dim_plot"]])
print(dim_plot)
dev.off()
