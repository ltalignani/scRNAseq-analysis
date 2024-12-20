log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# Définir un miroir CRAN si aucun n'est défini
#if (is.null(getOption("repos")) || getOption("repos")["CRAN"] == "@CRAN@") {
#  options(repos = c(CRAN = "https://cloud.r-project.org/"))
#}

# Fonction pour installer un package si non présent
#install_if_missing <- function(package) {
#  if (!require(package, character.only = TRUE)) {
#    install.packages(package)
#    library(package, character.only = TRUE)
#  }
#}


# Vérification et installation de BiocManager si nécessaire
#if (!requireNamespace("BiocManager", quietly = TRUE)) {
#  install.packages("BiocManager")
#}

# Vérification et installation de ComplexHeatmap si nécessaire
#if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
#  BiocManager::install("ComplexHeatmap")
#}


# Liste des packages à vérifier et installer si nécessaire
#packages <- c("scCustomize", "ggpubr", "hdf5r", "rliger")

# Boucle pour vérifier et installer chaque package
#for (pkg in packages) {
#  install_if_missing(pkg)
#}


library(dplyr)
library(scCustomize)
library(ggplot2)
# library(ggrepel)
library(ComplexHeatmap)

seurat_obj <- readRDS(snakemake@input[["celltype_seurat_object"]])

column_name <- snakemake@params[["column_name"]]
base_level <- snakemake@params[["base_level"]]
comparison <- snakemake@params[["comparison_variable"]]
cell_type <- snakemake@wildcards$celltype


cell_type_accessor_string <- paste0("seurat_obj$Cell_type")
group_accessor_string <- paste0("seurat_obj$", column_name)

seurat_obj$celltype_group <-
  paste(eval(parse(text = cell_type_accessor_string)),
    eval(parse(text = group_accessor_string)),
    sep = "_"
  )

Idents(seurat_obj) <- "celltype_group"


id_1 <- paste(cell_type, base_level, sep = "_")
id_2 <- paste(cell_type, comparison, sep = "_")


markers <- FindMarkers(seurat_obj, ident.1 = id_1, ident.2 = id_2)

# Calculer le heatmap pour les 25 premiers marqueurs
heatmapdf <- markers[1:min(25, nrow(markers)), ] # Limite à 25 ou au nombre total de lignes
row_ha <- rowAnnotation(
  "wt" = anno_barplot(heatmapdf$pct.1),
  "ko" = anno_barplot(heatmapdf$pct.2),
  width = unit(10, "cm")
)

ht0 <- Heatmap(heatmapdf$avg_log2FC,
  name = "Log2FC",
  cluster_rows = TRUE,
  row_labels = rownames(heatmapdf),
  right_annotation = row_ha,
  width = unit(1, "cm")
)

# Sauvegarder le heatmap dans un fichier PDF
pdf(file = snakemake@output[["log2fc_heatmap"]])
draw(ht0) # Utiliser `draw` pour rendre le heatmap dans le fichier PDF
dev.off()


write.csv(markers,
  file = snakemake@output[["all_markers"]], quote = FALSE
)
pos_markers <- FindMarkers(seurat_obj,
  ident.1 = id_1,
  ident.2 = id_2, only.pos = TRUE
)
pos_markers %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_markers

# Ajouter une colonne `gene` pour inclure les noms des gènes
top10_markers$gene <- rownames(top10_markers)

write.csv(top10_markers,
  file = snakemake@output[["top_10_markers"]], quote = FALSE
)


# heatmap <-
#   DoHeatmap(seurat_obj, features = top10_markers$gene) + NoLegend()

# pdf(file = snakemake@output[["heatmap"]])
# heatmap
# dev.off()

# volcano plot

markers$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
markers$diffexpressed[markers$avg_log2FC > 1 & markers$p_val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
markers$diffexpressed[markers$avg_log2FC < -1 & markers$p_val < 0.05] <- "DOWN"



markers$gene <- row.names(markers)
markers$delabel <- NA
markers$delabel[markers$diffexpressed != "NO"] <-
  markers$gene[markers$diffexpressed != "NO"]


volcano_plot <- ggplot(data = markers, aes(
  x = avg_log2FC, y = -log10(p_val),
  col = diffexpressed, label = delabel
)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = c("blue", "black", "red")) +
  labs(title = column_name)

pdf(file = snakemake@output[["volcano_plot"]])
volcano_plot
dev.off()

