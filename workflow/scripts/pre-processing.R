log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dplyr)
library(Seurat)
library(patchwork)
library(scCustomize)

model <- snakemake@params[["model"]]

# Lecture des matrices d'expression
expression_matrices <- Read10X_Multi_Directory(base_path = snakemake@params[["path"]])

# Lecture des échantillons
samples <- read.csv(snakemake@input[["samples"]], header = TRUE, sep = "\t")

# Fonction de filtrage
filter_invalid_cells_genes <- function(count_matrix) {
    # Identifiez les cellules valides
    valid_cells <- colSums(count_matrix) > 0
    valid_genes <- rowSums(count_matrix > 0) >= 3

    # Appliquer les deux filtres directement
    count_matrix <- count_matrix[valid_genes, valid_cells, drop = FALSE]

    # Vérifier les valeurs NA, NaN ou infinies
    if (any(is.na(count_matrix))) {
      stop("La matrice contient des NA après le filtrage.")
    }
    if (any(!is.finite(count_matrix))) {
      stop("La matrice contient des valeurs infinies après le filtrage.")
    }
    return(count_matrix)
}

# Fonction de prétraitement
do_preprocess <- function(obj) {
    obj <- PercentageFeatureSet(obj, pattern = "^mt-", col.name = "percent.mt")
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
    obj <- ScaleData(obj, vars.to.regress = "percent.mt")
    obj <- RunPCA(obj, features = VariableFeatures(object = obj))
    return(obj)
}

# Initialisation des listes pour stocker les objets Seurat
seurat_obj <- list()
intr_seurat_obj <- list()

# Boucle sur les échantillons
for (i in 1:length(samples[[1]])) {
  sample_name <- names(expression_matrices[i])

  # Filtrer les données pour éliminer les cellules/gènes invalides
  filtered_matrix <- filter_invalid_cells_genes(expression_matrices[[i]])

  # Loguer les dimensions après filtrage
  cat("Dimensions après filtrage pour", sample_name, ":",
      dim(filtered_matrix), "\n")

  # Créer un objet Seurat
  sam_seurat_objt <- CreateSeuratObject(counts = filtered_matrix, project = sample_name)

  # Vérification et suppression des cellules restantes avec 0 UMI
  sam_seurat_objt <- subset(sam_seurat_objt, subset = nFeature_RNA > 0)

  # Calcul des pourcentages mitochondriaux pour chaque cellule
  intr_seurat_obj[[sample_name]] <- PercentageFeatureSet(sam_seurat_objt,
                                                         pattern = "^mt-",
                                                         col.name = "percent.mt")

  # Vérification des anomalies avant SCTransform
  counts_matrix <- GetAssayData(sam_seurat_objt, layer = "counts")
  if (any(is.na(counts_matrix))) {
    stop("NA détecté dans les comptes avant SCTransform pour", sample_name)
  }

  # Normalisation et prétraitement
  seurat_obj[[sample_name]] <- do_preprocess(sam_seurat_objt)
}

# Sauvegarde des objets Seurat
saveRDS(intr_seurat_obj, file = snakemake@output[["intergrated_seurat_object"]])
saveRDS(seurat_obj, file = snakemake@output[["seurat_object"]])
