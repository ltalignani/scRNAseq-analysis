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
    valid_cells <- colSums(count_matrix) > 0
    valid_genes <- rowSums(count_matrix > 0) >= 3
    count_matrix <- count_matrix[valid_genes, valid_cells, drop = FALSE]

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

    # Filtrage des données pour éliminer cellules/gènes invalides
    filtered_matrix <- filter_invalid_cells_genes(expression_matrices[[i]])

    cat("Dimensions après filtrage pour", sample_name, ":", dim(filtered_matrix), "\n")

    # Création d'un objet Seurat
    sam_seurat_objt <- CreateSeuratObject(counts = filtered_matrix, project = sample_name)

    # Vérification et suppression des cellules restantes avec 0 UMI
    sam_seurat_objt <- subset(sam_seurat_objt, subset = nFeature_RNA > 0)

    # Ajout des pourcentages mitochondriaux
    sam_seurat_objt <- PercentageFeatureSet(sam_seurat_objt, pattern = "^mt-", col.name = "percent.mt")

    # Étape supplémentaire : filtrage basé sur des critères biologiques
    sam_seurat_objt <- subset(
        sam_seurat_objt,
        subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5
    )

    # Log des dimensions après filtrage biologique
    cat("Dimensions après filtrage biologique pour", sample_name, ":", dim(sam_seurat_objt), "\n")

    # Vérification des anomalies avant SCTransform
    counts_matrix <- GetAssayData(sam_seurat_objt, layer = "counts")
    if (any(is.na(counts_matrix))) {
        stop("NA détecté dans les comptes avant SCTransform pour", sample_name)
    }

    # Normalisation et prétraitement
    seurat_obj[[sample_name]] <- do_preprocess(sam_seurat_objt)

    # Stockage intermédiaire des objets Seurat
    intr_seurat_obj[[sample_name]] <- sam_seurat_objt
}

# Sauvegarde des objets Seurat
saveRDS(intr_seurat_obj, file = snakemake@output[["intergrated_seurat_object"]])
saveRDS(seurat_obj, file = snakemake@output[["seurat_object"]])
