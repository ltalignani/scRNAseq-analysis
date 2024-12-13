log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# Définir un dépôt CRAN par défaut
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Fonction pour installer un package si non présent
install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}

# Liste des packages à vérifier et installer si nécessaire
packages <- c("scCustomize", "ggpubr", "hdf5r", "liger")

# Boucle pour vérifier et installer chaque package
for (pkg in packages) {
  install_if_missing(pkg)
}


library(dplyr)
library(Seurat)
library(patchwork)
library(scCustomize)
library(tibble)

# Augmenter la taille maximale des objets exportés pour future
options(future.globals.maxSize = 2 * 1024^3) # 2 Go


model <- snakemake@params[["model"]]


samples <- read.csv(snakemake@input[["samples"]], header = TRUE, sep = ",")


integration_seurat_object <-
  readRDS(snakemake@input[["integration_seurat_object"]])

merged_seurat_object <- list()

do_seurat <- function(obj) {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)
  obj <- FindNeighbors(obj,
    dims = 1:snakemake@params[["dims"]], verbose = FALSE
  )
  obj <- FindClusters(obj,
    verbose = FALSE,
    resolution = snakemake@params[["resolution"]]
  )
  obj <- RunUMAP(obj,
    dims = 1:snakemake@params[["dims"]],
    verbose = FALSE, reduction = "pca",
    reduction.name = "umap.unintegrated"
  )
}

merge_seurat_objects <- function(x, y) {
  # Combine objects (modify for your specific merging logic)
  merged_seurat <- merge(x, y, add.cell.ids = c("x", "y"))
  return(merged_seurat)
}

# Merge all Seurat objects (one object per sample)
merged_seurat <- Reduce(merge_seurat_objects, integration_seurat_object)

# Reduce the number of layers: 4 -> 1
merged_seurat <- JoinLayers(merged_seurat)

samples <- rename(samples, orig.ident = sample)

groups <- snakemake@params[["groups"]]

accessor_string <- paste0("merged_seurat$", groups)

merged_seurat@meta.data <-
  merged_seurat@meta.data %>%
  rownames_to_column("cells") %>%
  full_join(samples, by = "orig.ident") %>%
  column_to_rownames("cells")
merged_seurat[["RNA"]] <-
  split(merged_seurat[["RNA"]], f = eval(parse(text = accessor_string)))

merged_seurat <- do_seurat(merged_seurat)

saveRDS(merged_seurat, file = snakemake@output[["merge_seurat_object"]])
