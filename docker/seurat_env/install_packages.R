# install_packages.R

# Définir le dépôt CRAN
options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_packages <- c(
  "dplyr", "SeuratObject", "Seurat", "scCustomize", "ggpubr",
  "hdf5r", "rliger", "ggplot2", "tidyr", "tibble",
  "stringr", "patchwork", "harmony", "devtools"
)

# Installer les packages
install.packages(cran_packages)

# Installer le package depuis github
devtools::install_github("immunogenomics/presto")

# Installer des packages Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Installer les packages Bioconductor avec versions exactes
bioc_packages <- list(
  "beachmat" = "2.18.0",
  "DelayedArray" = "0.28.0",
  "DelayedMatrixStats" = "1.24.0",
  "GenomicRanges" = "1.54.1",
  "GenomeInfoDb" = "1.38.1",
  "GenomeInfoDbData" = "1.2.11",
  "glmGamPoi" = "1.14.0",
  "HDF5Array" = "1.30.0",
  "IRanges" = "2.36.0",
  "MatrixGenerics" = "1.14.0",
  "rhdf5" = "2.46.1",
  "rhdf5filters" = "1.14.1",
  "Rhdf5lib" = "1.24.0",
  "S4Arrays" = "1.2.0",
  "S4Vectors" = "0.40.2",
  "SingleCellExperiment" = "1.24.0",
  "SparseArray" = "1.2.2",
  "sparseMatrixStats" = "1.14.0",
  "SummarizedExperiment" = "1.32.0",
  "XVector" = "0.42.0",
  "zlibbioc" = "1.48.0",
  "Biobase" = "2.62.0",
  "BiocGenerics" = "0.48.1"
)

# Installer chaque package Bioconductor avec la version exacte
for (pkg in names(bioc_packages)) {
  tryCatch(
    {
      BiocManager::install(pkg, version = bioc_packages[[pkg]], ask = FALSE)
    },
    error = function(e) {
      message(sprintf("Erreur lors de l'installation de %s version %s : %s", pkg, bioc_packages[[pkg]], e$message))
    }
  )
}
