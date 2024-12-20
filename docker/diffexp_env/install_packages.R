# install_packages.R
# Configurer les repos CRAN
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Liste des packages CRAN extraits de diff_exp.yaml
cran_packages <- c(
    "dplyr", "devtools", "ggplot2",
    "ggpubr", "ggrepel", "hdf5r",
    "rliger", "SeuratObject", "Seurat",
    "scCustomize", "sctransform", "tibble",
    "tidyr"
)

# Installer les packages CRAN
install.packages(cran_packages)

# Packages Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

bioc_packages <- c(
    "BiocGenerics", "DelayedArray", "MatrixGenerics",
    "SummarizedExperiment", "SingleCellExperiment", "SeuratObject",
    "scater", "scran", "sctransform", "ComplexHeatmap"
)

# Installer les packages Bioconductor
update.packages(ask = FALSE)
BiocManager::install(bioc_packages, ask = FALSE)

devtools::install_github("immunogenomics/presto")
