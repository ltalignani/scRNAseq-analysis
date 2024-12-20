# install_packages.R
# Configurer les repos CRAN
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Liste des packages CRAN extraits de enrichment.yaml
cran_packages <- c("ggplot2", "forcats", "cowplot", "dplyr", "aplot")

# Installer les packages CRAN
install.packages(cran_packages)

# Installer les packages Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

bioc_packages <- c(
  "AnnotationDbi", "org.Mm.eg.db", "enrichplot",
  "clusterProfiler", "ReactomePA"
)

update.packages(ask = FALSE)
BiocManager::install(bioc_packages, ask = FALSE)