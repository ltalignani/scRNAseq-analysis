# Snakemake workflow: SingleCell-10X-Genomics-RNA-seq

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/manuelphilip/Single-cell-analysis-Snakemake-workflow/workflows/Tests/badge.svg?branch=main)](https://github.com/manuelphilip/Single-cell-analysis-Snakemake-workflow/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for single-cell data analysis of 10X genomics data including cell type annotation, differential expression (marker gene identification), scRNA-seq integration


## Usage

# Description
 A Snakemake workflow designed to analyse single-cell data from cellranger count output. The workflow expects the cellranger count output which contains per sample bc_matrix (raw and filtered) under the `~/sample_name/outs` folder. For more information please refer [cellranger ouput](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-gex-overview)
 
The general steps are as follows:
All the steps are carried out using [Seurat](https://satijalab.org/seurat/articles/get_started_v5_new>)
 * Preprocessing 
 * Clustering and dimensional reduction
 * Marker identification 
 * Assigning cell types to clusters (Automate cell type assignment using marker genes)
 * Integrative analysis
 * Default differential expression tests across models 
 * Differential expression across samples within the same cell types
 * Gene Ontology analysis using [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
 * Pathway enrichment analysis using clusterProfiler
