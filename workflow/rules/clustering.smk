rule perform_clustering_and_dim_reduction:
    input:
        seurat_object="results/seurat/preprocessing/all.seurat_objt.rds",
    output:
        seurat_object="results/seurat/clustering/all.seurat_objt.rds",
    resources:
        cpus_per_task=20,
        mem_mb=94000,
        nodes=10,
        runtime=20,
    params:
        dims=config["clustering"]["dim"],
        resolution=config["clustering"]["resolution"],
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/seurat/clustering/all.seurat_object.log",
    script:
        "../scripts/clustr-dim-red.R"


rule plot_clustering_dim_reduction_plots:
    input:
        seurat_object="results/seurat/clustering/all.seurat_objt.rds",
    output:
        dim_plot=report(
            "results/plots/clustering/all.Dim-plot.pdf",
            caption="../report/clustering_dimplot.rst",
            category="Clustering",
            subcategory="global",
            labels={
                "plot": "Per sample Dim plot",
            },
        ),
    resources:
        cpus_per_task=20,
        mem_mb=94000,
        nodes=5,
    params:
        dims=config["clustering"]["dim"],
        resolution=config["clustering"]["resolution"],
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/seurat/clustering/all.clustering.log",
    script:
        "../scripts/plot-dimplot.R"


# Update 20250123: automatic subclustering
rule perform_subclustering:
    input:
        seurat_object="results/seurat/clustering/all.seurat_objt.rds",
    output:
        seurat_object="results/seurat/subclustering/all.seurat_objt.rds",
        subdim_plot="results/plots/subclustering/all.SubDim-plot.pdf",
    resources:
        cpus_per_task=20,
        mem_mb=94000,
        nodes=10,
        runtime=20,
    params:
        dims=config["clustering"]["dim"],  # Réutilise les paramètres de dimension
        resolution=config["clustering"]["resolution"],  # Réutilise la résolution
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/seurat/subclustering/all.subclustering.log",
    script:
        "../scripts/subclustering.R"
