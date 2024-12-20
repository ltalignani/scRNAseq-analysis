rule perform_preprocessing:
    input:
        samples="config/samples.tsv",
    output:
        seurat_object= "results/seurat/preprocessing/all.seurat_objt.rds",
        intergrated_seurat_object="results/seurat/preprocessing/all.seurat_integration_objt.rds",
    resources:
        cpus_per_task=20,
        runtime = 20,
        mem_mb=94000,
        nodes=10
    params:
        path=config["resources"]["path"],
        clustering=config["clustering"]["activate"]
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/seurat/preprocessing/all.seurat_object.log",
    script:
        "../scripts/pre--processing.R"


rule plot_preprocessing_plots:
    input:
        seurat_object="results/seurat/preprocessing/all.seurat_objt.rds",
    output:
        QC_vln_plot=report("results/plots/preprocessing/all.QC-Vln-plot.pdf",
            caption="../report/qc_volin_plot.rst",
            category="QC",
            subcategory="global",
            labels={
                    "plot": "pre-processing Vln plot",
                },
        ),
        QC_nCount_nFeatures=report("results/plots/preprocessing/all.QC-nCount_nFeatures.pdf",
            caption="../report/qc_nC_nF_plot.rst",
            category="QC",
            subcategory="global",
            labels={
                    "plot": "pre-processing nCount nFeatures",
                },
        ),
        variable_features=report("results/plots/preprocessing/all.Highly_variable_features-plot.pdf",
            caption="../report/variable_feature_plot.rst",
            category="QC",
            subcategory="global",
            labels={
                    "plot": "Highly variable features plot",
                },
        ),
        elbow_plot=report("results/plots/preprocessing/all.Elbow-plot.pdf",
            caption="../report/elbow_plot.rst",
            category="QC",
            subcategory="global",
            labels={
                    "plot": "Elbow plot",
                },        
        )
    resources:
        cpus_per_task=20,
        mem_mb=94000,
        nodes=5
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/plots/seurat-preprocessing/all.QC-Vln-plot.log",
    script:
        "../scripts/Plot_preprocessing_plots.R"


