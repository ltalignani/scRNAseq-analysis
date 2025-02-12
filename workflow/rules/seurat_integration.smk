rule perform_merge_seurat:
    input:
        integration_seurat_object="results/seurat/preprocessing/all.seurat_integration_objt.rds",
        samples="config/samples.tsv",
    output:
        merge_seurat_object="results/seurat/integration/merge/{model}.seurat_objt_before_integration.rds",
    resources:
        runtime = 20,
    params:
        model=get_model_samples,
        groups=lambda wc: config["diffexp"]["models"][wc.model][
            "column_name"
        ],
        base_level=lambda wc: config["diffexp"]["models"][wc.model][
            "base_level"
        ],
        dims=config["merged_sample_clustering"]["dim"],
        resolution=config["merged_sample_clustering"]["resolution"],
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/seurat/integration/merge/{model}.seurat_object_before_integration.log",
    script:
        "../scripts/seurat_before_integration.R"


rule plot_dim_plots_for_merged_seurat:
    input:
        seurat_object="results/seurat/integration/merge/{model}.seurat_objt_before_integration.rds",
    output:
        merge_dim_plot=report("results/plots/seurat/integration/merge/clustering/{model}.Dim-plot.pdf",
            caption="../report/merged_seurat_dimplot.rst",
            category="Sample cluster plot before integration",
            subcategory="Dimension reduction plots",
            labels={
                "model": "{model}", "plot": "Dimension reduction plot",
                },  
        ),       
    resources:
        cpus_per_task=20,
        mem_mb=94000,
        nodes=5
    params:
        column_name=lambda wc: config["diffexp"]["models"][wc.model][
            "column_name"
        ],
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/seurat/plots/seurat/integration/merge/clustering/{model}.Dim-plot.log",
    script:
        "../scripts/plot-merge-dimplot.R"

rule seurat_integration:
    input:
        merge_seurat_object="results/seurat/integration/merge/{model}.seurat_objt_before_integration.rds",
    output:
        intergrated_seurat_object="results/seurat/integration/{model}.seurat_objt_integration.rds",
        top_10_markers="results/tables/seurat/integration/markers/{model}.top-10-markers-per-cluster.tsv",
    resources:
        runtime = 20,
    params:
        model=get_model_samples,
        groups=lambda wc: config["diffexp"]["models"][wc.model][
            "column_name"
        ],
        base_level=lambda wc: config["diffexp"]["models"][wc.model][
            "base_level"
        ],
        dims=config["merged_sample_clustering"]["dim"],
        resolution=config["merged_sample_clustering"]["resolution"],
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/seurat/integration/{model}.seurat_objt_integration.log",
    script:
        "../scripts/seurat-integration.R"

rule assign_celltype_to_integrated_seurat:
    input:
        integrated_seurat_object="results/seurat/integration/{model}.seurat_objt_integration.rds",
        file_celltype=config["celltype_annotation"]["path"],
        top_10_celltype_markers="results/tables/seurat/integration/markers/{model}.top-10-markers-per-cluster.tsv",
    output:
        dim_plot=report("results/plots/seurat/integration/per_celltype/{model}.Dim-plot.pdf",
            caption="../report/celltype_assign_dimplot.rst",
            category="Dimension reduction plots of assigned celltypes",
            subcategory="Dimension reduction plots",
            labels={
                   "model": "{model}", "plot": "Dimension reduction plot",
                },  
        ),
        celltype_seurat_object="results/seurat/integration/celltype/{model}.celltype_integration_seurat_objt.rds",
        top_markers_plot="results/plots/seurat/integration/per_celltype/{model}.top-markers-plot.pdf",
    resources:
        cpus_per_task=20,
        mem_mb=94000,
        nodes=10
    conda:
        "../envs/seurat.yaml"
    log:
        "logs/seurat/integration/celltype/{model}.assign_celltype_to_integrated_seurat.log",
    script:
        "../scripts/plot__celltype__dimplot.R"

rule per_celltype_diffexp:
    input:
        celltype_seurat_object="results/seurat/integration/celltype/{model}.celltype_integration_seurat_objt.rds",
    output:
        all_markers=report("results/tables/seurat/integration/per_celltype/{model}.{celltype}.celltype-diff-exp-genes.tsv",
                caption="../report/celltype_all_markers.rst",
                category=" Per celltype differential expression ",
                subcategory="Differential expression tables",
                labels={
                   "model": "{model}", "celltype": "{celltype}", "table": "differentially expressed genes",
                }, 
        ),
        top_10_markers=report("results/tables/seurat/integration/per_celltype/{model}.{celltype}.top-10-celltype-markers.tsv",
                caption="../report/celltype_top_10_markers.rst",
                category=" Per celltype differential expression ",
                subcategory="Differential expression tables",
                labels={
                   "model": "{model}", "celltype": "{celltype}", "table": "top 10 markers",
                }, 
        
        ),       
        volcano_plot=report("results/plots/seurat/integration/per_celltype/volcano_plots/{model}.{celltype}.celltype-volcano_plot.pdf",
                caption="../report/celltype_volcano_plot.rst",
                category=" Per celltype differential expression ",
                subcategory="Per celltype volcano plots",
                labels={
                   "model": "{model}", "celltype": "{celltype}", "plot": "Volcano plots",
                },  
        ),
        log2fc_heatmap=report(
            "results/plots/seurat/integration/per_celltype/heatmaps/{model}.{celltype}.celltype-log2fc-heatmap.pdf",
            caption="../report/celltype_log2fc_heatmap.rst",
            category=" Per celltype differential expression ",
            subcategory="Log2FC Heatmaps",
            labels={
                "model": "{model}",
                "celltype": "{celltype}",
                "plot": "Log2FC heatmap",
            },
        ),
    resources:
        cpus_per_task=20,
        mem_mb=94000,
        nodes=10,
        runtime = 20,
    params:
        model=config["diffexp"]["models"],
        column_name=lambda wc: config["diffexp"]["models"][wc.model][
            "column_name"
        ],
        base_level=lambda wc: config["diffexp"]["models"][wc.model][
            "base_level"
        ],
        comparison_variable=lambda wc: config["diffexp"]["models"][wc.model][
            "comparison_variable"
        ],
        cell_type=config["cell_type_diff_exp"]["cell_type"],
    conda:
        "../envs/diff_exp.yaml"

    log:
        "logs/seurat/integration/per_celltype_diffexp/{model}.{celltype}.per_celltype_diffexp.log",
    script:
        "../scripts/celltype--diffexp.R"




# TO do plot, as this is not a mandatory output, need to fix the issue of failing

#rule plot_integration_plots:
#    input:
#        intergrated_seurat_object="results/seurat/integration/{model}.seurat_objt_integration.rds",
#    output:
#        dim_plot="results/plots/seurat_integration/integration/{model}.Dim-plot.pdf",
#    resources:
#        cpus_per_task=20,
#        mem_mb=94000,
#        nodes=10
#    conda:
#        "../envs/seurat.yaml"
#    log:
#        "logs/plots/seurat_intergration/integration/{model}.plot_integration_plots.log",
#    script:
#        "../scripts/plot-integration-plots.R"
