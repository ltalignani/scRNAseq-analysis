# Resultats produits par le pipeline

## Plots

### Celltype
  
- <Wt>.Dim-plot.pdf: Dimension reduction plot, contenant l'identification du type cellulaire de chaque cluster, de l'échantillon sauvage  
- <ko>.Dm-plot.pdf: Dimension reduction plot, contenant l'identification du type cellulaire de chaque cluster, de l'échantillon traité
  
  
### Clustering
  
- all.Dim.plot.pdf: Dimension reduction plot, contenant l'identification du numéro de chaque cluster, de l'échantillon sauvage et de l'échantillon traité, sur deux graphes différents dans le même fichier  
  
  
### Diffexp
  
- <Wt>.Heatmap-plot.pdf: Heatmap représentant l'expression différentielle des *top markers* au sein de chaque cluster, pour l'échantilonn sauvage  
- <ko>.Heatmap-plot.pdf: Heatmap représentant l'expression différentielle des *top markers* au sein de chaque cluster, pour l'échantilonn traité  
  
  
### Preprocessing
  
- all.Elbow-plot.pdf: ce graphe trace l'écart type de chaque composante principale pour une identification simple d'un "coude" sur le graphe. Ce coude permet de relever rapidement le seuil représentant le nombre de dimensions significatives, c'est à dire présentant une proportion de variance expliquée suffisante pour être considérée dans les analyses en aval.
- all.Highly_variable_features-plot.pdf:identifie les gènes outliers présentant une variation *cell-to-cell* élevée (hautement exprimés dans certaines cellules et faiblement exprimées dans d'autres) sur un plot de variabilité moyenne. Par défaut, ce graphe retourne 2000 gènes par dataset
- all.QC-Vln-plot.pdf: QC metrics, regroupant le nombre de gènes (nFeatures_RNA) par cellule (barcode), le nombre de molécules (nCount_RNA) par cellule, le pourcentage de gènes issus de mitochondries (^MT-) 

### Seurat / integration
  
  
#### per_celltype
  
- wt_vs_ko.Dim-plot.pdf: Dimensional reduction plot, présentant les clusters de cellules ainsi que le type cellulaire qui leur est associé. L'un des graphes les plus importants.  
- wt_vs_ko.top-markers-plot.pdf: table présentant les top markers par type cellulaire.  

##### heatmaps
  
- wt_vs_ko.<cell_type>.celltype-log2fc-heatmap.pdf: Heatmap du Log2 Fold-change des top markers par type cellulaire.

  
##### volcano_plots

- wt_vs_ko.<cell_type>.celltype-volcano_plot.pdf: 


#### enrichment (GO Terms)
  
0. classification des *GO Terms*
*GO Terms*:
- bp: Biological Process
- cc: Cellular Components
- mf: Molecular Function
<cell_type>: liste des types cellulaires identifiés lors de l'étape d'assignation

1. celltype: GOs par type cellulaire. Le nombre de fichier dépend du nombre de types cellulaires identifiés  
   - wt_vs_ko.<cell_type>.bar_plot_<bp|cc|mf>.pdf: visualisation des termes enrichis (*GO Terms*) sous forme de barres pour les BP|CC|MF. Il représente les scores d'enrichissement (*p-values* ajustées par des méthodes statistiques, permettant de contrôler le False Discovery Rate). Les barres en rouge classées en haut sont les plus significatives et les plus basses en bleu les moins significatives.  
   - wt_vs_ko.<cell_type>.dot_plot_<bp|cc|mf>.pdf: visualisation similaire au bar-plot mais incluant le score de taille de points, représentant le nombre de gènes impliqués dans chaque *GO Term*.  
   - wt_vs_ko.<cell_type>.tree_plot_<bp|cc|mf>.pdf: Clustering hiérarchique des *GO Terms*. Il repose sur les similitudes par paires de *GO Terms* calculés par la fonction pairwise_termsim(), qui utilise par défaut l'indice de similarité (JC) de Jaccard (ratio du nombre d'éléments communs à deux entités sur le nombre d'éléments uniques dans chaque entité).  
   - wt_vs_ko.<cell_type>.upset_plot_<bp|cc|mf>.pdf: visualisation des associations complexes entre les *GO Terms*. Les barres représentent le nombre de gènes impliqués dans l'interaction.  
   - wt_vs_ko.<cell_type>.enrichment_map_<bp|cc|mf>.pdf: La carte d'enrichissement organise les termes enrichis en un réseau avec des arêtes reliant des ensembles de gènes qui se chevauchent. De cette façon, les ensembles de gènes qui se chevauchent ont tendance à se regrouper, ce qui facilite l'identification du module fonctionnel.  


2. groups: *GO Terms* après intégration et fusion. 
   - wt_vs_ko.bar_plot_<bp|cc|mf>.pdf: visualisation des termes enrichis (*GO Terms*) sous forme de barres pour les BP|CC|MF. Il représente les scores d'enrichissement (*p-values* ajustées par des méthodes statistiques, permettant de contrôler le False Discovery Rate). Les barres en rouge classées en haut sont les plus significatives et les plus basses en bleu les moins significatives.  
   - wt_vs_ko.dot_plot_<bp|cc|mf>.pdf: visualisation similaire au bar-plot mail incluant le score de taille de points, représentant le nombre de gènes impliqués dans chaque *GO Term*.  
   - wt_vs_ko.tree_plot_<bp|cc|mf>.pdf: Clustering hiérarchique des *GO Terms*. Il repose sur les similitudes par paires de *GO Terms* calculés par la fonction pairwise_termsim(), qui utilise par défaut l'indice de similarité (JC) de Jaccard (ratio du nombre d'éléments communs à deux entités sur le nombre d'éléments uniques dans chaque entité).  
   - wt_vs_ko.upset_plot_<bp|cc|mf>.pdf: visualisation des associations complexes entre les *GO Terms*. Les barres représentent le nombre de gènes impliqués dans l'interaction.  
   - wt_vs_ko.enrichment_map_<bp|cc|mf>.pdf: La carte d'enrichissement organise les termes enrichis en un réseau avec des arêtes reliant des ensembles de gènes qui se chevauchent. De cette façon, les ensembles de gènes qui se chevauchent ont tendance à se regrouper, ce qui facilite l'identification du module fonctionnel.  
  

#### heatmaps
  
- wt_vs_ko.heatmap_plots.pdf
  

#### merge / clustering
  
- wt_vs_ko.Dim_plots.pdf: Dimensional reduction plot des échantillons wt et ko, après fusion et avant clustering.
  

#### pathway (Gene Set Enrichment pathway)
  
1. celltype: GOs par type cellulaire. Le nombre de fichier dépend du nombre de types cellulaires identifiés  
   - wt_vs_ko.<cell_type>.bar_plot_<bp|cc|mf>.pdf: visualisation des termes enrichis (*GSE Analysis*) sous forme de barres pour les BP|CC|MF. Il représente les scores d'enrichissement (*p-values* ajustées par des méthodes statistiques, permettant de contrôler le False Discovery Rate). Les barres en rouge classées en haut sont les plus significatives et les plus basses en bleu les moins significatives.  
   - wt_vs_ko.<cell_type>.dot_plot_<bp|cc|mf>.pdf: visualisation similaire au bar-plot mail incluant le score de taille de points, représentant le nombre de gènes impliqués dans chaque **GSE Analysis*.  

2. groups:*GSE Analysis* après intégration et fusion. 
   - wt_vs_ko.bar_plot_<bp|cc|mf>.pdf: visualisation des termes enrichis (*GSE Analysis*) sous forme de barres pour les BP|CC|MF. Il représente les scores d'enrichissement (*p-values* ajustées par des méthodes statistiques, permettant de contrôler le False Discovery Rate). Les barres en rouge classées en haut sont les plus significatives et les plus basses en bleu les moins significatives.  
   - wt_vs_ko.dot_plot_<bp|cc|mf>.pdf: visualisation similaire au bar-plot mail incluant le score de taille de points, représentant le nombre de gènes impliqués dans chaque *GSE Analysis*.  


-----------------------------------------------------------------------

## reports
  
- dag.pdf|png: Directed Acyclic Graph, représentant le déroulement du pipeline. Les bulles hachurées indiquent que la règle a bien été effectuée.  
- filegraph.pdf|png: graphe acyclique, représentant le déroulement du pipeline et indiquant quels sont les fichiers utilisés en input et produits en sortie, ainsi que leur emplacement.  
- files_summary.txt: fichier résumant le lien d'accès des fichiers produits, la date, la règle snakemake ayant produit le fichier, l'emplacement du fichier de log associé, le status du fichier et s'il a été mis à jour.  
- rulegraph.pdf|png: graphe acyclique, représentant le déroulement du pipeline.  

-----------------------------------------------------------------------

## seurat
  
Regroupe tous les objets Seurat, au format .rds, produits par les différentes étapes du pipeline:

### clustering

- all.seurat_obj.rds
  
  
### preprocessing
  
- all.seurat_integration_objt.rds
- all.seurat_obj.rds
  

### integration
  
- wt_vs_ko.seurat_obj_integration.rds

#### celltype
  
- wt_vs_ko.celltype_integration_seurat_obj.rds
  
#### merge

- wt_vs_ko.seurat_objt_before_integration.rds
  
-----------------------------------------------------------------------

## tables
  
Tables au format .tsv, produites par les différentes étapes du pipeline, servant à réaliser les graphes:

### diffexp
  
- <wt>.diff-exp-genes.tsv
- <ko>.diff-exp-genes.tsv
- <wt>.top-10-markers.tsv
- <ko>.top-10-markers.tsv

### seurat / integration
  
- wt_vs_ko.diff-exp-genes.tsv: GeneID,p_val,avg_log2FC,pct.1,pct.2,p_val_adj
- wt_vs_ko.top-10-markers.tsv: GeneID,p_val,avg_log2FC,pct.1,pct.2,p_val_adj
  

#### per_celltype
  
- wt_vs_ko.<cell_type>.diff-exp-genes.tsv: GeneID,p_val,avg_log2FC,pct.1,pct.2,p_val_adj
- wt_vs_ko.<cell_type>.top-10-markers.tsv: GeneID,p_val,avg_log2FC,pct.1,pct.2,p_val_adj
  

#### enrichment
  
1. celltype:
   - wt_vs_ko.<cell_type>.top_biological_process.tsv: table des processus biologiques, contenant l'ID (*GO Term accession number*),la Description, le GeneRatio,le BgRatio,la pvalue,la p.adjust,la qvalue,les geneID impliqués, et leur nombre (Count).  
   - wt_vs_ko.<cell_type>.top_cellular_components.tsv: table des composants cellulaires, contenant l'ID (*GO Term accession number*),la Description, le GeneRatio,le BgRatio,la pvalue,la p.adjust,la qvalue,les geneID impliqués, et leur nombre (Count). 
   - wt_vs_ko.<cell_type>.top_molecular_function.tsv: table des fonctions moléculaires, contenant l'ID (*GO Term accession number*),la Description, le GeneRatio,le BgRatio,la pvalue,la p.adjust,la qvalue,les geneID impliqués, et leur nombre (Count). 

2. groups:
   - wt_vs_ko.top_biological_process.tsv: table des processus biologiques, contenant l'ID (*GO Term accession number*),la Description, le GeneRatio,le BgRatio,la pvalue,la p.adjust,la qvalue,les geneID impliqués, et leur nombre (Count).  
   - wt_vs_ko.top_cellular_components.tsv: table des composants cellulaires, contenant l'ID (*GO Term accession number*),la Description, le GeneRatio,le BgRatio,la pvalue,la p.adjust,la qvalue,les geneID impliqués, et leur nombre (Count). 
   - wt_vs_ko.top_molecular_function.tsv: table des fonctions moléculaires, contenant l'ID (*GO Term accession number*),la Description, le GeneRatio,le BgRatio,la pvalue,la p.adjust,la qvalue,les geneID impliqués, et leur nombre (Count). 

#### markers
  
- wt_vs_ko.top-10-markers-per-cluster.tsv: table des top 10 markers par cluster, contenant la p_value (p_val), le Log2 Fold-change (avg_log2FC), la valeur de composante principale 1 (pct.1), la valeur de composante principale 2 (pct.2), la p-value ajustée (p_val_adj), le numéro du cluster attribué (cluster), et le nom du gène (gene)


#### pathway (Gene Set Enrichment Analysis)
  
1. celltype:
   - wt_vs_ko.<cell_type>.pathway_results.tsv: ID,Description,setSize,enrichmentScore,NES,pvalue,p.adjust,qvalue
2. groups:
   - wt_vs_ko.pathway_results.tsv: ID,Description,setSize,enrichmentScore,NES,pvalue,p.adjust,qvalue
