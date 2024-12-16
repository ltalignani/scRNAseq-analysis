# Initialisation des logs
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

print("Début de l'exécution du script")

# Chargement et mise à jour des packages nécessaires
print("Chargement des packages nécessaires")
required_packages <- c(
    "AnnotationDbi", "org.Mm.eg.db", "ggplot2", "forcats", "cowplot",
    "enrichplot", "clusterProfiler", "ReactomePA", "dplyr", "aplot"
)
for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
        print(paste("Installation du package :", pkg))
    }
    library(pkg, character.only = TRUE)
}

# Chargement des données d'expression
print("Chargement des données d'expression différentielle")
expression_file <- read.csv(snakemake@input[["diff_exp_genes"]], sep = ",", header = TRUE)

expression_file$gene <- expression_file$X

# Annotation des gènes
print("Annotation des gènes avec ENSEMBL et ENTREZID")
expression_file$ENSEMBL <- mapIds(org.Mm.eg.db,
    keys = expression_file$gene,
    column = "ENSEMBL",
    keytype = "SYMBOL",
    multiVals = "first"
)
expression_file$ENTREZID <- mapIds(org.Mm.eg.db,
    keys = expression_file$ENSEMBL,
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
)

# Paramètres
print("Chargement des paramètres depuis Snakemake")
gene_ontology_p_val <- snakemake@params[["gene_ontology_p_val"]]
gene_ontology_q_val <- snakemake@params[["gene_ontology_q_val"]]
pathway_significant <- snakemake@params[["pathway_significant"]]

# Filtrage des gènes significatifs
print("Filtrage des gènes significatifs")
ensembl_list <- expression_file %>%
    filter(abs(avg_log2FC) > 0.1 & p_val < 0.05) %>%
    select(ENSEMBL)

# Enrichissement GO
print("Enrichissement Gene Ontology (GO)")
enrich_go <- function(ont) {
    print(paste("Analyse d'enrichissement GO pour", ont))
    tryCatch(
        {
            enrichGO(
                gene = ensembl_list$ENSEMBL,
                OrgDb = org.Mm.eg.db,
                keyType = "ENSEMBL",
                ont = ont,
                pAdjustMethod = "BH",
                pvalueCutoff = gene_ontology_p_val,
                qvalueCutoff = gene_ontology_q_val
            )
        },
        error = function(e) {
            message(paste("Erreur lors de l'enrichissement GO pour", ont, ":", e$message))
            NULL
        }
    )
}

ego2_cc <- enrich_go("CC")
ego2_bp <- enrich_go("BP")
ego2_mf <- enrich_go("MF")

# Sauvegarde des tableaux de résultats GO
print("Sauvegarde des tableaux de résultats GO")
go_tables <- list(CC = ego2_cc, BP = ego2_bp, MF = ego2_mf)
lapply(names(go_tables), function(ont) {
    if (!is.null(go_tables[[ont]])) {
        write.csv(data.frame(go_tables[[ont]]), file = snakemake@output[[paste0("GO_", ont)]], quote = FALSE)
        print(paste("Résultats enregistrés pour GO", ont))
    }
})

# Analyse GSEA Pathway
print("Préparation des données pour GSEA Pathway")
gene_entrezid_fc <- expression_file %>%
    filter(p_val < 0.05) %>%
    select(ENTREZID, avg_log2FC)

genelist <- setNames(gene_entrezid_fc$avg_log2FC, as.character(gene_entrezid_fc$ENTREZID))
genelist <- sort(genelist, decreasing = TRUE)
genelist <- genelist[!duplicated(names(genelist))]

print("Analyse GSEA Pathway")
pathway <- tryCatch(
    {
        gsePathway(genelist,
            nPermSimple = 10000,
            pvalueCutoff = pathway_significant,
            pAdjustMethod = "BH",
            verbose = FALSE, organism = "mouse"
        )
    },
    error = function(e) {
        message("Erreur dans l'analyse GSEA Pathway :", e$message)
        NULL
    }
)

if (!is.null(pathway)) {
    pathway_table <- data.frame(pathway)
    write.csv(pathway_table, file = snakemake@output[["pathway"]], quote = FALSE)
    print("Résultats GSEA Pathway enregistrés")
} else {
    pathway_table <- NULL
    print("Aucun résultat significatif pour GSEA Pathway")
}

# Barplot pour GO : Cellulaire, Biologique, et Molecular Function
create_barplot <- function(result, output_path) {
    if (!is.null(result) && nrow(result) > 0) {
        print(paste("Création d'un Barplot pour", output_path))

        # Calcul du score qscore (log10 de l'ajustement de p)
        result <- mutate(result, qscore = -log10(p.adjust))

        # Création du barplot
        bar_plot <- barplot(result, x = "qscore") +
            theme(
                axis.text.x = element_text(size = 10), # Réduction de la taille du texte de l'axe X
                axis.text.y = element_text(size = 8), # Réduction de la taille du texte de l'axe Y
                plot.title = element_text(size = 14) # Taille du titre
            )

        # Sauvegarde du barplot en PDF
        pdf(file = output_path)
        print(bar_plot)
        dev.off()
        print(paste("Barplot enregistré pour", output_path))
    } else {
        # Si aucun résultat, créer un PDF vide avec un message
        pdf(file = output_path)
        plot.new()
        text(0.5, 0.5, "No significant results", cex = 1.5, font = 2)
        dev.off()
        print(paste("Aucun résultat pour Barplot", output_path))
    }
}

# Appel des fonctions pour chaque catégorie (CC, BP, MF)
create_barplot(ego2_cc, snakemake@output[["bar_plot_cc"]])
create_barplot(ego2_bp, snakemake@output[["bar_plot_bp"]])
create_barplot(ego2_mf, snakemake@output[["bar_plot_mf"]])


# Visualisation : Dotplots avec texte réduit sur l'axe des ordonnées
print("Création des visualisations : Dotplots")
create_plot <- function(result, output_path, plot_type = "dotplot", showCategory = 30) {
    print(paste("Création d'un", plot_type, "pour", output_path))
    if (!is.null(result) && nrow(result) > 0) {
        if (plot_type == "dotplot") {
            p <- dotplot(result, showCategory = min(showCategory, nrow(result))) +
                ggtitle(paste("Dotplot for", result@ontology)) +
                theme(
                    axis.text.y = element_text(size = 5), # Réduction de la taille du texte de l'axe Y
                    plot.title = element_text(size = 12) # Taille du titre
                )
        } else if (plot_type == "barplot") {
            result <- mutate(result, qscore = -log10(p.adjust))
            p <- barplot(result, x = "qscore") +
                theme(
                    axis.text.x = element_text(size = 10), # Réduction de la taille du texte de l'axe X
                    axis.text.y = element_text(size = 5), # Réduction de la taille du texte de l'axe Y
                    plot.title = element_text(size = 12) # Taille du titre
                )
        }
        pdf(file = output_path)
        print(p)
        dev.off()
        print(paste(plot_type, "enregistré pour", output_path))
    } else {
        pdf(file = output_path)
        plot.new()
        text(0.5, 0.5, "No significant results", cex = 1.5, font = 2)
        dev.off()
        print(paste("Aucun résultat pour", plot_type, output_path))
    }
}

# Appel de la fonction
create_plot(ego2_cc, snakemake@output[["dot_plot_cc"]], "dotplot")
create_plot(ego2_bp, snakemake@output[["dot_plot_bp"]], "dotplot")
create_plot(ego2_mf, snakemake@output[["dot_plot_mf"]], "dotplot")

# Treeplot
print("Création des Treeplots")
create_treeplot <- function(result, output_path, nCluster = 4) {
    if (!is.null(result) && nrow(result) > 0) {
        print(paste("Création d'un Treeplot pour", output_path))
        tryCatch({
            # Calcul des similarités entre termes
            result <- pairwise_termsim(result)
            
            # Ajuster dynamiquement le nombre de clusters si nécessaire
            if (nrow(result) < nCluster) {
                nCluster <- nrow(result)
                print(paste("Nombre de clusters ajusté à", nCluster, "en raison de résultats limités."))
            }

            # Générer les treeplots avec la méthode recommandée
            p1 <- treeplot(result, cluster.params = list(method = "complete"), nCluster = nCluster)
            p2 <- treeplot(result, cluster.params = list(method = "average"), nCluster = nCluster)

            # Combiner les deux treeplots dans un seul
            tree_plot <- aplot::plot_list(p1, p2, tag_levels = "A")

            # Sauvegarder le plot
            pdf(file = output_path, height = 20, width = 30)
            print(tree_plot)
            dev.off()
            print(paste("Treeplot enregistré pour", output_path))
        }, error = function(e) {
            # Gestion des erreurs
            print(paste("Erreur dans la création du Treeplot :", e$message))
            pdf(file = output_path)
            plot.new()
            text(0.5, 0.5, paste("Error:", e$message), cex = 1.5, font = 2)
            dev.off()
        })
    } else {
        # Cas où les résultats sont vides ou null
        pdf(file = output_path)
        plot.new()
        text(0.5, 0.5, "No significant results", cex = 1.5, font = 2)
        dev.off()
        print(paste("Aucun résultat pour Treeplot :", output_path))
    }
}

# Appel de la fonction avec les différents résultats enrichis
create_treeplot(ego2_cc, snakemake@output[["tree_plot_cc"]])
create_treeplot(ego2_bp, snakemake@output[["tree_plot_bp"]])
create_treeplot(ego2_mf, snakemake@output[["tree_plot_mf"]])

print("Création des Enrichment Maps avec matrice de similarité")

create_enrichment_map <- function(result, output_path, showCategory = 50) {
    if (!is.null(result) && nrow(result) > 0) {
        print(paste("Création d'une Enrichment Map pour", output_path))
        
        # Calcul de la matrice de similarité des termes
        result <- tryCatch(
            {
                pairwise_termsim(result)
            },
            error = function(e) {
                message("Erreur lors du calcul de la matrice de similarité :", e$message)
                return(NULL)
            }
        )

        # Vérification que le calcul de la matrice de similarité a réussi
        if (!is.null(result)) {
            # Gestion des erreurs spécifiques à la génération des cartes
            tryCatch(
                {
                    p1 <- emapplot(
                        result, 
                        cex.params = list(category_node = 1.5), 
                        layout.params = list(layout = "kk"), 
                        showCategory = showCategory
                    )
                    p2 <- emapplot(
                        result, 
                        cex.params = list(category_node = 1.5), 
                        layout.params = list(layout = "fr"), 
                        showCategory = showCategory
                    )
                    enrichment_map <- cowplot::plot_grid(p1, p2, ncol = 2, labels = LETTERS[1:2])
                    
                    # Sauvegarde du plot
                    pdf(file = output_path, height = 20, width = 30)
                    print(enrichment_map)
                    dev.off()
                    print(paste("Enrichment Map enregistrée pour", output_path))
                },
                error = function(e) {
                    message("Erreur lors de la création de l'Enrichment Map :", e$message)
                    
                    # Création d'un PDF vide en cas d'erreur
                    pdf(file = output_path)
                    plot.new()
                    text(0.5, 0.5, "Error generating Enrichment Map", cex = 1.5, font = 2)
                    dev.off()
                    print(paste("Erreur lors de la création de l'Enrichment Map :", output_path))
                }
            )
        } else {
            # Fichier PDF vide si la matrice de similarité n'a pas pu être générée
            pdf(file = output_path)
            plot.new()
            text(0.5, 0.5, "No similarity matrix generated", cex = 1.5, font = 2)
            dev.off()
            print(paste("Aucun résultat pour Enrichment Map", output_path))
        }
    } else {
        # Fichier PDF vide si aucun résultat significatif
        pdf(file = output_path)
        plot.new()
        text(0.5, 0.5, "No significant results", cex = 1.5, font = 2)
        dev.off()
        print(paste("Aucun résultat pour Enrichment Map", output_path))
    }
}

# Appel des fonctions pour chaque catégorie
create_enrichment_map(ego2_cc, snakemake@output[["enrichment_map_cc"]])
create_enrichment_map(ego2_bp, snakemake@output[["enrichment_map_bp"]])
create_enrichment_map(ego2_mf, snakemake@output[["enrichment_map_mf"]])


# Visualisation : UpSetPlot corrigée
print("Création des UpSetPlots")
create_upset_plot <- function(result, output_path) {
    if (!is.null(result) && nrow(result) > 0) {
        print(paste("Création d'un UpSetPlot pour", output_path))
        
        # Augmenter les marges pour éviter la troncature
        pdf(file = output_path, height = 10, width = 15)
        upset_plot <- upsetplot(result) +
            theme(
                text = element_text(size = 12),          # Taille globale du texte (améliorée)
                axis.text.x = element_text(size = 8),   # Taille des textes des catégories (axe X)
                axis.text.y = element_text(size = 6),   # Taille des textes des catégories (axe Y)
                axis.title.y = element_text(size = 12),  # Taille du titre de l'axe des ordonnées
                plot.margin = margin(1, 1, 1, 8, "cm")   # Augmentation des marges, notamment à gauche
            ) +
            labs(
                y = "Interaction Size",                 # Ajouter le label pour l'axe des ordonnées
                x = NULL                                # Supprime le label de l'axe des abscisses
            )
        
        # Affichage et sauvegarde
        print(upset_plot)
        dev.off()
        print(paste("UpSetPlot enregistré pour", output_path))
    } else {
        pdf(file = output_path)
        plot.new()
        text(0.5, 0.5, "No significant results", cex = 1.5, font = 2)
        dev.off()
        print(paste("Aucun résultat pour UpSetPlot", output_path))
    }
}

# Appel de la fonction
create_upset_plot(ego2_cc, snakemake@output[["upset_plot_cc"]])
create_upset_plot(ego2_bp, snakemake@output[["upset_plot_bp"]])
create_upset_plot(ego2_mf, snakemake@output[["upset_plot_mf"]])

# Pathway
if (nrow(pathway_table) > 30) {
    dotplot_pathway <- dotplot(pathway, showCategory = 30) +
        ggtitle("dotplot for Reactome Pathway")
    pdf(file = snakemake@output[["dotplot_pathway"]], height = 20, width = 20)
    dotplot_pathway
} else if (nrow(pathway_table) > 1) {
    dotplot_pathway <- dotplot(pathway, showCategory = nrow(pathway_table)) +
        ggtitle("dotplot for Reactome Pathway")
    pdf(file = snakemake@output[["dotplot_pathway"]], height = 20, width = 20)
    dotplot_pathway
} else {
    pdf(file = snakemake@output[["dotplot_pathway"]], height = 20, width = 20)
    dotplot_pathway <- "no significant term enriched"
    plot.new()
    text(.5, .5, dotplot_pathway, font = 2, cex = 1.5)
}


if (nrow(pathway_table) > 1) {
    arrange_pathway <- dplyr::arrange(pathway, abs(NES)) %>%
        dplyr::group_by(sign(NES))
    pdf(file = snakemake@output[["bar_plot_pathway"]])
    ggplot(arrange_pathway, aes(NES, fct_reorder(Description, NES), fill = qvalue), showCategory = 10) +
        geom_col(orientation = "y") +
        scale_fill_continuous(low = "red", high = "blue", guide = guide_colorbar(reverse = TRUE)) +
        theme_minimal() +
        ylab(NULL)
} else {
    pdf(file = snakemake@output[["bar_plot_pathway"]])
    txt <- "no significant term enriched"
    plot.new()
    text(.5, .5, txt, font = 2, cex = 1.5)
}

# Fermeture des logs
print("Fin de l'exécution du script")
sink()
sink(type = "message")
# Fermeture des logs
print("Fin de l'exécution du script")
sink()
sink(type = "message")

