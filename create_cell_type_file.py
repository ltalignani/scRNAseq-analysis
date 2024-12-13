import pandas as pd

# Remplacez le chemin ci-dessous par le chemin du fichier que vous avez téléchargé manuellement
input_file = "PanglaoDB_markers_27_Mar_2020.tsv"

# Charger le fichier dans un DataFrame en le décompressant
data = pd.read_csv(input_file, sep="\t")

# Filtrer pour ne garder que les lignes avec des gènes de Mus musculus
data = data[data['species'].str.contains('Mm')]

# Garder uniquement les colonnes d'intérêt : 'official gene symbol' et 'cell type'
data = data[['official gene symbol', 'cell type']]

# Renommer les colonnes pour correspondre au format souhaité
data.columns = ['Gene', 'Cell Type']

# Enregistrer le fichier en format CSV
output_file = "Mus_musculus_cell_types.csv"
data.to_csv(output_file, sep="\t", index=False)

print(f"Fichier '{output_file}' créé avec succès.")
