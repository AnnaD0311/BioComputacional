import pandas as pd

# Lê o GFF
features = pd.read_csv(
    "uniprotkb_s_cerevisiae.gff",
    delimiter="\t",
    comment="#",
    names=["Entry", "source", "type", "start", "end", "score", "strand", "phase", "attributes"],
    dtype=str,
    engine="python"
)

print("Total linhas no ficheiro:", len(features))

# Verifica tipos únicos
print("Tipos únicos:", features["type"].unique())

# Verifica se 'Modified residue' está presente
print("Modified residue count:", (features["type"] == "Modified residue").sum())

# Verifica se há 'phospho'
phospho_test = features[features["attributes"].str.contains("phospho", case=False, na=False)]
print("Total com 'phospho':", len(phospho_test))
print(phospho_test.head())
