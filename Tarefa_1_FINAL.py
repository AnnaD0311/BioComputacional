import pandas as pd

# === 1. Leitura dos ficheiros ===

features = pd.read_csv("uniprotkb_s_cerevisiae.gff",
                       sep='\t', comment='#', header=None, low_memory=False)
features.columns = ['Entry', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes', 'extra']

seq_df = pd.read_csv(
    "uniprotkb_s_cerevisiae.tsv",
    sep="\t", usecols=["Entry", "Sequence"]
)


# === 2. Filtragem das fosfoserinas ===

def filtrar_fosfoserinas(features):
    fosfo = features[
        (features["type"] == "Modified residue") &
        (features["attributes"].str.contains("Phosphoserine", case=False, na=False))
        ]

    fosfo_dict = {}
    for _, row in fosfo.iterrows():
        entry = row["Entry"]
        pos = int(row["start"])  # posições 1-based
        fosfo_dict.setdefault(entry, set()).add(pos)

    return fosfo_dict


fosfo_dict = filtrar_fosfoserinas(features)


# === 3. Análise das sequências e extração de serinas ===

def extrair_serinas(seq_df, fosfo_dict):
    records = []

    for _, row in seq_df.iterrows():
        entry = row["Entry"]
        seq = row["Sequence"]
        fosfo_pos = fosfo_dict.get(entry, set())

        for i, aa in enumerate(seq):
            if aa == 'S':
                pos_1based = i + 1
                known_P = pos_1based in fosfo_pos

                window = []
                for offset in range(-10, 11):
                    idx = i + offset
                    window.append(seq[idx] if 0 <= idx < len(seq) else 'X')

                records.append([entry, pos_1based, known_P] + window)

    return records


records = extrair_serinas(seq_df, fosfo_dict)


# === 4. Construção do DataFrame final ===

def construir_dataframe(records):
    columns = ["entry", "pos", "known P"] + [str(i) for i in range(-10, 11)]
    return pd.DataFrame(records, columns=columns)


df_serinas = construir_dataframe(records)

# === 5. Exportação para CSV ===

df_serinas.to_csv("serinas_fosforilaveis.csv", index=False)

# === 6. Estatísticas ===

total = len(df_serinas)
fosforilaveis = df_serinas["known P"].sum()
nao_fosforilaveis = total - fosforilaveis

print(f"Total de serinas analisadas: {total}")
print(f"Serinas fosforiláveis (anotadas como Phosphoserine): {fosforilaveis}")
print(f"Serinas não fosforiláveis: {nao_fosforilaveis}")

# === 7. Mostrar primeiras 9 linhas ===

print(df_serinas.iloc[:9])
