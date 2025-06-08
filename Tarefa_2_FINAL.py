import pandas as pd
import numpy as np

# === 1.1 Leitura dos ficheiros ===

features = pd.read_csv("uniprotkb_s_cerevisiae.gff",
                       sep='\t', comment='#', header=None, low_memory=False)
features.columns = ['Entry', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes', 'extra']

seq_df = pd.read_csv(
    "uniprotkb_s_cerevisiae.tsv",
    sep="\t", usecols=["Entry", "Sequence"]
)

# === 1.2 Filtragem das fosfoserinas ===
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


# === 1.3. Análise das sequências e extração de serinas ===

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


# === 1.4. Construção do DataFrame final ===

def construir_dataframe(records):
    columns = ["entry", "pos", "known P"] + [str(i) for i in range(-10, 11)]
    return pd.DataFrame(records, columns=columns)


df_serinas = construir_dataframe(records)


# === 2.1. Separar os dois grupos ===

df_knownP = df_serinas[df_serinas["known P"] == True]
df_notP = df_serinas[df_serinas["known P"] == False]

# === 2.2. Função para calcular log₂(f / f_global) ===

def calcular_log2_frequencias(df):
    posicoes = [str(i) for i in range(-10, 11) if i != 0]
    aminoacidos = list("ACDEFGHIKLMNPQRSTVWY")

    # Contar todas as ocorrências (globais)
    aa_counts_global = pd.Series(0, index=aminoacidos, dtype=int)
    for pos in posicoes:
        aa_counts_global += df[pos].value_counts().reindex(aminoacidos, fill_value=0)

    total_global = aa_counts_global.sum()
    f_global = aa_counts_global / total_global

    # Contar por posição
    log2_matrix = pd.DataFrame(index=aminoacidos, columns=posicoes, dtype=float)

    for pos in posicoes:
        counts = df[pos].value_counts().reindex(aminoacidos, fill_value=0)
        total_pos = counts.sum()
        f_pos = counts / total_pos
        log2 = np.log2(f_pos / f_global)
        log2_matrix[pos] = log2

    return log2_matrix

# === 2.3. Aplicar para ambos os grupos ===

log2_knownP = calcular_log2_frequencias(df_knownP)
log2_notP = calcular_log2_frequencias(df_notP)

# === 2.4. Estilização opcional para Jupyter Notebook ===

def style_table(styler):
    styler.format(precision=2)
    styler.background_gradient(axis=None, vmin=-3, vmax=3, cmap="seismic")
    return styler

# Exemplo para visualizar estilizado (em Jupyter Notebook):
# log2_knownP.style.pipe(style_table)
# log2_notP.style.pipe(style_table)



# Teste para PyCharm

#Esta é a tabela de valores de log₂(f / f_global) para as posições relativas em torno das serinas fosforiláveis (ou seja, serinas anotadas como Phosphoserine).
print(log2_knownP)
print(log2_notP)
log2_knownP.to_csv("log2_knownP.csv")
log2_notP.to_csv("log2_notP.csv")
