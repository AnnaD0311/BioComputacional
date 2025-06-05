import pandas as pd
import numpy as np

# === Leitura dos ficheiros ===

# Lê ficheiro de sequências (.tsv), apenas a coluna Entry e Sequence
sequences_df = pd.read_csv("uniprotkb_s_cerevisiae.tsv", sep="\t", usecols=["Entry", "Sequence"])

# Lê ficheiro GFF com anotações estruturais e funcionais
features = pd.read_csv(
    "uniprotkb_s_cerevisiae.gff",
    sep="\t",
    comment="#",
    names=["Entry", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
)

# === Filtra apenas as fosfoserinas ===
# Garante que a coluna 'attributes' é string (para evitar erros com .str)
features["attributes"] = features["attributes"].astype(str)

phospho = features[
    (features["type"] == "Modified residue") &
    features["attributes"].str.contains("Phosphoserine", case=False)
].copy()

phospho["start"] = phospho["start"].astype(int)

# === Tarefa 1: Criar janelas -10 a +10 em torno de cada serina ===

data = []
window_range = list(range(-10, 11))
window_range.remove(0)  # Posição 0 é a serina, que não vai para a tabela

for _, row in sequences_df.iterrows():
    entry = row["Entry"]
    seq = row["Sequence"]
    phospho_pos = phospho[phospho["Entry"] == entry]["start"].tolist()

    for i, aa in enumerate(seq):
        if aa != "S":
            continue

        window = []
        for offset in window_range:
            idx = i + offset
            if 0 <= idx < len(seq):
                window.append(seq[idx])
            else:
                window.append("X")

        data.append({
            "entry": entry,
            "pos": i + 1,  # Numeração 1-based
            "known_P": (i + 1) in phospho_pos,
            **{str(offset): aa for offset, aa in zip(window_range, window)}
        })

# DataFrame final com janelas e labels
final_df = pd.DataFrame(data)

# === Tarefa 2: Frequências e enriquecimento log2(f / f_global) ===

def compute_log2_enrichment(df, known=True):
    # Subconjunto: só serinas fosforiladas ou não-fosforiladas
    subset = df[df["known_P"] == known]
    pos_cols = [str(i) for i in range(-10, 11) if i != 0]
    aa20 = list("ACDEFGHIKLMNPQRSTVWY")

    # Frequência global (em todas as posições combinadas)
    all_aa = subset[pos_cols].values.flatten()
    global_counts = pd.Series(all_aa).value_counts()
    f_global = global_counts.reindex(aa20, fill_value=0) / global_counts.sum()

    # Frequências por posição
    freq_by_pos = pd.DataFrame(index=aa20, columns=pos_cols)

    for col in pos_cols:
        counts = subset[col].value_counts()
        if counts.sum() == 0:
            freq_by_pos[col] = [0 for aa in aa20]
        else:
            freq_by_pos[col] = [counts.get(aa, 0) / counts.sum() for aa in aa20]

    # log₂(f / f_global) com pequena constante para evitar log(0)
    enrichment = np.log2(freq_by_pos.divide(f_global + 1e-10, axis=0))

    return enrichment

# Calcula enriquecimentos
log2_knownP = compute_log2_enrichment(final_df, known=True)
log2_unknownP = compute_log2_enrichment(final_df, known=False)

# === Estilo para Jupyter Notebook (opcional) ===

def style_table(styler):
    styler.format(precision=2)
    styler.background_gradient(axis=None, vmin=-3, vmax=3, cmap="seismic")
    return styler

# Exemplo de visualização no notebook:
# log2_knownP.style.pipe(style_table)
# log2_unknownP.style.pipe(style_table)

# Também pode guardar para Excel/CSV se necessário:
# log2_knownP.to_csv("log2_knownP.csv")
# log2_unknownP.to_csv("log2_unknownP.csv")


#Para verificar o resultado é dar : 
# Mostrar as primeiras 5 linhas da tabela final com todas as serinas
#print("Primeiras serinas processadas:")
#print(final_df.head())

# Mostrar uma amostra do enriquecimento log2 para serinas fosforiladas
#print("\nLog2 enrichment (fosforiladas):")
#print(log2_knownP.round(2).iloc[:, :5])  # mostra só as primeiras 5 posições

# Mostrar uma amostra do enriquecimento log2 para serinas não fosforiladas
#print("\nLog2 enrichment (não fosforiladas):")
#print(log2_unknownP.round(2).iloc[:, :5])  # idem
 
 # Testar se existem linhas com "Phosphoserine"
phospho_test = features[features["attributes"].str.contains("Phosphoserine", case=False)]
print(f"Total de linhas com 'Phosphoserine': {len(phospho_test)}")
print(phospho_test.head())

