import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="whitegrid", context="talk")

features = pd.read_csv("uniprotkb_s_cerevisiae.gff", sep='\t', comment='#', header=None, low_memory=False)
features.columns = ['Entry', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes', 'extra']

seq_df = pd.read_csv("uniprotkb_s_cerevisiae.tsv", sep="\t", usecols=["Entry", "Sequence"])

def filtrar_fosfoserinas(features):
    fosfo = features[
        (features["type"] == "Modified residue") &
        (features["attributes"].str.contains("Phosphoserine", case=False, na=False))
    ]
    fosfo_dict = {}
    for _, row in fosfo.iterrows():
        entry = row["Entry"]
        pos = int(row["start"])
        fosfo_dict.setdefault(entry, set()).add(pos)
    return fosfo_dict

fosfo_dict = filtrar_fosfoserinas(features)

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
                window = [seq[i + offset] if 0 <= i + offset < len(seq) else 'X' for offset in range(-10, 11)]
                records.append([entry, pos_1based, known_P] + window)
    return records

records = extrair_serinas(seq_df, fosfo_dict)

def construir_dataframe(records):
    columns = ["entry", "pos", "known P"] + [str(i) for i in range(-10, 11)]
    return pd.DataFrame(records, columns=columns)

df_serinas = construir_dataframe(records)

df_knownP = df_serinas[df_serinas["known P"] == True]
df_notP = df_serinas[df_serinas["known P"] == False]

def calcular_log2_frequencias(df):
    posicoes = [str(i) for i in range(-10, 11) if i != 0]
    aminoacidos = list("ACDEFGHIKLMNPQRSTVWY")

    aa_counts_global = pd.Series(0, index=aminoacidos, dtype=int)
    for pos in posicoes:
        aa_counts_global += df[pos].value_counts().reindex(aminoacidos, fill_value=0)

    total_global = aa_counts_global.sum()
    f_global = aa_counts_global / total_global

    log2_matrix = pd.DataFrame(index=aminoacidos, columns=posicoes, dtype=float)
    for pos in posicoes:
        counts = df[pos].value_counts().reindex(aminoacidos, fill_value=0)
        f_pos = counts / counts.sum()
        log2_matrix[pos] = np.log2(f_pos / f_global)

    return log2_matrix

log2_knownP = calcular_log2_frequencias(df_knownP)
log2_notP = calcular_log2_frequencias(df_notP)

plt.figure(figsize=(20, 8))
sns.heatmap(log2_knownP, cmap="seismic", center=0, cbar_kws={"label": "log₂(f / f_global)"})
plt.title("Serinas Fosforiláveis - log₂(f / f_global)")
plt.xlabel("Posição relativa à Serina")
plt.ylabel("Aminoácido")
plt.show()

plt.figure(figsize=(20, 8))
sns.heatmap(log2_notP, cmap="seismic", center=0, cbar_kws={"label": "log₂(f / f_global)"})
plt.title("Serinas Não Fosforiláveis - log₂(f / f_global)")
plt.xlabel("Posição relativa à Serina")
plt.ylabel("Aminoácido")
plt.show()


import matplotlib.pyplot as plt
import seaborn as sns

def plotar_heatmap(df_log2, titulo):
    plt.figure(figsize=(14, 7))
    sns.heatmap(
        df_log2,
        cmap="seismic",
        center=0,
        vmin=-3,
        vmax=3,
        annot=True,
        fmt=".2f",
        linewidths=0.1,
        cbar_kws={'label': 'log₂(f / f_global)'}
    )
    plt.title(titulo, fontsize=14)
    plt.xlabel("Posição relativa à Serina")
    plt.ylabel("Aminoácido")
    plt.tight_layout()
    plt.show()

# Plot para serinas fosforiláveis
plotar_heatmap(log2_knownP, "Serinas Fosforiláveis - log₂(f / f_global)")

# Plot para serinas não fosforiláveis
plotar_heatmap(log2_notP, "Serinas Não Fosforiláveis - log₂(f / f_global)")


#Tarefa 3
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

# 1. Separar serinas fosforiláveis e não-fosforiláveis
df_knownP = df_serinas[df_serinas["known P"] == True]
df_notP = df_serinas[df_serinas["known P"] == False]

# 2. Reduzir aleatoriamente o número de serinas não-fosforiláveis
np.random.seed(10)
indexes2drop = np.random.choice(df_notP.index, 253000, replace=False)
df_notP_reduced = df_notP.drop(indexes2drop)

# 3. Juntar os dois subconjuntos para obter a tabela equilibrada
serine_table_balanced = pd.concat([df_knownP, df_notP_reduced])

# 4. Separar features (X) e labels (Y)
X = serine_table_balanced[[str(i) for i in range(-10, 11)]]
Y = serine_table_balanced["known P"]

# 5. Train-test split
X_train, X_test, Y_train, Y_test = train_test_split(
    X, Y, test_size=0.3, stratify=Y, random_state=42
)

# 6. Função para calcular log₂(f / f_global)
def calcular_log2_freq(df):
    aa_list = list("ACDEFGHIKLMNPQRSTVWY")
    posicoes = [str(i) for i in range(-10, 11) if i != 0]  # excluir a posição central 0

    # Frequência global (soma de todas as posições)
    aa_counts_global = pd.Series(0, index=aa_list, dtype=int)
    for pos in posicoes:
        aa_counts_global += df[pos].value_counts().reindex(aa_list, fill_value=0)
    f_global = aa_counts_global / aa_counts_global.sum()
    f_global = f_global.replace(0, 1e-10)  # evitar log2(0)

    # Frequência por posição e log₂(f / f_global)
    log2_matrix = pd.DataFrame(index=aa_list, columns=posicoes, dtype=float)
    for pos in posicoes:
        counts = df[pos].value_counts().reindex(aa_list, fill_value=0)
        f_pos = counts / counts.sum()
        f_pos = f_pos.replace(0, 1e-10)
        log2_matrix[pos] = np.log2(f_pos / f_global)
    return log2_matrix

# 7. Criar tabelas de log2 para cada classe
log2P = calcular_log2_freq(X_train[Y_train == True])
log2notP = calcular_log2_freq(X_train[Y_train == False])

# 8. Função para codificar X com base nas tabelas
def codificar_df(df, Y_vals, tabela_P, tabela_notP):
    codificados = []
    posicoes = [str(i) for i in range(-10, 11) if i != 0]
    aminoacidos_validos = set(tabela_P.index)

    for idx, row in df.iterrows():
        tabela = tabela_P if Y_vals.loc[idx] else tabela_notP
        linha = []
        for pos in posicoes:
            aa = row[pos]
            if aa in aminoacidos_validos:
                linha.append(tabela.loc[aa, pos])
            else:
                linha.append(0.0)  # fallback
        codificados.append(linha)
    
    return pd.DataFrame(codificados, index=df.index, columns=posicoes)

# 9. Codificar os conjuntos de treino e teste
X_train_encoded = codificar_df(X_train, Y_train, log2P, log2notP)
X_test_encoded = codificar_df(X_test, Y_test, log2P, log2notP)

# 10. Verificações finais
print("Shape treino codificado:", X_train_encoded.shape)
print("Shape teste codificado:", X_test_encoded.shape)
print("Distribuição treino:", Y_train.value_counts())
print("Distribuição teste:", Y_test.value_counts())