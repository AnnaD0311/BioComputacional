import pandas as pd
import numpy as np

# --- 1. Ler ficheiros da UniProt ---

# Lê apenas a coluna das sequências
seq_df = pd.read_csv("uniprotkb_s_cerevisiae.tsv", sep='\t', usecols=['Entry', 'Sequence'])

# Lê o ficheiro GFF com as anotações
features = pd.read_csv("uniprotkb_s_cerevisiae.gff", sep='\t', comment='#', header=None)
features.columns = ['Entry', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes', 'extra']

# Garantir que 'attributes' é string antes de usar .str
features['attributes'] = features['attributes'].astype(str)

# Filtrar apenas as fosfoserinas
phospho = features[(features['type'] == 'Modified residue') & features['attributes'].str.contains("Phosphoserine")]

# Converter posições para inteiros
phospho = phospho.copy()
phospho['start'] = phospho['start'].astype(int)

# Criar dicionário com entradas de fosfoserinas
phospho_map = phospho.groupby('Entry')['start'].apply(set).to_dict()

# --- Verificação opcional: garantir que as posições anotadas são serinas ---
for entry, positions in phospho_map.items():
    seq_row = seq_df[seq_df['Entry'] == entry]
    if seq_row.empty:
        print(f"Aviso: sequência para {entry} não encontrada no seq_df")
        continue
    seq = seq_row['Sequence'].values[0]
    for p in positions:
        if p <= 0 or p > len(seq):
            print(f"Aviso: posição {p} fora dos limites na proteína {entry}")
        elif seq[p - 1] != 'S':
            print(f"Aviso: posição {p} na proteína {entry} não é uma serina (é '{seq[p - 1]}')")

# --- 2. Extrair todas as serinas e contexto (-10, +10) ---

data = []


def extrair_janela(seq, pos, window_size=10):
    seq_len = len(seq)
    window = []
    for offset in range(-window_size, window_size + 1):
        idx = pos + offset
        if 0 <= idx < seq_len:
            window.append(seq[idx])
        else:
            window.append('X')  # fora dos limites
    return window


for entry, seq in zip(seq_df['Entry'], seq_df['Sequence']):
    phospho_pos = phospho_map.get(entry, set())  # posições anotadas como fosfoserina

    seq_array = np.array(list(seq))
    serina_indices = np.where(seq_array == 'S')[0]

    for idx in serina_indices:
        pos = idx + 1  # posição 1-based
        window = extrair_janela(seq, idx)
        known_p = pos in phospho_pos
        data.append([entry, pos, known_p] + window)

# --- 3. Montar DataFrame final ---

cols = ['entry', 'pos', 'known_P'] + list(range(-10, 11))
serine_df = pd.DataFrame(data, columns=cols)

# --- 4. Visualização e exportação ---

print(serine_df.head(10))
print(f"Total serinas: {len(serine_df)}")
print(f"Serinas fosforiláveis :{serine_df['known_P'].sum()}")
print(f"Serinas não fosforiláveis: {(~serine_df['known_P']).sum()}")

serine_df = pd.DataFrame(data, columns=cols)

# Garante que known_P é booleana
serine_df['known_P'] = serine_df['known_P'].astype(bool)

# <-- AQUI colocas a linha:
pos_cols = [i for i in range(-10, 11) if i != 0]


# análises da Tarefa 2
# Função para calcular log2(f / f_global) com prevenção de log2(0)
def compute_log2_ratios(df, pos_cols):
    aa_list = list("ACDEFGHIKLMNPQRSTVWY")

    # Frequência global
    all_counts = df[pos_cols].values.flatten()
    all_counts = all_counts[all_counts != 'X']  # excluir fora dos limites
    global_freq = pd.Series(all_counts).value_counts(normalize=True)

    # Garante que todos os 20 aa estão presentes
    for aa in aa_list:
        if aa not in global_freq:
            global_freq[aa] = 0.00001  # evita log2(0)

    global_freq = global_freq[aa_list]  # ordenado

    # Frequências por posição
    pos_freqs = {}
    for col in pos_cols:
        col_counts = df[col].value_counts(normalize=True)
        freqs = []
        for aa in aa_list:
            f = col_counts.get(aa, 0.00001)  # evita zero
            f_global = global_freq[aa]
            ratio = np.log2(f / f_global)
            freqs.append(ratio)
        pos_freqs[col] = freqs

    log2_df = pd.DataFrame(pos_freqs, index=aa_list)
    return log2_df


# Separar os dois conjuntos
df_known = serine_df[serine_df['known_P'] == True]
df_unknown = serine_df[serine_df['known_P'] == False]

# Calcular log2(f / f_global) para os dois grupos
log2_knownP = compute_log2_ratios(df_known, pos_cols)
log2_unknownP = compute_log2_ratios(df_unknown, pos_cols)


# Estilizar para Jupyter (opcional)
def style_table(styler):
    styler.format(precision=2)
    styler.background_gradient(axis=None, vmin=-3, vmax=3, cmap="seismic")
    return styler

# Se estiveres num Jupyter Notebook, descomenta:
# display(log2_knownP.style.pipe(style_table))
# display(log2_unknownP.style.pipe(style_table))

# Se quiseres salvar os resultados como CSV:
# log2_knownP.to_csv("log2_fosforiladas.csv")
# log2_unknownP.to_csv("log2_nao_fosforiladas.csv")