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

cols = ['entry', 'pos', 'known_P'] + [str(i) for i in range(-10, 11)]
serine_df = pd.DataFrame(data, columns=cols)

# --- 4. Visualização e exportação ---

print(serine_df.head(10))
print(f"Total serinas: {len(serine_df)}")
print(f"Serinas fosforiláveis :{serine_df['known_P'].sum()}")
print(f"Serinas não fosforiláveis: {(~serine_df['known_P']).sum()}")


