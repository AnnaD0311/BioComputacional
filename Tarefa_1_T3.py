import pandas as pd
import numpy as np
import re
from tqdm import tqdm

# --- 1. Ler ficheiros da UniProt ---

# Lê apenas a coluna das sequências
seq_df = pd.read_csv("uniprotkb_s_cerevisiae.tsv", sep='\t', usecols=['Entry', 'Sequence'])

# Lê o ficheiro GFF com as anotações
features = pd.read_csv("uniprotkb_s_cerevisiae.gff", sep='\t', comment='#', header=None)
features.columns = ['Entry', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes', 'extra']

# Filtrar apenas as fosfoserinas
phospho = features[(features['type'] == 'Modified residue') & features['attributes'].str.contains("Phosphoserine")]

# Converter posições para inteiros
phospho['start'] = phospho['start'].astype(int)

# Cria um dicionário com entradas de fosfoserinas
phospho_map = phospho.groupby('Entry')['start'].apply(set).to_dict()

# --- 2. Extrair todas as serinas e contexto (-10, +10) ---

data = []

# Percorrer cada sequência de proteína
for _, row in tqdm(seq_df.iterrows(), total=len(seq_df)):
    entry = row['Entry']
    seq = row['Sequence']
    phospho_pos = phospho_map.get(entry, set())  # Posições fosforiladas

    for i, aa in enumerate(seq):
        if aa == 'S':
            pos = i + 1  # posição biológica (começa em 1)

            # Extrair janela de -10 a +10
            window = []
            for offset in range(-10, 11):
                idx = i + offset
                if 0 <= idx < len(seq):
                    window.append(seq[idx])
                else:
                    window.append('X')  # fora dos limites -> X

            # Verifica se a posição está anotada como fosfoserina
            known_p = pos in phospho_pos

            data.append([entry, pos, known_p] + window)

# --- 3. Montar DataFrame final ---

cols = ['entry', 'pos', 'known_P'] + [str(i) for i in range(-10, 11)]
serine_df = pd.DataFrame(data, columns=cols)

# --- 4. Visualização e exportação inicial ---

print(serine_df.head(10))  # exemplo
print(f"Total serinas: {len(serine_df)}")
print(f"Serinas fosforiladas: {serine_df['known_P'].sum()}")
print(f"Serinas não fosforiladas: {(~serine_df['known_P']).sum()}")

# Salvar se quiser
serine_df.to_csv("serinas_contexto.csv", index=False)
