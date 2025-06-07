import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Carregar a tabela da Tarefa 1 já criada (ex: "serinas_df.csv")
df = pd.read_csv("serinas_df.csv")  # Certifica-te de que guardaste a Tarefa 1 neste ficheiro

# Lista de posições -10 a +10 (exceto a posição 0, onde está a serina)
positions = [str(i) for i in range(-10, 11) if i != 0]
aminoacids = list('ACDEFGHIKLMNPQRSTVWY')  # Os 20 aminoácidos padrão

# Separar as serinas fosforiláveis e não fosforiláveis
fosforilaveis = df[df['known_P'] == True]
nao_fosforilaveis = df[df['known_P'] == False]

# Se necessário, corta aleatoriamente as serinas não fosforiláveis para balancear
nao_fosforilaveis = nao_fosforilaveis.sample(n=len(fosforilaveis), random_state=42)

# Função para calcular frequência por posição
def calcular_frequencias(df_subset):
    freq = pd.DataFrame(0, index=aminoacids, columns=positions)

    for pos in positions:
        valores = df_subset[pos].value_counts()
        for aa in aminoacids:
            freq.loc[aa, pos] = valores.get(aa, 0)

    # Converter para frequências relativas
    freq = freq.div(freq.sum(axis=0), axis=1)
    return freq

# Frequências por posição
freq_fosforilaveis = calcular_frequencias(fosforilaveis)
freq_nao_fosforilaveis = calcular_frequencias(nao_fosforilaveis)

# Frequência global dos aminoácidos (em todas as janelas e posições, exceto serina)
def calcular_frequencia_global(df_subset):
    todos = df_subset[positions].values.flatten()
    total = len(todos)
    f_global = pd.Series({aa: np.sum(todos == aa) / total for aa in aminoacids})
    return f_global

f_global_fosforilaveis = calcular_frequencia_global(fosforilaveis)
f_global_nao_fosforilaveis = calcular_frequencia_global(nao_fosforilaveis)

# Calcular log2(f / f_global) para cada posição
def calcular_log2(freq_por_posicao, f_global):
    return np.log2(freq_por_posicao.div(f_global, axis=0))

log2_fosforilaveis = calcular_log2(freq_fosforilaveis, f_global_fosforilaveis)
log2_nao_fosforilaveis = calcular_log2(freq_nao_fosforilaveis, f_global_nao_fosforilaveis)

# Mostrar as duas tabelas (log2 ratios)
print("Log2(f / f_global) para serinas fosforiláveis:")
print(log2_fosforilaveis.round(2))

print("\nLog2(f / f_global) para serinas não fosforiláveis:")
print(log2_nao_fosforilaveis.round(2))

# Se estiveres num Jupyter Notebook, aplica estilos:
def style_table(styler):
    styler.format(precision=2)
    styler.background_gradient(axis=None, vmin=-3, vmax=3, cmap="seismic")
    return styler

# Se quiseres guardar as tabelas:
log2_fosforilaveis.to_csv("log2_fosforilaveis.csv")
log2_nao_fosforilaveis.to_csv("log2_nao_fosforilaveis.csv")
