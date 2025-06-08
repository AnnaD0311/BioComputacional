#Passo 1: Ler os ficheiros
import pandas as pd

# Lê o ficheiro .tsv com as sequências das proteínas
seqs = pd.read_csv("uniprotkb_s_cerevisiae.tsv", sep="\t", usecols=["Entry", "Sequence"])

# Lê o ficheiro .gff
features = pd.read_csv("uniprotkb_s_cerevisiae.gff", sep="\t", comment="#", header=None)
features.columns = ["Entry", "source", "type", "start", "end", "score", "strand", "phase", "attributes", "extra"]



# Passo 2: Extrair apenas as fosfoserinas do GFF
fosfoserinas = features[
    (features["type"] == "Modified residue") &
    (features["attributes"].str.contains("Phosphoserine"))
]



#  Passo 3: Criar o DataFrame com as serinas e janelas -10 a +10
def extrair_janela(seq, pos, tamanho=10):
    janela = []
    for offset in range(-tamanho, tamanho + 1):
        if offset == 0:
            continue  # Ignorar a posição da serina
        idx = pos + offset
        if 0 <= idx < len(seq):
            janela.append(seq[idx])
        else:
            janela.append("X")
    return janela

# Montagem do DataFrame com todas as serinas de todas as proteínas:
linhas = []

for _, row in seqs.iterrows():
    entry = row["Entry"]
    seq = row["Sequence"]
    for i, aa in enumerate(seq):
        if aa == "S":  # Serina
            pos_uniprot = i + 1  # Uniprot começa em 1
            janela = extrair_janela(seq, i)
            knownP = ((fosfoserinas["Entry"] == entry) & (fosfoserinas["start"] == pos_uniprot)).any()
            linhas.append([entry, pos_uniprot, knownP] + janela)

colunas = ["entry", "pos", "known P"] + [str(i) for i in range(-10, 0)] + [str(i) for i in range(1, 11)]
serinas_df = pd.DataFrame(linhas, columns=colunas)



#  Passo 4: Criar dados de treino balanceados
# Separar serinas fosforiladas e não fosforiladas
positivas = serinas_df[serinas_df["known P"] == True]  # Serinas com known P = True (fosforiladas)
negativas = serinas_df[serinas_df["known P"] == False]  # Serinas com known P = False (não fosforiladas)

# Amostragem aleatória para balancear os negativos com as positivas
negativas_amostradas = negativas.sample(n=len(positivas), random_state=42)  # Amostra o mesmo número de negativos que positivos

# Combina as serinas fosforiladas e as serinas negativas amostradas
dados_treino = pd.concat([positivas, negativas_amostradas])

# Embaralhar os dados (baralhar as linhas do DataFrame)
dados_treino = dados_treino.sample(frac=1, random_state=42)  # frac=1 significa "embaralhar todas as linhas"
