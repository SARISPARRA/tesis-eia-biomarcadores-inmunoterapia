import os
import pandas as pd
import numpy as np
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests

# === 1. Cargar los archivos ===
# Ajusta estos nombres según los archivos que tengas

path_cancer = r"C:\Users\Personal\Documents\Saris\Universidad\Semestre 9\TG\ReadCounts\EXCEL"
path_noresponden = r"C:\Users\Personal\Documents\Saris\Universidad\Semestre 9\TG\Carolina"
path_responden = r"C:\Users\Personal\Documents\Saris\Universidad\Semestre 9\TG\Carolina"
path_sanos = r"C:\Users\Personal\Documents\Saris\Universidad\Semestre 9\TG\Archivos GTEX"

# Cargar cada grupo (ajusta los nombres exactos de archivo)
df_cancer_norm = pd.read_csv(os.path.join(path_cancer, "_normalized_counts_cancer.csv"), index_col=0)
df_noresponden_norm = pd.read_csv(os.path.join(path_noresponden, "_normalized_counts_NoResponden.csv"), index_col=0)
df_responden_norm = pd.read_csv(os.path.join(path_responden, "_normalized_counts_Responden.csv"), index_col=0)
df_sanos_norm = pd.read_csv(os.path.join(path_sanos, "_normalized_counts_sanos.csv"), index_col=0)

# --- LIMPIEZA DE ÍNDICES ---
def clean_index(df):
    df.index = df.index.str.strip()  # elimina espacios
    df = df[~df.index.duplicated(keep='first')]  # elimina genes duplicados
    return df

df_sanos_norm = clean_index(df_sanos_norm)
df_cancer_norm = clean_index(df_cancer_norm)
df_responden_norm = clean_index(df_responden_norm)
df_noresponden_norm = clean_index(df_noresponden_norm)

# --- UNIFICAR POR GENES COMUNES ---
common_genes = set(df_sanos_norm.index) & set(df_cancer_norm.index) & set(df_responden_norm.index) & set(df_noresponden_norm.index)

# Filtrar todos para que tengan los mismos genes (en el mismo orden)
df_sanos_norm = df_sanos_norm.loc[sorted(common_genes)]
df_cancer_norm = df_cancer_norm.loc[sorted(common_genes)]
df_responden_norm = df_responden_norm.loc[sorted(common_genes)]
df_noresponden_norm = df_noresponden_norm.loc[sorted(common_genes)]

# === UNIR TODAS LAS MUESTRAS EN UN SOLO DATAFRAME ===
all_data = pd.concat([
    df_sanos_norm,
    df_cancer_norm,
    df_responden_norm,
    df_noresponden_norm
], axis=1)

# === CREAR LISTA DE GRUPOS ===
groups = (
    ['Sano'] * df_sanos_norm.shape[1] +
    ['Cancer'] * df_cancer_norm.shape[1] +
    ['Responden'] * df_responden_norm.shape[1] +
    ['NoResponden'] * df_noresponden_norm.shape[1]
)

# === PRUEBA DE KRUSKAL–WALLIS ===
pvals = []

for gene in all_data.index:
    try:
        values = [
            all_data.loc[gene, np.array(groups) == g].values
            for g in ['Sano', 'Cancer', 'Responden', 'NoResponden']
        ]
        stat, p = kruskal(*values)
    except Exception:
        p = np.nan
    pvals.append(p)

# === AJUSTE POR FDR ===
pvals = np.array(pvals)
padj = multipletests(pvals, method='fdr_bh')[1]

# === GUARDAR RESULTADOS ===
result_df = pd.DataFrame({
    'Geneid': all_data.index,
    'p_value': pvals,
    'padj': padj
})
result_df.to_csv("kruskal_wallis_resultados.csv", index=False)

print("✅ Análisis Kruskal–Wallis terminado.")
print("📁 Resultados guardados en: kruskal_wallis_resultados.csv")
