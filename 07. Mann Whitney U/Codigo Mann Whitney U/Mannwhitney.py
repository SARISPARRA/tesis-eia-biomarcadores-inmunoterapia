import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

# Definir las rutas
path_cancer = r"C:\Users\Personal\Documents\Saris\Universidad\Semestre 9\TG\ReadCounts\EXCEL"
path_noresponden = r"C:\Users\Personal\Documents\Saris\Universidad\Semestre 9\TG\Carolina"
path_responden = r"C:\Users\Personal\Documents\Saris\Universidad\Semestre 9\TG\Carolina"
path_sanos = r"C:\Users\Personal\Documents\Saris\Universidad\Semestre 9\TG\Archivos GTEX"

# Cargar archivos
df_cancer = pd.read_csv(os.path.join(path_cancer, "_normalized_counts_cancer.csv"), index_col=0)
df_noresponden = pd.read_csv(os.path.join(path_noresponden, "_normalized_counts_NoResponden.csv"), index_col=0)
df_responden = pd.read_csv(os.path.join(path_responden, "_normalized_counts_Responden.csv"), index_col=0)
df_sanos = pd.read_csv(os.path.join(path_sanos, "_normalized_counts_sanos.csv"), index_col=0)

# Genes de interés
genes_of_interest = [
    "SFTPC", "DPYSL5", "SLC6A12", "CSRP3", "KLRC1", "PKD1L2", "FAT4", "CXCL9",
    "ARL4C", "SST", "CPNE8", "CTAGE15", "RNASEH2A", "PAPLN", "BICC1", "ZNF649",
    "TVP23B", "TRMT1L", "BIRC2", "QRSL1"
]

# Filtrar genes de interés
def filter_genes(df, genes):
    return df[df.index.isin(genes)]

df_cancer = filter_genes(df_cancer, genes_of_interest)
df_noresponden = filter_genes(df_noresponden, genes_of_interest)
df_responden = filter_genes(df_responden, genes_of_interest)
df_sanos = filter_genes(df_sanos, genes_of_interest)

# -------------------------
# Crear DataFrame largo para visualización
# -------------------------
def melt_group(df, group_name):
    df_melted = df.T.melt(var_name="Gene", value_name="RPM")
    df_melted["Group"] = group_name
    return df_melted

df_long = pd.concat([
    melt_group(df_cancer, "Cancer"),
    melt_group(df_noresponden, "NoResponden"),
    melt_group(df_responden, "Responden"),
    melt_group(df_sanos, "Sanos")
])

# -------------------------
# Mann-Whitney U tests
# -------------------------
p_values = []
for gene in genes_of_interest:
    cancer_data = df_cancer.loc[gene].values
    noresponden_data = df_noresponden.loc[gene].values
    responden_data = df_responden.loc[gene].values
    sanos_data = df_sanos.loc[gene].values

    # Comparaciones entre pares
    _, p1 = mannwhitneyu(cancer_data, noresponden_data, alternative='two-sided')
    _, p2 = mannwhitneyu(cancer_data, responden_data, alternative='two-sided')
    _, p3 = mannwhitneyu(cancer_data, sanos_data, alternative='two-sided')
    _, p4 = mannwhitneyu(noresponden_data, responden_data, alternative='two-sided')
    _, p5 = mannwhitneyu(noresponden_data, sanos_data, alternative='two-sided')
    _, p6 = mannwhitneyu(responden_data, sanos_data, alternative='two-sided')

    p_values.append({
        'Gene': gene,
        'Cancer vs NoResponden': p1,
        'Cancer vs Responden': p2,
        'Cancer vs Sanos': p3,
        'NoResponden vs Responden': p4,
        'NoResponden vs Sanos': p5,
        'Responden vs Sanos': p6
    })

df_pvalues = pd.DataFrame(p_values)
df_pvalues.to_csv("p_values_results.csv", index=False)

# -------------------------
# Crear boxplots mejorados
# -------------------------
sns.set(style="whitegrid")

num_genes = len(genes_of_interest)
ncols = 5
nrows = int(np.ceil(num_genes / ncols))

fig, axes = plt.subplots(nrows, ncols, figsize=(20, 16))
axes = axes.flatten()

for idx, gene in enumerate(genes_of_interest):
    ax = axes[idx]
    subset = df_long[df_long["Gene"] == gene]
    
    sns.boxplot(x="Group", y="RPM", data=subset, ax=ax, palette="Set2", showfliers=False)
    sns.stripplot(x="Group", y="RPM", data=subset, ax=ax, color='black', size=2, jitter=True, alpha=0.5)
    
    ax.set_title(gene, fontsize=10)
    ax.set_xlabel("")
    ax.set_ylabel("RPM", fontsize=8)
    ax.tick_params(axis='x', labelrotation=45, labelsize=8)
    ax.tick_params(axis='y', labelsize=8)

# Eliminar ejes vacíos si hay
for j in range(idx + 1, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.savefig("boxplots_genes_mejorados.png", dpi=300)
plt.close()

