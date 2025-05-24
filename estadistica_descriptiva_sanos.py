import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# Lista de archivos CSV
files = [
    'gene_reads_adipose_subcutaneous.csv', 'gene_reads_bladder.csv',
    'gene_reads_esophagus_mucosa.csv', 'gene_reads_kidney_cortex.csv',
    'gene_reads_lung.csv', 'gene_reads_minor_salivary_gland.csv',
    'gene_reads_skin_not_sun_exposed_suprapubic.csv', 'gene_reads_skin_sun_exposed_lower_leg.csv'
]

# Ruta base
base_path = r"C:\Users\Personal\Documents\Saris\Universidad\Semestre 9\TG\Archivos GTEX"

# DataFrames finales
suma_totales = pd.DataFrame(columns=['Columna', 'Suma Total', 'Archivo'])
promedios_transcriptos = pd.DataFrame(columns=['Archivo', 'Transcripto', 'Promedio'])

# Procesar cada archivo
for file in files:
    file_path = os.path.join(base_path, file)

    # Leer el archivo CSV
    df = pd.read_csv(file_path, sep=None, engine="python")

    # Extraer los nombres de los transcriptos (columna 3)
    transcriptos = df.iloc[:, 2]

    # Seleccionar columnas desde la cuarta en adelante
    df_selected = df.iloc[:, 3:].apply(pd.to_numeric, errors='coerce')

    # Calcular suma por columna
    suma_df = pd.DataFrame({'Columna': df_selected.columns, 'Suma Total': df_selected.sum()})
    suma_df['Archivo'] = file
    suma_totales = pd.concat([suma_totales, suma_df], ignore_index=True)

    # Calcular promedio por transcripto
    promedios = df_selected.mean(axis=1)
    df_resultados = pd.DataFrame({'Transcripto': transcriptos, 'Promedio': promedios, 'Archivo': file})
    promedios_transcriptos = pd.concat([promedios_transcriptos, df_resultados], ignore_index=True)

# Ajuste del paso para mostrar menos etiquetas en el eje X
step_suma = max(1, len(suma_totales) // 20)  # Mostrar como máximo 20 etiquetas
step_promedio = max(1, len(promedios_transcriptos) // 20)  # Mostrar como máximo 20 etiquetas

# Gráfica 1: Suma total de lecturas por columna
plt.figure(figsize=(14, 6))
plt.scatter(range(len(suma_totales)), suma_totales['Suma Total'], marker='o', color='b', alpha=1, s=20)
plt.title('Suma total de lecturas por paciente')
plt.xlabel('Columna')
plt.ylabel('Suma Total de Lecturas')

# Mostrar solo algunas etiquetas en el eje X
plt.xticks(
    ticks=range(0, len(suma_totales), step_suma),  # Elegir etiquetas con intervalo
    labels=suma_totales['Columna'][::step_suma],
    rotation=90, fontsize=8
)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.subplots_adjust(bottom=0.25)  # Ajuste para etiquetas largas
plt.show()

# Gráfica 2: Promedio de lecturas por transcripto
plt.figure(figsize=(14, 6))
plt.scatter(range(len(promedios_transcriptos)), promedios_transcriptos['Promedio'], marker='o', color='r', alpha=1, s=20)
plt.title('Promedio de lecturas por transcripto')
plt.xlabel('Transcripto')
plt.ylabel('Promedio de lecturas')

# Mostrar solo algunas etiquetas en el eje X
plt.xticks(
    ticks=range(0, len(promedios_transcriptos), step_promedio),
    labels=promedios_transcriptos['Transcripto'][::step_promedio],
    rotation=90, fontsize=8
)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.subplots_adjust(bottom=0.25)  # Ajuste para etiquetas largas
plt.show()

# Crear histograma de suma total de pacientes
plt.figure(figsize=(12, 6))
plt.hist(suma_totales['Suma Total'], bins=20, edgecolor='black', alpha=0.7)
plt.title('Histograma del número de lecturas por paciente')
plt.xlabel('Total de lecturas')
plt.ylabel('Frecuencia absoluta')  
plt.grid(axis='y', alpha=0.75)
plt.axvline(suma_totales['Suma Total'].mean(), color='b', linestyle='dashed', linewidth=1, label='Media')
plt.gca().xaxis.set_major_formatter(mticker.ScalarFormatter())
plt.gca().yaxis.set_major_formatter(mticker.ScalarFormatter())
plt.ticklabel_format(style='plain', axis='both')
plt.legend()
plt.show()


# Histograma con transformación log2
plt.figure(figsize=(12, 6))
plt.hist(np.log2(1 + promedios_transcriptos['Promedio']), bins=20, edgecolor='black', alpha=0.7)
plt.title('Histograma de promedio de lecturas por transcripto (Transformación log₂)')
plt.xlabel('log₂(1 + Promedio por fila)')
plt.ylabel('Frecuencia absoluta')
plt.grid(axis='y', alpha=0.75)
plt.axvline(np.log2(1 + promedios_transcriptos['Promedio']).mean(), color='r', linestyle='dashed', linewidth=1, label='Media')
plt.legend()
plt.show()

# Función para calcular estadísticas descriptivas
def calcular_estadisticas(df, columna):
    return {
        'Promedio': df[columna].mean(),
        'Desviación_Estandar': df[columna].std(),
        'Mínimo': df[columna].min(),
        'Máximo': df[columna].max(),
        'Percentil_25': df[columna].quantile(0.25),
        'Percentil_50': df[columna].median(),
        'Percentil_75': df[columna].quantile(0.75)
    }

# Calcular estadísticas descriptivas
estadisticas_suma_df = pd.DataFrame(calcular_estadisticas(suma_totales, 'Suma Total'), index=[0])
estadisticas_promedios_df = pd.DataFrame(calcular_estadisticas(promedios_transcriptos, 'Promedio'), index=[0])

# Mostrar resultados
print("\nEstadísticas Descriptivas - Suma total por paciente:")
print(estadisticas_suma_df)
print("\nEstadísticas Descriptivas - Promedio de transcriptos:")
print(estadisticas_promedios_df)
