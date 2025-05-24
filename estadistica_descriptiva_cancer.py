import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Ruta base
base_path = r"C:\Users\Personal\Documents\Saris\Universidad\Semestre 9\TG\ReadCounts\EXCEL"

# Lista de archivos
files = [
    'readCounts_Wei_etal.xlsx', 'readCounts_Santosh_etal.xlsx', 
    'readCounts_Przybyl_etal.xlsx', 'readCounts_Li_etal.xlsx', 
    'readCounts_Le_etal.xlsx', 'readCounts_Alldredge_etal.xlsx',
    'readCounts_Zhang_etal.xlsx','readCounts_Barrasetal.xlsx'
]

# DataFrames finales
suma_totales = pd.DataFrame(columns=['Columna', 'Suma Total', 'Archivo'])
promedios_filas = pd.DataFrame(columns=['Archivo', 'Transcripto', 'Promedio'])

# Procesar cada archivo
for file in files:
    file_path = os.path.join(base_path, file)
    
    # Leer archivo Excel
    df = pd.read_excel(file_path)
    df_selected = df.iloc[:, 6:].apply(pd.to_numeric, errors='coerce')  # Seleccionar columnas desde la séptima en adelante
    transcriptos = df.iloc[:, 0]  # Extraer los nombres de los transcriptos
    
    # Calcular suma por columna
    suma_df = pd.DataFrame({'Columna': df_selected.columns, 'Suma Total': df_selected.sum()})
    suma_df['Archivo'] = file
    suma_totales = pd.concat([suma_totales, suma_df], ignore_index=True)
    
    # Calcular promedio por transcripto
    promedios = df_selected.mean(axis=1)
    df_resultados = pd.DataFrame({'Transcripto': transcriptos, 'Promedio': promedios, 'Archivo': file})
    promedios_filas = pd.concat([promedios_filas, df_resultados], ignore_index=True)
    
    

# Ajuste del paso para mostrar menos etiquetas en el eje X
step_suma = max(1, len(suma_totales) // 20)  # Mostrar como máximo 20 etiquetas
step_promedio = max(1, len(promedios_filas) // 20)  # Mostrar como máximo 20 etiquetas

# Gráfica 1: Suma total de lecturas por columna
plt.figure(figsize=(14, 6))
plt.scatter(range(len(suma_totales)), suma_totales['Suma Total'], marker='o', color='b', alpha=1, s=20)
plt.title('Suma total de lecturas por paciente')
plt.xlabel('Columna')
plt.ylabel('Suma Total de Lecturas')
plt.xticks(
    ticks=range(0, len(suma_totales), step_suma),
    labels=suma_totales['Columna'][::step_suma],
    rotation=90, fontsize=8
)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.subplots_adjust(bottom=0.25)
plt.show()

# Gráfica 2: Promedio de lecturas por transcripto
plt.figure(figsize=(14, 6))
plt.scatter(range(len(promedios_filas)), promedios_filas['Promedio'], marker='o', color='r', alpha=1, s=20)
plt.title('Promedio de lecturas por transcripto (sin filtrar)')
plt.xlabel('Transcripto')
plt.ylabel('Promedio de lecturas')
plt.xticks(
    ticks=range(0, len(promedios_filas), step_promedio),
    labels=promedios_filas['Transcripto'][::step_promedio],
    rotation=90, fontsize=8
)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.subplots_adjust(bottom=0.25)
plt.show()

# Filtrar los promedios menores o iguales a 25000
promedios_filtrados = promedios_filas[promedios_filas['Promedio'] <= 25000]

# Gráfica 3: Promedio de lecturas por transcripto (Filtrado <= 25000)
plt.figure(figsize=(12, 6))
plt.scatter(range(len(promedios_filtrados)), promedios_filtrados['Promedio'], marker='o', color='g', alpha=1, s=20)
plt.title('Promedio de lecturas por transcripto (Filtrado <= 25000)')
plt.xlabel('Transcripto')
plt.ylabel('Promedio de lecturas')
plt.xticks(
    ticks=range(0, len(promedios_filtrados), step_promedio),
    labels=promedios_filtrados['Transcripto'][::step_promedio],
    rotation=90, fontsize=8
)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.subplots_adjust(bottom=0.25)
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
estadisticas_promedios_df = pd.DataFrame(calcular_estadisticas(promedios_filas, 'Promedio'), index=[0])
estadisticas_promedios_filtrados_df = pd.DataFrame(calcular_estadisticas(promedios_filtrados, 'Promedio'), index=[0])

# Mostrar resultados
print("\nEstadísticas Descriptivas - Suma total por paciente:")
print(estadisticas_suma_df)
print("\nEstadísticas Descriptivas - Promedio de transcriptos sin filtro (Total):")
print(estadisticas_promedios_df)
print("\nEstadísticas Descriptivas - Promedio de transcriptos con filtro (≤ 25000) (Total):")
print(estadisticas_promedios_filtrados_df)

# Histograma suma total paciente
plt.figure(figsize=(12, 6))
plt.hist(suma_totales['Suma Total'], bins=20, edgecolor='black', alpha=0.7)
plt.title('Histograma del número de lecturas por paciente')
plt.xlabel('Total de lecturas')
plt.ylabel('Frecuencia absoluta')
plt.grid(axis='y', alpha=0.75)
plt.axvline(suma_totales['Suma Total'].mean(), color='b', linestyle='dashed', linewidth=1, label='Media')
plt.legend()
plt.show()

# Histograma con transformación log2
plt.figure(figsize=(12, 6))
plt.hist(np.log2(1 + promedios_filtrados['Promedio']), bins=20, edgecolor='black', alpha=0.7)
plt.title('Histograma de promedio de lecturas por transcripto (Transformación log₂)')
plt.xlabel('log₂(1 + Promedio por fila)')
plt.ylabel('Frecuencia absoluta')
plt.grid(axis='y', alpha=0.75)
plt.axvline(np.log2(1 + promedios_filtrados['Promedio']).mean(), color='g', linestyle='dashed', linewidth=1, label='Media')
plt.legend()
plt.show()
