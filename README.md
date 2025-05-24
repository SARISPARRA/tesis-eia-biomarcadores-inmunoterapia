# 🧪 Análisis de Biomarcadores en Inmunoterapia PD-1/PD-L1

Este repositorio contiene todo el flujo de trabajo bioinformático desarrollado como parte del trabajo de grado **"Análisis comparativo de la expresión de biomarcadores asociados con la respuesta a la inmunoterapia PD-1/PD-L1 usando herramientas bioinformáticas"**.

El objetivo principal es identificar genes diferencialmente expresados que puedan funcionar como **biomarcadores predictivos** de respuesta a inmunoterapia en pacientes con cáncer tratados con inhibidores de PD-1 y PD-L1, usando datos públicos de RNA-seq y herramientas como DESeq2, GSEA, Kruskal-Wallis y Mann-Whitney U.

---

## 🗂 Estructura del Repositorio

| Carpeta | Descripción |
|--------|-------------|
| `01. Metadata` | Archivos que documentan la organización, clasificación y origen de las muestras (GTEx, GEO, ENA, etc.). |
| `02. Cuantificacion_expresion_rnaseq` | Scripts y resultados del procesamiento de secuencias RNA-seq: limpieza con Cutadapt, alineación con STAR y cuantificación con FeatureCounts. |
| `03. Estadistica descriptiva` | Análisis descriptivo de la expresión génica: histogramas, transformaciones log₂(1+x), filtros de lecturas, y visualizaciones. |
| `04. Análisis de expresión diferencial` | Comparaciones entre cáncer, sanos, respondedores y no respondedores. Incluye resultados de DESeq2 y tablas de genes significativos. |
| `05. GSEA` | Análisis de enriquecimiento funcional (GSEA) con base de datos KEGG, diferenciando entre perfiles tumorales y respuesta inmunológica. |
| `06. Kruskal Wallis` | Evaluación de expresión génica entre los cuatro grupos (sano, cáncer, responde, no responde) con una prueba no paramétrica robusta. |
| `07. Mann Whitney U` | Comparaciones entre pares de grupos usando la prueba U de Mann-Whitney, centrado en los 20 genes de interés identificados como biomarcadores. |

---

## 📌 Metodología General

1. **Obtención de datos:** RNA-seq de 3,369 muestras sanas (GTEx) y 231 tumorales (GEO/ENA), clasificadas según respuesta a inmunoterapia.
2. **Procesamiento:** Cutadapt + FastQC/MultiQC → STAR → FeatureCounts.
3. **Análisis:**
   - Estadística descriptiva.
   - Análisis de expresión diferencial (DESeq2).
   - GSEA (KEGG).
   - Pruebas no paramétricas: Kruskal-Wallis y Mann-Whitney U.
4. **Validación de genes candidatos:** 20 genes seleccionados por Carolina Castaño fueron evaluados como posibles **biomarcadores predictivos**.

---

## 🔬 Resultados Clave

- Genes como `DPYSL5`, `BIRC2`, `SFTPC`, `PAPLN`, y `CXCL9` mostraron diferencias altamente significativas en expresión entre grupos.
- Estos genes podrían ayudar a predecir la **respuesta positiva a inmunoterapia**, diferenciando pacientes que sí responden.
- Los análisis integrados (AED, Kruskal-Wallis, Mann-Whitney) validan su utilidad potencial como **biomarcadores clínicos**.

---

## 🎯 Objetivo Final

Este repositorio busca apoyar la validación del modelo de biomarcadores propuesto en la **tesis doctoral de Carolina Castaño**, demostrando cómo el análisis computacional puede facilitar la identificación de perfiles genómicos con relevancia clínica en inmunoterapia contra el cáncer.

---

## 📎 Requisitos

- Python 3.10 / R ≥ 4.2
- Cutadapt ≥ 4.9
- STAR ≥ 2.7.11b
- FastQC, MultiQC, DESeq2
- pandas, matplotlib, seaborn, scipy


