library(DESeq2)
library(ggplot2)

# 📂 **Definir las rutas de los archivos**
ruta_base_cancer <- "C:/Users/Personal/Documents/Saris/Universidad/Semestre 9/TG/ReadCounts/EXCEL"
archivo_cancer <- file.path(ruta_base_cancer, "fusionado_cancer.csv")

ruta_base_sanos <- "C:/Users/Personal/Documents/Saris/Universidad/Semestre 9/TG/Archivos GTEX"
archivo_sanos <- file.path(ruta_base_sanos, "fusionado_sanos.csv")

# 📌 **Leer archivos CSV**
counts_raw_cancer <- read.csv(archivo_cancer, sep = ";", check.names = FALSE)
counts_raw_sanos <- read.csv(archivo_sanos, sep = ";", check.names = FALSE)

# 🔍 **Extraer datos de expresión y nombres de genes**
counts_cancer <- counts_raw_cancer[, 3:ncol(counts_raw_cancer)]  # Desde la 3ra columna
genes_cancer <- counts_raw_cancer[, 1]  # Nombres de genes

counts_sanos <- counts_raw_sanos[, 4:ncol(counts_raw_sanos)]  # Desde la 4ta columna
genes_sanos <- counts_raw_sanos[, 3]  # Nombres de genes

# 🔍 **Convertir datos a numéricos**
counts_cancer <- as.matrix(suppressWarnings(sapply(counts_cancer, as.numeric)))
counts_sanos <- as.matrix(suppressWarnings(sapply(counts_sanos, as.numeric)))

# 🔍 **Reemplazar NA por 0**
counts_cancer[is.na(counts_cancer)] <- 0
counts_sanos[is.na(counts_sanos)] <- 0

# 🔍 **Filtrar genes con expresión total menor a 5**
genes_filtrados_cancer <- rowSums(counts_cancer) >= 5
genes_filtrados_sanos <- rowSums(counts_sanos) >= 5

counts_cancer_filtrado <- counts_cancer[genes_filtrados_cancer, , drop = FALSE]
genes_cancer_filtrado <- genes_cancer[genes_filtrados_cancer]

counts_sanos_filtrado <- counts_sanos[genes_filtrados_sanos, , drop = FALSE]
genes_sanos_filtrado <- genes_sanos[genes_filtrados_sanos]

# 🔍 **Encontrar genes comunes entre cáncer y sanos**
genes_comunes <- intersect(genes_cancer_filtrado, genes_sanos_filtrado)

# 📌 **Asegurar que ambas matrices tengan las mismas filas**
counts_cancer_final <- counts_cancer_filtrado[match(genes_comunes, genes_cancer_filtrado), , drop = FALSE]
counts_sanos_final <- counts_sanos_filtrado[match(genes_comunes, genes_sanos_filtrado), , drop = FALSE]

# 🔍 **Verificar que las filas coincidan**
stopifnot(nrow(counts_cancer_final) == nrow(counts_sanos_final))

# 📌 **Unir los datos**
counts_combined <- cbind(counts_sanos_final, counts_cancer_final)
genes_finales <- genes_comunes

# 📌 **Crear metadatos (colData)**
colData <- data.frame(
  sample = colnames(counts_combined),
  cohort = c(rep("Sanos", ncol(counts_sanos_final)), rep("Cancer", ncol(counts_cancer_final)))
)
colData$cohort <- as.factor(colData$cohort)

# 📌 **Crear objeto DESeqDataSet**
dds <- DESeqDataSetFromMatrix(countData = counts_combined, colData = colData, design = ~ cohort)

# 📌 **Ejecutar DESeq para análisis diferencial**
dds <- DESeq(dds)

# 📌 **Obtener resultados del análisis diferencial**
res <- results(dds, alpha = 0.05, lfcThreshold = 0.25)

# 📌 **Guardar resultados en CSV**
res_df <- as.data.frame(res)
res_df <- cbind(Geneid = genes_finales, res_df)  # Agregar nombres de genes

archivo_resultados <- file.path(ruta_base_cancer, "differential_expression_results.csv")
write.csv(res_df, archivo_resultados, row.names = FALSE)

print(paste("✅ Análisis de expresión diferencial completado y guardado en:", archivo_resultados))

