library(DESeq2)

# 📂 **Definir rutas de los archivos**
ruta_base_cancer <- "C:/Users/Personal/Documents/Saris/Universidad/Semestre 9/TG/ReadCounts/EXCEL"
archivo_cancer <- file.path(ruta_base_cancer, "fusionado_cancer.csv")

ruta_base_responden <- "C:/Users/Personal/Documents/Saris/Universidad/Semestre 9/TG/Carolina"
archivo_responden <- file.path(ruta_base_responden, "Matriz_Responden.csv")

# 📌 **Leer archivos CSV**
counts_raw_cancer <- read.csv(archivo_cancer, sep = ";", check.names = FALSE)
counts_raw_responden <- read.csv(archivo_responden, sep = ";", check.names = FALSE)

# 🔍 **Extraer datos de expresión y nombres de genes**
counts_cancer <- counts_raw_cancer[, 3:ncol(counts_raw_cancer)] 
genes_cancer <- counts_raw_cancer[, 1]  # Nombres de genes

counts_responden <- counts_raw_responden[, 3:ncol(counts_raw_responden)]  # Desde la 3ra columna
genes_responden <- counts_raw_responden[, 1]  # Nombres de genes

# 🔍 **Convertir datos a numéricos**
counts_cancer <- as.matrix(sapply(counts_cancer, as.numeric))
counts_responden <- as.matrix(sapply(counts_responden, as.numeric))

# 🔍 **Reemplazar NA por 0**
counts_cancer[is.na(counts_cancer)] <- 0
counts_responden[is.na(counts_responden)] <- 0

# 🔍 **Filtrar genes con expresión baja**
genes_filtrados_cancer <- rowSums(counts_cancer) >= 10
genes_filtrados_responden <- rowSums(counts_responden) >= 10

counts_cancer_filtrado <- counts_cancer[genes_filtrados_cancer, , drop = FALSE]
genes_cancer_filtrado <- genes_cancer[genes_filtrados_cancer]

counts_responden_filtrado <- counts_responden[genes_filtrados_responden, , drop = FALSE]
genes_responden_filtrado <- genes_responden[genes_filtrados_responden]

# 🔍 **Encontrar genes comunes**
genes_comunes <- intersect(genes_cancer_filtrado, genes_responden_filtrado)

# 📌 **Asegurar que ambas matrices tengan las mismas filas**
counts_cancer_final <- counts_cancer_filtrado[match(genes_comunes, genes_cancer_filtrado), , drop = FALSE]
counts_responden_final <- counts_responden_filtrado[match(genes_comunes, genes_responden_filtrado), , drop = FALSE]

# 🔍 **Verificar que las filas coincidan**
stopifnot(nrow(counts_cancer_final) == nrow(counts_responden_final))

# 📌 **Unir los datos**
counts_combined <- cbind(counts_cancer_final, counts_responden_final)
genes_finales <- genes_comunes

# 📌 **Crear metadatos (colData)**
colData <- data.frame(
  sample = colnames(counts_combined),
  cohort = c(rep("cancer", ncol(counts_cancer_final)), rep("responden", ncol(counts_responden_final)))
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

archivo_resultados <- file.path(ruta_base_responden, "differential_expression_responden_vs_cancer.csv")
write.csv(res_df, archivo_resultados, row.names = FALSE)

print(paste("✅ Análisis de expresión diferencial completado y guardado en:", archivo_resultados))
