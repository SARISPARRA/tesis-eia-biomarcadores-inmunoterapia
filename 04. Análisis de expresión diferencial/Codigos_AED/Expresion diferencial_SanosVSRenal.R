library(DESeq2)

# 📂 **Definir rutas de los archivos**
ruta_base_sanos <- "C:/Users/Personal/Documents/Saris/Universidad/Semestre 9/TG/Archivos GTEX"
archivo_sanos <- file.path(ruta_base_sanos, "fusionado_sanos.csv")

ruta_base_Utero_Ovarios <- "C:/Users/Personal/Documents/Saris/Universidad/Semestre 9/TG/ReadCounts/CANCER"
archivo_Utero_Ovarios <- file.path(ruta_base_Utero_Ovarios, "Utero_Ovarios.csv")

# 📌 **Leer archivos CSV**
counts_raw_sanos <- read.csv(archivo_sanos, sep = ";", check.names = FALSE)
counts_raw_Utero_Ovarios <- read.csv(archivo_Utero_Ovarios, sep = ";", check.names = FALSE)

# 🔍 **Extraer datos de expresión y nombres de genes**
counts_sanos <- counts_raw_sanos[, 4:ncol(counts_raw_sanos)]  # Desde la 4ta columna
genes_sanos <- counts_raw_sanos[, 3]  # Nombres de genes

counts_Utero_Ovarios <- counts_raw_Utero_Ovarios[, 3:ncol(counts_raw_Utero_Ovarios)]  # Desde la 3ra columna
genes_Utero_Ovarios <- counts_raw_Utero_Ovarios[, 1]  # Nombres de genes

# 🔍 **Convertir datos a numéricos**
counts_sanos <- as.matrix(sapply(counts_sanos, as.numeric))
counts_Utero_Ovarios <- as.matrix(sapply(counts_Utero_Ovarios, as.numeric))

# 🔍 **Reemplazar NA por 0**
counts_sanos[is.na(counts_sanos)] <- 0
counts_Utero_Ovarios[is.na(counts_Utero_Ovarios)] <- 0

# 🔍 **Filtrar genes con expresión baja**
genes_filtrados_sanos <- rowSums(counts_sanos) >= 10
genes_filtrados_Utero_Ovarios <- rowSums(counts_Utero_Ovarios) >= 10

counts_sanos_filtrado <- counts_sanos[genes_filtrados_sanos, , drop = FALSE]
genes_sanos_filtrado <- genes_sanos[genes_filtrados_sanos]

counts_Utero_Ovarios_filtrado <- counts_Utero_Ovarios[genes_filtrados_Utero_Ovarios, , drop = FALSE]
genes_Utero_Ovarios_filtrado <- genes_Utero_Ovarios[genes_filtrados_Utero_Ovarios]

# 🔍 **Encontrar genes comunes**
genes_comunes <- intersect(genes_sanos_filtrado, genes_Utero_Ovarios_filtrado)

# 📌 **Asegurar que ambas matrices tengan las mismas filas**
counts_sanos_final <- counts_sanos_filtrado[match(genes_comunes, genes_sanos_filtrado), , drop = FALSE]
counts_Utero_Ovarios_final <- counts_Utero_Ovarios_filtrado[match(genes_comunes, genes_Utero_Ovarios_filtrado), , drop = FALSE]

# 🔍 **Verificar que las filas coincidan**
stopifnot(nrow(counts_sanos_final) == nrow(counts_Utero_Ovarios_final))

# 📌 **Unir los datos**
counts_combined <- cbind(counts_sanos_final, counts_Utero_Ovarios_final)
genes_finales <- genes_comunes

# 📌 **Crear metadatos (colData)**
colData <- data.frame(
  sample = colnames(counts_combined),
  cohort = c(rep("Sanos", ncol(counts_sanos_final)), rep("Utero_Ovarios", ncol(counts_Utero_Ovarios_final)))
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

archivo_resultados <- file.path(ruta_base_Utero_Ovarios, "differential_expression_Utero_Ovarios_vs_sanos.csv")
write.csv(res_df, archivo_resultados, row.names = FALSE)

print(paste("✅ Análisis de expresión diferencial completado y guardado en:", archivo_resultados))
