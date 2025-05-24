# Instalar paquetes si no están disponibles
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")

# Cargar librerías
library(ggplot2)
library(ggrepel)

# Definir umbrales de significancia
pvalue_cutoff <- 0.05  
lfc_cutoff <- 0.25    

# Cargar los resultados del análisis de expresión diferencial
ruta_base <- "C:/Users/Personal/Documents/Saris/Universidad/Semestre 9/TG/Carolina"
archivo_resultados <- file.path(ruta_base, "differential_expression_noresponden_vs_cancer.csv")
results_df <- read.csv(archivo_resultados)

# Reemplazar valores de padj = 0
results_df$padj[results_df$padj == 0] <- 1e-10  

# Filtrar valores NA en padj
results_df <- na.omit(results_df)

# Definir límite superior eje Y
y_lim <- min(max(-log10(results_df$padj), na.rm = TRUE) * 1.2, 100)

# Clasificar genes
results_df$Estado <- "No significativo"
results_df$Estado[results_df$padj < pvalue_cutoff & results_df$log2FoldChange > lfc_cutoff] <- "Sobreexpresados"
results_df$Estado[results_df$padj < pvalue_cutoff & results_df$log2FoldChange < -lfc_cutoff] <- "Subexpresados"

# Convertir Estado en factor
results_df$Estado <- factor(results_df$Estado, levels = c("Subexpresados", "No significativo", "Sobreexpresados"))

# 🔍 Seleccionar los 15 genes más sobreexpresados y subexpresados
top_sobreexpresados <- results_df[results_df$Estado == "Sobreexpresados", ]
top_sobreexpresados <- top_sobreexpresados[order(-top_sobreexpresados$log2FoldChange), ][1:15, ]

top_subexpresados <- results_df[results_df$Estado == "Subexpresados", ]
top_subexpresados <- top_subexpresados[order(top_subexpresados$log2FoldChange), ][1:15, ]

# Combinar los genes seleccionados
top_genes <- rbind(top_sobreexpresados, top_subexpresados)


# 🎨 Generar Volcano Plot
volcan_plot <- ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj), color = Estado)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_text_repel(data = top_genes, aes(label = Geneid), 
                  size = 3, 
                  box.padding = 0.5, 
                  segment.color = "gray60",
                  max.overlaps = 50) +
  scale_color_manual(values = c("red", "black", "blue")) +
  labs(title = "Expresión Diferencial: noresponden vs cancer",
       x = "Cambio en Log2 de Expresión (Log2 Fold Change)",
       y = "-Log10 del p-valor ajustado (-Log10 padj)") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "gray") +
  ylim(0, y_lim)

# 📈 Calcular estadísticas adicionales

# Total de genes analizados
total_genes <- nrow(results_df)

# Genes upregulated
genes_up <- sum(results_df$log2FoldChange > lfc_cutoff & results_df$padj < pvalue_cutoff)
porcentaje_up <- round((genes_up / total_genes) * 100, 2)

# Genes downregulated
genes_down <- sum(results_df$log2FoldChange < -lfc_cutoff & results_df$padj < pvalue_cutoff)
porcentaje_down <- round((genes_down / total_genes) * 100, 2)

# Genes con baja expresión (<5 reads) — suponiendo que tienes una columna llamada 'baseMean' con el promedio de lecturas
if ("baseMean" %in% colnames(results_df)) {
  genes_baja_expresion <- sum(results_df$baseMean < 5)
  porcentaje_baja_expresion <- round((genes_baja_expresion / total_genes) * 100, 2)
} else {
  genes_baja_expresion <- NA
  porcentaje_baja_expresion <- NA
}

# 📝 Imprimir resumen
cat("\nResumen de análisis:\n")
cat("Total de genes analizados: ", total_genes, "\n")
cat("Genes upregulated (LFC > ", lfc_cutoff, "): ", genes_up, " (", porcentaje_up, "%)\n", sep = "")
cat("Genes downregulated (LFC < -", lfc_cutoff, "): ", genes_down, " (", porcentaje_down, "%)\n", sep = "")

if (!is.na(genes_baja_expresion)) {
  cat("Genes con baja expresión (<5 reads): ", genes_baja_expresion, " (", porcentaje_baja_expresion, "%)\n", sep = "")
} else {
  cat("⚠️ Advertencia: No se encontró la columna 'baseMean', no se pudo calcular baja expresión.\n")
}


# 🔹 Guardar gráfico
ruta_guardado <- "C:/Users/Personal/Documents/Saris/Universidad/Semestre 9/TG/figuras/Top30_noresponden_vs_cancers_.png"
ggsave(filename = ruta_guardado, plot = volcan_plot, width = 8, height = 6, dpi = 300)

cat("\n✅ Gráfico guardado en:", ruta_guardado, "\n")

