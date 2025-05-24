
# Cargar librerías
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(readr)
library(dplyr)
library(ggplot2)

# Define la ruta
ruta_base_cancer <- "C:/Users/Personal/Documents/Saris/Universidad/Semestre 9/TG/Carolina"
archivo <- file.path(ruta_base_cancer, "differential_expression_responden_vs_cancer.csv")

# Cargar archivo
df <- read_csv(archivo)

# Filtra genes con log2FoldChange y Geneid no NA
df <- df %>%
  filter(!is.na(log2FoldChange)) %>%
  filter(!is.na(Geneid))

# Convierte a vector nombrado
gene_vector <- df$log2FoldChange
names(gene_vector) <- df$Geneid

# Ordena el vector
gene_vector <- sort(gene_vector, decreasing = TRUE)

# Mapea SYMBOL a ENTREZID
genes_entrez <- bitr(names(gene_vector), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Elimina genes que no se mapearon
gene_vector <- gene_vector[genes_entrez$SYMBOL]
names(gene_vector) <- genes_entrez$ENTREZID

# Ejecuta GSEA
gsea_kegg <- gseKEGG(
  geneList = gene_vector,
  organism = "hsa",
  pvalueCutoff = 0.1,
  verbose = FALSE
)

# Visualiza los resultados
dotplot(gsea_kegg, showCategory = 10, title = "Análisis GSEA de genes diferencialmente expresados: respondedores a PD/PDL1 vs tejido cáncer")

# Crear el gráfico
p <- dotplot(
  gsea_kegg,
  showCategory = 10,
  title = "Análisis GSEA de genes diferencialmente expresados: respondedores a PD/PDL1 vs tejido cáncer"
)

# Guardar el gráfico como PNG
ggsave(
  filename = "gsea_responden_vs_cancer.png",
  plot = p,
  path = "C:/Users/Personal/Documents/Saris/Universidad/Semestre 9/TG/figuras",
  width = 10,
  height = 8,
  dpi = 300
)
