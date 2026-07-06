# ==============================================================================
# Trabajo Fin de Máster
# Máster Universitario en Bioinformática - UNIR
#
# Título: Análisis multiómico de expresión génica y splicing alternativo en 
#         carcinoma renal de células claras (KIRC).
#
# Script: 03_ensembl_translation_quality_control.R
# Correspondiente al apartado: 3.3.1. Traducción de identificadores de Ensembl
#                              y 3.3.2. Control de calidad.
# Descripción: Traducción de identificadores Ensembl a símbolos génicos,
#              agregación de transcritos y control de calidad de la matriz de
#              expresión génica.
# ==============================================================================

# ==============================================================================
# 1. CARGA DE LIBRERÍAS
# ==============================================================================

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

if (!require("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")
if (!require("data.table", quietly = TRUE)) install.packages("data.table")

library(data.table)
library(biomaRt)

# ==============================================================================
# 2. CARGA DE DATOS
# ==============================================================================

archivo_entrada <- "TFM/procesados/matriz_genes_final.rds"

if (file.exists(archivo_entrada)) {
  
  matriz_genes_final <- readRDS(archivo_entrada)
  setDT(matriz_genes_final)
  
  message("Matriz de expresión génica cargada correctamente.")
  
} else {
  
  stop("Error: No se encuentra 'TFM/procesados/matriz_genes_final.rds'. Ejecuta previamente el Script 01.")
  
}

# ==============================================================================
# 3. TRADUCCIÓN DE IDENTIFICADORES ENSEMBL
# ==============================================================================

message("Conectando con Ensembl para la traducción de identificadores...")

# Conexión con la base de datos Ensembl
mart <- useMart(
  "ensembl",
  dataset = "hsapiens_gene_ensembl"
)

# Traducción de identificadores Ensembl a símbolos génicos
genes_map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = matriz_genes_final$gene_id,
  mart = mart
)

setDT(genes_map)
setnames(genes_map, c("gene_id", "symbol"))

# ==============================================================================
# 4. CONTROL DE CALIDAD DE LA MATRIZ DE EXPRESIÓN
# ==============================================================================

message("Aplicando el control de calidad de la matriz de expresión...")

# Integración de los símbolos génicos
matriz_counts_preparada <- merge(
  matriz_genes_final,
  genes_map,
  by = "gene_id"
)

# Reorganización de la matriz
matriz_counts_preparada[, gene_id := NULL]
setcolorder(matriz_counts_preparada, "symbol")

# Eliminación de genes sin anotación
matriz_counts_preparada <- matriz_counts_preparada[
  symbol != "" & !is.na(symbol),
]

# Agregación de conteos correspondientes al mismo símbolo génico
matriz_counts_definitiva <- matriz_counts_preparada[
  ,
  lapply(.SD, sum),
  by = symbol
]

# Eliminación de genes sin expresión en la cohorte
genes_activos <- rowSums(
  matriz_counts_definitiva[, -1, with = FALSE],
  na.rm = TRUE
) > 0

matriz_counts_definitiva <- matriz_counts_definitiva[
  genes_activos,
]

# ==============================================================================
# 5. EXPORTACIÓN DE RESULTADOS
# ==============================================================================

if (!dir.exists("TFM/procesados")) {
  dir.create("TFM/procesados", recursive = TRUE)
}

saveRDS(
  matriz_counts_definitiva,
  "TFM/procesados/matriz_counts_definitiva.rds"
)

message("Proceso completado correctamente.")
message(paste("Genes activos:", nrow(matriz_counts_definitiva)))
message(paste("Pacientes sincronizados:", ncol(matriz_counts_definitiva) - 1))
