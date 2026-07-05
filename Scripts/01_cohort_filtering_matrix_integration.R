# ==============================================================================
# Trabajo Fin de Máster
# Máster Universitario en Bioinformática - UNIR
#
# Título: Análisis del splicing alternativo asociado a la progresión tumoral
#          en carcinoma renal de células claras mediante RNA-seq
#
# Script: 01_cohort_filtering_matrix_integration.R
# Correspondiente al apartado: 3.2. Filtrado de la cohorte clínica y unificación
#                              de matrices.
# Descripción: Filtrado de la cohorte, procesamiento de las matrices de expresión
#              génica y splicing alternativo, e integración de los datos clínicos.
# ==============================================================================

# ==============================================================================
# 1. CONFIGURACIÓN DE RUTAS
# ==============================================================================

# Configuración de las rutas de trabajo
ruta <- "C:/Users/PC/Downloads/TFM/DATOS"
ruta_metadatos <- paste0(ruta, "/TCGA_KIRC")
ruta_gene_counts <- paste0(ruta, "/gdc")
ruta_psi <- "C:/Users/PC/Downloads/TFM/DATOS/PSI_download_KIRC/PSI_download_KIRC.txt"

# ==============================================================================
# 2. CARGA DE LIBRERÍAS
# ==============================================================================

if (!require("data.table")) install.packages("data.table")
library(data.table)

# ==============================================================================
# 3. CONSTRUCCIÓN DE LA MATRIZ DE EXPRESIÓN GÉNICA
# ==============================================================================

message("Construcción de la matriz de expresión génica...")

# Lectura de los metadatos de las muestras
metadata <- read.table(
  file.path(ruta_metadatos, "gdc_sample_sheet.2026-04-07.tsv"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

# Localización de los archivos de conteos génicos
ficheros_en_disco <- list.files(
  path = ruta_gene_counts,
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

# Asociación entre los archivos descargados y los metadatos
tabla_rutas <- data.frame(
  ruta_completa = ficheros_en_disco,
  file_name_disco = basename(ficheros_en_disco)
)

metadata_final <- merge(
  metadata,
  tabla_rutas,
  by.x = "File Name",
  by.y = "file_name_disco"
)

message(paste("Muestras identificadas:", nrow(metadata_final)))

# Función para la lectura de los archivos de conteos génicos
leer_y_limpiar <- function(i) {
  
  ruta <- metadata_final$ruta_completa[i]
  nombre <- metadata_final$`Sample ID`[i]
  
  datos <- fread(
    ruta,
    skip = 1,
    select = c(1, 4)
  )
  
  setnames(datos, c("gene_id", nombre))
  
  # Conservación exclusiva de identificadores Ensembl
  datos <- datos[grep("ENSG", gene_id)]
  
  return(datos)
}

message("Leyendo archivos de conteos génicos...")

lista_de_tablas <- lapply(
  1:nrow(metadata_final),
  leer_y_limpiar
)

message("Integrando las muestras en una única matriz...")

matriz_genes <- Reduce(
  function(x, y) merge(x, y, by = "gene_id", all = TRUE),
  lista_de_tablas
)

# Eliminación de la versión de los identificadores Ensembl
matriz_genes$gene_id <- gsub("\\..*", "", matriz_genes$gene_id)

# ==============================================================================
# 4. PROCESAMIENTO DE LA MATRIZ DE SPLICING ALTERNATIVO
# ==============================================================================

message("Procesamiento de la matriz de splicing alternativo...")

# Lectura de la matriz de valores PSI
datos_psi_raw <- fread(
  ruta_psi,
  header = TRUE,
  quote = "",
  fill = TRUE,
  check.names = FALSE
)

# Extracción de la información clínica
clinica_psi <- datos_psi_raw[1:13, ]

# Extracción de la matriz cuantitativa de valores PSI
matriz_psi <- datos_psi_raw[14:nrow(datos_psi_raw), ]

# Homogeneización de los nombres de las muestras
colnames(matriz_psi) <- gsub("_", "-", colnames(matriz_psi))

identificacion <- matriz_psi[, 1:10]
valores_numericos <- matriz_psi[, 11:ncol(matriz_psi)]

# Conversión de los valores PSI a formato numérico
valores_numericos <- as.data.frame(
  lapply(valores_numericos, as.numeric)
)

matriz_psi_final <- cbind(
  identificacion,
  valores_numericos
)

# ==============================================================================
# 5. INTEGRACIÓN DE MATRICES Y FILTRADO DE LA COHORTE
# ==============================================================================

message("Integración de las matrices y selección de pacientes comunes...")

# Homogeneización de los identificadores de las muestras
colnames(matriz_psi_final) <- gsub("\\.", "-", colnames(matriz_psi_final))
colnames(clinica_psi) <- gsub("_", "-", colnames(clinica_psi))

nombres_genes <- colnames(matriz_genes)
nombres_genes_cortos <- gsub("-(01|11)[A-Z]$", "", nombres_genes)

colnames(matriz_genes) <- nombres_genes_cortos
colnames(matriz_genes)[1] <- "gene_id"

# Identificación de los pacientes comunes
pacientes_genes <- colnames(matriz_genes)[-1]
pacientes_psi <- colnames(matriz_psi_final)[11:ncol(matriz_psi_final)]

pacientes_comunes <- intersect(
  pacientes_genes,
  pacientes_psi
)

message(paste("Pacientes comunes identificados:", length(pacientes_comunes)))

# Filtrado de las matrices
matriz_genes_final <- matriz_genes[
  ,
  c("gene_id", pacientes_comunes),
  with = FALSE
]

matriz_psi_final <- matriz_psi_final[
  ,
  c(colnames(matriz_psi_final)[1:10], pacientes_comunes),
  with = FALSE
]

# Construcción de la matriz clínica
clinica_final <- as.data.frame(
  t(clinica_psi[, ..pacientes_comunes])
)

# Asignación de las variables clínicas
colnames(clinica_final) <- clinica_psi$symbol
rownames(clinica_final) <- pacientes_comunes

clinica_final$Age <- as.numeric(
  as.character(clinica_final$Age)
)

# ==============================================================================
# 6. EXPORTACIÓN DE RESULTADOS
# ==============================================================================

# Creación del directorio de salida si no existe
if (!dir.exists("TFM/procesados")) {
  dir.create("TFM/procesados", recursive = TRUE)
}

saveRDS(
  matriz_genes_final,
  "TFM/procesados/matriz_genes_final.rds"
)

saveRDS(
  matriz_psi_final,
  "TFM/procesados/matriz_psi_final.rds"
)

saveRDS(
  clinica_final,
  "TFM/procesados/clinica_final.rds"
)

message("Proceso completado correctamente.")

message(
  paste(
    "Matriz de expresión génica:",
    ncol(matriz_genes_final) - 1,
    "muestras"
  )
)

message(
  paste(
    "Matriz de splicing alternativo:",
    ncol(matriz_psi_final) - 10,
    "muestras"
  )
)

message(
  paste(
    "Matriz clínica:",
    nrow(clinica_final),
    "pacientes"
  )
)