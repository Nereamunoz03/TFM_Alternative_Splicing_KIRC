# ==============================================================================
# Trabajo Fin de Máster
# Máster Universitario en Bioinformática - UNIR
#
# Título: Análisis del splicing alternativo asociado a la progresión tumoral
#          en carcinoma renal de células claras mediante RNA-seq
#
# Script: 11_tumor_status_association.R
# Correspondiente al apartado: 3.6. Asociación de los genes candidatos con el
#                              estado tumoral.
# Descripción: Comparación de los niveles de expresión de los genes candidatos
#              multiómicos entre pacientes con y sin evidencia de tumor mediante
#              análisis no paramétrico y generación de una tabla resumen.
# ==============================================================================

# ==============================================================================
# 1. CARGA DE LIBRERÍAS
# ==============================================================================

if (!require("gtsummary", quietly = TRUE)) install.packages("gtsummary")
if (!require("gt", quietly = TRUE)) install.packages("gt")
if (!require("data.table", quietly = TRUE)) install.packages("data.table")
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!require("reshape2", quietly = TRUE)) install.packages("reshape2")

library(data.table)
library(ggplot2)
library(reshape2)
library(gtsummary)
library(gt)

# ==============================================================================
# 2. CARGA DE DATOS
# ==============================================================================

message("Cargando datos clínicos y de expresión génica...")

ruta_clinica <- "TFM/procesados/clinica_filtrada_estadios.rds"
ruta_counts  <- "TFM/procesados/matriz_counts_definitiva.rds"

if (file.exists(ruta_clinica) & file.exists(ruta_counts)) {
  
  clinica_filtrada_estadios <- readRDS(ruta_clinica)
  matriz_counts_definitiva  <- readRDS(ruta_counts)
  
  setDT(matriz_counts_definitiva)
  
  clinica_dt <- as.data.table(
    clinica_filtrada_estadios,
    keep.rownames = "Sample"
  )
  
  message("Datos cargados correctamente.")
  
} else {
  
  stop("Error: No se encuentran los archivos necesarios en 'TFM/procesados/'. Ejecuta previamente los Scripts 02 y 03.")
  
}

genes_candidatos <- c(
  "SAA2",
  "APOLD1",
  "UGT1A10",
  "MAPT",
  "DNASE1L3",
  "SCGB3A2"
)

# ==============================================================================
# 3. PREPARACIÓN DE LOS DATOS
# ==============================================================================

message("Preparando la cohorte para el análisis...")

clinica_tumor <- clinica_dt[
  Tumor_Status %in% c("TUMOR FREE", "WITH TUMOR"),
]

counts_candidatos <- matriz_counts_definitiva[
  symbol %in% genes_candidatos
]

df_genes <- as.data.frame(
  t(counts_candidatos[, -1, with = FALSE])
)

colnames(df_genes) <- counts_candidatos$symbol
df_genes$Sample <- rownames(df_genes)

setDT(df_genes)

df_analisis <- merge(
  clinica_tumor[, .(Sample, Tumor_Status)],
  df_genes,
  by = "Sample"
)

# ==============================================================================
# 4. TRANSFORMACIÓN LOGARÍTMICA
# ==============================================================================

message("Aplicando transformación logarítmica...")

for (gen in genes_candidatos) {
  
  df_analisis[[gen]] <- log2(
    df_analisis[[gen]] + 1
  )
  
}

# ==============================================================================
# 5. ANÁLISIS ESTADÍSTICO
# ==============================================================================

message("Generando la tabla resumen...")

df_tabla_clinica <- df_analisis[
  ,
  .(
    Tumor_Status,
    SAA2,
    APOLD1,
    UGT1A10,
    MAPT,
    DNASE1L3,
    SCGB3A2
  )
]

tabla_tumor_status <- df_tabla_clinica %>%
  tbl_summary(
    by = Tumor_Status,
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})"
    ),
    digits = list(
      all_continuous() ~ 2
    ),
    label = list(
      SAA2     ~ "SAA2 (Log2 Counts)",
      APOLD1   ~ "APOLD1 (Log2 Counts)",
      UGT1A10  ~ "UGT1A10 (Log2 Counts)",
      MAPT     ~ "MAPT (Log2 Counts)",
      DNASE1L3 ~ "DNASE1L3 (Log2 Counts)",
      SCGB3A2  ~ "SCGB3A2 (Log2 Counts)"
    )
  ) %>%
  add_p(
    test = all_continuous() ~ "wilcox.test",
    pvalue_fun = function(x) style_pvalue(x, digits = 3)
  ) %>%
  bold_labels() %>%
  bold_p(t = 0.05) %>%
  modify_spanning_header(
    all_stat_cols() ~ "**Estado clínico del paciente (TCGA-KIRC)**"
  ) %>%
  modify_header(
    label = "**Gen candidato**"
  )

# ==============================================================================
# 6. EXPORTACIÓN DE RESULTADOS
# ==============================================================================

if (!dir.exists("TFM/procesados")) {
  dir.create("TFM/procesados", recursive = TRUE)
}

gtsave(
  as_gt(tabla_tumor_status),
  "TFM/procesados/tabla_tumor_status.html"
)

message("Proceso completado correctamente.")