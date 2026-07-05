# ==============================================================================
# Trabajo Fin de Máster
# Máster Universitario en Bioinformática - UNIR
#
# Título: Análisis del splicing alternativo asociado a la progresión tumoral
#          en carcinoma renal de células claras mediante RNA-seq
#
# Script: ANEXO_A_Muestras_TCGA.R
# Correspondiente al Anexo A del repositorio.
# Descripción: Generación del listado completo de muestras TCGA-KIRC incluidas
#              en la cohorte final del estudio.
# ==============================================================================

# ==============================================================================
# 1. CARGA DE LIBRERÍAS
# ==============================================================================

if (!require("openxlsx", quietly = TRUE)) install.packages("openxlsx")

library(openxlsx)

# ==============================================================================
# 2. CARGA DE DATOS
# ==============================================================================

message("Cargando cohorte clínica...")

ruta_clinica <- "TFM/procesados/clinica_final.rds"

if (file.exists(ruta_clinica)) {
  
  clinica_final <- readRDS(ruta_clinica)
  
} else {
  
  stop("Error: No se encuentra 'clinica_final.rds'. Ejecuta previamente el Script 01.")
  
}

# ==============================================================================
# 3. CONSTRUCCIÓN DE LA TABLA DEL ANEXO
# ==============================================================================

codigos_tcga <- rownames(clinica_final)

anexo_A <- data.frame(
  Nº = seq_along(codigos_tcga),
  `Código TCGA` = codigos_tcga,
  Estadio = clinica_final$Pathologic_Tumor_Stage
)

# ==============================================================================
# 4. GENERACIÓN DEL LIBRO DE EXCEL
# ==============================================================================

if (!dir.exists("TFM/procesados")) {
  
  dir.create("TFM/procesados", recursive = TRUE)
  
}

wb <- createWorkbook()

addWorksheet(
  wb,
  "Muestras TCGA"
)

writeData(
  wb,
  sheet = 1,
  x = anexo_A
)

# ==============================================================================
# 5. FORMATO DEL DOCUMENTO
# ==============================================================================

headerStyle <- createStyle(
  textDecoration = "bold",
  fgFill = "#EDEDED",
  halign = "center",
  valign = "center",
  border = "Bottom"
)

addStyle(
  wb,
  sheet = 1,
  style = headerStyle,
  rows = 1,
  cols = 1:ncol(anexo_A),
  gridExpand = TRUE
)

addFilter(
  wb,
  sheet = 1,
  row = 1,
  cols = 1:ncol(anexo_A)
)

freezePane(
  wb,
  sheet = 1,
  firstRow = TRUE
)

setColWidths(
  wb,
  sheet = 1,
  cols = 1:ncol(anexo_A),
  widths = "auto"
)

# ==============================================================================
# 6. EXPORTACIÓN
# ==============================================================================

saveWorkbook(
  wb,
  "TFM/procesados/ANEXO_A_Muestras_TCGA.xlsx",
  overwrite = TRUE
)

# ==============================================================================