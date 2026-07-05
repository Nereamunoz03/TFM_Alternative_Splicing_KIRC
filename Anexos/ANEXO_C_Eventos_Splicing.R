# ==============================================================================
# Trabajo Fin de Máster
# Máster Universitario en Bioinformática - UNIR
#
# Título: Análisis del splicing alternativo asociado a la progresión tumoral
#          en carcinoma renal de células claras mediante RNA-seq
#
# Script: ANEXO_C_Eventos_Splicing.R
# Correspondiente al Anexo C del repositorio.
# Descripción: Exportación de la lista completa de eventos de splicing
#              diferencial identificados en ambos contrastes clínicos.
# ==============================================================================

# ==============================================================================
# 1. CARGA DE LIBRERÍAS
# ==============================================================================

if (!require("openxlsx", quietly = TRUE)) install.packages("openxlsx")

library(openxlsx)

# ==============================================================================
# 2. CARGA DE DATOS
# ==============================================================================

message("Cargando resultados de splicing diferencial...")

ruta_ext <- "TFM/procesados/ds_extremos.rds"
ruta_prog <- "TFM/procesados/ds_progresivo.rds"

if (file.exists(ruta_ext) && file.exists(ruta_prog)) {
  
  ds_extremos <- readRDS(ruta_ext)
  ds_progresivo <- readRDS(ruta_prog)
  
} else {
  
  stop("Error: No se encuentran los archivos 'ds_extremos.rds' y/o 'ds_progresivo.rds'. Ejecuta previamente el Script 08.")
  
}

# ==============================================================================
# 3. PREPARACIÓN DE LAS TABLAS
# ==============================================================================

ds_extremos_filtrado <- subset(
  ds_extremos,
  abs(Delta) > 0.1
)

ds_progresivo_filtrado <- subset(
  ds_progresivo,
  abs(Delta) > 0.1
)

ds_extremos_filtrado <- ds_extremos_filtrado[, c(
  "Gen",
  "ID_Evento",
  "Tipo",
  "Delta"
)]

ds_progresivo_filtrado <- ds_progresivo_filtrado[, c(
  "Gen",
  "ID_Evento",
  "Tipo",
  "Delta"
)]

colnames(ds_extremos_filtrado) <- c(
  "Gen (Symbol)",
  "ID de Evento",
  "Tipo de Evento (Anotación TCGA)",
  "Magnitud del Cambio (ΔΨ)"
)

colnames(ds_progresivo_filtrado) <- c(
  "Gen (Symbol)",
  "ID de Evento",
  "Tipo de Evento (Anotación TCGA)",
  "Magnitud del Cambio (ΔΨ)"
)

# ==============================================================================
# 4. ANOTACIÓN DE LOS TIPOS DE EVENTO
# ==============================================================================

traducir_evento <- function(x) {
  
  recode <- c(
    "ES" = "Exon Skipping (ES)",
    "AP" = "Alternate Promoter (AP)",
    "AT" = "Alternate Terminator (AT)",
    "RI" = "Retained Intron (RI)",
    "AA" = "Alternate Acceptor (AA)",
    "AD" = "Alternate Donor (AD)",
    "ME" = "Mutually Exclusive Exons (ME)"
  )
  
  recode[x]
  
}

ds_extremos_filtrado$`Tipo de Evento (Anotación TCGA)` <-
  traducir_evento(
    ds_extremos_filtrado$`Tipo de Evento (Anotación TCGA)`
  )

ds_progresivo_filtrado$`Tipo de Evento (Anotación TCGA)` <-
  traducir_evento(
    ds_progresivo_filtrado$`Tipo de Evento (Anotación TCGA)`
  )

ds_extremos_filtrado$Regulación <- ifelse(
  ds_extremos_filtrado$`Magnitud del Cambio (ΔΨ)` > 0,
  "Inclusión",
  "Exclusión"
)

ds_progresivo_filtrado$Regulación <- ifelse(
  ds_progresivo_filtrado$`Magnitud del Cambio (ΔΨ)` > 0,
  "Inclusión",
  "Exclusión"
)

# ==============================================================================
# 5. GENERACIÓN DEL LIBRO DE EXCEL
# ==============================================================================

if (!dir.exists("TFM/procesados")) {
  
  dir.create("TFM/procesados", recursive = TRUE)
  
}

wb <- createWorkbook()

addWorksheet(
  wb,
  "Stage I vs Stage IV"
)

writeData(
  wb,
  sheet = 1,
  x = ds_extremos_filtrado
)

addWorksheet(
  wb,
  "Early vs Late"
)

writeData(
  wb,
  sheet = 2,
  x = ds_progresivo_filtrado
)

# ==============================================================================
# 6. FORMATO DEL DOCUMENTO
# ==============================================================================

headerStyle <- createStyle(
  textDecoration = "bold",
  fgFill = "#EDEDED",
  halign = "center",
  valign = "center",
  border = "Bottom"
)

for (i in 1:2) {
  
  datos <- if (i == 1) ds_extremos_filtrado else ds_progresivo_filtrado
  
  addStyle(
    wb,
    sheet = i,
    style = headerStyle,
    rows = 1,
    cols = 1:ncol(datos),
    gridExpand = TRUE
  )
  
  addFilter(
    wb,
    sheet = i,
    row = 1,
    cols = 1:ncol(datos)
  )
  
  freezePane(
    wb,
    sheet = i,
    firstRow = TRUE
  )
  
  setColWidths(
    wb,
    sheet = i,
    cols = 1:ncol(datos),
    widths = "auto"
  )
  
}

# ==============================================================================
# 7. EXPORTACIÓN
# ==============================================================================

saveWorkbook(
  wb,
  "TFM/procesados/ANEXO_C_Eventos_Splicing.xlsx",
  overwrite = TRUE
)

# ==============================================================================