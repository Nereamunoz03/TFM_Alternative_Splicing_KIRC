# ==============================================================================
# Trabajo Fin de Máster
# Máster Universitario en Bioinformática - UNIR
#
# Título: Análisis del splicing alternativo asociado a la progresión tumoral
#          en carcinoma renal de células claras mediante RNA-seq
#
# Script: ANEXO_B_DEGs.R
# Correspondiente al Anexo B del repositorio.
# Descripción: Exportación de la lista completa de genes diferencialmente
#              expresados (DEGs) obtenidos en ambos contrastes clínicos.
# ==============================================================================

# ==============================================================================
# 1. CARGA DE LIBRERÍAS
# ==============================================================================

if (!require("openxlsx", quietly = TRUE)) install.packages("openxlsx")

library(openxlsx)

# ==============================================================================
# 2. CARGA DE DATOS
# ==============================================================================

message("Cargando resultados de expresión génica diferencial...")

ruta_ext <- "TFM/procesados/dge_extremos.rds"
ruta_prog <- "TFM/procesados/dge_progresivo.rds"

if (file.exists(ruta_ext) && file.exists(ruta_prog)) {
  
  dge_extremos <- readRDS(ruta_ext)
  dge_progresivo <- readRDS(ruta_prog)
  
} else {
  
  stop("Error: No se encuentran los archivos 'dge_extremos.rds' y/o 'dge_progresivo.rds'. Ejecuta previamente el Script 04.")
  
}

# ==============================================================================
# 3. PREPARACIÓN DE LAS TABLAS
# ==============================================================================

deg_extremos <- subset(
  dge_extremos,
  adj.P.Val < 0.05 & abs(logFC) > 1
)

deg_progresivo <- subset(
  dge_progresivo,
  adj.P.Val < 0.05 & abs(logFC) > 1
)

deg_extremos <- deg_extremos[, c("symbol", "logFC", "adj.P.Val")]
deg_progresivo <- deg_progresivo[, c("symbol", "logFC", "adj.P.Val")]

colnames(deg_extremos) <- c(
  "Gen (Symbol)",
  "Magnitud del Cambio (Log2FC)",
  "Significación estadística (p-adj)"
)

colnames(deg_progresivo) <- c(
  "Gen (Symbol)",
  "Magnitud del Cambio (Log2FC)",
  "Significación estadística (p-adj)"
)

deg_extremos$Regulación <- ifelse(
  deg_extremos$`Magnitud del Cambio (Log2FC)` > 0,
  "Sobreexpresado",
  "Subexpresado"
)

deg_progresivo$Regulación <- ifelse(
  deg_progresivo$`Magnitud del Cambio (Log2FC)` > 0,
  "Sobreexpresado",
  "Subexpresado"
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
  "Stage I vs Stage IV"
)

writeData(
  wb,
  sheet = 1,
  x = deg_extremos
)

addWorksheet(
  wb,
  "Early vs Late"
)

writeData(
  wb,
  sheet = 2,
  x = deg_progresivo
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

for (i in 1:2) {
  
  datos <- if (i == 1) deg_extremos else deg_progresivo
  
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
# 6. EXPORTACIÓN
# ==============================================================================

saveWorkbook(
  wb,
  "TFM/procesados/ANEXO_B_DEGs.xlsx",
  overwrite = TRUE
)

# ==============================================================================