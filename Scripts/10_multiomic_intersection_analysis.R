# ==============================================================================
# Trabajo Fin de Máster
# Máster Universitario en Bioinformática - UNIR
#
# Título: Análisis del splicing alternativo asociado a la progresión tumoral
#          en carcinoma renal de células claras mediante RNA-seq
#
# Script: 10_multiomic_intersection_analysis.R
# Correspondiente al apartado: 3.5. Integración multiómica.
# Descripción: Integración de los resultados de expresión génica diferencial
#              y splicing diferencial para identificar biomarcadores
#              candidatos compartidos entre ambos niveles moleculares.
# ==============================================================================

# ==============================================================================
# 1. CARGA DE LIBRERÍAS
# ==============================================================================

if (!require("VennDiagram", quietly = TRUE)) install.packages("VennDiagram")
if (!require("data.table", quietly = TRUE)) install.packages("data.table")
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(data.table)
library(VennDiagram)
library(ggplot2)
library(grid)

options(venn.print.log = FALSE)

# ==============================================================================
# 2. CARGA DE DATOS
# ==============================================================================

message("Cargando resultados de expresión génica y splicing diferencial...")

ruta_dge_ext  <- "TFM/procesados/dge_extremos.rds"
ruta_dge_prog <- "TFM/procesados/dge_progresivo.rds"
ruta_ds_ext   <- "TFM/procesados/ds_extremos.rds"
ruta_ds_prog  <- "TFM/procesados/ds_progresivo.rds"

if (all(file.exists(c(
  ruta_dge_ext,
  ruta_dge_prog,
  ruta_ds_ext,
  ruta_ds_prog
)))) {
  
  dge_extremos   <- readRDS(ruta_dge_ext)
  dge_progresivo <- readRDS(ruta_dge_prog)
  ds_extremos    <- readRDS(ruta_ds_ext)
  ds_progresivo  <- readRDS(ruta_ds_prog)
  
  setDT(dge_extremos)
  setDT(dge_progresivo)
  setDT(ds_extremos)
  setDT(ds_progresivo)
  
  message("Datos cargados correctamente.")
  
} else {
  
  stop("Error: No se encuentran todos los archivos necesarios en 'TFM/procesados/'. Ejecuta previamente los Scripts 04 y 08.")
  
}

if (!dir.exists("TFM/procesados")) {
  dir.create("TFM/procesados", recursive = TRUE)
}

if (!dir.exists("TFM/graficos")) {
  dir.create("TFM/graficos", recursive = TRUE)
}

# ==============================================================================
# 3. INTEGRACIÓN MULTIÓMICA: CONTRASTE EXTREMO
# ==============================================================================

message("Procesando el contraste Stage I vs Stage IV...")

# Identificación de genes candidatos en cada nivel molecular
genes_ds_ext  <- unique(ds_extremos$Gen)
genes_dge_ext <- unique(
  dge_extremos[
    adj.P.Val < 0.05 & abs(logFC) > 1,
    symbol
  ]
)

# Intersección entre expresión diferencial y splicing diferencial
master_extremos <- intersect(
  genes_ds_ext,
  genes_dge_ext
)

# Diagrama de Venn
venn_ext <- venn.diagram(
  x = list(
    Splicing = genes_ds_ext,
    Expresion = genes_dge_ext
  ),
  category.names = c(
    "Splicing Alternativo (DS)",
    "Expresión Diferencial (DGE)"
  ),
  filename = NULL,
  fill = c("#4f8c9d", "#d95f02"),
  alpha = c(0.6, 0.6),
  cat.col = c("#4f8c9d", "#d95f02"),
  cat.fontface = "bold",
  cat.dist = c(0.08, 0.08),
  margin = 0.1,
  fontface = "bold"
)

png(
  "TFM/graficos/Venn_Integracion_Extremos.png",
  width = 1800,
  height = 1800,
  res = 300
)

grid.draw(venn_ext)

dev.off()

# Integración de la información molecular
tabla_dge_ext <- dge_extremos[
  symbol %in% master_extremos,
  .(symbol, logFC, adj.P.Val)
]

info_ds_ext <- ds_extremos[
  ,
  .(Gen, ID_Evento, Nombre_Biologico, Delta)
]

setnames(
  info_ds_ext,
  "Gen",
  "symbol"
)

master_extremos_tabla <- merge(
  tabla_dge_ext,
  info_ds_ext,
  by = "symbol",
  allow.cartesian = TRUE
)

master_extremos_tabla <- master_extremos_tabla[
  order(adj.P.Val)
]

saveRDS(
  master_extremos,
  "TFM/procesados/master_extremos.rds"
)

write.csv(
  master_extremos_tabla,
  "TFM/procesados/master_extremos_tabla.csv",
  row.names = FALSE
)

# ==============================================================================
# 4. INTEGRACIÓN MULTIÓMICA: CONTRASTE PROGRESIVO
# ==============================================================================

message("Procesando el contraste Early vs Late...")

# Identificación de genes candidatos en cada nivel molecular
genes_ds_prog <- unique(ds_progresivo$Gen)

genes_dge_prog <- unique(
  dge_progresivo[
    adj.P.Val < 0.05 & abs(logFC) > 1,
    symbol
  ]
)

# Intersección entre expresión diferencial y splicing diferencial
master_progresivo <- intersect(
  genes_ds_prog,
  genes_dge_prog
)

# Diagrama de Venn
venn_prog <- venn.diagram(
  x = list(
    Splicing = genes_ds_prog,
    Expresion = genes_dge_prog
  ),
  category.names = c(
    "Splicing Alternativo (DS)",
    "Expresión Diferencial (DGE)"
  ),
  filename = NULL,
  fill = c("#4daf4a", "#984ea3"),
  alpha = c(0.6, 0.6),
  cat.col = c("#4daf4a", "#984ea3"),
  cat.fontface = "bold",
  cat.dist = c(0.08, 0.08),
  margin = 0.1,
  fontface = "bold"
)

png(
  "TFM/graficos/Venn_Integracion_Progresivo.png",
  width = 1800,
  height = 1800,
  res = 300
)

grid.draw(venn_prog)

dev.off()

# Integración de la información molecular
tabla_dge_prog <- dge_progresivo[
  symbol %in% master_progresivo,
  .(symbol, logFC, adj.P.Val)
]

info_ds_prog <- ds_progresivo[
  ,
  .(Gen, ID_Evento, Nombre_Biologico, Delta)
]

setnames(
  info_ds_prog,
  "Gen",
  "symbol"
)

master_progresivo_tabla <- merge(
  tabla_dge_prog,
  info_ds_prog,
  by = "symbol",
  allow.cartesian = TRUE
)

master_progresivo_tabla <- master_progresivo_tabla[
  order(adj.P.Val)
]

saveRDS(
  master_progresivo,
  "TFM/procesados/master_progresivo.rds"
)

write.csv(
  master_progresivo_tabla,
  "TFM/procesados/master_progresivo_tabla.csv",
  row.names = FALSE
)

# ==============================================================================
# 5. RESUMEN DEL ANÁLISIS
# ==============================================================================

message("Resumen de la integración multiómica:")

message(
  "Contraste Stage I vs Stage IV - Genes con splicing diferencial: ",
  length(genes_ds_ext)
)

message(
  "Contraste Stage I vs Stage IV - Genes con expresión diferencial: ",
  length(genes_dge_ext)
)

message(
  "Contraste Stage I vs Stage IV - Genes compartidos: ",
  length(master_extremos)
)

message(
  "Contraste Early vs Late - Genes con splicing diferencial: ",
  length(genes_ds_prog)
)

message(
  "Contraste Early vs Late - Genes con expresión diferencial: ",
  length(genes_dge_prog)
)

message(
  "Contraste Early vs Late - Genes compartidos: ",
  length(master_progresivo)
)


message("Proceso completado correctamente.")