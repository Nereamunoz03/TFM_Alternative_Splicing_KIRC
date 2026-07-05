# ==============================================================================
# Trabajo Fin de Máster
# Máster Universitario en Bioinformática - UNIR
#
# Título: Análisis del splicing alternativo asociado a la progresión tumoral
#          en carcinoma renal de células claras mediante RNA-seq
#
# Script: 08_splicing_event_annotation.R
# Correspondiente al apartado: 3.4.3. Anotación funcional y distribución de los
#                              tipos de splicing.
# Descripción: Anotación biológica de los eventos de splicing diferencial
#              mediante la nomenclatura TCGA-SpliceSeq y preparación de los
#              resultados para su representación gráfica.
# ==============================================================================

# ==============================================================================
# 1. CARGA DE LIBRERÍAS
# ==============================================================================

if (!require("data.table", quietly = TRUE)) install.packages("data.table")
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(data.table)
library(ggplot2)

# ==============================================================================
# 2. CARGA DE DATOS
# ==============================================================================

message("Cargando resultados de splicing diferencial...")

ruta_ext  <- "TFM/procesados/ds_extremos.rds"
ruta_prog <- "TFM/procesados/ds_progresivo.rds"

if (file.exists(ruta_ext) & file.exists(ruta_prog)) {
  
  ds_extremos <- readRDS(ruta_ext)
  ds_progresivo <- readRDS(ruta_prog)
  
  setDT(ds_extremos)
  setDT(ds_progresivo)
  
  message("Datos cargados correctamente.")
  
} else {
  
  stop("Error: No se encuentran los archivos 'ds_extremos.rds' o 'ds_progresivo.rds'. Ejecuta previamente el Script 06.")
  
}

# ==============================================================================
# 3. DICCIONARIO DE ANOTACIÓN DE EVENTOS DE SPLICING
# ==============================================================================

diccionario_splicing <- data.frame(
  Tipo = c("ES", "AP", "AT", "RI", "AA", "AD", "ME"),
  Nombre_Biologico = c(
    "Exon Skipping (ES)",
    "Alternate Promoter (AP)",
    "Alternate Terminator (AT)",
    "Intron Retention (RI)",
    "Alternate Acceptor (AA)",
    "Alternate Donor (AD)",
    "Mutually Exclusive Exons (ME)"
  )
)

setDT(diccionario_splicing)

# ==============================================================================
# 4. ANOTACIÓN DE LOS EVENTOS DE SPLICING
# ==============================================================================

message("Anotando los eventos de splicing diferencial...")

if ("Nombre_Biologico" %in% colnames(ds_extremos)) {
  ds_extremos[, Nombre_Biologico := NULL]
}

if ("Nombre_Biologico" %in% colnames(ds_progresivo)) {
  ds_progresivo[, Nombre_Biologico := NULL]
}

ds_extremos <- merge(
  ds_extremos,
  diccionario_splicing,
  by = "Tipo",
  all.x = TRUE
)

ds_progresivo <- merge(
  ds_progresivo,
  diccionario_splicing,
  by = "Tipo",
  all.x = TRUE
)

ds_extremos[
  is.na(Nombre_Biologico),
  Nombre_Biologico := "Otros eventos"
]

ds_progresivo[
  is.na(Nombre_Biologico),
  Nombre_Biologico := "Otros eventos"
]

ds_extremos <- ds_extremos[
  order(abs(Delta), decreasing = TRUE),
]

ds_progresivo <- ds_progresivo[
  order(abs(Delta), decreasing = TRUE),
]

columnas_orden <- c(
  "ID_Evento",
  "Gen",
  "Tipo",
  "Nombre_Biologico",
  "Delta"
)

resto_cols_ext <- setdiff(
  colnames(ds_extremos),
  columnas_orden
)

resto_cols_prog <- setdiff(
  colnames(ds_progresivo),
  columnas_orden
)

setcolorder(
  ds_extremos,
  c(columnas_orden, resto_cols_ext)
)

setcolorder(
  ds_progresivo,
  c(columnas_orden, resto_cols_prog)
)

# ==============================================================================
# 5. EXPORTACIÓN DE RESULTADOS
# ==============================================================================

saveRDS(
  ds_extremos,
  "TFM/procesados/ds_extremos.rds"
)

saveRDS(
  ds_progresivo,
  "TFM/procesados/ds_progresivo.rds"
)

message("Resultados anotados y almacenados correctamente.")


# ==============================================================================
# 6. GENERACIÓN DE LOS GRÁFICOS DE DISTRIBUCIÓN
# ==============================================================================

plot_distribucion_splicing <- function(df, tipo_contraste) {
  
  # Recuento de eventos por tipo de splicing
  df_counts <- df[, .N, by = Nombre_Biologico]
  
  # Número total de eventos diferencialmente regulados
  gran_total <- sum(df_counts$N)
  
  texto_subtitulo <- paste0(
    "Contraste ",
    tipo_contraste,
    " (KIRC) | Total de eventos significativos = ",
    gran_total
  )
  
  # Representación gráfica
  p <- ggplot(
    df_counts,
    aes(
      x = reorder(Nombre_Biologico, N),
      y = N,
      fill = Nombre_Biologico
    )
  ) +
    geom_bar(
      stat = "identity",
      alpha = 0.85,
      color = "black",
      width = 0.65
    ) +
    geom_text(
      aes(label = N),
      hjust = -0.3,
      color = "black",
      size = 3.5,
      fontface = "bold"
    ) +
    coord_flip() +
    theme_minimal() +
    scale_fill_viridis_d(
      option = "plasma",
      direction = -1
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.12))
    ) +
    labs(
      title = "Distribución de los tipos de eventos de splicing alternativo",
      subtitle = texto_subtitulo,
      x = "",
      y = "Número de eventos diferenciales (|ΔΨ| > 0.10)"
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10, face = "italic", color = "grey25"),
      axis.text = element_text(size = 9, color = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    )
  
  return(p)
  
}

# ==============================================================================
# 7. EXPORTACIÓN DE LAS FIGURAS
# ==============================================================================

message("Generando gráficos de distribución de eventos de splicing...")

if (!dir.exists("TFM/graficos")) {
  dir.create("TFM/graficos", recursive = TRUE)
}

grafico_ext <- plot_distribucion_splicing(
  ds_extremos,
  "Extremo: Stage I vs Stage IV"
)

grafico_prog <- plot_distribucion_splicing(
  ds_progresivo,
  "Progresivo: Early vs Late"
)

ggsave(
  "TFM/graficos/Distribucion_Splicing_Extremos.png",
  grafico_ext,
  width = 8.5,
  height = 4.5,
  dpi = 300
)

ggsave(
  "TFM/graficos/Distribucion_Splicing_Progresivo.png",
  grafico_prog,
  width = 8.5,
  height = 4.5,
  dpi = 300
)

message("Proceso completado correctamente.")