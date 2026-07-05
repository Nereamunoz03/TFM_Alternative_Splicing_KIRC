# ==============================================================================
# Trabajo Fin de Máster
# Máster Universitario en Bioinformática - UNIR
#
# Título: Análisis del splicing alternativo asociado a la progresión tumoral
#          en carcinoma renal de células claras mediante RNA-seq
#
# Script: 07_differential_splicing_visualization.R
# Correspondiente al apartado: 3.4.4. Visualización de los eventos candidatos.
# Descripción: Generación de Volcano Plots y mapas de calor para la
#              representación de los principales eventos de splicing
#              diferencial identificados en ambos contrastes clínicos.
# ==============================================================================

# ==============================================================================
# 1. CARGA DE LIBRERÍAS
# ==============================================================================

if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!require("ggrepel", quietly = TRUE)) install.packages("ggrepel")
if (!require("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!require("data.table", quietly = TRUE)) install.packages("data.table")

library(ggplot2)
library(ggrepel)
library(pheatmap)
library(data.table)

# ==============================================================================
# 2. CARGA DE DATOS
# ==============================================================================

message("Cargando resultados de splicing diferencial...")

rutas <- c(
  "TFM/procesados/matriz_psi_final.rds",
  "TFM/procesados/clinica_filtrada_estadios.rds",
  "TFM/procesados/ds_extremos.rds",
  "TFM/procesados/ds_progresivo.rds"
)

if (all(file.exists(rutas))) {
  
  matriz_psi_final <- readRDS(rutas[1])
  clinica_filtrada_estadios <- readRDS(rutas[2])
  ds_extremos <- readRDS(rutas[3])
  ds_progresivo <- readRDS(rutas[4])
  
  setDT(matriz_psi_final)
  
  message("Datos cargados correctamente.")
  
} else {
  
  stop("Error: No se encuentran todos los archivos necesarios en 'TFM/procesados/'. Ejecuta previamente los scripts anteriores.")
  
}

if (!dir.exists("TFM/graficos")) {
  dir.create("TFM/graficos", recursive = TRUE)
}

genes_candidatos <- c(
  "GNAS",
  "GRAMD1C",
  "MYH10",
  "APOLD1",
  "GLUL"
)

# ==============================================================================
# 3. VISUALIZACIÓN DEL CONTRASTE EXTREMO
# ==============================================================================

ds_extremos$color <- "Estable"
ds_extremos$color[ds_extremos$Delta > 0.1] <- "Sube en Stage IV"
ds_extremos$color[ds_extremos$Delta < -0.1] <- "Baja en Stage IV"

ds_extremos$etiqueta <- ifelse(
  ds_extremos$Gen %in% genes_candidatos &
    ds_extremos$color != "Estable",
  ds_extremos$Gen,
  ""
)

v_ext <- ggplot(
  ds_extremos,
  aes(
    x = Delta,
    y = abs(Delta) * 10,
    color = color
  )
) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_text_repel(
    aes(label = etiqueta),
    color = "black",
    size = 3.5,
    max.overlaps = 15,
    fontface = "bold",
    box.padding = 0.5
  ) +
  scale_color_manual(
    values = c(
      "Baja en Stage IV" = "royalblue",
      "Estable" = "grey85",
      "Sube en Stage IV" = "firebrick"
    )
  ) +
  theme_minimal() +
  geom_vline(
    xintercept = c(-0.1, 0.1),
    linetype = "dashed",
    alpha = 0.5
  ) +
  labs(
    title = "Volcano Plot de Splicing (Stage I vs Stage IV)",
    x = "Delta PSI (ΔΨ)",
    y = "Magnitud del efecto (|ΔΨ| × 10)",
    color = "Estado"
  )

ggsave(
  "TFM/graficos/Volcano_Splicing_Extremos.png",
  v_ext,
  width = 7,
  height = 5,
  dpi = 300
)

pacientes_S1 <- intersect(
  rownames(clinica_filtrada_estadios)[
    clinica_filtrada_estadios$Pathologic_Tumor_Stage == "Stage I"
  ],
  colnames(matriz_psi_final)
)

pacientes_S4 <- intersect(
  rownames(clinica_filtrada_estadios)[
    clinica_filtrada_estadios$Pathologic_Tumor_Stage == "Stage IV"
  ],
  colnames(matriz_psi_final)
)

top_eventos_ext <- head(ds_extremos$ID_Evento, 20)

matriz_plot_ext <- matriz_psi_final[
  matriz_psi_final$`as-id` %in% top_eventos_ext,
]

columnas_ext <- c(pacientes_S1, pacientes_S4)

matriz_normal_ext <- as.data.frame(matriz_plot_ext)

plot_data_ext <- as.matrix(
  matriz_normal_ext[, columnas_ext]
)

rownames(plot_data_ext) <- paste(
  matriz_normal_ext$symbol,
  matriz_normal_ext$`as-id`,
  sep = "_"
)

annotation_ext <- data.frame(
  Stage = factor(
    c(
      rep("Stage I", length(pacientes_S1)),
      rep("Stage IV", length(pacientes_S4))
    )
  )
)

rownames(annotation_ext) <- columnas_ext

png(
  "TFM/graficos/Heatmap_Splicing_Extremos.png",
  width = 2400,
  height = 1800,
  res = 300
)

pheatmap(
  plot_data_ext,
  annotation_col = annotation_ext,
  show_colnames = FALSE,
  main = "Top 20 eventos de splicing diferencial (Stage I vs Stage IV)",
  color = colorRampPalette(
    c("navy", "white", "firebrick3")
  )(50),
  scale = "row",
  na_col = "grey90"
)

dev.off()

# ==============================================================================
# 4. VISUALIZACIÓN DEL CONTRASTE PROGRESIVO
# ==============================================================================

ds_progresivo$color <- "Estable"
ds_progresivo$color[ds_progresivo$Delta > 0.1] <- "Sube en Late"
ds_progresivo$color[ds_progresivo$Delta < -0.1] <- "Baja en Late"

ds_progresivo$etiqueta <- ifelse(
  ds_progresivo$Gen %in% genes_candidatos &
    ds_progresivo$color != "Estable",
  ds_progresivo$Gen,
  ""
)

v_prog <- ggplot(
  ds_progresivo,
  aes(
    x = Delta,
    y = abs(Delta) * 10,
    color = color
  )
) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_text_repel(
    aes(label = etiqueta),
    color = "black",
    size = 3.5,
    max.overlaps = 15,
    fontface = "bold",
    box.padding = 0.5
  ) +
  scale_color_manual(
    values = c(
      "Baja en Late" = "royalblue",
      "Estable" = "grey85",
      "Sube en Late" = "firebrick"
    )
  ) +
  theme_minimal() +
  geom_vline(
    xintercept = c(-0.1, 0.1),
    linetype = "dashed",
    alpha = 0.5
  ) +
  labs(
    title = "Volcano Plot de Splicing (Early vs Late)",
    x = "Delta PSI (ΔΨ)",
    y = "Magnitud del efecto (|ΔΨ| × 10)",
    color = "Estado"
  )

ggsave(
  "TFM/graficos/Volcano_Splicing_Progresivo.png",
  v_prog,
  width = 7,
  height = 5,
  dpi = 300
)

pacientes_early <- intersect(
  rownames(clinica_filtrada_estadios)[
    clinica_filtrada_estadios$Pathologic_Tumor_Stage %in% c("Stage I", "Stage II")
  ],
  colnames(matriz_psi_final)
)

pacientes_late <- intersect(
  rownames(clinica_filtrada_estadios)[
    clinica_filtrada_estadios$Pathologic_Tumor_Stage %in% c("Stage III", "Stage IV")
  ],
  colnames(matriz_psi_final)
)

top_eventos_prog <- head(ds_progresivo$ID_Evento, 20)

matriz_plot_prog <- matriz_psi_final[
  matriz_psi_final$`as-id` %in% top_eventos_prog,
]

columnas_prog <- c(pacientes_early, pacientes_late)

matriz_normal_prog <- as.data.frame(matriz_plot_prog)

plot_data_prog <- as.matrix(
  matriz_normal_prog[, columnas_prog]
)

rownames(plot_data_prog) <- paste(
  matriz_normal_prog$symbol,
  matriz_normal_prog$`as-id`,
  sep = "_"
)

annotation_prog <- data.frame(
  Progression = factor(
    c(
      rep("Early", length(pacientes_early)),
      rep("Late", length(pacientes_late))
    )
  )
)

rownames(annotation_prog) <- columnas_prog

png(
  "TFM/graficos/Heatmap_Splicing_Progresivo.png",
  width = 2400,
  height = 1800,
  res = 300
)

pheatmap(
  plot_data_prog,
  annotation_col = annotation_prog,
  show_colnames = FALSE,
  main = "Top 20 eventos de splicing diferencial (Early vs Late)",
  color = colorRampPalette(
    c("navy", "white", "firebrick3")
  )(50),
  scale = "row",
  na_col = "grey90"
)

dev.off()

message("Proceso completado correctamente.")