# ==============================================================================
# Trabajo Fin de Máster
# Máster Universitario en Bioinformática - UNIR
#
# Título: Análisis multiómico de expresión génica y splicing alternativo en 
#         carcinoma renal de células claras (KIRC).
#
# Script: 06_differential_splicing_analysis.R
# Correspondiente al apartado: 3.4.1. Análisis exploratorio y cálculo del 
#                                     splicing alternativo diferencial.
# Descripción: Análisis exploratorio de los perfiles de splicing alternativo
#              mediante PCA y preparación de los datos para el análisis de
#              splicing diferencial.
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

message("Cargando matriz de splicing y metadatos clínicos...")

ruta_psi     <- "TFM/procesados/matriz_psi_final.rds"
ruta_clinica <- "TFM/procesados/clinica_filtrada_estadios.rds"

if (file.exists(ruta_psi) & file.exists(ruta_clinica)) {
  
  matriz_psi_final <- readRDS(ruta_psi)
  clinica_filtrada_estadios <- readRDS(ruta_clinica)
  
  setDT(matriz_psi_final)
  
  message("Datos cargados correctamente.")
  
} else {
  
  stop("Error: No se encuentran los archivos necesarios en 'TFM/procesados/'. Ejecuta previamente los Scripts 01 y 02.")
  
}

# ==============================================================================
# 3. ANÁLISIS EXPLORATORIO MEDIANTE PCA
# ==============================================================================

# Definición de los grupos de progresión
clinica_filtrada_estadios$Progression <- ifelse(
  clinica_filtrada_estadios$Pathologic_Tumor_Stage %in% c("Stage I", "Stage II"),
  "Early",
  "Late"
)

# Extracción de la matriz PSI
psi_mat <- as.matrix(
  matriz_psi_final[, 11:ncol(matriz_psi_final), with = FALSE]
)

colnames(psi_mat) <- colnames(matriz_psi_final)[11:ncol(matriz_psi_final)]

# Sincronización con la cohorte clínica
pacientes_validos <- rownames(clinica_filtrada_estadios)

psi_mat <- psi_mat[, pacientes_validos]

clinica_pca <- clinica_filtrada_estadios[pacientes_validos, ]

# Eliminación de eventos con más de un 20 % de valores ausentes
keep <- rowMeans(is.na(psi_mat)) < 0.20
psi_mat <- psi_mat[keep, ]

# Imputación mediante la media de cada evento
for (i in seq_len(nrow(psi_mat))) {
  
  psi_mat[i, is.na(psi_mat[i, ])] <- mean(
    psi_mat[i, ],
    na.rm = TRUE
  )
  
}

# Selección de los eventos con mayor variabilidad
vars <- apply(
  psi_mat,
  1,
  var
)

top_events <- order(
  vars,
  decreasing = TRUE
)[1:2000]

psi_mat_top <- psi_mat[top_events, ]

# Cálculo del PCA
pca_psi <- prcomp(
  t(psi_mat_top),
  scale. = TRUE
)

# Varianza explicada
var_exp_psi <- round(
  100 * (pca_psi$sdev^2 / sum(pca_psi$sdev^2)),
  1
)

# Preparación de los datos para la representación
pca_psi_df <- data.frame(
  PC1 = pca_psi$x[, 1],
  PC2 = pca_psi$x[, 2],
  Progression = clinica_pca$Progression
)

# Representación gráfica
pca_psi_plot <- ggplot(
  pca_psi_df,
  aes(
    x = PC1,
    y = PC2,
    color = Progression
  )
) +
  geom_point(
    size = 3,
    alpha = 0.6
  ) +
  scale_color_manual(
    values = c(
      "Early" = "#4daf4a",
      "Late" = "#984ea3"
    )
  ) +
  theme_bw() +
  labs(
    title = "Análisis de componentes principales basado en los perfiles de splicing alternativo de la cohorte TCGA-KIRC",
    x = paste0("PC1 (", var_exp_psi[1], "%)"),
    y = paste0("PC2 (", var_exp_psi[2], "%)"),
    color = "Progresión"
  )

ggsave(
  "TFM/graficos/PCA_PSI_KIRC.png",
  pca_psi_plot,
  width = 8,
  height = 6,
  dpi = 300
)

message("PCA generado correctamente.")


# ==============================================================================
# 4. ANÁLISIS DE SPLICING DIFERENCIAL
# ==============================================================================

message("Ejecutando el contraste progresivo (Early vs Late)...")

# Definición de los grupos clínicos
pacientes_early <- rownames(clinica_filtrada_estadios)[
  clinica_filtrada_estadios$Pathologic_Tumor_Stage %in%
    c("Stage I", "Stage II")
]

pacientes_late <- rownames(clinica_filtrada_estadios)[
  clinica_filtrada_estadios$Pathologic_Tumor_Stage %in%
    c("Stage III", "Stage IV")
]

# Cálculo del PSI medio por grupo
media_early <- rowMeans(
  matriz_psi_final[, ..pacientes_early],
  na.rm = TRUE
)

media_late <- rowMeans(
  matriz_psi_final[, ..pacientes_late],
  na.rm = TRUE
)

# Construcción de la tabla de resultados
ds_progresivo <- data.frame(
  ID_Evento   = matriz_psi_final$`as-id`,
  Gen         = matriz_psi_final$symbol,
  Tipo        = matriz_psi_final$`splice-type`,
  Media_Early = media_early,
  Media_Late  = media_late,
  Delta       = media_late - media_early
)

# Selección de los eventos con mayor cambio
ds_progresivo <- ds_progresivo[
  abs(ds_progresivo$Delta) > 0.1 &
    !is.na(ds_progresivo$Delta),
]

ds_progresivo <- ds_progresivo[
  order(abs(ds_progresivo$Delta), decreasing = TRUE),
]

# ==============================================================================
# 5. CONTRASTE EXTREMO
# ==============================================================================

message("Ejecutando el contraste extremo (Stage I vs Stage IV)...")

# Definición de los grupos clínicos
pacientes_S1 <- rownames(clinica_filtrada_estadios)[
  clinica_filtrada_estadios$Pathologic_Tumor_Stage == "Stage I"
]

pacientes_S4 <- rownames(clinica_filtrada_estadios)[
  clinica_filtrada_estadios$Pathologic_Tumor_Stage == "Stage IV"
]

# Cálculo del PSI medio por grupo
media_S1 <- rowMeans(
  matriz_psi_final[, ..pacientes_S1],
  na.rm = TRUE
)

media_S4 <- rowMeans(
  matriz_psi_final[, ..pacientes_S4],
  na.rm = TRUE
)

# Construcción de la tabla de resultados
ds_extremos <- data.frame(
  ID_Evento = matriz_psi_final$`as-id`,
  Gen       = matriz_psi_final$symbol,
  Tipo      = matriz_psi_final$`splice-type`,
  Media_S1  = media_S1,
  Media_S4  = media_S4,
  Delta     = media_S4 - media_S1
)

# Selección de los eventos con mayor cambio
ds_extremos <- ds_extremos[
  abs(ds_extremos$Delta) > 0.1 &
    !is.na(ds_extremos$Delta),
]

ds_extremos <- ds_extremos[
  order(abs(ds_extremos$Delta), decreasing = TRUE),
]

# ==============================================================================
# 6. EXPORTACIÓN DE RESULTADOS
# ==============================================================================

if (!dir.exists("TFM/procesados")) {
  dir.create("TFM/procesados", recursive = TRUE)
}

saveRDS(
  ds_progresivo,
  "TFM/procesados/ds_progresivo.rds"
)

saveRDS(
  ds_extremos,
  "TFM/procesados/ds_extremos.rds"
)

message("Proceso completado correctamente.")
message(paste(
  "Eventos de splicing diferencial (Early vs Late):",
  nrow(ds_progresivo)
))

message(paste(
  "Eventos de splicing diferencial (Stage I vs Stage IV):",
  nrow(ds_extremos)
))
