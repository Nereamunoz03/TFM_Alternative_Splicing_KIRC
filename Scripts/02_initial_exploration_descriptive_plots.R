# ==============================================================================
# Trabajo Fin de Máster
# Máster Universitario en Bioinformática - UNIR
#
# Título: Análisis del splicing alternativo asociado a la progresión tumoral
#          en carcinoma renal de células claras mediante RNA-seq
#
# Script: 02_initial_exploration_descriptive_plots.R
# Correspondiente al apartado: 3.2.1. Análisis exploratorio inicial y
#                              caracterización de la cohorte.
# Descripción: Exploración descriptiva de la cohorte TCGA-KIRC, generación de
#              gráficos demográficos y filtrado de registros clínicos incompletos.
# ==============================================================================

# ==============================================================================
# 1. CARGA DE LIBRERÍAS Y DATOS
# ==============================================================================

library(ggplot2)

# Carga independiente de la matriz clínica
if (!file.exists("TFM/procesados/clinica_final.rds")) {
  stop("Error: No se encuentra el archivo 'clinica_final.rds'. Ejecuta primero el Script 01.")
}

clinica <- readRDS("TFM/procesados/clinica_final.rds")

# Creación del directorio de salida
if (!dir.exists("TFM/graficos/descriptivos")) {
  dir.create("TFM/graficos/descriptivos", recursive = TRUE)
}

# Ordenación de los niveles del estadio patológico
clinica$Pathologic_Tumor_Stage <- factor(
  clinica$Pathologic_Tumor_Stage,
  levels = c("Stage I", "Stage II", "Stage III", "Stage IV", "null")
)

# ==============================================================================
# 2. RESUMEN DESCRIPTIVO DE LA COHORTE
# ==============================================================================

message("Resumen descriptivo de la cohorte:")

print("Distribución de Edad:")
print(summary(clinica$Age))

print("Distribución de Género:")
print(table(clinica$Gender))

print("Distribución de Estadio Patológico:")
print(table(clinica$Pathologic_Tumor_Stage))

# ==============================================================================
# 3. GENERACIÓN DE GRÁFICOS DESCRIPTIVOS
# ==============================================================================

message("Generando gráficos descriptivos...")

# Distribución por sexo
ggplot(clinica, aes(x = Gender, fill = Gender)) +
  geom_bar() +
  theme_minimal() +
  labs(
    title = "Distribución por Sexo en la Cohorte KIRC",
    subtitle = paste("Total:", nrow(clinica), "pacientes"),
    x = "Género",
    y = "Número de Pacientes"
  ) +
  scale_fill_manual(values = c(
    "FEMALE" = "#f8766d",
    "MALE" = "#00bfc4"
  ))

ggsave(
  "TFM/graficos/descriptivos/01_distribucion_genero.png",
  width = 8,
  height = 6,
  dpi = 300
)

# Distribución por edad
ggplot(clinica, aes(x = Age)) +
  geom_histogram(
    binwidth = 5,
    fill = "steelblue",
    color = "white"
  ) +
  theme_minimal() +
  labs(
    title = "Distribución de Edad en la Cohorte KIRC",
    x = "Edad (años)",
    y = "Frecuencia"
  )

ggsave(
  "TFM/graficos/descriptivos/02_distribucion_edad.png",
  width = 8,
  height = 6,
  dpi = 300
)

# Distribución inicial por estadio patológico
ggplot(clinica, aes(x = Pathologic_Tumor_Stage, fill = Pathologic_Tumor_Stage)) +
  geom_bar() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
    title = "Distribución por Estadio Patológico (Original)",
    subtitle = paste("Total inicial:", nrow(clinica), "pacientes"),
    x = "Estadio Patológico",
    y = "Número de Pacientes"
  ) +
  scale_fill_manual(values = c(
    "Stage I"   = "#fee0d2",
    "Stage II"  = "#fc9272",
    "Stage III" = "#de2d26",
    "Stage IV"  = "#a50f15",
    "null"      = "#969696"
  ))

ggsave(
  "TFM/graficos/descriptivos/03_estadios_originales.png",
  width = 8,
  height = 6,
  dpi = 300
)

# ==============================================================================
# 4. DEPURACIÓN DE LA VARIABLE DE ESTADIO PATOLÓGICO
# ==============================================================================

message("Depurando la variable de estadio patológico...")

# Eliminación de registros con estadio patológico no disponible
clinica_filtrada <- subset(
  clinica,
  Pathologic_Tumor_Stage != "null"
)

# Reordenación de los niveles del estadio patológico
clinica_filtrada$Pathologic_Tumor_Stage <- factor(
  clinica_filtrada$Pathologic_Tumor_Stage,
  levels = c("Stage I", "Stage II", "Stage III", "Stage IV")
)

# ==============================================================================
# 5. VISUALIZACIÓN DE LA COHORTE DEPURADA
# ==============================================================================

# Distribución final por estadio patológico
ggplot(clinica_filtrada, aes(x = Pathologic_Tumor_Stage, fill = Pathologic_Tumor_Stage)) +
  geom_bar() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
    title = "Distribución por Estadio Patológico (Curado y Ordenado)",
    subtitle = paste("Total válido:", nrow(clinica_filtrada), "pacientes"),
    x = "Estadio Patológico",
    y = "Número de Pacientes"
  ) +
  scale_fill_manual(values = c(
    "Stage I"   = "#fee0d2",
    "Stage II"  = "#fc9272",
    "Stage III" = "#de2d26",
    "Stage IV"  = "#a50f15"
  ))

ggsave(
  "TFM/graficos/descriptivos/04_estadios_curados.png",
  width = 8,
  height = 6,
  dpi = 300
)

# ==============================================================================
# 6. EXPORTACIÓN DE RESULTADOS
# ==============================================================================

# Exportación de la matriz clínica depurada
saveRDS(
  clinica_filtrada,
  "TFM/procesados/clinica_filtrada_estadios.rds"
)

message("Proceso completado correctamente.")
message(paste("Pacientes iniciales:", nrow(clinica)))
message(paste("Pacientes incluidos en el análisis:", nrow(clinica_filtrada)))
