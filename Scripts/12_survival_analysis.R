# ==============================================================================
# Trabajo Fin de Máster
# Máster Universitario en Bioinformática - UNIR
#
# Título: Análisis del splicing alternativo asociado a la progresión tumoral
#          en carcinoma renal de células claras mediante RNA-seq
#
# Script: 12_survival_analysis.R
# Correspondiente al apartado: 3.7. Análisis de supervivencia clínica.
# Descripción: Preparación de la información clínica de supervivencia mediante
#              la integración de los datos de seguimiento de TCGA para su
#              posterior análisis de Kaplan-Meier.
# ==============================================================================

# ==============================================================================
# 1. CARGA DE LIBRERÍAS
# ==============================================================================

if (!require("survival", quietly = TRUE)) install.packages("survival")
if (!require("survminer", quietly = TRUE)) install.packages("survminer")
if (!require("data.table", quietly = TRUE)) install.packages("data.table")
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(data.table)
library(survival)
library(survminer)
library(ggplot2)

# ==============================================================================
# 2. CARGA DE DATOS
# ==============================================================================

message("Cargando matriz de expresión génica...")

ruta_counts <- "TFM/procesados/matriz_counts_definitiva.rds"

if (file.exists(ruta_counts)) {
  
  matriz_counts_definitiva <- readRDS(ruta_counts)
  setDT(matriz_counts_definitiva)
  
} else {
  
  stop("Error: No se encuentra 'matriz_counts_definitiva.rds' en 'TFM/procesados/'. Ejecuta previamente el Script 03.")
  
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
# 3. CARGA Y PREPARACIÓN DE LA INFORMACIÓN CLÍNICA
# ==============================================================================

message("Cargando información clínica de TCGA...")

ruta_tar_gz <- "C:/Users/PC/Downloads/TFM/DATOS/TCGA_KIRC/clinical.cart.2026-04-07.tar.gz"

if (!file.exists(ruta_tar_gz)) {
  
  stop("Error: No se encuentra el archivo clínico especificado.")
  
}

temp_dir <- tempdir()

untar(
  ruta_tar_gz,
  files = "clinical.tsv",
  exdir = temp_dir
)

ruta_descomprimida <- file.path(
  temp_dir,
  "clinical.tsv"
)

if (!file.exists(ruta_descomprimida)) {
  
  stop("Error: No se pudo extraer el archivo 'clinical.tsv'.")
  
}

clinica_raw <- fread(
  ruta_descomprimida,
  header = TRUE,
  sep = "\t",
  na.strings = "'--'"
)

clinica_raw[
  ,
  Sample_Limpio := substr(demographic.submitter_id, 1, 12)
]

clinica_recortada <- clinica_raw[
  ,
  .(
    Sample = Sample_Limpio,
    Vital_Status = demographic.vital_status,
    Days_Death = as.numeric(demographic.days_to_death),
    Days_FollowUp = as.numeric(diagnoses.days_to_last_follow_up)
  )
]

clinica_clean <- clinica_recortada[
  ,
  .(
    Vital_Status = first(Vital_Status),
    Days_Death = max(Days_Death, na.rm = TRUE),
    Days_FollowUp = max(Days_FollowUp, na.rm = TRUE)
  ),
  by = Sample
]

clinica_clean[
  is.infinite(Days_Death),
  Days_Death := NA
]

clinica_clean[
  is.infinite(Days_FollowUp),
  Days_FollowUp := NA
]

# ==============================================================================
# 4. PREPARACIÓN DE LOS DATOS DE SUPERVIVENCIA
# ==============================================================================

message("Calculando tiempos de seguimiento...")

clinica_clean[
  ,
  Evento := ifelse(Vital_Status == "Dead", 1, 0)
]

clinica_clean[
  ,
  Tiempo_Dias := ifelse(
    Evento == 1,
    Days_Death,
    Days_FollowUp
  )
]

clinica_clean[
  ,
  Tiempo_Meses := Tiempo_Dias / 30.43
]

clinica_surv_definitiva <- clinica_clean[
  !is.na(Tiempo_Meses) & Tiempo_Meses > 0,
  .(
    Sample,
    Evento,
    Tiempo_Meses
  )
]


# ==============================================================================
# 5. INTEGRACIÓN DE LOS DATOS DE EXPRESIÓN Y SUPERVIVENCIA
# ==============================================================================

message("Integrando expresión génica y datos de supervivencia...")

counts_candidatos <- matriz_counts_definitiva[
  symbol %in% genes_candidatos
]

df_genes <- as.data.frame(
  t(counts_candidatos[, -1, with = FALSE])
)

colnames(df_genes) <- counts_candidatos$symbol

df_genes$Sample <- substr(
  rownames(df_genes),
  1,
  12
)

setDT(df_genes)

for (gen in genes_candidatos) {
  
  df_genes[[gen]] <- log2(df_genes[[gen]] + 1)
  
}

df_surv_final <- merge(
  clinica_surv_definitiva,
  df_genes,
  by = "Sample"
)

message(
  "Pacientes disponibles para el análisis: ",
  nrow(df_surv_final)
)

# ==============================================================================
# 6. ANÁLISIS DE SUPERVIVENCIA DE KAPLAN-MEIER
# ==============================================================================

if (!dir.exists("TFM/graficos/supervivencia")) {
  
  dir.create(
    "TFM/graficos/supervivencia",
    recursive = TRUE
  )
  
}

message("Generando curvas de supervivencia...")

for (gen in genes_candidatos) {
  
  dt_bucle <- copy(df_surv_final)
  
  mediana_gen <- median(
    dt_bucle[[gen]],
    na.rm = TRUE
  )
  
  dt_bucle[
    ,
    Grupo_Expresion := ifelse(
      get(gen) >= mediana_gen,
      "Alta Expresión",
      "Baja Expresión"
    )
  ]
  
  fit <- survfit(
    Surv(Tiempo_Meses, Evento) ~ Grupo_Expresion,
    data = dt_bucle
  )
  
  km_plot <- ggsurvplot(
    fit,
    data = dt_bucle,
    pval = TRUE,
    pval.method = TRUE,
    conf.int = FALSE,
    risk.table = TRUE,
    risk.table.col = "strata",
    palette = c("#d73027", "#4575b4"),
    legend.labs = c("Alta Expresión", "Baja Expresión"),
    legend.title = paste("Nivel de", gen),
    title = paste(
      "Supervivencia Global (Kaplan-Meier):",
      gen
    ),
    xlab = "Tiempo transcurrido (meses post-diagnóstico)",
    ylab = "Probabilidad de supervivencia global",
    ggtheme = theme_minimal()
  )
  
  nombre_archivo <- paste0(
    "TFM/graficos/supervivencia/KM_Real_",
    gen,
    ".png"
  )
  
  png(
    nombre_archivo,
    width = 2400,
    height = 2000,
    res = 300
  )
  
  print(km_plot)
  
  dev.off()
  
  message(
    "Curva generada para el gen: ",
    gen
  )
  
}
