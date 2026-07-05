# ==============================================================================
# Trabajo Fin de Máster
# Máster Universitario en Bioinformática - UNIR
#
# Título: Análisis del splicing alternativo asociado a la progresión tumoral
#          en carcinoma renal de células claras mediante RNA-seq
#
# Script: 05_gene_expression_visualization.R
# Correspondiente al apartado: 3.3.4. Visualización de los perfiles de
#                              transcripción.
# Descripción: Generación de Volcano Plots y preparación de los perfiles de
#              expresión génica de los genes candidatos para su visualización.
# ==============================================================================

# ==============================================================================
# 1. CARGA DE LIBRERÍAS
# ==============================================================================

if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!require("ggrepel", quietly = TRUE)) install.packages("ggrepel")
if (!require("data.table", quietly = TRUE)) install.packages("data.table")
if (!require("reshape2", quietly = TRUE)) install.packages("reshape2")

library(ggplot2)
library(ggrepel)
library(data.table)
library(reshape2)

# ==============================================================================
# 2. CARGA DE DATOS
# ==============================================================================

message("Cargando resultados de expresión génica y datos clínicos...")

ruta_ext     <- "TFM/procesados/dge_extremos.rds"
ruta_prog    <- "TFM/procesados/dge_progresivo.rds"
ruta_counts  <- "TFM/procesados/matriz_counts_definitiva.rds"
ruta_clinica <- "TFM/procesados/clinica_filtrada_estadios.rds"

if (all(file.exists(c(ruta_ext, ruta_prog, ruta_counts, ruta_clinica)))) {
  
  dge_extremos             <- readRDS(ruta_ext)
  dge_progresivo           <- readRDS(ruta_prog)
  matriz_counts_definitiva <- readRDS(ruta_counts)
  clinica_filtrada_estadios <- readRDS(ruta_clinica)
  
  setDT(matriz_counts_definitiva)
  
  message("Datos cargados correctamente.")
  
} else {
  
  stop("Error: No se encuentran todos los archivos necesarios en 'TFM/procesados/'. Ejecuta previamente los scripts correspondientes.")
  
}

# Creación del directorio de salida
if (!dir.exists("TFM/graficos")) {
  dir.create("TFM/graficos", recursive = TRUE)
}

# ==============================================================================
# 3. GENERACIÓN DE VOLCANO PLOTS
# ==============================================================================

plot_volcano_safe <- function(df, titulo) {
  
  # Eliminación de registros sin valor de significación ajustado
  df <- df[!is.na(df$adj.P.Val), ]
  
  # Clasificación de genes según los criterios de significación
  df$diffexpressed <- "NO"
  df$diffexpressed[df$logFC > 1 & df$adj.P.Val < 0.05] <- "UP"
  df$diffexpressed[df$logFC < -1 & df$adj.P.Val < 0.05] <- "DOWN"
  
  # Etiquetado de los genes candidatos
  genes_interes <- c("SAA1", "PITX2", "FOXP3", "PAEP")
  
  df$label <- ifelse(
    df$symbol %in% genes_interes,
    df$symbol,
    NA
  )
  
  # Reducción del número de puntos representados
  df_light <- df[df$adj.P.Val < 0.5, ]
  
  p <- ggplot(
    df_light,
    aes(
      x = logFC,
      y = -log10(adj.P.Val),
      color = diffexpressed
    )
  ) +
    geom_point(
      alpha = 0.4,
      size = 1
    ) +
    scale_color_manual(
      values = c(
        "DOWN" = "royalblue",
        "NO" = "lightgrey",
        "UP" = "firebrick"
      )
    ) +
    geom_text_repel(
      aes(label = label),
      color = "black",
      size = 4,
      max.overlaps = Inf,
      fontface = "bold",
      box.padding = 0.5,
      point.padding = 0.3,
      force = 2
    ) +
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed",
      alpha = 0.5
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      alpha = 0.5
    ) +
    theme_classic() +
    labs(
      title = titulo,
      subtitle = "FDR < 0.05 | |Log2FC| > 1",
      x = "Log2 Fold Change",
      y = "-Log10 P-ajustado"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "right"
    )
  
  return(p)
  
}

message("Generando Volcano Plots...")

volcano_ext <- plot_volcano_safe(
  dge_extremos,
  "Expresión Génica Diferencial (Stage I vs Stage IV)"
)

volcano_prog <- plot_volcano_safe(
  dge_progresivo,
  "Expresión Génica Diferencial (Early vs Late)"
)

ggsave(
  "TFM/graficos/Volcano_Extremos.png",
  volcano_ext,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  "TFM/graficos/Volcano_Progresivo.png",
  volcano_prog,
  width = 7,
  height = 5,
  dpi = 300
)


# ==============================================================================
# 4. PREPARACIÓN DE LOS DATOS PARA LA VISUALIZACIÓN
# ==============================================================================

genes_box <- c("SAA1", "PITX2", "FOXP3", "PAEP")

message("Preparando los perfiles de expresión de los genes candidatos...")

df_genes <- as.data.frame(
  t(
    matriz_counts_definitiva[
      symbol %in% genes_box,
      -1,
      with = FALSE
    ]
  )
)

colnames(df_genes) <- matriz_counts_definitiva[
  symbol %in% genes_box
]$symbol

df_genes$Sample <- rownames(df_genes)

clinica_filtrada_estadios$Progression <- ifelse(
  clinica_filtrada_estadios$Pathologic_Tumor_Stage %in% c("Stage I", "Stage II"),
  "Early",
  "Late"
)

df_plot <- merge(
  df_genes,
  clinica_filtrada_estadios[
    ,
    c("Pathologic_Tumor_Stage", "Progression"),
    drop = FALSE
  ],
  by.x = "Sample",
  by.y = "row.names"
)

df_melt <- melt(
  df_plot,
  id.vars = c(
    "Sample",
    "Pathologic_Tumor_Stage",
    "Progression"
  )
)

# ==============================================================================
# 5. GENERACIÓN DE BOXPLOTS
# ==============================================================================

boxplot_estadios <- ggplot(
  df_melt,
  aes(
    x = Pathologic_Tumor_Stage,
    y = log2(value + 1),
    fill = Pathologic_Tumor_Stage
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7,
    color = "black"
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.15,
    size = 0.4,
    color = "darkgrey"
  ) +
  facet_wrap(
    ~variable,
    scales = "free_y",
    nrow = 2
  ) +
  scale_fill_brewer(palette = "Reds") +
  theme_bw() +
  labs(
    title = "Perfiles de expresión de los genes candidatos por estadio clínico",
    y = "Nivel de expresión Log2(Counts + 1)",
    x = "Estadio patológico"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none",
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold")
  )

boxplot_progresion <- ggplot(
  df_melt,
  aes(
    x = Progression,
    y = log2(value + 1),
    fill = Progression
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7,
    color = "black"
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.15,
    size = 0.4,
    color = "darkgrey"
  ) +
  facet_wrap(
    ~variable,
    scales = "free_y",
    nrow = 2
  ) +
  scale_fill_manual(
    values = c(
      "Early" = "#4daf4a",
      "Late" = "#984ea3"
    )
  ) +
  theme_bw() +
  labs(
    title = "Perfiles de expresión de los genes candidatos por grupo de progresión",
    y = "Nivel de expresión Log2(Counts + 1)",
    x = "Grupo de progresión"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    legend.position = "none",
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold")
  )

# ==============================================================================
# 6. EXPORTACIÓN DE RESULTADOS
# ==============================================================================

ggsave(
  "TFM/graficos/Boxplots_Genes_Clave_Estadios.png",
  boxplot_estadios,
  width = 8,
  height = 7,
  dpi = 300
)

ggsave(
  "TFM/graficos/Boxplots_Genes_Clave_Progresion.png",
  boxplot_progresion,
  width = 8,
  height = 7,
  dpi = 300
)

message("Proceso completado correctamente.")