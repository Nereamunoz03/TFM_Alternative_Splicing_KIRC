# ==============================================================================
# Trabajo Fin de Máster
# Máster Universitario en Bioinformática - UNIR
#
# Título: Análisis del splicing alternativo asociado a la progresión tumoral
#          en carcinoma renal de células claras mediante RNA-seq
#
# Script: 09_functional_enrichment_analysis.R
# Correspondiente al apartado: 3.4.4. Análisis de enriquecimiento funcional.
# Descripción: Análisis de enriquecimiento funcional mediante Gene Ontology
#              (Procesos Biológicos) y KEGG Pathways a partir de los genes
#              asociados a los eventos de splicing diferencial.
# ==============================================================================

# ==============================================================================
# 1. CARGA DE LIBRERÍAS
# ==============================================================================

if (!require("enrichR", quietly = TRUE)) install.packages("enrichR")
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!require("data.table", quietly = TRUE)) install.packages("data.table")

library(enrichR)
library(ggplot2)
library(data.table)

# ==============================================================================
# 2. CARGA DE DATOS
# ==============================================================================

message("Cargando resultados de splicing diferencial...")

ruta_ext <- "TFM/procesados/ds_extremos.rds"
ruta_prog <- "TFM/procesados/ds_progresivo.rds"

if (file.exists(ruta_ext) & file.exists(ruta_prog)) {
  
  ds_extremos <- readRDS(ruta_ext)
  ds_progresivo <- readRDS(ruta_prog)
  
  message("Datos cargados correctamente.")
  
} else {
  
  stop("Error: No se encuentran los archivos 'ds_extremos.rds' o 'ds_progresivo.rds'. Ejecuta previamente el Script 08.")
  
}

if (!dir.exists("TFM/procesados")) {
  dir.create("TFM/procesados", recursive = TRUE)
}

if (!dir.exists("TFM/graficos")) {
  dir.create("TFM/graficos", recursive = TRUE)
}

# Bases de datos utilizadas en Enrichr
dbs <- c(
  "GO_Biological_Process_2023",
  "KEGG_2021_Human"
)

# ==============================================================================
# 3. FUNCIÓN DE ENRIQUECIMIENTO FUNCIONAL
# ==============================================================================

ejecutar_enriquecimiento_tfm <- function(df_splicing, etiqueta_contraste) {
  
  genes_tfm <- unique(df_splicing$Gen)
  genes_tfm <- genes_tfm[!is.na(genes_tfm) & genes_tfm != ""]
  
  message(
    paste(
      "\nAnalizando el contraste",
      etiqueta_contraste,
      "con",
      length(genes_tfm),
      "genes."
    )
  )
  
  enriquecimiento <- enrichr(
    genes_tfm,
    dbs
  )
  
  #--------------------------------------------------------------------------
  # Gene Ontology
  #--------------------------------------------------------------------------
  
  res_go <- enriquecimiento[["GO_Biological_Process_2023"]]
  
  res_go_sig <- res_go[
    res_go$Adjusted.P.value < 0.05,
  ]
  
  if (nrow(res_go_sig) == 0) {
    
    message(
      paste(
        "No se identificaron términos GO con FDR < 0.05 para",
        etiqueta_contraste,
        ". Se empleará P < 0.05."
      )
    )
    
    res_go_sig <- res_go[
      res_go$P.value < 0.05,
    ]
    
  }
  
  top_go <- head(
    res_go_sig[
      order(res_go_sig$Combined.Score, decreasing = TRUE),
    ],
    10
  )
  
  if (nrow(top_go) > 0) {
    
    top_go$Term <- gsub(
      " \\(GO:.*\\)",
      "",
      top_go$Term
    )
    
    g_go <- ggplot(
      top_go,
      aes(
        x = reorder(Term, Combined.Score),
        y = Combined.Score
      )
    ) +
      geom_bar(
        stat = "identity",
        fill = "indianred3",
        color = "black",
        alpha = 0.85,
        width = 0.7
      ) +
      coord_flip() +
      theme_minimal() +
      labs(
        title = "Procesos biológicos enriquecidos (Gene Ontology)",
        subtitle = paste("Contraste:", etiqueta_contraste),
        x = "",
        y = "Combined Score"
      ) +
      theme(
        axis.text.y = element_text(size = 9, color = "black"),
        plot.title = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank()
      )
    
    ggsave(
      paste0(
        "TFM/graficos/Enriquecimiento_GO_",
        gsub(" ", "_", etiqueta_contraste),
        ".png"
      ),
      g_go,
      width = 10,
      height = 5,
      dpi = 300
    )
    
    write.csv(
      res_go_sig,
      paste0(
        "TFM/procesados/tabla_enriquecimiento_GO_",
        gsub(" ", "_", etiqueta_contraste),
        ".csv"
      ),
      row.names = FALSE
    )
    
  }
  
  #--------------------------------------------------------------------------
  # KEGG Pathways
  #--------------------------------------------------------------------------
  
  res_kegg <- enriquecimiento[["KEGG_2021_Human"]]
  
  res_kegg_sig <- res_kegg[
    res_kegg$Adjusted.P.value < 0.05,
  ]
  
  if (nrow(res_kegg_sig) == 0) {
    
    message(
      paste(
        "No se identificaron rutas KEGG con FDR < 0.05 para",
        etiqueta_contraste,
        ". Se empleará P < 0.05."
      )
    )
    
    res_kegg_sig <- res_kegg[
      res_kegg$P.value < 0.05,
    ]
    
  }
  
  top_kegg <- head(
    res_kegg_sig[
      order(res_kegg_sig$Combined.Score, decreasing = TRUE),
    ],
    10
  )
  
  if (nrow(top_kegg) > 0) {
    
    top_kegg$Term <- gsub(
      " - Homo sapiens \\(human\\)",
      "",
      top_kegg$Term
    )
    
    g_kegg <- ggplot(
      top_kegg,
      aes(
        x = reorder(Term, Combined.Score),
        y = Combined.Score
      )
    ) +
      geom_bar(
        stat = "identity",
        fill = "dodgerblue4",
        color = "black",
        alpha = 0.85,
        width = 0.7
      ) +
      coord_flip() +
      theme_minimal() +
      labs(
        title = "Rutas enriquecidas (KEGG)",
        subtitle = paste("Contraste:", etiqueta_contraste),
        x = "",
        y = "Combined Score"
      ) +
      theme(
        axis.text.y = element_text(size = 9, color = "black"),
        plot.title = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank()
      )
    
    ggsave(
      paste0(
        "TFM/graficos/Enriquecimiento_KEGG_",
        gsub(" ", "_", etiqueta_contraste),
        ".png"
      ),
      g_kegg,
      width = 10,
      height = 5,
      dpi = 300
    )
    
    write.csv(
      res_kegg_sig,
      paste0(
        "TFM/procesados/tabla_enriquecimiento_KEGG_",
        gsub(" ", "_", etiqueta_contraste),
        ".csv"
      ),
      row.names = FALSE
    )
    
  }
  
}

# ==============================================================================
# 4. EJECUCIÓN DEL ANÁLISIS
# ==============================================================================

message("Iniciando el análisis de enriquecimiento funcional...")

ejecutar_enriquecimiento_tfm(
  ds_extremos,
  "Stage I vs Stage IV"
)

ejecutar_enriquecimiento_tfm(
  ds_progresivo,
  "Early vs Late"
)

message("Proceso completado correctamente.")