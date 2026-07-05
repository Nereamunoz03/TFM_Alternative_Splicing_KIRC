# ==============================================================================
# Trabajo Fin de Máster
# Máster Universitario en Bioinformática - UNIR
#
# Título: Análisis del splicing alternativo asociado a la progresión tumoral
#          en carcinoma renal de células claras mediante RNA-seq
#
# Script: 04_differential_gene_expression_limma_voom.R
# Correspondiente al apartado: 3.3.3. Normalización de datos, análisis
#                              exploratorio transcriptómico y modelado
#                              estadístico mediante Limma-Voom.
# Descripción: Normalización de la matriz de expresión génica, análisis
#              exploratorio mediante PCA y análisis de expresión génica
#              diferencial utilizando el método Limma-Voom.
# ==============================================================================

# ==============================================================================
# 1. CARGA DE LIBRERÍAS
# ==============================================================================

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("edgeR", quietly = TRUE)) BiocManager::install("edgeR")
if (!require("limma", quietly = TRUE)) BiocManager::install("limma")

library(data.table)
library(edgeR)
library(limma)

# ==============================================================================
# 2. CARGA Y PREPARACIÓN DE LOS DATOS
# ==============================================================================

message("Cargando datos de expresión génica y metadatos clínicos...")

ruta_counts  <- "TFM/procesados/matriz_counts_definitiva.rds"
ruta_clinica <- "TFM/procesados/clinica_filtrada_estadios.rds"

if (file.exists(ruta_counts) & file.exists(ruta_clinica)) {
  
  matriz_counts_definitiva <- readRDS(ruta_counts)
  clinica_filtrada_estadios <- readRDS(ruta_clinica)
  
  setDT(matriz_counts_definitiva)
  
  message("Datos cargados correctamente.")
  
} else {
  
  stop("Error: No se encuentran los archivos necesarios en 'TFM/procesados/'. Ejecuta previamente los Scripts 02 y 03.")
  
}

# Conversión de la matriz de conteos a formato numérico
counts_mat <- as.matrix(
  matriz_counts_definitiva[, -1, with = FALSE]
)

rownames(counts_mat) <- matriz_counts_definitiva$symbol

# Sincronización de los datos clínicos con la matriz de expresión
clinica_modelado <- clinica_filtrada_estadios[
  colnames(counts_mat),
]

# Definición de la variable de progresión tumoral
clinica_modelado$Progression <- ifelse(
  clinica_modelado$Pathologic_Tumor_Stage %in% c("Stage I", "Stage II"),
  "Early",
  "Late"
)

# ==============================================================================
# 3. ANÁLISIS EXPLORATORIO MEDIANTE PCA
# ==============================================================================

if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(ggplot2)

# Normalización de los datos mediante TMM
dge_pca <- DGEList(counts = counts_mat)
dge_pca <- calcNormFactors(dge_pca)

design_pca <- model.matrix(
  ~1,
  data = clinica_modelado
)

# Transformación Voom
v_pca <- voom(
  dge_pca,
  design_pca,
  plot = FALSE
)

# Selección de los genes con mayor variabilidad
vars <- apply(
  v_pca$E,
  1,
  var
)

top_genes <- names(
  sort(vars, decreasing = TRUE)
)[1:2000]

# Cálculo del análisis de componentes principales
pca <- prcomp(
  t(v_pca$E[top_genes, ]),
  scale. = TRUE
)

# Varianza explicada
var_exp <- round(
  100 * (pca$sdev^2 / sum(pca$sdev^2)),
  1
)

# Preparación de los datos para la representación gráfica
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Stage = clinica_modelado$Pathologic_Tumor_Stage,
  Progression = clinica_modelado$Progression
)

# Creación del directorio de salida
if (!dir.exists("TFM/graficos")) {
  dir.create("TFM/graficos", recursive = TRUE)
}

# Representación del PCA
pca_plot <- ggplot(
  pca_df,
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
    title = "Análisis de componentes principales basado en los datos de expresión génica de la cohorte TCGA-KIRC",
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)"),
    color = "Progresión"
  )

ggsave(
  "TFM/graficos/PCA_TCGA_KIRC.png",
  pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

message("PCA generado correctamente.")


# ==============================================================================
# 4. ANÁLISIS DE EXPRESIÓN GÉNICA DIFERENCIAL
# ==============================================================================

ejecutar_dge_voom <- function(matriz, vector_grupos, niveles) {
  
  # Selección de las muestras correspondientes al contraste
  idx <- vector_grupos %in% niveles
  sub_matriz <- matriz[, idx]
  
  # Adaptación de los nombres de los grupos para su utilización en limma
  sub_grupos_raw <- gsub(" ", "_", vector_grupos[idx])
  niveles_validos <- gsub(" ", "_", niveles)
  
  sub_grupos <- factor(
    sub_grupos_raw,
    levels = niveles_validos
  )
  
  # Normalización mediante TMM
  dge <- DGEList(counts = sub_matriz)
  dge <- calcNormFactors(dge)
  
  # Construcción de la matriz de diseño
  design <- model.matrix(~0 + sub_grupos)
  colnames(design) <- niveles_validos
  
  # Transformación Voom
  v <- voom(
    dge,
    design,
    plot = FALSE
  )
  
  # Ajuste del modelo lineal
  fit <- lmFit(v, design)
  
  # Definición del contraste
  contraste_formula <- paste0(
    niveles_validos[2],
    " - ",
    niveles_validos[1]
  )
  
  cont_matrix <- makeContrasts(
    contrasts = contraste_formula,
    levels = design
  )
  
  fit2 <- contrasts.fit(
    fit,
    cont_matrix
  )
  
  fit2 <- eBayes(fit2)
  
  # Obtención de la tabla completa de resultados
  res <- topTable(
    fit2,
    coef = 1,
    number = Inf,
    sort.by = "P"
  )
  
  res$symbol <- rownames(res)
  
  return(res)
  
}

# ==============================================================================
# 5. EJECUCIÓN DE LOS CONTRASTES
# ==============================================================================

message("Ejecutando el contraste extremo (Stage I vs Stage IV)...")

dge_extremos <- ejecutar_dge_voom(
  counts_mat,
  clinica_modelado$Pathologic_Tumor_Stage,
  c("Stage I", "Stage IV")
)

message("Ejecutando el contraste progresivo (Early vs Late)...")

dge_progresivo <- ejecutar_dge_voom(
  counts_mat,
  clinica_modelado$Progression,
  c("Early", "Late")
)

# ==============================================================================
# 6. EXPORTACIÓN DE RESULTADOS
# ==============================================================================

if (!dir.exists("TFM/procesados")) {
  dir.create("TFM/procesados", recursive = TRUE)
}

saveRDS(
  dge_extremos,
  "TFM/procesados/dge_extremos.rds"
)

saveRDS(
  dge_progresivo,
  "TFM/procesados/dge_progresivo.rds"
)

# Resumen de genes diferencialmente expresados
hits_ext <- dge_extremos[
  dge_extremos$adj.P.Val < 0.05 &
    abs(dge_extremos$logFC) > 1,
]

hits_pro <- dge_progresivo[
  dge_progresivo$adj.P.Val < 0.05 &
    abs(dge_progresivo$logFC) > 1,
]

message("Proceso completado correctamente.")
message(paste("Genes diferencialmente expresados (Stage I vs Stage IV):", nrow(hits_ext)))
message(paste("Genes diferencialmente expresados (Early vs Late):", nrow(hits_pro)))

# Resumen de las principales dianas de interés
dianas_secundarias <- dge_extremos[
  dge_extremos$symbol %in% c("FOXP3", "PAEP"),
  c("symbol", "logFC", "adj.P.Val")
]

if (nrow(dianas_secundarias) > 0) {
  
  message("Resumen de las dianas FOXP3 y PAEP:")
  print(dianas_secundarias)
  
} else {
  
  message("Las dianas FOXP3 y PAEP no se localizaron en los resultados del contraste extremo.")
  
}
