# Código fuente

Esta carpeta contiene el código fuente completo desarrollado durante el análisis bioinformático del Trabajo Fin de Máster.

El pipeline fue implementado íntegramente en el lenguaje de programación **R** y se estructuró en **12 scripts independientes**, organizados de forma secuencial para facilitar la reproducibilidad del flujo de trabajo. Cada script puede ejecutarse de forma independiente, cargando automáticamente los archivos intermedios (`.rds`) generados en pasos previos. Este diseño permite reproducir cada etapa del análisis sin necesidad de ejecutar nuevamente todo el pipeline.

## Organización de los scripts

| Script | Descripción |
|---------|-------------|
| **Script 01** | Filtrado de la cohorte, procesamiento de las matrices de expresión génica y splicing alternativo, e integración de los datos moleculares y clínicos. |
| **Script 02** | Exploración descriptiva de la cohorte clínica, generación de gráficos descriptivos y depuración de los estadios tumorales. |
| **Script 03** | Traducción de identificadores Ensembl a símbolos génicos y control de calidad de la matriz de expresión génica. |
| **Script 04** | Normalización mediante Limma-Voom, análisis exploratorio mediante PCA y análisis de expresión génica diferencial (DGE). |
| **Script 05** | Visualización de los resultados de expresión génica diferencial mediante *volcano plots* y *boxplots*. |
| **Script 06** | Análisis de splicing alternativo diferencial (DS) y análisis exploratorio mediante PCA basado en valores PSI. |
| **Script 07** | Visualización de los eventos de splicing diferencial mediante *volcano plots* y mapas de calor. |
| **Script 08** | Anotación funcional de los eventos de splicing y análisis de la distribución de los distintos tipos de eventos. |
| **Script 09** | Análisis de enriquecimiento funcional mediante Gene Ontology (GO) y KEGG Pathways. |
| **Script 10** | Integración multiómica entre expresión génica diferencial y splicing alternativo mediante diagramas de Venn e identificación de genes candidatos. |
| **Script 11** | Evaluación de la asociación entre los genes candidatos y el estado tumoral mediante análisis estadístico comparativo. |
| **Script 12** | Análisis de supervivencia global mediante curvas de Kaplan–Meier para los genes candidatos identificados. |

## Requisitos

Los scripts fueron desarrollados utilizando **R** y diversas librerías de **CRAN** y **Bioconductor**, entre las que se incluyen:

- BiocManager
- data.table
- ggplot2
- ggrepel
- reshape2
- edgeR
- limma
- biomaRt
- survival
- survminer
- pheatmap
- enrichR
- VennDiagram
- gtsummary
- gt
- grid

El orden de ejecución de los scripts corresponde al flujo de trabajo descrito en el apartado de Material y Métodos de la memoria.
