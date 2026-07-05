# Código fuente

Esta carpeta contiene el código fuente completo desarrollado durante el análisis bioinformático del Trabajo Fin de Máster.

El pipeline se implementó íntegramente en el lenguaje de programación **R** y se estructuró en **12 scripts independientes**, organizados de forma secuencial para facilitar la reproducibilidad del flujo de trabajo.

## Organización de los scripts

| Script | Descripción |
|---------|-------------|
| Script 01 | Filtrado inicial de la cohorte y unificación de matrices. |
| Script 02 | Caracterización clínica y demográfica de la cohorte. |
| Script 03 | Traducción de identificadores Ensembl y control de calidad de la matriz de expresión génica. |
| Script 04 | Normalización mediante Limma-Voom y análisis de expresión génica diferencial (DGE). |
| Script 05 | Visualización de los resultados de expresión génica diferencial. |
| Script 06 | Análisis de splicing alternativo diferencial (DS). |
| Script 07 | Visualización de los eventos de splicing alternativo. |
| Script 08 | Anotación funcional y clasificación de los tipos de eventos de splicing. |
| Script 09 | Análisis de enriquecimiento funcional (GO y KEGG). |
| Script 10 | Integración multiómica de expresión génica y splicing alternativo. |
| Script 11 | Análisis de supervivencia mediante curvas de Kaplan-Meier. |
| Script 12 | Generación de las figuras y tablas finales del estudio. |

## Requisitos

Los scripts fueron desarrollados utilizando **R** y diversas librerías del proyecto **Bioconductor**, entre ellas:

- edgeR
- limma
- biomaRt
- data.table
- ggplot2
- pheatmap
- enrichR
- openxlsx

El orden de ejecución de los scripts corresponde al flujo de trabajo descrito en el apartado de Material y Métodos de la memoria.
