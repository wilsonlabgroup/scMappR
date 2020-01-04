# Single Cell MappeR (scMappR)

Dustin Sokolowski

## Description

  We built an R package called scMappR, which contains statistical tools to interrogate the cell-type specificity of any gene list given a matrix of cell-types and genes associated with those cell-types (a signature matrix).  We further processed re-aligned scRNA-seq data from panglaoDB (cite) and provided signature matrices across hundreds of tissues from mouse and human samples (described below). This R package has a number of uses including the (i) processing scRNA-seq count data and automated cell-type naming using Seurat V3 and enrichment of CellMarker and Panglao databases (process_dgTMatrix_lists()), (ii) tissue-by-cell-type gene set enrichment (cellmarker_enrich()), (iii) cell-type specific enrichment of a gene list within a particular tissue (tissue_scMappR_custom/tissue_scMappR_internal), (iv) weighted cell-type specific reranking of a list of differentially expressed genes (scMappR_and_pathway_analysis()).

## Installation
1. Github (Development Version)



devtools::install_github("DustinSokolowski/scMappR")

library(scMappR)

2. CRAN (Stable Release)

## Data Download
Link to data used in scMappR: https://github.com/DustinSokolowski/scMappR_Data

Currently, to run scMappR locally, please download all .rda files in this data download repository. In many of the functions, the "rda_path" argument can be changed to wherever you would like to download these files to. It however assumes "~/scMappR/data".

Internally, if these functions do not detect these rda files, they will temporarily download them with the downloader R package; however, these rda files must already be downloaded to use the examples that are not automatically run.

## Reference
< Insert publication here > 
This is another test line
