# Single Cell MappeR (scMappR)

Dustin Sokolowski

## Description
We built an R package called scMappR, which contains statistical tools to interrogate the cell-type specificity of any gene list given a matrix of cell-types and genes associated with those cell-types (a signature matrix).  We further processed re-aligned scRNA-seq data from panglaoDB (cite) and provided signature matrices across hundreds of tissues from mouse and human samples (described below). This R package has a number of uses including the (i) processing scRNA-seq count data and automated cell-type naming using Seurat V3 and enrichment of CellMarker and Panglao databases (process_dgTMatrix_lists()), (ii) tissue-by-cell-type gene set enrichment (tissue_cell_fishers()),  (iii) weighted cell-type specific reranking of a list of differentially expressed genes (scMappR_and_pathway_analysis()).

## Installation
1. Github (Development Version)

2. CRAN (Stable Release)

## Data Download
< link to data repository >

## Reference
< Insert publication here > 
