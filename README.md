
# single-cell mapper (scMappR)

### Dustin Sokolowski: dustin-dot-sokolowski-at-sickkids-dot-ca

### Date: 07/01/2020


## Description

*Description of scMappR R package adapted from pre-print*. The single cell mapper (scMappR) R package contains a suite of bioinformatic tools that provide experimentally relevant cell-type specific information to a list of differentially expressed genes (DEG). The function `scMappR_and_pathway_analysis` reranks DEGs to generate cell-type specificity scores called cell-weighted fold-changes. Users input a list of DEGs, normalized counts, and a signature matrix into this function. scMappR then re-weights bulk DEGs by cell-type specific expression from the signature matrix, cell-type proportions from RNA-seq deconvolution and the ratio of cell-type proportions between the two conditions to account for changes in cell-type proportion. With cwFold-changes calculated, scMappR uses two approaches to utilize cwFold-changes to complete cell-type specific pathway analysis. The `process_dgTMatrix_lists` function in the scMappR package contains an automated scRNA-seq processing pipeline where users input scRNA-seq count data, which is made compatible for scMappR and other R packages that analyze scRNA-seq data. We further used this to store hundreds up regularly updating signature matrices. The functions `tissue_by_celltype_enrichment`, `tissue_scMappR_internal`, and `tissue_scMappR_custom` combine these consistently processed scRNAseq count data with gene-set enrichment tools to allow for cell-type marker enrichment of a generic gene list (e.g. GWAS hits). Reference: Sokolowski,D.J., Faykoo-Martinez,M., Erdman,L., Hou,H., Chan,C., Zhu,H., Holmes,M.M., Goldenberg,A. and Wilson,M.D. (2020) Single-cell mapper (scMappR): using scRNA-seq to infer cell-type specificities of differentially expressed genes. BioRxiv, 10.1101/2020.08.24.265298.

## Installation

Currently, there is only  a development version. scMappR relies on the following dependencies which should be downlaoded/updated with scMappR automatically. Please ensure that these packages are not open when installing scMappR. 

  * ggplot2 - CRAN
  * pheatmap - CRAN
  * graphics - CRAN
  * Seurat - CRAN
  * GSVA - Bioconductor
  * stats - CRAN
  * utils - CRAN
  * downloader - CRAN
  * pcaMethods - Bioconductor
  * grDevices - CRAN
  * gProfileR - CRAN
  * limSolve - CRAN 
  
Install GSVA and pcaMethods from bioconductor first, as `devtools::install_github()` will automatically install CRAN dependencies. 

1. Github (Development Version)


```{r install_developter, eval=FALSE}

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

if (!requireNamespace("pcaMethods", quietly = TRUE))
BiocManager::install("pcaMethods")

if (!requireNamespace("GSVA", quietly = TRUE))
BiocManager::install("GSVA")

devtools::install_github("wilsonlabgroup/scMappR")



```


2. CRAN (Stable Release)


```{r install_cran, eval=FALSE}

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

if (!requireNamespace("pcaMethods", quietly = TRUE))
BiocManager::install("pcaMethods")

if (!requireNamespace("GSVA", quietly = TRUE))
BiocManager::install("GSVA")

install.packages("scMappR")

```

## Data Download

Link to data used in scMappR: https://github.com/wilsonlabgroup/scMappR_Data

To run scMappR locally, please download all .rda files in this data download repository. In many of the functions, the "rda_path" argument can be changed to wherever you would like to download these files to. It however assumes "~/scMappR/data".

If these functions do not detect these rda files, they will temporarily download them with the downloader R package; however, these rda files must already be downloaded to use the examples that are not automatically run.

### Download signature matrices from scMappR_Data repo.

```{r tissue_scMappR_internal, eval=FALSE }

signature_matrix_download <- get_signature_matrices("all") # get all matrices and celltype labels


```


## Primary Functionalities of scMappR.

Below describes the primary ways scMappR can contextualize gene lists and process data. It is strongly reccomended to set `toSave = TRUE` in functions and, when appropraite, `internet = TRUE`. Otherwise scMappR will not print files/directories and many of the results will not be printed.

* `tissue_scMappR_custom()`: This function visualizes signature matrix, clusters subsetted genes, completes enrichment of individual cell-types and co-enrichment. 
* `tissue_scMappR_internal()`: This function loops through every signature matrix in a particular tissue and generates heatmaps, cell-type preferences, and co-enrichment.
* `tissue_by_celltype_enrichment()`: This function completes a fishers exact test of an input list of genes against one of the two curated tissue by cell-type marker datasets from scMappR.
* `scMappR_and_pathway_analysis()`: This function generates scMappR Transformed Variables (STVs), visualies them in a heatmap, and completes pathway enrichment of STVs and bulk gene list.
* `process_dgTMatrix_lists()`: This function takes a list of count matrices, processes them, calls cell-types, and genreates signature matrices.

### Cell-type markers in a list of genes.

#### Internal signature matrix.

Input a list of human or mouse gene symbols as well as a tissue of interest to intentify and visualize enriched cell-types and cell-types co-enriched by the same genes. This process is repeated for every signature matrix (i.e. scRNA-seq study) present for the inputted tissue.

```{r tissue_scMappR_internal, eval=FALSE }


data(Preoptic_Area)
Signature <- POA_example$POA_Rank_signature
rowname <- get_gene_symbol(Signature) 
rownames(Signature) <- rowname$rowname
genes <- rownames(Signature)[1:60]
 rda_path1 = "~/Documents/scMappR/data" 
internal <- tissue_scMappR_internal(genes,"mouse",output_directory = "scMappR_Test",
tissue = "hypothalamus",rda_path = rda_path1, toSave = TRUE)

```

#### Provided signature matrix.

Complete the same process but with a signature matrix and gene list provided by the user. Here, gene symbols do not have to be human or mouse, the symbols in the list must match the signature.

```{r tissue_scMappR_custom, eval=FALSE}

data(POA_example)
Signature <- POA_example$POA_Rank_signature
 rowname <- get_gene_symbol(Signature)
 rownames(Signature) <- rowname$rowname
 genes <- rownames(Signature)[1:200]

 internal <- tissue_scMappR_custom(genes,Signature,output_directory = "scMappR_Test_custom", toSave = TRUE)

```

#### Non tissue-specific

While computing custom signature matrices, there are cell-type markers across tissues and studies. These markers are stored in a gmt file format and can be used for cell-type enrichment. This may be useful if there is not a particular tissue in mind. This function inputs a list of human or mouse gene symbols.

```{r tissue_by_celltype_enrichment, eval=FALSE }


data(POA_example)
POA_generes <- POA_example$POA_generes
POA_OR_signature <- POA_example$POA_OR_signature
POA_Rank_signature <- POA_example$POA_Rank_signature
Signature <- POA_Rank_signature
rowname <- get_gene_symbol(Signature)
rownames(Signature) <- rowname$rowname
genes <- rownames(Signature)[1:100]

enriched <- tissue_by_celltype_enrichment(gene_list = genes, species = "mouse",p_thresh = 0.05, isect_size = 3)


```



### Scaling and visualizing Differentially Expressed Genes from bulk RNA-seq data 

This function requires a normalized count matrix from RNA-seq, a signature matrix such that there are fewer cell-types than there are samples, and a list of differentially expressed genes (with log2Fold-Change and adjusted P-value). 

Then, it will calculate scMappR Transformed Values (STVs) before estimating if there are changes in cell-type proportion between samples. If `toSave = TRUE` then STV's will be visualized. Additionally, if `internet = TRUE`, scMappR will iteratively re-order DEGs based on their STV's and complete pathway analysis with g:ProfileR or gprofiler2.

```{r scMappR_and_pathway_analysis, eval=FALSE}

data(PBMC_scMappR) 
bulk_DE_cors <- PBMC_example$bulk_DE_cors 
bulk_normalized <- PBMC_example$bulk_normalized 
odds_ratio_in <- PBMC_example$odds_ratio_in 
case_grep <- "_female" 
control_grep <- "_male" 
max_proportion_change <- 10 
theSpecies <- "human" 

toOut <- scMappR_and_pathway_analysis(bulk_normalized, odds_ratio_in, 
                                      bulk_DE_cors, case_grep = case_grep,
                                      control_grep = control_grep, rda_path = "", 
                                      max_proportion_change = 10, print_plots = TRUE, 
                                      plot_names = "scMappR_vignette_", theSpecies = "human", 
                                      output_directory = "scMappR_vignette_",
                                      sig_matrix_size = 3000, up_and_downregulated = TRUE, 
                                      internet = TRUE, toSave = TRUE)

```

### Generating a signature matrix and processing scRNA-seq count data

A matrix, dgTMatrix, or list of these matrices are inputted where the rows are genes and the columns are indiviudal cells. The gene names must be human or mouse gene symbols. If the dataset is being processed from PanglaoDB, then the rownames are GeneSymbol-ENSMBL, here, set `species_name = -9`. 

This function returns the signature matrix and cell-type labels. If `toSave = TRUE` signature matrices, all cell-type markers, the average expression of genes from each cell-type, and cell-type labels from gsva are stored as files in the working directory. Additionally, if `saveSCObject = TRUE`, then the Seurat object is also saved in the working directory.

```{r process_dgTMatrix_lists, eval=FALSE}

data(sm)
toProcess <- list(example = sm)
tst1 <- process_dgTMatrix_lists(toProcess, name = "testPropcess", species_name = -9, naming_preference = "eye",
rda_path = "~/scMappR/data", panglao_set = "TRUE", toSave = TRUE, saveSCObject = TRUE)


```

