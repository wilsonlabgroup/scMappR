#' Preoptic_Area
#' 
#' Data for tissue_scMappR_custom, tissue_scMappR_internal, generes_to_heatmap.
#'
#' A named list called POA_example (pre-optic area example) containing three objects, POA_generes: a list of truncated dataframes containing summary statistics for each cell-type marker, POA_OR_signature a truncated signature matrix of odds ratio's for cell-types in the POA, and POA_Rank_signature a truncated signature matrix of -log10(Padj) for cell-type markers in the POA.
#'
#' @rdname POA_example
#' @name POA_example
#' 
#' @usage data(POA_example)
#'
#' @format A list containing three objects: summary statistics of cell-type markers, a signature matrix of odds ratios, and a signature matrix of ranks.
#' \describe{
#'   \item{POA_generes}{ A list of 27 data frames containing (up to 30) cell-type markers. Each element of the list is a dataframe where rows are genes, and columns are p-value, log2FC, percentage of cells expressing gene in cell-type, percentage of cells expressing gene in other cell-types, and FDR adjusted p-value. }
#'   \item{POA_OR_signature}{ A 266 x 27 matrix where rows are genes, columns are cell-types and matrix is filled with the odds-ratio that a gene is in each cell-type.}
#'   \item{POA_Rank_signature}{A 266 x 27 matrix of  matrix where rows are genes, columns are cell-types and matrix is filled with the rank := -log10(P_fdr) that a gene is in each cell-type. }
#' }
#' @examples 
#' data(POA_example)
#' @keywords datasets
#' 
NULL
