#' single_cell_process
#' 
#' Example data for processing scRNA-seq count data with Seurat.
#'
#' A dgCMatrix object containing count data for scRNA-seq processing.
#'
#' @usage data(single_cell_process)
#'
#' @format A 752 x 236 matrix of class dgCMatrix where rows are genes and columns are cells. Data matrix is filled with counts detected from scRNAseq
#' \describe{
#'   \item{TCTCTAACACAGGCCT}{ Barcode of one of the sequenced cells present. Each column is the count from a scRNA-seq dataset reprocessed by PanglaoDB.}
#' }
#' @examples 
#' data(single_cell_process)
#' @keywords datasets
#' 
NULL