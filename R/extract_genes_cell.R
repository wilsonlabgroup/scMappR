#' Extract Markers
#'
#' Extracting cell-type markers from a signature matrix.
#'
#' This function takes a signature matrix and 
#' extracts cell-type markers above a p-value or fold-change threshold.
#'
#'
#' @rdname extract_genes_cell
#' @name extract_genes_cell
#'
#' @param geneHeat The heatmap of ranks from your scRNA-seq dataset with your genes subsetted
#' @param cellTypes The cell-types that you're interested in extracting. They need to be colnames (not-case sensitive)
#' @param val How associated a gene is with a particualr cell type to include in your list - default is slightly
#' @param isMax If you are taking the single best CT marker (T/F) -- TRUE not reccomended
#' @param isPvalue If the signature matrix is raw p-value (T/F) -- TRUE not reccomended 
#' 
#'
#' @return \code{extract_genes_cell} A list of genes above the threshold for each sample. \cr
#'
#' @import matrixStats
#' @import DeconRNASeq
#' @import S4Vectors
#' @import ggplot2
#' @import gplots
#' @import graphics
#' @import Seurat
#' @import GSVA
#' @import stats
#' @import utils
#'
#' @examples
#' 
#' # load in signature matrices
#' load("data/Preoptic_region_example.rda")
#' # data(Preoptic_region_example)
#' Signature <- POA_Rank_signature
#'  RowName <- get_gene_symbol(Signature)
#'  rownames(Signature) <-RowName$rowname
#'  # extract genes with a -log10(Padj > 1)
#'  Signat <- extract_genes_cell(Signature)
#'  
#' @export

extract_genes_cell <- function(geneHeat, cellTypes = "ALL", val = 1, isMax = F, isPvalue = FALSE) {

# This function takes a signature matrix and extracts cell-type markers above a p-value or fold-change threshold. 

  #Args:   
    # geneHeat: the heatmap of ranks from your scRNA-seq dataset with your genes subsetted
    # cellTypes: The cell-types that you're interested in extracting. They need to be colnames (not-case sensitive)
    # val: how associated a gene is with a particualr cell type to include in your list - default is slightly above 1
    # isMax true or false, if you want to sort genes into the CT marker that represents theme best
    # isPvalue: if the signature matrix is simply a raw P-value: not reccomended
  #Returns: a list of genes above the threshold for each sample.

colnames(geneHeat) <- toupper(colnames(geneHeat))
cellTypes <- toupper(cellTypes)

geneHeat <- geneHeat[rowSums(geneHeat) > 1,] # extract genes with any CT specificity
if(class(geneHeat) == "numeric" | class(geneHeat) == "character" ) {
  geneHeat <- data.frame(t(geneHeat))
}  

genes_extracted <- list()
if(cellTypes == "ALL") { # if you extract all cell-type just make it your output
  cellTypes <- colnames(geneHeat)
}

if(isMax == TRUE) {  #take the single top gene
  # this way, you take the cell type that marks each gene the best.
  # no gene is double counted into multiple cell types.
  warning("Extracting the single most expressed marker in each cell-type. Not reccomeneded.")
  cn <- colnames(geneHeat)[max.col(geneHeat,"first")] # the cell type with the most signficant
  names(cn) <- rownames(geneHeat)
  for(i in 1:length(colnames(geneHeat))) {
    genes_extracted[[i]] <- names(cn)[which(cn == colnames(geneHeat)[i])]
    names(genes_extracted)[i] <- colnames(geneHeat)[i]
  }
  return(genes_extracted) 
} else { # take top genes past a threshold
  # a gene can be counted into multiple cell types and you're just saying that a 
  # gene is in it if it's p-value passes a certain threshold -- This is probably
  # what fits the premise of scMappR slightly better
  if(isPvalue == TRUE) {
    val = -1*log10(val)
  } else {
    val = val
  }
  for(i in 1:length(colnames(geneHeat))) {
    genes_extracted[[i]] <- rownames(geneHeat)[geneHeat[,i] > val]
    names(genes_extracted)[i] <- names(genes_extracted)[i] <- colnames(geneHeat)[i]
    
  }
  return(genes_extracted)
}
}
