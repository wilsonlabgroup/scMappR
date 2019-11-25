#' Extract Top Markers
#'
#'  Internal -- Extracts strongect cell-type markers from a Seurat object
#'
#' Internal, this function runs through a list of outputs from FindMarkers objects in Seurat
#' at will take genes past a padj and FC threshold. Then it extracts the topNum number of genes
#' if you have not used the FindMarkers function, then a list of summary statistics with 
#' fld change designated by avg_logFC and p-val by p_val_adj
#'
#' @rdname topgenes_extract
#' @name topgenes_extract
#'
#' @param generes A list of cell-tpe markers with fold-changes and P-vlaues (FindMarkers output in Seurat)
#' @param padj The p-value (FDR) cutoff
#' @param FC The fold-change cutoff
#' @param topNum The number of genes to extract
#'
#' @return \code{topgenes_extract} Returns a numeric vector. \cr
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
#' @import downloader
#'
#' @examples
#' 
#'  
#'  # load generes object
#'  load("data/Preoptic_region_example.rda")
#'  data(Preoptic_region_example)
#'  topGenes <- topgenes_extract(POA_generes)
#' 
NULL
#' @rdname topgenes_extract
#' @export


topgenes_extract <- function(generes,  padj = 0.05, FC = 1.5, topNum = 30) {
  # Internal, this function runs through a list of outputs from FindMarkers objects in Seurat
  # at will take genes past a padj and FC threshold. Then it extracts the topNum number of genes
  # if you have not used the FindMarkers function, then a list of summary statistics with 
  # fld change designated by avg_logFC and p-val by p_val_adj
  
  # Args:
  # generes: a list of cell-type markers as the output of the FindMarkers function 
  # padj: the p-value (FDR) cutoff
  # FC: The fold-change
  # topNum: the number of genes to extract
  
  # Returns:
  # a list of gene symbols showing the top most "topNum" cell-markers in each cell-type (i.e. for every element in generes)
  topGenes <- list()
  for(i in 1:length(generes)) {
    #extract all genes that are above the cutoff
    genes <- rownames(generes[[i]])[generes[[i]]$avg_logFC > log2(FC) & generes[[i]]$p_val_adj < padj]
    if(length(genes) > topNum) {
      #if there are more than the cutoff (traditionally 30) number of DEGs then only take those
      genes <- genes[1:topNum]
    }
    topGenes[[i]] <- genes
  }
  names(topGenes) <- names(generes)
  return(topGenes)  
}

