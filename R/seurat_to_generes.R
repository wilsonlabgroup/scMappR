#' Identify all cell-type markers 
#'
#' Takes processed Seurat matrix and identifies cell-type markers with FindMarkers
#'
#' Internal: This function runs the FindMarkers function from seurat in a loop, will use the seurat v2 or seurat v3 object after identifying which seurat object is inputted. 
#' It then takes the output of the FindMarkers and puts it in a list, returning it.
#' 
#' @rdname seurat_to_generes
#' @name seurat_to_generes
#'
#' @param pbmc Processed seurat object.
#' 
#'  
#'
#' @return \code{seurat_to_generes} A list of genes where their over-representation in the i'th cell-type is computed. Each element contains the gene name, adjusted p-value, and the log2FC of each gene being present in that cell-type. \cr
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
#' load("~/scMappR/data/cell_process_example.rda")
#' #data(cell_process_example)
#' toProcess <- list(example = sm)
#' tst1 <- process_from_count(toProcess, "testProcess")
#' cellnames <- gsva_cellIdentify(tst1, "mouse", "brain", "~/scMappR/data")
#' generes <- seurat_to_generes(tst1)
#' 
#' @export


seurat_to_generes <- function(pbmc){
  # Internal: This function runs the FindMarkers function from seurat in a loop, will use the seurat v2 or seurat v3 object 
  # after identifying which seurat object is inputted. It then takes the output of the FindMarkers and puts it in a list, returning it
  # Args:
  # pbmc -- a processed suerat object
  # Returns:
  # A list of genes where their over-representation in the i'th cell-type is computed. Each element contains the gene name, adjusted p-value, and the log2FC of each gene being present in that cell-type.
  id <- try(pbmc@ident, silent = T)
  # try for seurat v2
  
  generes <- list()
  count <-1
  if(class(id) != "try-error") {
    # if it's V2
    for(i in sort(unique(pbmc@ident))) {
      
      de_genes <- try(Seurat::FindMarkers(pbmc, ident.1 = i, test.use = "wilcox"))
      if(class(de_genes) == "try-error" | class(de_genes) == "NULL" ) {
        
        next
      }
      
      generes[[count]] <- de_genes
      names(generes)[count] <- i
      count <- count+1
      print("DE")
      print(i)
    }
    
    
  }
  for(i in sort(unique(pbmc@active.ident))) {
    # if object is Seurat V3
    de_genes <- try(Seurat::FindMarkers(pbmc, ident.1 = i, test.use = "wilcox"))
    if(class(de_genes) == "try-error" | class(de_genes) == "NULL" ) {
      
      next
    }
    
    generes[[count]] <- de_genes
    names(generes)[count] <- i
    count <- count+1
    print("DE")
    print(i)
  }
  return(generes)
}

