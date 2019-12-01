#' Internal -- get gene symbol from Panglao.db matrix
#'
#' internal -- removes ensembl signature appended to signature matrix from panglao and figure out species by type of ensembl is appended to gene names.
#'
#' Internal: This function runs the FindMarkers function from seurat in a loop, will use the seurat v2 or seurat v3 object after identifying which seurat object is inputted. 
#' It then takes the output of the FindMarkers and puts it in a list, returning it.
#' 
#' @rdname get_gene_symbol
#' @name get_gene_symbol
#'
#' @param wilcoxon_rank_mat_t Matrix where row names are "GeneSymbol-Ensembl" (human or mouse)
#' 
#'  
#'
#' @return \code{seurat_to_generes} A list containing the gene-symbols only as well as if the species is mouse or human. \cr
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
#'  # load signature
#' data(Preoptic_Area) 
#' POA_OR_signature <- POA_example$POA_OR_signature
#' symbols <- get_gene_symbol(POA_OR_signature)
#' 
NULL
#' @rdname get_gene_symbol
#' @export


get_gene_symbol <- function(wilcoxon_rank_mat_t) {
  # internal -- removes ensembl signature appended to signature matrix from panglao
  
  # figure out species by type of ensembl is appended to gene names
  # Args: 
  # wilcoxon_rank_mat_t: often a signature matrix but can be any dataframe with rownames such as geneSymbol-Ensembl
  # Returns:
  # A list containing the gene-symbols only as well as if the species is mouse or human
  the_human <- length(grep("ENSG00", rownames(wilcoxon_rank_mat_t)))
  the_mouse <- length(grep("ENSMUSG00", rownames(wilcoxon_rank_mat_t)))
  
  #say which species
  if(the_human >= the_mouse) {
    theSpecies <- "human"
  }
  if(the_mouse > the_human) {
    theSpecies <- "mouse"
  }
  sm <- wilcoxon_rank_mat_t
  
  # make sure that mitochondrial genes are flagged
  RN = rownames(sm)
  # take off numbers after decimal in ensembl names
  RN_1 <- sub('(.*)[.](.*)','\\1',RN)
  
  # remove 16 or 19 characters if human ensembl or mouse ensembl names are appended, respectively
  if(theSpecies == "human") {
    
    RN_2 <- substr(RN_1, 1, nchar(RN_1) - 16)
  }
  if(theSpecies == "mouse" ) {
    
    RN_2 <- substr(RN_1, 1, nchar(RN_1) - 19)
  }
  
  return(list(rowname = RN_2, species = theSpecies))
}

