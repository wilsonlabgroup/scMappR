#' Internal -- get gene symbol from Panglao.db matrix.
#'
#' Internal -- removes Ensembl signature appended to signature matrix from Panglao and figure out species by pre-fix Ensembl of the Ensembl ID that is appended to gene names.
#'
#' Internal: This function removes the ENGMUS/ENGS tag from Panglao created gene names (symbol-ENGS).
#' From the ENSG/ENSMUS, this function determines if the species is mouse/human and returns the gene symbols.  
#' 
#' @rdname get_gene_symbol
#' @name get_gene_symbol
#'
#' @param wilcoxon_rank_mat_t Matrix where row names are "GeneSymbol-Ensembl" (human or mouse).
#' 
#'  
#'
#' @return List with the following elements:
#' \item{rowname}{Genes in the signature matrix excluding the ensemble name.}
#' \item{species}{"mouse" or "human" depending on appended ensembl symbols.}
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_text theme coord_flip labs element_text
#' @importFrom pheatmap pheatmap
#' @importFrom graphics barplot plot
#' @importFrom Seurat AverageExpression CreateSeuratObject PercentageFeatureSet SCTransform SelectIntegrationFeatures PrepSCTIntegration FindIntegrationAnchors IntegrateData DefaultAssay RunPCA RunUMAP FindNeighbors FindClusters ScaleData FindMarkers
#' @importFrom GSVA gsva
#' @importFrom stats fisher.test median p.adjust reorder t.test sd var complete.cases
#' @importFrom utils combn read.table write.table head tail
#' @importFrom downloader download
#' @importFrom grDevices pdf dev.off colorRampPalette
#' @importFrom gprofiler2 gost
#' @importFrom gProfileR gprofiler
#' @importFrom pcaMethods prep pca R2cum
#' @importFrom limSolve lsei
#' @importFrom pbapply pblapply
#' @importFrom ADAPTS estCellPercent
#'
#' @examples
#' 
#'  # load signature
#' data(POA_example) 
#' POA_OR_signature <- POA_example$POA_OR_signature
#' symbols <- get_gene_symbol(POA_OR_signature)
#' 
#' 
#' @export
#' 
get_gene_symbol <- function(wilcoxon_rank_mat_t) {
  # internal -- removes ensembl signature appended to signature matrix from panglao
  
  # figure out species by type of ensembl is appended to gene names
  # Args: 
  # wilcoxon_rank_mat_t: often a signature matrix but can be any dataframe with rownames such as geneSymbol-Ensembl
  # Returns:
  # A list containing the gene-symbols only as well as if the species is mouse or human
  
  if(is.null(rownames(wilcoxon_rank_mat_t))) {
    stop("wilcoxon_rank_mat_t must be an object with row names. Currently, rownames == NULL")
  }
  
  the_human <- length(grep("ENSG00", rownames(wilcoxon_rank_mat_t)))
  the_mouse <- length(grep("ENSMUSG00", rownames(wilcoxon_rank_mat_t)))

  if((the_mouse == 0 & the_human == 0)[1]) {
    stop("Matrix is not internal and this function should not be required. gene symbols of genes, bulk matrix, and signature matrix should already match.")
  }
  
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

