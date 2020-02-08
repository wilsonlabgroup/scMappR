#' Extract Top Markers
#'
#' Internal -- Extracts strongest cell-type markers from a Seurat object.
#'
#' Internal, this function runs through a list of outputs from FindMarkers objects in Seurat
#' and will extract genes past a padj and fold-change threshold. Then it extracts the topNum number of genes.
#' if you have not used the FindMarkers function, then a list of summary statistics with 
#' fold-change designated by avg_logFC and p-val by p_val_adj.
#'
#' @rdname topgenes_extract
#' @name topgenes_extract
#'
#' @param generes A list of cell-tpe markers with fold-changes and p-vlaues (FindMarkers output in Seurat).
#' @param padj The p-value (FDR) cutoff.
#' @param FC The fold-change cutoff.
#' @param topNum The number of genes to extract.
#'
#' @return \code{topgenes_extract} Returns a list of character vectors with the top (topNum) of gene markers for each cell-type. \cr
#' 
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_text theme coord_flip labs element_text
#' @importFrom gplots heatmap.2
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
#'
#' @examples 
#'  
#'  \donttest{
#'  # load generes object
#'  data(POA_example)
#'  topGenes <- topgenes_extract(POA_example$POA_generes)
#' }
#' 
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
  
  if(class(generes) != "list") {
    stop("generes must be of class list.")
  }
  if(all(is.numeric(padj), is.numeric(FC), is.numeric(topNum))[1] == FALSE) {
    stop("padj, FC, and topNum must all be of class numeric.")
  }
  
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

