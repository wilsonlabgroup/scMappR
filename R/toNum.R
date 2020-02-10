#' To Numeric.
#'
#' This function checks if your vector is not a character and if it is, then converts it to a numeric.
#'
#'
#' @rdname toNum
#' @name toNum
#'
#' @param x A character, factor, or numeric vector.
#'
#' @return \code{toNum} Returns a numeric vector. \cr
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
#'  
#'  # vector of factors
#' fact <- factor(c("1", "2", "3", "4"))
#' # convert to numeric
#' num <- toNum(fact)
#'  
#' @export

toNum <- function(x) {
  # this function checks if your vector is not a numeric and will then convert it to a numeric
  # args
  # x = vector that is character, factor, or numeric
  # Returns
  # Vector as a numeric 
  if(!(class(x)[1] %in% c("character", "factor", "numeric"))) {
    stop("x must be in class character, factor, or numeric")
  }
  
  
  if(is.character(x)) return(as.numeric(x))
  
  if(is.factor(x)) return(as.numeric(levels(x))[x])
  if(is.numeric(x)) return(x)
  if(class(x)[1] == "integer") return(as.numeric(x))
} 


