#' To Character.
#'
#' This function checks if your vector is not a character and if not, will  convert it to a character.
#'
#'
#' @rdname tochr
#' @name tochr
#'
#' @param x A character, factor or numeric vector.
#'
#' @return \code{tochr} Returns a character vector. \cr
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
#'  fact <- factor(c("a", "b", "c", "d"))
#'  # convert to character
#'  char <- tochr(fact)
#'  
#' @export

tochr <- function(x) {
  # this function checks if your vector is not a character and will then convert it to a character
  
  # Args.
	# x = vector that is character, factor, or numeric
  # Returns
	# Vector as a character 
  if(!(class(x)[1] %in% c("character", "factor", "numeric","integer"))) {
    stop("x must be in class character, factor, or numeric")
  }
  
  if(is.character(x)) return(x)
  if(is.factor(x)) return(as.character(levels(x))[x])
  if(is.numeric(x)) return(as.character(x))
  if(class(x)[1] == "integer") return(as.character(x))
} 

