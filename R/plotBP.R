#' Plot g:ProfileR Barplot
#' 
#' Make a barplot of the top biological factors enriched by g:ProfileR
#'
#' This function takes a g:ProfleR output and prints the top "top_bp" most significantly
#' enriched P-values before plotting the rank of their P-values#' 
#'
#'
#' @rdname plotBP
#' @name plotBP
#'
#' @param ordered_back_all output of the g:ProfileR function
#' @param top_bp The number of pathways you want to plot
#' 
#' @return \code{plotBP} A barplot of the number of "top_bp" pathways, ranked by -log10(Pfdr). \cr
#' 
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_text theme coord_flip labs element_text
#' @importFrom gplots heatmap.2
#' @importFrom graphics barplot
#' @importFrom Seurat AverageExpression CreateSeuratObject PercentageFeatureSet SCTransform SelectIntegrationFeatures PrepSCTIntegration FindIntegrationAnchors IntegrateData DefaultAssay RunPCA RunUMAP FindNeighbors FindClusters ScaleData FindMarkers
#' @importFrom GSVA gsva
#' @importFrom stats fisher.test median p.adjust reorder t.test
#' @importFrom utils combn
#' @importFrom downloader download
#' @importFrom grDevices pdf dev.off colorRampPalette
#' @importFrom gProfileR gprofiler
#' @importFrom pcaMethods prep pca R2cum
#' @importFrom limSolve lsei
#'
#' @examples 
#' \donttest{
#' data("Preoptic_Area")
#'  POA_generes <- POA_example$POA_generes
#'  POA_OR_signature <- POA_example$POA_OR_signature
#'  POA_Rank_signature <- POA_example$POA_Rank_signature
#' Signature <- as.data.frame(POA_Rank_signature)
#' rowname <- get_gene_symbol(Signature)
#' rownames(Signature) <- rowname$rowname
#' BP <- gProfileR::gprofiler(rowname$rowname, "mmusculus", src_filter = c("GO:BP", "KEGG", "REAC"), 
#'                            max_set_size = 2000, exclude_iea = TRUE)
#' TF <- gProfileR::gprofiler(rowname$rowname, "mmusculus", src_filter = c("TF"),
#'                           max_set_size = 5000, exclude_iea = TRUE)
#' bp <- plotBP(BP)
#' tf <- make_TF_barplot(TF)
#' 
#'  }
#' @export
#' 
plotBP <- function(ordered_back_all, top_bp = 10) {
  # This function takes a g:ProfleR output and prints the top "top_bp" most significantly
  # enriched P-values before plotting the rank of their P-values
  # Args:
  # ordered_back_all: output of the g:ProfileR function
  # The number of pathways you want to plot
  # Returns:
  # A barplot of the number of "top_bp" pathways, ranked by -log10(Pfdr)
  
  ordered_back_all$term.name <- tolower(ordered_back_all$term.name)
  ordered_back_all$log10 <- -1*log10(ordered_back_all$p.value) # ranks of g:profileR
  if(nrow(ordered_back_all) > top_bp) {
    ordered_back_all <- ordered_back_all[1:top_bp,]
  }
  term.name <- ordered_back_all$term.name
  log10 <- ordered_back_all$log10
  # Plot the barplot and set the size of the text of each pathway to fit 
  g <- ggplot2::ggplot(ordered_back_all, ggplot2::aes(x = stats::reorder(term.name, log10), y = log10)) + ggplot2::geom_bar(stat = "identity", fill = "turquoise") + ggplot2::coord_flip() +  ggplot2::labs(y = "-log10(Padj)", x = "Gene Ontology") 
  y <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(face=NULL, color="black", 
                                            size=12, angle=35),
                 axis.text.y = ggplot2::element_text(face=NULL, color="black", 
                                            size=12, angle=35), 
                 axis.title= ggplot2::element_text(size=16, color = "black"))
  print(y)
  return(y)
}
