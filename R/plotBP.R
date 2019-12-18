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
#' \donttest{
#' data(POA_example)
#'  POA_generes <- POA_example$POA_generes
#'  POA_OR_signature <- POA_example$POA_OR_signature
#'  POA_Rank_signature <- POA_example$POA_Rank_signature
#' Signature <- as.data.frame(POA_Rank_signature)
#' rowname <- get_gene_symbol(Signature)
#' rownames(Signature) <- rowname$rowname
#' ordered_back_all <- gprofiler2::gost(query = rowname$rowname[1:100], organism = "mmusculus",
#'  ordered_query = TRUE, significant = TRUE, exclude_iea = FALSE, multi_query = FALSE,
#'   measure_underrepresentation = FALSE, evcodes = FALSE, user_threshold = 0.05,
#'    correction_method = "fdr",   numeric_ns = "", sources = c("GO:BP", "KEGG", "REAC"))  
#' ordered_back_all <- ordered_back_all$result
#' ordered_back_all <- ordered_back_all[ordered_back_all$term_size > 15
#'  & ordered_back_all$term_size < 2000 & ordered_back_all$intersection_size > 2,]
#' ordered_back_all_tf <- gprofiler2::gost(query = rowname$rowname[1:150], organism = "mmusculus",
#'  ordered_query = TRUE, significant = TRUE, exclude_iea = FALSE, multi_query = FALSE,
#'   measure_underrepresentation = FALSE, evcodes = FALSE, user_threshold = 0.05,
#'    correction_method = "fdr",  numeric_ns = "", sources = c("TF"))  
#' ordered_back_all_tf <- ordered_back_all_tf$result
#' ordered_back_all_tf <- ordered_back_all_tf[ordered_back_all_tf$term_size > 15 
#' & ordered_back_all_tf$term_size < 5000 & ordered_back_all_tf$intersection_size > 2,]
#' TF = ordered_back_all_tf
#' BP <- ordered_back_all
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
  
  if(class(ordered_back_all) != "data.frame" & class(ordered_back_all) != "matrix") {
    stop("ordered_back_all_tf must be of class data.frame or matrix")
  }
  if(class(ordered_back_all) == "matrix") {
    warning("converting ordered_back_all_tf matrix to dataframe") 
    ordered_back_all <- as.data.frame(ordered_back_all) 
    
  }
  if(!("term_name" %in% colnames(ordered_back_all)) | !("p_value" %in% colnames(ordered_back_all)) ) {
    stop("ordered_back_all must contain two columns, term_name and p_vale")
  }
  if(class(top_bp) != "numeric") {
    stop("top_bp must be of class numeric.")
  }
  
  
  if(nrow(ordered_back_all) == 0) {
    g <- ggplot2::ggplot() + ggplot2::geom_bar(stat = "identity", fill = "mediumpurple") + ggplot2::coord_flip() +  ggplot2::labs(y = "-log10(Padj)", x = "Gene Ontology") 
    y <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(face=NULL, color="black", 
                                                                size=12, angle=35),
                            axis.text.y = ggplot2::element_text(face=NULL, color="black", 
                                                                size=12, angle=35), 
                            axis.title=ggplot2::element_text(size=16, color = "black"))
    print(y)
    return(y)
    
  }
  ordered_back_all$p_value <- toNum(ordered_back_all$p_value)
  ordered_back_all$term_name <- tolower(tochr(ordered_back_all$term_name))
  ordered_back_all$log10 <- -1*log10(ordered_back_all$p_value) # ranks of g:profileR
  ordered_back_all <- ordered_back_all[order(ordered_back_all$log10, decreasing = TRUE),]
  if(nrow(ordered_back_all) > top_bp) {
    ordered_back_all <- ordered_back_all[1:top_bp,]
  }
  
  term_name <- ordered_back_all$term_name
  log10 <- ordered_back_all$log10
  # Plot the barplot and set the size of the text of each pathway to fit 
  g <- ggplot2::ggplot(ordered_back_all, ggplot2::aes(x = stats::reorder(term_name, log10), y = log10)) + ggplot2::geom_bar(stat = "identity", fill = "turquoise") + ggplot2::coord_flip() +  ggplot2::labs(y = "-log10(Padj)", x = "Gene Ontology") 
  y <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(face=NULL, color="black", 
                                            size=12, angle=35),
                 axis.text.y = ggplot2::element_text(face=NULL, color="black", 
                                            size=12, angle=35), 
                 axis.title= ggplot2::element_text(size=16, color = "black"))
  print(y)
  return(y)
}
