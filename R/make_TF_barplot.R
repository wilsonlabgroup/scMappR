#' Plot gprofileR Barplot (TF)
#' 
#' Make a barplot of the top transcription factors enriched by gprofileR.
#'
#' This function takes a gprofileR output and prints the top "top_tfs" most significantly
#' enriched p-values before plotting the rank of their p-values.
#'
#'
#' @rdname make_TF_barplot
#' @name make_TF_barplot
#'
#' @param ordered_back_all_tf Output of the gprofileR function.
#' @param top_tf The number of pathways you want to plot.
#' 
#' @return \code{make_TF_barplot} A barplot of the number of "top_tf" tf names (not motifs), ranked by -log10(Pfdr). \cr
#' 
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_text theme coord_flip labs element_text geom_bar theme_classic xlab ylab scale_fill_manual element_line
#' @importFrom pheatmap pheatmap
#' @importFrom graphics barplot plot
#' @importFrom Seurat AverageExpression CreateSeuratObject PercentageFeatureSet SCTransform SelectIntegrationFeatures PrepSCTIntegration FindIntegrationAnchors IntegrateData DefaultAssay RunPCA RunUMAP FindNeighbors FindClusters ScaleData FindMarkers
#' @importFrom GSVA gsva
#' @importFrom stats fisher.test median p.adjust reorder t.test sd var complete.cases ks.test dist shapiro.test mad
#' @importFrom utils combn read.table write.table head tail
#' @importFrom downloader download
#' @importFrom grDevices pdf dev.off colorRampPalette
#' @importFrom gprofiler2 gost
#' @importFrom gProfileR gprofiler
#' @importFrom pcaMethods prep pca R2cum
#' @importFrom limSolve lsei
#' @importFrom pbapply pblapply
#' @importFrom ADAPTS estCellPercent
#' @importFrom reshape melt
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
#' ordered_query = TRUE, significant = TRUE, exclude_iea = FALSE, multi_query = FALSE, 
#' measure_underrepresentation = FALSE, evcodes = FALSE, user_threshold = 0.05, 
#' correction_method = "fdr",   numeric_ns = "", sources = c("GO:BP", "KEGG", "REAC"))  
#' ordered_back_all <- ordered_back_all$result
#' ordered_back_all <- ordered_back_all[ordered_back_all$term_size > 15 &
#'  ordered_back_all$term_size < 2000 & ordered_back_all$intersection_size > 2,]
#' ordered_back_all_tf <- gprofiler2::gost(query = rowname$rowname[1:150], organism = "mmusculus",
#'  ordered_query = TRUE, significant = TRUE, exclude_iea = FALSE, multi_query = FALSE,
#'  measure_underrepresentation = FALSE, evcodes = FALSE, user_threshold = 0.05,
#'   correction_method = "fdr",  numeric_ns = "", sources = c("TF"))  
#' ordered_back_all_tf <- ordered_back_all_tf$result
#' ordered_back_all_tf <- ordered_back_all_tf[ordered_back_all_tf$term_size > 15
#'  & ordered_back_all_tf$term_size < 5000 & ordered_back_all_tf$intersection_size > 2,]
#' TF = ordered_back_all_tf
#' BP <- ordered_back_all
#' bp <- plotBP(BP)
#' tf <- make_TF_barplot(TF)
#'  }
#' @export
#' 
make_TF_barplot <- function(ordered_back_all_tf, top_tf = 5) {
  # This function takes the TF output and prints the top "top_tf" most enriched transcription factors
  # Args:
  # ordered_back_all_tf = transcription factor enrichment from gprofileR
  # top_tf = the number of TFs being plotted
  # Returns:
  # The top "top_5" TF names, ordered by -log10(Pfdr)
  
  ordered_back_all_tf_class <- class(ordered_back_all_tf)[1] %in% c("data.frame", "matrix")
  if(ordered_back_all_tf_class[1] == FALSE) {
    stop("ordered_back_all_tf must be of class data.frame or matrix")
  }
  if(is.matrix(ordered_back_all_tf)) {
    warning("converting ordered_back_all_tf matrix to dataframe") 
    ordered_back_all_tf <- as.data.frame(ordered_back_all_tf) 
    
  }
  term_name_p_val <- ("term_name" %in% colnames(ordered_back_all_tf))  & ("p_value" %in% colnames(ordered_back_all_tf))
  if(term_name_p_val[1] == FALSE ) {
    stop("ordered_back_all_tf must contain two columns, term_name and p_value")
  }
  
  if(!is.numeric(top_tf)) {
    stop("top_tf must be of class numeric.")
  }
  
  if(nrow(ordered_back_all_tf) == 0) {
    g <- ggplot2::ggplot() + ggplot2::geom_bar(stat = "identity", fill = "mediumpurple") + ggplot2::coord_flip() +  ggplot2::labs(y = "-log10(Padj)", x = "TF Motif") 
    y <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(face=NULL, color="black", 
                                                                size=12, angle=35),
                            axis.text.y = ggplot2::element_text(face=NULL, color="black", 
                                                                size=12, angle=35), 
                            axis.title=ggplot2::element_text(size=16, color = "black"))
    print(y)
    return(y)
    
  }
  ordered_back_all_tf$p_value <- toNum(ordered_back_all_tf$p_value)
  ordered_back_all_tf$term_name <- tochr(ordered_back_all_tf$term_name)
  take1 <- function(x) return(x[1]) # take the first element of a list
  sp <- strsplit(tochr(ordered_back_all_tf$term_name), ";") # split the ane of the TF output
  tfs <- unlist(lapply(sp, take1))
  tfs <- gsub("Factor:","",gsub("-","", tochr(tfs))) # remove extra text
  ordered_back_all_tf$tf <- tochr(tfs)
  nodup <- ordered_back_all_tf[!duplicated(tochr(tfs)),] # keep the most signficant TF motif
  
  ndup_1_10 <- nodup[order(nodup$p_value),]
  if(nrow(ndup_1_10) > top_tf) { # take the top TF numberof factors
    ndup_1_10 <- ndup_1_10[1:top_tf,]
  }
  ndup_1_10$log10 <- -1*toNum(log10(ndup_1_10$p_value)) # make ranks
  # ggplot barplot
  log10 <- ndup_1_10$log10
  tf <- ndup_1_10$tf
  g <- ggplot2::ggplot(ndup_1_10, ggplot2::aes(x = stats::reorder(tf, log10), y = log10)) + ggplot2::geom_bar(stat = "identity", fill = "mediumpurple") + ggplot2::coord_flip() +  ggplot2::labs(y = "-log10(Padj)", x = "TF Motif") 
  y <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(face=NULL, color="black", 
                                                              size=12, angle=0),
                          axis.text.y = ggplot2::element_text(face=NULL, color="black", 
                                                              size=12, angle=0), 
                          axis.title=ggplot2::element_text(size=16, color = "black"))
  y <- y + ggplot2::theme_classic()
  print(y)
  return(y)
}