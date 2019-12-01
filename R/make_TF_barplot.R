#' Plot g:ProfileR Barplot (TF)
#' 
#' Make a barplot of the top transcription factors enriched by g:ProfileR
#'
#' This function takes a g:ProfleR output and prints the top "top_tfs" most significantly
#' enriched P-values before plotting the rank of their P-values
#'
#'
#' @rdname make_TF_barplot
#' @name make_TF_barplot
#'
#' @param ordered_back_all output of the g:ProfileR function
#' @param top_tf The number of pathways you want to plot
#' 
#' @return \code{make_TF_barplot} A barplot of the number of "top_tf" tf names (not motifs), ranked by -log10(Pfdr). \cr
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
#' @import grDevices
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
#'                            max_set_size = 2000, exclude_iea = T)
#' TF <- gProfileR::gprofiler(rowname$rowname, "mmusculus", src_filter = c("TF"),
#'                           max_set_size = 5000, exclude_iea = T)
#' bp <- plotBP(BP)
#' tf <- make_TF_barplot(TF)
#'  }
#' @export
#' 
make_TF_barplot <- function(ordered_back_all_tf, top_tf = 5) {
  # This function takes the TF output and prints the top "top_tf" most enriched transcription factors
  # Args:
  # ordered_back_all_tf = transcription factor enrichment from g:ProfileR
  # top_tf = the number of TFs being plotted
  # Returns:
  # The top "top_5" TF names, ordered by -log10(Pfdr)
  
  take1 <- function(x) return(x[1]) # take the first element of a list
  sp <- strsplit(ordered_back_all_tf$term.name, ";") # split the ane of the TF output
  tfs <- unlist(lapply(sp, take1))
  tfs <- gsub("Factor:","",gsub("-","", tochr(tfs))) # remove extra text
  ordered_back_all_tf$tf <- tochr(tfs)
  nodup <- ordered_back_all_tf[!duplicated(tochr(tfs)),] # keep the most signficant TF motif
  ndup_1_10 <- nodup[order(nodup$p.value),]
  if(nrow(ndup_1_10) > top_tf) { # take the top TF numberof factors
    ndup_1_10 <- ndup_1_10[1:top_tf,]
  }
  ndup_1_10$log10 <- -1*log10(ndup_1_10$p.value) # make ranks
  # ggplot barplot
  log10 <- ndup_1_10$log10
  tf <- ndup_1_10$tf
  g <- ggplot2::ggplot(ndup_1_10, ggplot2::aes(x = stats::reorder(tf, log10), y = log10)) + ggplot2::geom_bar(stat = "identity", fill = "mediumpurple") + ggplot2::coord_flip() +  ggplot2::labs(y = "-log10(Padj)", x = "TF Motif") 
  y <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(face=NULL, color="black", 
                                            size=12, angle=35),
                 axis.text.y = ggplot2::element_text(face=NULL, color="black", 
                                            size=12, angle=35), 
                 axis.title=ggplot2::element_text(size=16, color = "black"))
  print(y)
  return(y)
}
