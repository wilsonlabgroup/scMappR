#' compare_deconvolution_methods
#'
#' This function calculates cell-type proprotions of an inputted bulk sample using DeconRNA-seq, WGCNA, and DCQ methods. Outputted cell-type proportions are then compared.
#'
#'
#' @rdname compare_deconvolution_methods
#' @name compare_deconvolution_methods
#'
#' @param count_file Normalized (CPM, TPM, RPKM) RNA-seq count matrix where rows are gene symbols and columns are individuals. Either the object tself of the path of a .tsv file.
#' @param signature_matrix Signature matrix (odds ratios) of cell-type specificity of genes. Either the object itself or a pathway to a .RData file containing an object named "wilcoxon_rank_mat_or" - generally internal.
#' @param print_plot print the barplot of estimated cell-type propotions from each method into the R console (logical: TRUE/FALSE)
#' @param order_celltype Specify the order that cell-type are placed on the barplot. NULL = alphabetical, otherwise a character vector of cell-type labels (i.e. column names of the signature matrix).
#'
#'
#' @return List with the following elements:
#' \item{cellWeighted_Foldchange}{data frame of cellweightedFold changes for each gene.}
#' \item{cellType_Proportions}{data frame of cell-type proportions from DeconRNA-seq.}
#' \item{leave_one_out_proportions}{data frame of average cell-type proportions for case and control when gene is removed.}
#' \item{processed_signature_matrix}{signature matrix used in final analysis.}
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
#'  
#' data(PBMC_example)
#' norm_counts <- PBMC_example$bulk_normalized
#' signature <- PBMC_example$odds_ratio_in
#' tst <- compare_deconvolution_methods(norm_counts, signature, FALSE, order_celltype = c("I_mono", "C_mono", "CD8_CM", "CD8_TE", "B_SM", "B_NSM", "B_naive"))
#' }
#'  
#' @export
#' 


compare_deconvolution_methods <- function(count_file, signature_matrix, print_plot = FALSE, order_celltype = NULL) {
  norm_counts_i <- count_file
  wilcox_or_signature <- signature_matrix
  
  # Identify clel-type proporitons with all genes
  Decon <- DeconRNAseq_CRAN(as.data.frame(norm_counts_i),as.data.frame(wilcox_or_signature))$out.all
  WGCNA <- t(ADAPTS::estCellPercent(wilcox_or_signature,norm_counts_i, method = "proportionsInAdmixture") / 100)
  WGCNA <- WGCNA[,colnames(WGCNA) != "others"]
  DCQ <- t(ADAPTS::estCellPercent(wilcox_or_signature,norm_counts_i, method = "DCQ") / 100)
  DCQ <-DCQ[,colnames(DCQ) != "others"]
  
  # store them
  proportions <- list(DeconRNAseq = Decon, WGCNA = WGCNA, DCQ = DCQ)
  
  # get average proportions per cell-type
  Decon_Mean <- colMeans(Decon)
  WGCNA_Mean <- colMeans(WGCNA)
  DCQ_Mean <- colMeans(DCQ)
  
  # store them 
  l <- list(DeconRNAseq_avgExpression = Decon_Mean, WGCNA_avgExpression = WGCNA_Mean, DCQ_avgExpression = DCQ_Mean)
  l1 <- l
  names(l1) <- c("DeconRNAseq", "WGCNA", "DCQ")
  
  # Perform KS test to see if there are differences in proportions
  Decon_WGCNA_distribution <- stats::ks.test(Decon_Mean, WGCNA_Mean)  
  Decon_DCQ_distribution <- stats::ks.test(Decon_Mean, DCQ_Mean)  
  DCQ_WGCNA_distribution <- stats::ks.test(WGCNA_Mean, DCQ_Mean) 
  
  Decon_WGCNA_cor <- stats::cor.test(Decon_Mean, WGCNA_Mean)  
  Decon_DCQ_cor <- stats::cor.test(Decon_Mean, DCQ_Mean)  
  DCQ_WGCNA_cor <- stats::cor.test(WGCNA_Mean, DCQ_Mean) 
  
  # Store it
  lD <- list(DeconRNAseq_WGCNA_ksTest = Decon_WGCNA_distribution,
             DeconRNAseq_DCQ_ksTest = Decon_DCQ_distribution,
             DCQ_WGCNA_ksTest = DCQ_WGCNA_distribution)
  
  lC <- list(DeconRNAseq_WGCNA_Pearson = Decon_WGCNA_cor,
             DeconRNAseq_DCQ_Pearson = Decon_DCQ_cor,
             DCQ_WGCNA_Pearson = DCQ_WGCNA_cor)
  

  
  # Plot est. cell-type proporitons in bar plots.
  forBarplot <- reshape::melt(t(do.call("rbind", l1)))
  colnames(forBarplot) <- c("CellType", "DeconvolutionMethod", "avgProportion")
  if(!is.null(order_celltype)) {
    forBarplot$CellType <- factor(forBarplot$CellType,levels = order_celltype)
    
  }
  p <- ggplot2::ggplot(forBarplot, ggplot2::aes(x = CellType, y = avgProportion, fill = DeconvolutionMethod, group=DeconvolutionMethod)) + 
    ggplot2::geom_bar(position = "dodge", stat = "identity") + ggplot2::theme_classic() + ggplot2::xlab("cell-type") + ggplot2::ylab("average cell-type proportion") +
    ggplot2::scale_fill_manual("DeconvolutionMethod", values = c("DCQ" = "black", "WGCNA" = "grey", "DeconRNAseq" = "darkblue")) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1), panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))
  if(print_plot) {
    print(p)
  }
  return(list(KS_test_avg_expression = lD,Pearson_cor_avg_espression = lC, cellType_proportions = proportions,average_celltype_proportions = l, barplot_to_print = p))
}
