#' Measure cell-type specificity of cell-weighted Fold-changes
#' 
#' This function normalizes cwFold-changes by each gene to help visualize the cell-type specificity of DEGs. It then tests if a cell-type has a large change in correlation from bulk DEGs. Finally, it identifies genes that may be specific to each cell-type.
#'
#' cwFold-changes and re-normalized and re-processed to interrogate cell-type specificity at the level of the cell-type and at the level of the gene.
#' At the level of the cell-type, cwFold-changes are correlated to bulk DEGs. The difference in rank between bulk DEGs and cwFold-changes are also compared.
#' At the level of the gene, cwFold-changes are re-normalized so that each gene sums to 1. Normalization of their distributions are tested with a Shapiro test. 
#' Then, outlier cell-types for each gene are measured by testing for `sd_cutoff`'s mad or sd's greater than the median or mean depending on if the cwFold-change is non-normally or normally distributed respectively.
#' Cell-types considered outliers are then further filtered so their normalized cwFold-changes are greater than the cell-type proportions of that gene and `gene_cutoff` if the user sets it.
#' 
#' @rdname cwFoldChange_evaluate
#' @name cwFoldChange_evaluate
#' 
#' @param cwFC A matrix or data frame of cell-weighted fold-changes of DEGs. Rows are DEGs and columns are cell-types. 
#' @param celltype_prop A matrix or data frame of cell-type proportions. Rows are different cell-types and columns are different samples. These cell-type proportions can come from any source (not just scMappR).
#' @param DEG_list An object with the first column as gene symbols within the bulk dataset (doesn't have to be in signature matrix), second column is the adjusted p-value, and the third the log2FC path to a .tsv file containing this info is also acceptable.
#' @param gene_cutoff Additional cut-off of normalized cwFold-change to see if a gene is cut-off.
#' @param sd_cutoff Number of standard deviations or median absolute deviations to calculate outliers.
#' 
#' @return List with the following elements:
#' \item{gene_level_investigation}{data frame of genes showing the Euclidian distances between cwFold-change and null vector as well as if cwFold-changes are distributed.}
#' \item{celltype_level_investigation}{data frame of Spearman's and Pearson's correlation between bulk DEGs and cwFold-changes.}
#' \item{cwFoldchange_vs_bulk_rank_change}{data frame of the change in rank of DEG between the bulk fold-change and cwFold-change.}
#' \item{cwFoldChange_normalized}{cwFold-change normalized such that each gene sums to 1.}
#' \item{cwFoldchange_gene_assigned}{List of cell-types where genes are designated to cell-type specific differential expression.}
#' \item{cwFoldchange_gene_flagged_FP}{Mapped cwFoldchanges that are flagged as false-positives. These are genes that are driven by the reciprical ratio of cell-type proportions between case and control. These genes may be DE in a non-cell-type specific manner but are falsely assigned to cell-types with very large differences in proportion between condition.} 
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
#' 
#' data(PBMC_example)
#' bulk_DE_cors <- PBMC_example$bulk_DE_cors
#' bulk_normalized <- PBMC_example$bulk_normalized
#' odds_ratio_in <- PBMC_example$odds_ratio_in
#' case_grep <- "_female"
#' control_grep <- "_male"
#' max_proportion_change <- 10
#' print_plots <- FALSE
#' theSpecies <- "human"
#' toOut <- scMappR_and_pathway_analysis(bulk_normalized, odds_ratio_in, 
#'                                       bulk_DE_cors, case_grep = case_grep,
#'                                       control_grep = control_grep, rda_path = "", 
#'                                       max_proportion_change = 10, print_plots = TRUE, 
#'                                        plot_names = "tst1", theSpecies = "human", 
#'                                        output_directory = "tester",
#'                                        sig_matrix_size = 3000, up_and_downregulated = FALSE, 
#'                                        internet = FALSE)
#'                                        
#' cwFC1 <- toOut$cellWeighted_Foldchange
#' prop1 <- toOut$cellType_Proportions
#' DE <- bulk_DE_cors
#' eval_test <- cwFoldChange_evaluate(cwFC1, prop1, DE)
#' 
#' @export
#' 

cwFoldChange_evaluate <- function(cwFC, celltype_prop, DEG_list, gene_cutoff = NULL, sd_cutoff = 3) {

  # Pre-processing DEG list
  DEG_list_class <- class(DEG_list)[1] %in% c("character", "data.frame", "matrix") 
  

  if(DEG_list_class[1] == FALSE) {
    stop("DEG_list must be of class character, data.frame, or matrix.")
  }
  if(is.character(DEG_list)) {
    DEGs <- utils::read.table(DEG_list, header = FALSE, as.is = TRUE, sep = "\t")
  } else {
    DEGs <- as.data.frame(DEG_list)
  }
  
  colnames(DEGs) <- c("gene_name", "padj", "log2fc")
  DEGs$gene_name <- tochr(DEGs$gene_name)
  DEGs$padj <- toNum(DEGs$padj)
  DEGs$log2fc <- toNum(DEGs$log2fc)
  rownames(DEGs) <- DEGs$gene_name
  DEGs$FoldChange <- 2^DEGs$log2fc * sign(DEGs$log2fc)
  
  # Checking cell-type proportions and cwFC's are inputted correctly.
  cwFC_class <- class(cwFC)[1] %in% c("matrix", "data.frame")
  celltype_prop_class <- class(celltype_prop)[1] %in% c("matrix", "data.frame")
  allGood <- all(cwFC_class, celltype_prop_class)
  if(!allGood) stop("Both celltype_prop and cwFC_class must be of class 'matrix' or 'data.frame'.")

  # Overlapping genes and cell-types.
  inters <- intersect(rownames(cwFC), rownames(DEGs))
  intDEG <- any(length(inters) < nrow(DEGs), length(inters) < nrow(cwFC))
  if(length(inters) == 0 ) stop("none of your cwFC's overlap with DEGs. Check if genes are the same (case sensitive, ensembl vs gene symbol etc.)")
  if(intDEG) warning("Not all genes overlap between DEGs and cwFoldChanges. Keeping overlapping genes but this may be an issue with inputted data.")
  if(intDEG) message("Not all genes overlap between DEGs and cwFoldChanges. Keeping overlapping genes but this may be an issue with inputted data.")
  cwFC <- cwFC[inters,]
  DEGs <- DEGs[inters,]
  
  # Overlapping cell-types in proportions and cwFold-changes.
  intersProp <- intersect(colnames(cwFC), colnames(celltype_prop))
  intersPropT <- any(length(intersProp) < ncol(cwFC), length(intersProp) < ncol(celltype_prop))
  if(length(intersProp) == 0 ) stop("none of your cell-types in cwFC's overlap with the cell-types in the cell-type proportions. Check if genes are the same (case sensitive).")
  if(intersPropT) warning("Not all cell-types overlap between the cell-weighted fold-changes and cell-types in the proportions.")
  if(intersPropT) message("Not all cell-types overlap between the cell-weighted fold-changes and cell-types in the proportions.")
  cwFC <- cwFC[,intersProp]
  celltype_prop <- celltype_prop[,intersProp]
  
  # Average cell-type proporitons 
  props <- colMeans(celltype_prop)
  cwFC_rownorm <- function(x) {
    y <- (x / sum(x)) * sign(x)
    return(y)
  }
  
  # normalize cwFC's to suum to 1.
  cwFC_norm <- t(apply(X = cwFC, 1, cwFC_rownorm))
  
  # distance from normalized cwFold-changes to 1 
  oneRow <- rep(1, ncol(cwFC_norm))
  distsNorm <- function(x) {
    y <- stats::dist(t(cbind(x,oneRow)))[[1]]
    
    return(y)
  }
  
  
  distVec <- apply(cwFC_norm, 1,distsNorm)
  
  # Test if normalized cwFold-changes fit a normal distribution
  shapiro_split <- function(x) {
    y <- stats::shapiro.test(x)
    w <- unname(y$statistic)
    p <- unname(y$p.value)
    return(c(w,p))
  }
  
  norm <-  t(apply(cwFC_norm,1,shapiro_split))
  colnames(norm) <- c("W_shapiro", "p_shapiro")
  
  
  
  # testing the normal distribution of cell-types.
  df <- as.data.frame(cbind(distVec, norm), stringasFactors = FALSE)  
  df$p_shapiro_adjusted <- stats::p.adjust(df$p_shapiro, method = "fdr")
  colnames(df) = c("Distance_from_1", "W_shapirotest", "p_shapirotest",  "p_adj_shapirotest")
  
  # Testing if cell-types are highly correlated between cwFold-changes of a cell-type and bulk DEGs.
  column_tests <- function(x) {
    spear <- stats::cor.test(x, DEGs$FoldChange, method = "spearman")
    rho <- unname(spear$estimate)
    p_spear <- unname(spear$p.value)
    pear <- stats::cor.test(x, DEGs$FoldChange, method = "pearson")
    R_squared <- unname(pear$estimate)
    p_pear <- unname(pear$p.value)
    return(c(rho, p_spear, R_squared, p_pear))
  }
  df_col <- apply(cwFC, 2, column_tests)
  rownames(df_col) <- c("rho", "p_spear", "R_squared", "p_pear")
  ###### Looking at the rank order change in DEG
  rank_change_list <- list()
  abscwFC <- abs(cwFC)
  DEG_names <- rownames(abscwFC)
  # look at the difference in rank between cwFold-change and bulk DEG.
  for (i in 1:ncol(cwFC)) {
    abscwFC <- abscwFC[order(abscwFC[,i], decreasing = TRUE), ]
    
    rank_CT <- 1:length(rownames(abscwFC))
    names(rank_CT) <- rownames(abscwFC)
    DEGs$abs <- abs(DEGs$log2fc)
    DEGs1 <- DEGs[order(DEGs$abs, decreasing = TRUE),]
    
    rank_DEGs <- 1:length(DEGs1$gene_name)
    names(rank_CT) <- rownames(abscwFC)
    names(rank_DEGs) <- rownames(DEGs1)
    
    rank_CT <- rank_CT[ DEG_names]
    rank_DEGs <- rank_DEGs[ DEG_names]
    rank_change_list[[i]] <- rank_CT 
  }
  rank_change_df <- as.data.frame(do.call("cbind", rank_change_list), stringasFactors= FALSE)
  colnames(rank_change_df) <- colnames(cwFC)
  rank_change_df$bulk <-  rank_DEGs
  gene_level_investigation <- df
  celltype_level_investigation <- df_col
  cwFoldchange_vs_bulk_rank_change <- rank_change_df
  cwFoldchange_vs_bulk_rank_change <- cwFoldchange_vs_bulk_rank_change[order(cwFoldchange_vs_bulk_rank_change$bulk, decreasing = FALSE),]
  cwFoldChange_normalized <- cwFC_norm
  abs_cwFC_norm <- abs(cwFC_norm)
  
  # Get normal and non-normal DEGs.
  normal_genes <- rownames(df)[df$p_adj_shapirotest >= 0.05]
  abnormal_genes <- rownames(df)[df$p_adj_shapirotest < 0.05]
  #########

  # outliers for non-normal DEGs  
  getOutliersMedian <- function(outlierName) {
   
    outlierTest <- abs_cwFC_norm[outlierName,]
    
    lower_bound <- stats::median(outlierTest) - sd_cutoff * stats::mad(outlierTest)
    upper_bound <- stats::median(outlierTest) + sd_cutoff * stats::mad(outlierTest)
    
    outlier_ind_median <- which( outlierTest > upper_bound & outlierTest > props)
    
    return( outlier_median = outlier_ind_median)
  }
  
  # outliers for normal DEGs
  getOutliersMean <- function(outlierName) {
    outlierTest <- abs_cwFC_norm[outlierName,]
    lower_bound <- mean(outlierTest) - sd_cutoff * stats::sd(outlierTest)
    upper_bound <- mean(outlierTest) + sd_cutoff * stats::sd(outlierTest)
    
    outlier_ind_mean <- which( outlierTest > upper_bound & outlierTest > props)
    
    
    return(outlier_mean = outlier_ind_mean)
  }
  
  # get and combine genes that fit a normal or non-normal distribution.
  if(length(normal_genes) > 0 & length(abnormal_genes) > 0 ) {
    meanOutlier <- lapply(normal_genes, getOutliersMean)
    medianOutlier <- lapply(abnormal_genes, getOutliersMedian)
    names(meanOutlier) <- normal_genes
    names(medianOutlier) <- abnormal_genes
    outlier <- c(meanOutlier, medianOutlier)
  }
  # combine things if there is an outlier or not
  if(length(normal_genes) > 0 & length(abnormal_genes) == 0 ) {
    message("All cwFold-changes follow a normal distribution within the gene. This suggests low cell-type specificity")
    meanOutlier <- lapply(normal_genes, getOutliersMean)
    names(meanOutlier) <- normal_genes
    outlier <- c(meanOutlier)
  }
  if(length(normal_genes) == 0 & length(abnormal_genes) > 0 ) {
    message("All DEGs are not normally distributed, suggesting very consistent cell-type specificity across genes.")
    meanOutlier <- lapply(normal_genes, getOutliersMean)
    names(meanOutlier) <- normal_genes
    outlier <- c(meanOutlier)
  }
  
  outlier <- outlier[lengths(outlier) >0 ]
  if(length(outlier) == 0) {
    positiveList <- vector("list", length(names(props)))
    names(positiveList) <- names(props)
    cwFoldchange_gene_assigned <- cwFoldchange_gene_flagged_FP <- list()
  } else {
    if(is.null(gene_cutoff)) {
      gene_cutoff <- 1/ncol(cwFC_norm)
    }
    # Get the cell-types that have DEGs whose DE may be cell-type specific.
    positiveList <- list()
    for(p in names(outlier)) {
      val <- outlier[[p]]
      for(z in names(val)) {
        
        fc <- cwFoldChange_normalized[p,z]
        
        if(abs(fc) > gene_cutoff) {
          names(fc) <- p
          positiveList[[z]] <- c(positiveList[[z]], fc)
        }
        
        
      }
    }
    cwFoldchange_gene_assigned <- list()
    cwFoldchange_gene_flagged_FP <- list()
    for(p in names(positiveList)) {
      g1 <- positiveList[[p]]
      g1Abs <- abs(g1)
      nonPerfect <- g1[g1Abs > 0 & g1Abs < 0.999]
      Perfect <- g1[ g1Abs > 0.999 ]
      nonP_sorted <- sort(table(nonPerfect),decreasing=TRUE)
      if(!unname(is.na(nonP_sorted[2]))) {
        if(unname(nonP_sorted[1]) > (unname(nonP_sorted[2])*5)) {
          nonPerfect1 <- nonPerfect[round(nonPerfect,5) == round(as.numeric(names(nonP_sorted)[1]),5)]
          nonPerfect2 <- nonPerfect[round(nonPerfect,5) != round(as.numeric(names(nonP_sorted)[1]),5)] 
          if(length(Perfect) > 0) {
            nonPerfect2 <- c(Perfect, nonPerfect2)
          }
          cwFoldchange_gene_assigned[[p]] <- nonPerfect2
          cwFoldchange_gene_flagged_FP[[p]] <- nonPerfect1
        } else {
          cwFoldchange_gene_assigned[[p]] <- g1
          cwFoldchange_gene_flagged_FP[[p]] <- numeric(0)
        }
        
      }
      
    }
    
    
    
  }

  return(list(gene_level_investigation = gene_level_investigation, celltype_level_investigation = celltype_level_investigation, cwFoldchange_vs_bulk_rank_change= cwFoldchange_vs_bulk_rank_change, cwFoldChange_normalized = cwFoldChange_normalized, cwFoldchange_gene_assigned = cwFoldchange_gene_assigned, cwFoldchange_gene_flagged_FP=cwFoldchange_gene_flagged_FP)) 
}
