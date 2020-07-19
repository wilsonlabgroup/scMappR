#' Pathway enrichment for cwFold-changes
#' 
#' This function runs through each list of cell weighted Fold changes (cwFold-changes) and completes both pathway and transcription factor (TF) enrichment.
#'
#' This function takes a matrix of cellWeighted_Foldchange and a species (human, mouse, or a character directly compatible with g:ProfileR).
#' Before completing pathway analysis with g:ProfileR. Enriched pathways are stored in a list and returned.
#'
#' @rdname gProfiler_cellWeighted_Foldchange
#' @name gProfiler_cellWeighted_Foldchange
#'
#' @param cellWeighted_Foldchange_matrix Matrix of cell weighted Fold changes from the deconvolute_and_contextualize functions.
#' @param species Human, mouse, or a charcter that is compatible with gProfileR.
#' @param background A list of background genes to test against.
#' @param gene_cut The top number of genes in pathway analysis.
#' @param newGprofiler Using gProfileR or gprofiler2, (T/F).
#' 
#' @return List with the following elements:
#' \item{BP}{gprofiler enrichment of biological pathways for each cell-type}
#' \item{TF}{gprofiler enrichment of transcription factors for eachc cell-type.}
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
#'
#' @examples 
#' \donttest{
#' 
#' data(PBMC_example)
#' 
#' bulk_DE_cors <- PBMC_example$bulk_DE_cors
#' bulk_normalized <- PBMC_example$bulk_normalized
#' odds_ratio_in <- PBMC_example$odds_ratio_in
#' 
#' case_grep <- "_female"
#' control_grep <- "_male"
#' max_proportion_change <- 10
#' print_plots <- FALSE
#' theSpecies <- "human"
#' norm <- deconvolute_and_contextualize(bulk_normalized, odds_ratio_in,
#'                                       bulk_DE_cors, case_grep = case_grep, 
#'                                       control_grep = control_grep,
#'                                        max_proportion_change = max_proportion_change,
#'                                        print_plots = print_plots,
#'                                        theSpecies = theSpecies)
#' background = rownames(bulk_normalized)
#' STVs <- gProfiler_cellWeighted_Foldchange(norm$cellWeighted_Foldchange, theSpecies,
#'  background, gene_cut = -9, newGprofiler = FALSE)
#' 
#'  }
#'  
#' @export
#' 
gProfiler_cellWeighted_Foldchange <- function(cellWeighted_Foldchange_matrix, species , background , gene_cut, newGprofiler ) {
  
  # Args: 
  # cellWeighted_Foldchange_matrix: matrix of cell weighted fold-changes from the deconvolute_and_contextualize functions
  # species: human, mouse, or a charcter that is compatible with gProfileR
  # background: a list of background genes to test against
  # Returns:
  # A List of significantly enriched pathways and TFs (correction_method = FDR, hier_sorting = moderate), for every cell-type
  
  
  if(!is.data.frame(cellWeighted_Foldchange_matrix)) {
    stop("cellWeighted_Foldchange_matrix much be a class == data.frame of scMappR Transformed Values from the deconvolute_and_contextualize function.")
  }
  if(!is.character(background)) {
    stop("background must be class == character of gene symbols as a background for gprofiler")
  }
  
  if(!is.numeric(gene_cut)) {
    stop("gene_cut must be a positive integer giving the numer of top genes for gprofiler or -9 if not used.")
  }
  if(!(is.logical(newGprofiler))) {
    stop("newGprofiler must be logical (TRUE/FALSE)")
  }
  
  if(!(species %in% c("human", "mouse"))) {
    if(species != -9) {
      stop("species is not 'human' 'mouse' or '-9' (case sensitive), please try again with this filled.")
    }
  }
  
  
  if(species == "human") {
    message("Assuming species = Human")
    theSpecies <- "hsapiens"
    
  }
  if(species == "mouse") {
    message("Assuming species = Mouse")
    theSpecies <- "mmusculus"
    
    
  }
  if(!(species %in% c("human", "mouse"))) {
    warning("We are assuming that you inputted 'species' as something directly compatible with g:ProfileR") 
    theSpecies <- species
    
  }
  
  gProfiler_internal <- function(x, theSpecies1 = theSpecies, background_genes = background, cutgenes = gene_cut, NewGprofiler = newGprofiler) {
    #Internal -- This function takes a signle STV and computes gProfiler BP and TF enrichment with the list re-ordered for each cell-type
    # Args:
    # X: a numeric index giving the cell-type in the STV matrix
    # theSpecies: the species that will be used for gProfiler
    # background_genes: the detected genes in your RNA-seq experiment (what you want to use as a background)
    # Returns:
    # significantly enriched BPs or TFs for at least one cell-type
    
    
    message(paste0("Re-ordering by absolute value of STVs on cell-type ", paste(colnames(cellWeighted_Foldchange_matrix)[x])))
    cellWeighted_Foldchange_matrix1 <- cellWeighted_Foldchange_matrix[order(abs(cellWeighted_Foldchange_matrix[,x]), decreasing = TRUE) ,]
    
    cellWeighted_Foldchange_matrix1 <- rownames(cellWeighted_Foldchange_matrix1)[abs(cellWeighted_Foldchange_matrix1[,x]) > 1e-10]
    if(cutgenes != -9) {
      message(paste0("Taking the top ", cutgenes, " genes for pathway analysis."))
      cellWeighted_Foldchange_matrix1 <- cellWeighted_Foldchange_matrix1[1:cutgenes]
    }
    if(NewGprofiler == FALSE) {
    ordered_back_all <- gProfileR::gprofiler(cellWeighted_Foldchange_matrix1, theSpecies1, ordered_query = TRUE, min_set_size = 3, max_set_size = 2000, src_filter = c("GO:BP", "REAC", "KEGG"),custom_bg = background_genes, correction_method = "fdr", min_isect_size = 3, hier_filtering = "moderate")
    ordered_back_all_tf <- gProfileR::gprofiler(cellWeighted_Foldchange_matrix1, theSpecies1, ordered_query = TRUE, min_set_size = 3, max_set_size = 5000, src_filter = c("TF"),custom_bg = background_genes, correction_method = "fdr", min_isect_size = 3, hier_filtering = "moderate")
    ordered_back_all$term_name <- ordered_back_all$term.name
    ordered_back_all$p_value <- ordered_back_all$p.value
    ordered_back_all_tf$term_name <- ordered_back_all_tf$term.name
    ordered_back_all_tf$p_value <- ordered_back_all_tf$p.value
    } else {
    ordered_back_all <- gprofiler2::gost(query = cellWeighted_Foldchange_matrix1, organism = theSpecies1, ordered_query = TRUE, significant = TRUE, exclude_iea = FALSE, multi_query = FALSE, measure_underrepresentation = FALSE, evcodes = FALSE, user_threshold = 0.05, correction_method = "fdr",  custom_bg =background_genes, numeric_ns = "", sources = c("GO:BP", "KEGG", "REAC"))  
    if(is.null(ordered_back_all)) { #if nothing is significant, gost returns null. making compatible with plotBP and make_TF_barplot
      ordered_back_all <- as.data.frame(matrix(c(1," "),nrow = 1))
      colnames(ordered_back_all) <- c("p_value", "term_name")
      ordered_back_all <- ordered_back_all[-1,]
    } else {
      ordered_back_all <- ordered_back_all$result
      ordered_back_all <- ordered_back_all[ordered_back_all$term_size > 15 & ordered_back_all$term_size < 2000 & ordered_back_all$intersection_size > 2,]
    }
    ordered_back_all_tf <- gprofiler2::gost(query = cellWeighted_Foldchange_matrix1, organism = theSpecies1, ordered_query = TRUE, significant = TRUE, exclude_iea = FALSE, multi_query = FALSE, measure_underrepresentation = FALSE, evcodes = FALSE, user_threshold = 0.05, correction_method = "fdr",  custom_bg =background_genes, numeric_ns = "", sources = c("TF"))  
    if(is.null(ordered_back_all_tf)) {
      ordered_back_all_tf <- as.data.frame(matrix(c(1," "),nrow = 1))
      colnames(ordered_back_all_tf) <- c("p_value", "term_name")
      ordered_back_all_tf <- ordered_back_all_tf[-1,]
      
    } else {
      ordered_back_all_tf <- ordered_back_all_tf$result
      ordered_back_all_tf <- ordered_back_all_tf[ordered_back_all_tf$term_size > 15 & ordered_back_all_tf$term_size < 5000 & ordered_back_all_tf$intersection_size > 2,]
    }
    }
    
    return(list(BPs = ordered_back_all, TFs = ordered_back_all_tf))
    
  }
  message(paste(colnames(cellWeighted_Foldchange_matrix), " "))
  message(theSpecies)
  paths <- lapply(1:ncol(cellWeighted_Foldchange_matrix), gProfiler_internal)
  
  BPs <- TFs <- list()
  # Converting BPs and TFs for every cell-type into a list
  for(i in 1:ncol(cellWeighted_Foldchange_matrix)) {
    BPs[[i]] <- paths[[i]]$BPs
    TFs[[i]] <- paths[[i]]$TFs
  }
  names(BPs) <- names(TFs) <- colnames(cellWeighted_Foldchange_matrix)
  return(list(BP = BPs, TF = TFs))
  
}
